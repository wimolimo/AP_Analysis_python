"""Temperature-dependence traces, interval means, and kinetic relationship plots.

All datetimes in this module are UTC. CLOUD18 timestamps are already UTC and
are therefore kept as timezone-naive ``datetime`` values throughout.
"""

from __future__ import annotations

from datetime import datetime
import re
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.odr import Model, ODR, RealData
from scipy.stats import t as student_t

from . import config as cfg
from .loading import load_data, read_trace_timeseries
from .models import Channel
from .plotting import _draw_stages
from .stages import load_stages

TimeRange = tuple[datetime, datetime]


def _validate_time_ranges(ranges: Sequence[TimeRange]) -> list[TimeRange]:
    """Validate timezone-naive UTC intervals and return them as a list."""
    validated: list[TimeRange] = []
    for start, stop in ranges:
        if not isinstance(start, datetime) or not isinstance(stop, datetime):
            raise TypeError("Time ranges must contain (datetime, datetime) pairs.")
        if start.tzinfo is not None or stop.tzinfo is not None:
            raise ValueError("Use timezone-naive datetime objects; CLOUD18 times are UTC.")
        if stop < start:
            raise ValueError("A time interval has end < start.")
        validated.append((start, stop))
    return validated


def _chemical_label(formula: str) -> str:
    """Return a Matplotlib label with chemical subscripts and ionic charge."""
    def format_part(part: str) -> str:
        charge = part[-1] if part.endswith(("+", "-")) else ""
        if charge:
            part = part[:-1]
        formatted = re.sub(r"([A-Z][a-z]?)(\d+)", r"\1$_{\2}$", part)
        return formatted + (rf"$^{{{charge}}}$" if charge else "")

    parts = str(formula).strip().split(".", 1)
    label = format_part(parts[0])
    return label if len(parts) == 1 else label + r"$\cdot$" + format_part(parts[1])


def _style_axis(ax: plt.Axes) -> None:
    """Apply the shared publication style to one Matplotlib axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(direction="out")


def _interval_label(index: int) -> str:
    return rf"$t_{{{int(index)}}}$"


def _shade_intervals(
    ax: plt.Axes, ranges: Sequence[TimeRange], label_intervals: bool = False,
) -> None:
    """Shade averaging intervals and optionally label them with matching tags."""
    for index, (start, stop) in enumerate(ranges, start=1):
        ax.axvspan(
            start,
            stop,
            color=cfg.IMPORTANT_RANGE_COLOR,
            alpha=cfg.IMPORTANT_RANGE_ALPHA,
            linewidth=0,
            zorder=0,
        )
        if label_intervals:
            ax.text(
                start + (stop - start) / 2, 0.85, _interval_label(index),
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=cfg.PAPER_LABEL_SIZE,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1),
                clip_on=True, zorder=5,
            )


def calculate_interval_means(
    data: pd.DataFrame,
    compounds: Sequence[str],
    time_col: str,
    important_time_ranges: Sequence[TimeRange],
) -> pd.DataFrame:
    """Calculate one AP-trace mean per compound and averaging interval.

    Returns one row per interval with ``range_index``, ``start``, ``end``,
    ``n_points`` and, for each requested compound, its mean, number of valid
    observations and standard error of the mean (``<compound>_sem``).
    """
    ranges = _validate_time_ranges(important_time_ranges)
    if not ranges:
        raise ValueError("important_time_ranges must not be empty.")

    time_values = pd.to_datetime(data[time_col], errors="coerce")
    rows: list[dict] = []

    for index, (start, stop) in enumerate(ranges, start=1):
        selected = data.loc[time_values.between(start, stop), compounds]
        row = {
            "range_index": index,
            "start": start,
            "end": stop,
            "n_points": len(selected),
        }
        for compound in compounds:
            values = pd.to_numeric(selected[compound], errors="coerce").dropna()
            row[compound] = values.mean()
            row[f"{compound}_n_points"] = int(values.size)
            row[f"{compound}_sem"] = (
                float(values.std(ddof=1) / np.sqrt(values.size))
                if values.size >= 2
                else np.nan
            )
        if selected.empty:
            print(f"Warning: T{index} contains no AP trace points: {start} – {stop}")
        rows.append(row)

    return pd.DataFrame(rows)


def calculate_channel_interval_means(
    channel: Channel,
    important_time_ranges: Sequence[TimeRange],
    column_name: Optional[str] = None,
) -> pd.DataFrame:
    """Calculate interval means for a standard CLOUD18 ``Channel``.

    ``Channel.time`` is stored as Unix milliseconds by the general loader.
    The same UTC T1, T2, ... intervals used for the AP traces are applied
    directly to those timestamps.

    Parameters
    ----------
    channel:
        Loaded time series, for example the HORUS ``HO2_ppt`` channel.
    important_time_ranges:
        UTC averaging intervals T1, T2, ... .
    column_name:
        Name of the returned mean column. Defaults to ``channel.name``.

    Returns
    -------
    pandas.DataFrame
        Columns ``range_index``, ``<column_name>``, ``<column_name>_n_points``
        and ``<column_name>_sem``. The latter is the standard error of the
        interval mean.
    """
    ranges = _validate_time_ranges(important_time_ranges)
    if not ranges:
        raise ValueError("important_time_ranges must not be empty.")

    name = column_name or channel.name
    times = pd.to_datetime(channel.time, unit="ms", errors="coerce")
    values = pd.to_numeric(pd.Series(channel.values), errors="coerce").to_numpy()

    rows: list[dict] = []
    for index, (start, stop) in enumerate(ranges, start=1):
        mask = np.asarray((times >= start) & (times <= stop))
        selected = values[mask]
        finite = selected[np.isfinite(selected)]

        rows.append({
            "range_index": index,
            name: float(np.mean(finite)) if finite.size else np.nan,
            f"{name}_n_points": int(finite.size),
            f"{name}_sem": (
                float(np.std(finite, ddof=1) / np.sqrt(finite.size))
                if finite.size >= 2
                else np.nan
            ),
        })

        if finite.size == 0:
            print(f"Warning: T{index} contains no valid {name} points: {start} – {stop}")

    return pd.DataFrame(rows)


def add_reaction_partner_interval_means(
    means: pd.DataFrame,
    horus_path: str,
    time_range: TimeRange,
    important_time_ranges: Sequence[TimeRange],
    precursor_compound: str,
    reaction_partner: str,
) -> tuple[pd.DataFrame, str]:
    """Add reaction-partner means and the precursor × partner kinetic quantity.

    Currently ``reaction_partner="HO2"`` is implemented. The selection is kept
    explicit so additional partners such as O3 can later be added without
    changing the public ``temp_dep`` interface.

    Parameters
    ----------
    means:
        AP interval-mean table.
    horus_path:
        Path to the HORUS HOx/ROx file used when the partner is HO2.
    time_range:
        Overall UTC loading interval.
    important_time_ranges:
        UTC averaging intervals T1, T2, ... .
    precursor_compound:
        AP precursor column to multiply by the partner mean.
    reaction_partner:
        Chemical reaction partner. Currently only ``"HO2"`` is supported.

    Returns
    -------
    means_with_partner, kinetic_x_column
        Updated interval-mean table and the name of the derived kinetic
        x-column.
    """
    partner = reaction_partner.strip().upper()

    if partner == "HO2":
        partner_channel = "HO2_ppt"
        partner_label = "HO2"
        partner_data = load_data(horus_path, time_range)
    else:
        raise NotImplementedError(
            f"Reaction partner {reaction_partner!r} is not implemented yet. "
            "Currently supported: 'HO2'."
        )

    if partner_channel not in partner_data.channels:
        available = ", ".join(sorted(partner_data.channels))
        raise KeyError(
            f"Channel {partner_channel!r} for reaction partner {partner_label} "
            f"was not found. Available channels: {available}"
        )

    partner_means = calculate_channel_interval_means(
        partner_data[partner_channel],
        important_time_ranges,
        column_name=partner_channel,
    )

    out = means.merge(
        partner_means,
        on="range_index",
        how="left",
        validate="one_to_one",
    )
    kinetic_x_column = f"{precursor_compound} * {partner_channel}"
    out[kinetic_x_column] = out[precursor_compound] * out[partner_channel]
    kinetic_x_sem_column = f"{kinetic_x_column}_sem"

    # First-order uncertainty propagation for a product x = a*b. This form is
    # also well behaved if either mean happens to be zero. It assumes that the
    # interval-mean errors of precursor and reaction partner are independent.
    out[kinetic_x_sem_column] = np.sqrt(
        (out[partner_channel] * out[f"{precursor_compound}_sem"]) ** 2
        + (out[precursor_compound] * out[f"{partner_channel}_sem"]) ** 2
    )

    return out, kinetic_x_column


def plot_kinetic_relationship(
    grouped_means: pd.DataFrame,
    x_column: str,
    y_compound: str,
    x_label: Optional[str] = None,
    show: bool = True,
    savepath: Optional[str] = None,
    label_intervals: bool = False,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Plot interval means with errors and fit zero-intercept ODR models.

    Parameters
    ----------
    grouped_means:
        Combined interval-mean table. It must contain a ``temperature`` column,
        the kinetic x quantity and the product compound.
    x_column:
        Column containing precursor × reaction-partner values.
    y_compound:
        Product compound plotted on the y-axis.
    x_label:
        Optional formatted x-axis label.
    show:
        Show the figure with Matplotlib.
    savepath:
        Optional output path for the figure.
    label_intervals:
        Label points using their original range_index, matching trace intervals.
        Labels use the temperature series color.

    Returns
    -------
    figure, axis, regression_results
        The combined kinetic plot and a table containing the zero-intercept
        ODR slope, its standard error, goodness-of-fit values and the number of
        points for each temperature.
    """
    x_error_column = f"{x_column}_sem"
    y_error_column = f"{y_compound}_sem"
    required = [
        "temperature",
        x_column,
        y_compound,
        x_error_column,
        y_error_column,
    ]
    if label_intervals:
        required.append("range_index")
    missing = [column for column in required if column not in grouped_means.columns]
    if missing:
        raise ValueError("Missing kinetic columns: " + ", ".join(missing))

    plt.rcParams.update(cfg.PAPER_RC_PARAMS)
    fig, ax = plt.subplots(figsize=cfg.KINETIC_FIGSIZE)
    regression_rows: list[dict] = []

    temperatures = list(pd.unique(grouped_means["temperature"]))
    colors = cfg.PAPER_LINE_COLORS

    for index, temperature in enumerate(temperatures):
        valid = grouped_means.loc[
            grouped_means["temperature"] == temperature,
            [x_column, y_compound, x_error_column, y_error_column]
            + (["range_index"] if label_intervals else []),
        ].replace([np.inf, -np.inf], np.nan).dropna()

        # ODR requires strictly positive standard deviations. A zero SEM means
        # that the interval provides no usable estimate of measurement spread.
        valid = valid.loc[
            (valid[x_error_column] > 0) & (valid[y_error_column] > 0)
        ]

        if valid.empty:
            print(f"Warning: no valid kinetic points for {temperature}.")
            continue

        color = colors[index % len(colors)]
        x = valid[x_column].to_numpy(dtype=float)
        y = valid[y_compound].to_numpy(dtype=float)
        x_error = valid[x_error_column].to_numpy(dtype=float)
        y_error = valid[y_error_column].to_numpy(dtype=float)

        ax.errorbar(
            x,
            y,
            xerr=x_error,
            yerr=y_error,
            fmt="o",
            markersize=np.sqrt(cfg.KINETIC_MARKER_SIZE),
            color=color,
            markeredgewidth=0,
            elinewidth=cfg.PAPER_LINE_WIDTH,
            capsize=2,
            label=str(temperature),
            zorder=3,
        )

        if label_intervals:
            for xi, yi, interval_index in zip(x, y, valid["range_index"]):
                ax.annotate(
                    _interval_label(interval_index), (xi, yi),
                    xytext=(5, 5), textcoords="offset points",
                    color=color, fontsize=cfg.PAPER_LABEL_SIZE,
                    ha="left", va="bottom", zorder=5,
                )

        if len(valid) < 2 or np.allclose(x, 0.0):
            print(
                f"Warning: at least two points and a non-zero x value are "
                f"required for an ODR fit at {temperature}."
            )
            continue

        def zero_intercept_model(beta: np.ndarray, x_values: np.ndarray) -> np.ndarray:
            return beta[0] * x_values

        initial_slope = float(np.dot(x, y) / np.dot(x, x))
        odr_result = ODR(
            RealData(x, y, sx=x_error, sy=y_error),
            Model(zero_intercept_model),
            beta0=[initial_slope],
        ).run()

        slope = float(odr_result.beta[0])
        slope_stderr = float(odr_result.sd_beta[0])
        x_fit = np.linspace(0.0, x.max() * 1.05, 200)
        y_fit = slope * x_fit

        ax.plot(
            x_fit,
            y_fit,
            color=color,
            linewidth=cfg.PAPER_LINE_WIDTH,
            zorder=2,
        )

        fitted_y = slope * x
        ss_res = float(np.sum((y - fitted_y) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        centered_r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        uncentered_ss_tot = float(np.sum(y ** 2))
        uncentered_r_squared = (
            1.0 - ss_res / uncentered_ss_tot
            if uncentered_ss_tot > 0
            else np.nan
        )
        degrees_of_freedom = len(valid) - 1
        p_value = (
            float(2.0 * student_t.sf(abs(slope / slope_stderr), degrees_of_freedom))
            if degrees_of_freedom > 0 and slope_stderr > 0
            else np.nan
        )

        regression_rows.append({
            "temperature": temperature,
            "slope": slope,
            "intercept": 0.0,
            "r_squared": centered_r_squared,
            "r_squared_uncentered": uncentered_r_squared,
            "p_value": p_value,
            "slope_stderr": slope_stderr,
            "residual_variance": float(odr_result.res_var),
            "sum_squared_weighted_residuals": float(odr_result.sum_square),
            "odr_stop_reason": "; ".join(odr_result.stopreason),
            "n_points": len(valid),
        })

    ax.set_xlabel(x_label or x_column)
    ax.set_ylabel(_chemical_label(y_compound) + " (ppt)")
    ax.legend(title="Temperature")
    _style_axis(ax)
    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
    if show:
        plt.show()

    return fig, ax, pd.DataFrame(regression_rows)


def _prepare_trace_compounds(
    compounds: Sequence[str],
    adduct: str,
) -> list[str]:
    """Return exact AP trace names, appending the default adduct if necessary."""
    normalized_adduct = (
        adduct if not adduct or adduct.startswith(".") else "." + adduct
    )
    return [
        name if "." in name else f"{name}{normalized_adduct}"
        for name in (str(compound).strip() for compound in compounds)
    ]


def _plot_temperature_trace(
    data: pd.DataFrame,
    trace_compounds: Sequence[str],
    time_col: str,
    time_range: TimeRange,
    important_time_ranges: Sequence[TimeRange],
    temperature: str,
    overlay: bool,
    show: bool,
    savepath: Optional[str],
    stages: Optional[str] = "onlytype",
    label_intervals: bool = False,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Create the time-trace figure for one temperature condition."""
    plot_time = pd.to_datetime(data[time_col], errors="coerce")

    if overlay:
        fig, ax = plt.subplots(figsize=cfg.TRACE_FIGSIZE)
        axes = [ax]
        _shade_intervals(ax, important_time_ranges, label_intervals=label_intervals)
        for compound in trace_compounds:
            ax.plot(
                plot_time,
                data[compound],
                label=_chemical_label(compound),
                linewidth=cfg.PAPER_LINE_WIDTH,
            )
        ax.set_ylabel("Signal (ppt)")
        ax.set_xlabel("Time (UTC)")
        ax.set_yscale("log")
        ax.legend(handlelength=1.5, borderaxespad=0.2)
        _style_axis(ax)
    else:
        fig, axes_array = plt.subplots(
            len(trace_compounds),
            1,
            sharex=True,
            squeeze=False,
            figsize=(cfg.TRACE_FIGSIZE[0], max(1.8 * len(trace_compounds), 2.5)),
        )
        axes = list(axes_array[:, 0])
        for ax, compound in zip(axes, trace_compounds):
            _shade_intervals(ax, important_time_ranges, label_intervals=label_intervals)
            ax.plot(plot_time, data[compound], linewidth=cfg.PAPER_LINE_WIDTH)
            ax.set_ylabel("Signal (ppt)")
            ax.set_yscale("log")
            ax.text(
                0.01,
                0.95,
                _chemical_label(compound),
                transform=ax.transAxes,
                ha="left",
                va="top",
            )
            _style_axis(ax)
        axes[-1].set_xlabel("Time (UTC)")

    for ax in axes:
        ax.set_xlim(*time_range)
        if stages is not None:
            _draw_stages(ax, stages, xlims=time_range)

    # Keep the temperature visible in the upper-right corner.
    axes[0].text(
        0.99,
        0.97,
        str(temperature),
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=cfg.PAPER_LABEL_SIZE,
    )

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes


def temp_dep(
    temperature_runs: Sequence[dict],
    compounds: Sequence[str],
    horus_path: str,
    precursor_compound: str,
    product_compound: str,
    reaction_partner: str = "HO2",
    adduct: str = ".NH3H+",
    show: bool = True,
    overlay: bool = True,
    kinetic_savepath: Optional[str] = None,
    stages_path: Optional[str] = None,
    stages: Optional[str] = "onlytype",
    label_intervals: bool = False,
) -> tuple[
    list[tuple[plt.Figure, list[plt.Axes]]],
    pd.DataFrame,
    pd.DataFrame,
    plt.Figure,
    plt.Axes,
]:
    """Analyse several temperatures and compare their kinetic slopes.

    Each entry in ``temperature_runs`` describes one temperature condition::

        {
            "temperature": "-50 °C",
            "trace_path": "...csv",
            "time_range": (start, stop),
            "important_time_ranges": [(start1, stop1), ...],
            "trace_savepath": "output/...png",  # optional
        }

    For every condition, AP compounds and the selected reaction partner are
    averaged over that condition's important intervals. The kinetic quantity is

        x = mean(precursor) * mean(reaction partner)
        y = mean(product)

    A separate time-trace figure is produced for each temperature. All kinetic
    points are then combined into one figure, including x/y standard-error bars
    and a separate zero-intercept orthogonal-distance regression (``y = k*x``)
    for each temperature.

    All datetimes are interpreted as timezone-naive UTC.

    Parameters
    ----------
    temperature_runs:
        Sequence of dictionaries defining temperature label, AP trace path,
        overall UTC time range, important averaging intervals and optionally a
        time-trace save path.
    compounds:
        AP compounds to load and display.
    horus_path:
        Path to the HORUS file. Used for ``reaction_partner="HO2"``.
    precursor_compound:
        AP precursor multiplied by the reaction-partner mean.
    product_compound:
        AP product plotted on the kinetic y-axis.
    reaction_partner:
        Reaction partner multiplied with the precursor. Currently ``"HO2"``.
    adduct:
        Default adduct appended to formulas without an explicit adduct.
    show:
        Show generated Matplotlib figures.
    overlay:
        Overlay AP compounds in each time-trace figure when True.
    kinetic_savepath:
        Optional output path for the combined kinetic plot.
    stages_path:
        Optional CLOUD stages.txt path, loaded once for all trace figures.
        If omitted, use stages previously loaded with load_stages.
    stages:
        Stage-label mode for trace figures; "onlytype" shows the stage type.
        Labels stay at a fixed height in the axes when zooming vertically.
        Set to None to hide stage markers and labels.
    label_intervals:
        Add matching t_i tags to shaded intervals and kinetic points. Numbering
        starts at 1 for each temperature, following important_time_ranges order.
        Kinetic labels match each temperature series color. Skipped points do
        not renumber the remaining tags. Default is False.

    Returns
    -------
    trace_figures, means, regression_results, kinetic_fig, kinetic_ax
        ``means`` contains all interval means from all temperatures, including
        the ``temperature`` column. ``regression_results`` contains one row per
        successfully fitted temperature.
    """
    if not temperature_runs:
        raise ValueError("temperature_runs must contain at least one condition.")
    if not compounds:
        raise ValueError("compounds must contain at least one compound.")

    plt.rcParams.update(cfg.PAPER_RC_PARAMS)
    trace_compounds = _prepare_trace_compounds(compounds, adduct)
    normalized_adduct = (
        adduct if not adduct or adduct.startswith(".") else "." + adduct
    )
    precursor_name = (
        precursor_compound
        if "." in precursor_compound
        else f"{precursor_compound}{normalized_adduct}"
    )
    product_name = (
        product_compound
        if "." in product_compound
        else f"{product_compound}{normalized_adduct}"
    )

    if precursor_name not in trace_compounds:
        raise ValueError(f"Precursor {precursor_name!r} is not in compounds.")
    if product_name not in trace_compounds:
        raise ValueError(f"Product {product_name!r} is not in compounds.")

    if stages_path is not None and stages is not None:
        load_stages(stages_path)

    all_means: list[pd.DataFrame] = []
    trace_figures: list[tuple[plt.Figure, list[plt.Axes]]] = []
    kinetic_x_column: Optional[str] = None

    for run in temperature_runs:
        required = {
            "temperature",
            "trace_path",
            "time_range",
            "important_time_ranges",
        }
        missing = required.difference(run)
        if missing:
            raise ValueError(
                "Temperature run is missing: " + ", ".join(sorted(missing))
            )

        temperature = str(run["temperature"])
        time_range = _validate_time_ranges([run["time_range"]])[0]
        intervals = _validate_time_ranges(run["important_time_ranges"])
        if not intervals:
            raise ValueError(
                f"important_time_ranges must not be empty for {temperature}."
            )

        data = read_trace_timeseries(
            run["trace_path"],
            trace_compounds,
            time_range,
        )
        if data.empty:
            raise RuntimeError(
                f"No AP trace data found for {temperature} inside "
                f"{time_range[0]} to {time_range[1]} UTC."
            )

        time_col = next(
            (column for column in data.columns if column.lower() == "time"),
            None,
        )
        if time_col is None:
            raise RuntimeError(
                f"Loaded AP trace data for {temperature} contains no Time column."
            )

        means = calculate_interval_means(
            data,
            trace_compounds,
            time_col,
            intervals,
        )
        means, current_x_column = add_reaction_partner_interval_means(
            means=means,
            horus_path=horus_path,
            time_range=time_range,
            important_time_ranges=intervals,
            precursor_compound=precursor_name,
            reaction_partner=reaction_partner,
        )
        means.insert(0, "temperature", temperature)
        all_means.append(means)

        if kinetic_x_column is None:
            kinetic_x_column = current_x_column
        elif kinetic_x_column != current_x_column:
            raise RuntimeError("Inconsistent kinetic x-column between runs.")

        trace_figures.append(
            _plot_temperature_trace(
                data=data,
                trace_compounds=trace_compounds,
                time_col=time_col,
                time_range=time_range,
                important_time_ranges=intervals,
                temperature=temperature,
                overlay=overlay,
                show=show,
                savepath=run.get("trace_savepath"),
                stages=stages,
                label_intervals=label_intervals,
            )
        )

    combined_means = pd.concat(all_means, ignore_index=True)

    partner_label = reaction_partner.strip().upper()
    partner_plot_label = r"HO$_2$" if partner_label == "HO2" else partner_label
    x_label = (
        _chemical_label(precursor_name)
        + rf" $\times$ {partner_plot_label} (ppt$^2$)"
    )

    kinetic_fig, kinetic_ax, regression_results = plot_kinetic_relationship(
        grouped_means=combined_means,
        x_column=kinetic_x_column,
        y_compound=product_name,
        x_label=x_label,
        show=show,
        savepath=kinetic_savepath,
        label_intervals=label_intervals,
    )

    return (
        trace_figures,
        combined_means,
        regression_results,
        kinetic_fig,
        kinetic_ax,
    )
