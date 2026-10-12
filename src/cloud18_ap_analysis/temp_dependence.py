"""Temperature-dependence traces, interval means, and kinetic relationship plots.

All datetimes in this module are UTC. CLOUD18 timestamps are already UTC and
are therefore kept as timezone-naive ``datetime`` values throughout.
"""

from __future__ import annotations

from datetime import datetime
import re
from typing import Optional, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.odr import Model, ODR, RealData
from scipy.stats import t as student_t

from . import config as cfg
from .loading import load_data, read_trace_timeseries
from .models import Channel
from .plotting import _draw_stages
from .stages import load_stages
from .units import (convert_channel_concentration, concentration_label,
                    concentration_yscale, validate_concentration_unit, STANDARD_PRESSURE_PA)

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
    concentration_unit: str = "native",
    temperature_c=None,
    pressure_pa=STANDARD_PRESSURE_PA,
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

    partner_trace = convert_channel_concentration(
        partner_data[partner_channel], "pptv", concentration_unit,
        temperature_c=temperature_c, pressure_pa=pressure_pa,
    )
    partner_means = calculate_channel_interval_means(
        partner_trace,
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
    color_by_ap: bool = False,
    concentration_unit: Optional[str] = None,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Plot interval means with errors and fit zero-intercept ODR models.

    The legend reports each slope k plus/minus its ODR standard error (1 SE).
    Both x and y interval-mean uncertainties are used in the original-space fit.

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
    color_by_ap:
        Color marker fills by AP_pptv using one shared colorbar. Temperature
        remains identified by marker shapes, outlines, and fit-line colors.
        Missing AP means are gray and do not change the fitted data.
    concentration_unit:
        Unit of the supplied means (no conversion here). If None, read their
        concentration_unit attribute, falling back to native pptv.

    Returns
    -------
    figure, axis, regression_results
        The combined kinetic plot and a table containing the zero-intercept
        ODR slope, its standard error, goodness-of-fit values and the number of
        points for each temperature.
    """
    concentration_unit = concentration_unit or grouped_means.attrs.get("concentration_unit", "native")
    validate_concentration_unit(concentration_unit)
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
    if color_by_ap:
        required.append("AP_pptv")
    missing = [column for column in required if column not in grouped_means.columns]
    if missing:
        raise ValueError("Missing kinetic columns: " + ", ".join(missing))

    plt.rcParams.update(cfg.PAPER_RC_PARAMS)
    fig, ax = plt.subplots(figsize=cfg.KINETIC_FIGSIZE)
    regression_rows: list[dict] = []

    temperatures = list(pd.unique(grouped_means["temperature"]))
    colors = cfg.PAPER_LINE_COLORS

    ap_norm = None
    ap_cmap = plt.get_cmap("plasma").copy()
    ap_cmap.set_bad("0.6")
    missing_ap = False
    if color_by_ap:
        kinetic_columns = [x_column, y_compound, x_error_column, y_error_column]
        drawable = np.isfinite(grouped_means[kinetic_columns]).all(axis=1)
        drawable &= (grouped_means[x_error_column] > 0) & (grouped_means[y_error_column] > 0)
        concentrations = grouped_means.loc[drawable, "AP_pptv"].to_numpy(dtype=float)
        finite_ap = concentrations[np.isfinite(concentrations)]
        missing_ap = bool(np.any(~np.isfinite(concentrations)))
        if finite_ap.size:
            lo, hi = float(finite_ap.min()), float(finite_ap.max())
            if lo == hi:
                margin = max(abs(lo) * 0.05, 0.01)
                lo, hi = lo - margin, hi + margin
            ap_norm = Normalize(lo, hi)
        else:
            print("Warning: no finite AP interval means for kinetic points; showing gray points.")

    for index, temperature in enumerate(temperatures):
        valid = grouped_means.loc[
            grouped_means["temperature"] == temperature,
            [x_column, y_compound, x_error_column, y_error_column]
            + (["range_index"] if label_intervals else [])
            + (["AP_pptv"] if color_by_ap else []),
        ].replace([np.inf, -np.inf], np.nan).dropna(
            subset=[x_column, y_compound, x_error_column, y_error_column]
            + (["range_index"] if label_intervals else [])
        )

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
            fmt="none" if color_by_ap else "o",
            markersize=np.sqrt(cfg.KINETIC_MARKER_SIZE),
            color=color,
            markeredgewidth=0,
            elinewidth=cfg.PAPER_LINE_WIDTH,
            capsize=2,
            label="_nolegend_" if color_by_ap else str(temperature),
            zorder=3,
        )

        if color_by_ap:
            marker = ("o", "s", "^", "D", "v", "P", "X")[index % 7]
            if ap_norm is not None:
                ax.scatter(
                    x, y, c=np.ma.masked_invalid(valid["AP_pptv"].to_numpy(dtype=float)),
                    cmap=ap_cmap, norm=ap_norm, plotnonfinite=True,
                    marker=marker, s=cfg.KINETIC_MARKER_SIZE,
                    edgecolors=color, linewidths=0.8, zorder=4,
                )
            else:
                ax.scatter(x, y, color="0.6", marker=marker,
                           s=cfg.KINETIC_MARKER_SIZE, edgecolors=color, zorder=4)
            ax.plot([], [], marker=marker, color=color, markerfacecolor="white",
                    label=str(temperature), linewidth=cfg.PAPER_LINE_WIDTH)

        if label_intervals:
            for xi, yi, interval_index in zip(x, y, valid["range_index"]):
                ax.annotate(
                    _interval_label(interval_index), (xi, yi),
                    xytext=(5, 5), textcoords="offset points",
                    color=color, fontsize=cfg.PAPER_LABEL_SIZE,
                    ha="left", va="bottom", zorder=5,
                )

        if len(valid) < 2 or not np.any(x != 0.0):
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
        if concentration_unit == "molecules/cm3":
            positive_x = x[x > 0]
            x_fit = (np.geomspace(positive_x.min() / 1.05, positive_x.max() * 1.05, 200)
                     if positive_x.size else np.array([]))
        else:
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
            "slope_unit": "cm3/molecule" if concentration_unit == "molecules/cm3" else "1/pptv",
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

    axis_scale = concentration_yscale(concentration_unit)
    ax.set_xscale(axis_scale)
    ax.set_yscale(axis_scale)
    ax.set_xlabel(x_label or x_column)
    ax.set_ylabel(_chemical_label(y_compound) + " [" + concentration_label(concentration_unit) + "]")
    if color_by_ap and ap_norm is not None:
        fig.colorbar(ScalarMappable(norm=ap_norm, cmap=ap_cmap), ax=ax,
                     label="Mean chamber AP [pptv]")
    # Extend the temperature legend with the fitted slope and its ODR standard
    # error. Use a shared exponent so values and uncertainties are easy to read.
    fit_labels = {}
    k_unit = (r"\mathrm{cm^3\,molecule^{-1}}" if concentration_unit == "molecules/cm3"
              else r"\mathrm{pptv^{-1}}")
    for result in regression_rows:
        k, se = result["slope"], result["slope_stderr"]
        magnitude = max(abs(k), abs(se))
        exponent = int(np.floor(np.log10(magnitude))) if np.isfinite(magnitude) and magnitude > 0 else 0
        scale = 10.0 ** exponent
        fit_labels[str(result["temperature"])] = (
            str(result["temperature"]) + "\n"
            + rf"$k=({k / scale:.3g} \pm {se / scale:.2g})"
            + rf"\times 10^{{{exponent}}}\ {k_unit}$"
        )
    handles, labels = ax.get_legend_handles_labels()
    labels = [fit_labels.get(label, label) for label in labels]
    if color_by_ap and missing_ap:
        handles.append(Line2D([], [], marker="o", linestyle="none", color="0.6"))
        labels.append("AP unavailable")
    ax.legend(handles, labels, title=r"Temperature; $k \pm 1$ SE")
    _style_axis(ax)
    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

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
    concentration_unit: str = "native",
    yscale: Optional[str] = None,
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
        ax.set_ylabel("Concentration [" + concentration_label(concentration_unit) + "]")
        ax.set_xlabel("Time (UTC)")
        ax.set_yscale(concentration_yscale(concentration_unit, yscale))
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
            ax.set_ylabel("Concentration [" + concentration_label(concentration_unit) + "]")
            ax.set_yscale(concentration_yscale(concentration_unit, yscale))
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
    else:
        plt.close(fig)

    return fig, axes


def _nominal_temperature_c(run):
    """Resolve a numeric plot temperature without affecting conversion inputs."""
    if "plot_temperature_c" in run:
        value = float(run["plot_temperature_c"])
    elif "temperature_c" in run:
        value = float(run["temperature_c"])
    else:
        match = re.fullmatch(r"\s*([+-]?\d+(?:\.\d+)?)\s*(?:°?\s*C)?\s*",
                             str(run["temperature"]))
        if match is None:
            return np.nan
        value = float(match.group(1))
    if not np.isfinite(value) or value <= -273.15:
        raise ValueError("Plot temperature must be finite and above absolute zero (Celsius)")
    return value


def plot_k_temperature(regression_results, concentration_unit="native", *,
                       temperature_order=None, show=True, savepath=None, title=None):
    """Plot fitted concentration slopes versus nominal Celsius temperature.

    Uses temperature_c, slope and slope_stderr columns; error bars are 1 SE.
    No temperature-dependence model is fitted. Missing/nonfinite results are
    skipped with a message. Axes are linear so zero/negative slopes are visible.
    """
    validate_concentration_unit(concentration_unit)
    fig, ax = plt.subplots(figsize=cfg.KINETIC_FIGSIZE)
    order = list(temperature_order or regression_results.get("temperature", []))
    count = 0
    for _, row in regression_results.iterrows():
        values = [row["temperature_c"], row["slope"], row["slope_stderr"]]
        if not np.all(np.isfinite(values)) or row["slope_stderr"] < 0:
            print(f"Skipping k(T) point for {row['temperature']}: missing temperature or finite slope/error.")
            continue
        index = order.index(row["temperature"])
        ax.errorbar(row["temperature_c"], row["slope"], yerr=row["slope_stderr"],
                    fmt="o", color=cfg.PAPER_LINE_COLORS[index % len(cfg.PAPER_LINE_COLORS)],
                    markersize=np.sqrt(cfg.KINETIC_MARKER_SIZE), capsize=3,
                    label=str(row["temperature"]))
        count += 1
    if count:
        ax.legend(title="Temperature")
    else:
        print("No usable slopes with numeric temperatures for the k(T) plot.")
        ax.text(0.5, 0.5, "No usable k(T) data", transform=ax.transAxes,
                ha="center", va="center")
    ax.set_xlabel(r"Temperature [$^\circ$C]")
    unit = r"cm$^3$ molecule$^{-1}$" if concentration_unit == "molecules/cm3" else r"pptv$^{-1}$"
    ax.set_ylabel(f"Concentration slope k [{unit}]")
    ax.set_title(title or "Concentration slope versus temperature")
    ax.text(0.02, 0.98, "Error bars: ±1 SE", transform=ax.transAxes,
            ha="left", va="top", fontsize=cfg.PAPER_LABEL_SIZE)
    _style_axis(ax)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig, ax


def _normalize_show(show):
    """Validate the public display selector before loading any data."""
    available = {"traces", "kinetics", "k(T)"}
    if isinstance(show, bool):
        return available if show else set()
    if isinstance(show, str):
        names = [show]
    elif isinstance(show, (list, tuple)) and all(isinstance(name, str) for name in show):
        names = show
    else:
        raise TypeError('show must be True, False, a plot name, or a list/tuple of plot names')
    unknown = set(names) - available
    if unknown:
        raise ValueError(f"Unknown show selection: {sorted(unknown)}. Choose 'traces', 'kinetics', or 'k(T)'.")
    return set(names)


def temp_dep(
    temperature_runs: Sequence[dict],
    compounds: Sequence[str],
    horus_path: str,
    precursor_compound: str,
    product_compound: str,
    reaction_partner: str = "HO2",
    adduct: str = ".NH3H+",
    show: bool | str | Sequence[str] = True,
    overlay: bool = True,
    kinetic_savepath: Optional[str] = None,
    stages_path: Optional[str] = None,
    stages: Optional[str] = "onlytype",
    label_intervals: bool = False,
    color_by_ap: bool = False,
    ap_path: Optional[str] = None,
    concentration_unit: str = "native",
    pressure_pa=STANDARD_PRESSURE_PA,
    yscale: Optional[str] = None,
    k_temperature_savepath: Optional[str] = None,
) -> tuple[
    list[tuple[plt.Figure, list[plt.Axes]]],
    pd.DataFrame,
    pd.DataFrame,
    plt.Figure,
    plt.Axes,
    plt.Figure,
    plt.Axes,
]:
    """Plot temperature-specific traces and fit their interval-mean kinetics.

    For each temperature run, shade the requested averaging intervals and
    calculate one kinetic point per interval:

    ``x = mean(precursor) * mean(reaction partner)``
    ``y = mean(product)``

    Fit ``y = k*x`` separately for each temperature using zero-intercept ODR
    with both x and y uncertainties. The kinetic legend displays each slope
    k and its standard error (1 SE). All datetimes must be timezone-naive UTC,
    e.g. ``datetime(2025, 9, 25, 10, 0)``; interval endpoints are inclusive.

    Parameters
    ----------
    temperature_runs : sequence of dict
        Required, nonempty list/tuple of run definitions. One trace figure is
        created per dictionary, in the supplied order. Dictionary fields:

        * ``temperature`` (required, str): display label, e.g. ``"-50 C"``.
          Used to group kinetic fits; use a distinct label for each condition.
          This label is NOT used as the numerical conversion temperature.
        * ``plot_temperature_c`` (optional, float): nominal Celsius temperature
          for the k(T) x-coordinate, e.g. -50.0. Does not affect concentration
          conversion. Defaults to temperature_c, then a numeric Celsius label
          such as "-50 °C". Unresolvable labels are skipped with a message.
        * ``trace_path`` (required, str or Path): calibrated compound trace
          CSV readable by ``read_trace_timeseries``; values are assumed pptv.
        * ``time_range`` (required, tuple of datetime): ``(start, stop)`` for
          loading and displaying that run; stop must be >= start.
        * ``important_time_ranges`` (required, nonempty sequence of tuples):
          ``[(start1, stop1), (start2, stop2), ...]`` for averaging. Choose
          intervals inside ``time_range``. Their input order defines t_1,
          t_2, etc.; numbering restarts for each temperature.
        * ``trace_savepath`` (optional, str or Path): figure filename, e.g.
          ``"output/traces_m50.png"``; omit or use None to avoid saving.
        * ``temperature_path`` (optional, str or Path): measurement file
          readable by ``load_data`` with a temperature channel in Celsius.
          Used only for ``concentration_unit="molecules/cm3"``.
        * ``temperature_channel`` (optional, str): temperature column name
          in that file; default ``"TE_Calib_3"``.
        * ``temperature_c`` (optional, float): constant Celsius temperature,
          e.g. ``-50.0``, as an alternative to measured temperature.
          In molecules/cm3 mode, supply either ``temperature_path`` or
          ``temperature_c``. If both are supplied, the file takes precedence.

    compounds : sequence of str
        Required, nonempty names of trace compounds to load and plot, e.g.
        ``["C10H15O4.NH3H+", "C10H16O4.NH3H+"]``. Bare formulas such as
        ``"C10H15O4"`` receive the ``adduct`` suffix. Both the precursor and
        product must occur in this list after adduct expansion.
    horus_path : str or Path
        Required HORUS measurement filename, readable by ``load_data``.
        Must contain ``HO2_ppt`` for the currently supported HO2 partner;
        those measurements are interpreted as pptv.
    precursor_compound : str
        Required precursor name, e.g. ``"C10H15O4.NH3H+"`` or ``"C10H15O4"``.
        Its interval mean is multiplied by the partner mean to form x.
    product_compound : str
        Required product name, e.g. ``"C10H16O4.NH3H+"`` or ``"C10H16O4"``.
        Its interval mean forms y. Uses the same adduct rules as the precursor.
    reaction_partner : str, default "HO2"
        Currently only ``"HO2"`` is implemented (case-insensitive; surrounding
        whitespace is ignored). Other partners raise NotImplementedError.
    adduct : str, default ".NH3H+"
        Suffix appended to compound names without a dot. Examples:
        ``".NH3H+"``, ``"NH3H+"``, ``".H+"``, or ``""`` for no suffix.
        Explicitly adducted names remain unchanged. The resulting names must
        match columns in the trace files.
    show : bool, str, or list/tuple of str, default True
        Select which plot windows open:

        * ``True``: show all trace, kinetic, and k(T) figures.
        * ``False`` or ``[]``: show no windows.
        * ``"traces"``: show only the time-trace figures.
        * ``"kinetics"``: show only the combined kinetic figure.
        * ``"k(T)"``: show only slopes versus temperature, with 1 SE bars
          and no temperature-dependence fit.
        * ``["kinetics", "k(T)"]``: combine any selections (tuples work too).

        Names are case-sensitive; unknown names raise ValueError before data
        loading. Selection order does not change the plotting order. All
        calculations, figure creation, return values and requested file saves
        still happen. Unselected figures are closed from pyplot's window
        registry but remain usable through the returned Figure objects.
        Showing a window may block until it closes.
    overlay : bool, default True
        True overlays all compounds on one axis per temperature. False gives
        each compound a separate subplot sharing the time axis. Does not
        change the combined kinetic plot or averaging calculations.
    kinetic_savepath : str or Path or None, default None
        Filename for the combined kinetic figure; None means do not save it.
        The extension selects the format, e.g. .png or .pdf. Create parent
        directories first, for this path and any ``trace_savepath`` values.
    stages_path : str or Path or None, default None
        CLOUD stages.txt file to load once for all trace figures. None reuses
        stages already loaded with ``load_stages``. Ignored when stages=None.
    stages : str or None, default "onlytype"
        Stage-label content on trace plots:

        * ``"onlytype"``: stage type only.
        * ``"stage"``: stage number/name only.
        * ``"type"``: stage number/name followed by its type.
        * ``"description"``: description, falling back to stage name.
        * ``"comments"``: comments, falling back to stage name.
        * ``"details"``: multiline log details, falling back to stage name.
        * ``"full"``: stage name, type, and description combined.
        * ``None``: hide stage markers and labels.

        Labels keep a fixed vertical position when zooming the y-axis.
        Without loaded stage data, no markers appear.
    label_intervals : bool, default False
        True adds matching t_i tags to shaded averaging intervals and kinetic
        points. False hides those tags but keeps interval shading. Omitted
        kinetic points do not renumber the remaining tags.
    color_by_ap : bool, default False
        False colors kinetic points by temperature. True colors their fills
        by mean chamber alpha-pinene concentration for the same interval,
        with one shared pptv colorbar. Temperature remains identifiable by
        marker shapes, outlines, and fit-line colors. Requires ``ap_path``.
        Missing AP means appear gray and do not alter the kinetic fit.
    ap_path : str or Path or None, default None
        Calibrated VOCUS file containing channel ``"AP"`` in ppbv. Required
        only when color_by_ap=True; otherwise ignored. AP is multiplied by
        1000 to display pptv, independently of ``concentration_unit``.
    concentration_unit : {"native", "molecules/cm3"}, default "native"
        ``"native"`` keeps compound and HO2 concentrations in pptv.
        ``"molecules/cm3"`` converts each measurement using temperature,
        pressure, the ideal gas law, and Avogadro's constant BEFORE averaging.
        Each run then needs measured or constant temperature (see above).
        Source channel names remain identifiers, even if they contain _ppt.
        AP colorbar values always stay in pptv.
    pressure_pa : float or Channel, default 101325.0
        Chamber pressure in Pa, e.g. ``101325.0`` for one standard atmosphere.
        A Channel supplies time-dependent pressure: values in Pa, timestamps
        in Unix milliseconds. Used only in molecules/cm3 mode. Pressure and
        absolute temperature must be positive. Measured conditions are
        linearly interpolated to each concentration timestamp without
        extrapolation; missing coverage gives NaN, and no coverage raises an
        error. Temperature/pressure uncertainties are not propagated.
    yscale : str or None, default None
        Trace y-axis scale only. None selects ``"log"`` for molecules/cm3 or
        ``"linear"`` for native units. Use ``"linear"`` or ``"log"`` to
        override; other Matplotlib scale names such as ``"symlog"`` are
        also forwarded. Nonpositive values cannot be shown on log axes.
        Kinetic x/y axes independently use log in molecules/cm3 mode and
        linear in native mode. Fits always use the original concentration
        values, not their logarithms.

    k_temperature_savepath : str or Path or None, default None
        Optional filename for the k(T) plot. Saved regardless of show selection;
        its parent directory must exist. None disables saving this figure.

    Returns
    -------
    trace_figures : list of (Figure, list of Axes)
        One figure and its trace axes per input run.
    means : pandas.DataFrame
        Interval means for all runs, with temperature, range_index, start/end,
        counts, compound/partner means and SEMs, plus the kinetic x product
        and its propagated SEM. AP_pptv, AP_pptv_sem and AP_pptv_n_points are
        included when color_by_ap=True. ``means.attrs["concentration_unit"]``
        records the chosen mode; original column names are preserved.
    regression_results : pandas.DataFrame
        One row per fitted temperature: slope, slope_stderr (1 SE), fixed
        zero intercept, slope_unit, centered/uncentered R-squared, p-value,
        ODR diagnostics and n_points. Slope units are cm3/molecule or 1/pptv;
        no per-second factor is implied by this concentration relationship.
    kinetic_fig : matplotlib.figure.Figure
        Combined kinetic figure, with AP colorbar when requested.
    kinetic_ax : matplotlib.axes.Axes
        Main kinetic axes, for further customization.
    k_temperature_fig : matplotlib.figure.Figure
        Slopes versus nominal Celsius temperature, without a fit.
    k_temperature_ax : matplotlib.axes.Axes
        Linear k(T) axes. The two new return values follow the original five.

    Notes
    -----
    SEMs describe variability of interval means, not full instrumental or
    calibration uncertainty. The x-product uncertainty assumes independent
    precursor and partner errors. Points with missing/nonfinite values or
    nonpositive x/y SEMs are excluded. A fit needs at least two usable points
    and an x-array that is not entirely zero.

    Examples
    --------
    Use your real input filenames; create the output directory before saving::

        from datetime import datetime
        from cloud18_ap_analysis import temp_dep

        runs = [{
            "temperature": "-50 C",
            "trace_path": "traces_m50.csv",
            "time_range": (datetime(2025, 9, 25), datetime(2025, 9, 27)),
            "important_time_ranges": [
                (datetime(2025, 9, 25, 10), datetime(2025, 9, 25, 20)),
                (datetime(2025, 9, 26, 2), datetime(2025, 9, 26, 6)),
            ],
            "temperature_path": "temp_dew_frost_-50C.csv",
            "temperature_channel": "TE_Calib_3",
        }]
        traces, means, fits, fig, ax, kt_fig, kt_ax = temp_dep(
            temperature_runs=runs,
            compounds=["C10H15O4", "C10H16O4"],
            horus_path="HORUS.txt",
            precursor_compound="C10H15O4",
            product_compound="C10H16O4",
            concentration_unit="molecules/cm3",  # "native" restores pptv
            pressure_pa=101325.0,
            stages_path="stages.txt",
            stages="onlytype",
            label_intervals=True,
            color_by_ap=True,
            ap_path="AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv",
            overlay=True,
            yscale=None,
            show=False,
        )
    """
    shown_plots = _normalize_show(show)
    validate_concentration_unit(concentration_unit)
    if color_by_ap and ap_path is None:
        raise ValueError("ap_path is required when color_by_ap=True.")
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

    plot_temperatures = {}
    for run in temperature_runs:
        label = str(run["temperature"])
        value = _nominal_temperature_c(run)
        if label in plot_temperatures and not np.isclose(plot_temperatures[label], value, equal_nan=True):
            raise ValueError(f"Conflicting plot temperatures for {label}")
        plot_temperatures[label] = value

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

        temperature_c = None
        if concentration_unit == "molecules/cm3":
            if "temperature_path" in run:
                temperature_c = load_data(run["temperature_path"])[
                    run.get("temperature_channel", "TE_Calib_3")
                ]
            elif "temperature_c" in run:
                temperature_c = run["temperature_c"]
            else:
                raise ValueError("Each run needs temperature_path or temperature_c for molecules/cm3")
            data = data.copy()
            timestamps = pd.to_datetime(data[time_col], errors="coerce")
            time_ms = timestamps.to_numpy(dtype="datetime64[ms]").astype("int64").astype(float)
            time_ms[timestamps.isna().to_numpy()] = np.nan
            for compound in trace_compounds:
                data[compound] = convert_channel_concentration(
                    Channel(compound, time_ms, data[compound].to_numpy()), "pptv", "molecules/cm3",
                    temperature_c=temperature_c, pressure_pa=pressure_pa,
                ).values

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
            concentration_unit=concentration_unit,
            temperature_c=temperature_c,
            pressure_pa=pressure_pa,
        )
        if color_by_ap:
            ap_data = load_data(ap_path, time_range)
            if "AP" not in ap_data.channels:
                raise KeyError("AP channel not found in " + str(ap_path))
            ap_means = calculate_channel_interval_means(
                convert_channel_concentration(ap_data["AP"], "ppbv", "pptv"),
                intervals, column_name="AP_pptv",
            )
            means = means.merge(ap_means, on="range_index", how="left",
                                validate="one_to_one")
        means.insert(0, "temperature", temperature)
        all_means.append(means)

        if kinetic_x_column is None:
            kinetic_x_column = current_x_column
        elif kinetic_x_column != current_x_column:
            raise RuntimeError("Inconsistent kinetic x-column between runs.")

        with plt.ioff():
            trace_figures.append(
                _plot_temperature_trace(
                    data=data,
                    trace_compounds=trace_compounds,
                    time_col=time_col,
                    time_range=time_range,
                    important_time_ranges=intervals,
                    temperature=temperature,
                    overlay=overlay,
                    show="traces" in shown_plots,
                    savepath=run.get("trace_savepath"),
                    stages=stages,
                    label_intervals=label_intervals,
                    concentration_unit=concentration_unit,
                    yscale=yscale,
                )
            )

    combined_means = pd.concat(all_means, ignore_index=True)
    combined_means.attrs["concentration_unit"] = concentration_unit
    combined_means.attrs["AP_pptv_unit"] = "pptv"

    partner_label = reaction_partner.strip().upper()
    partner_plot_label = r"HO$_2$" if partner_label == "HO2" else partner_label
    x_label = (
        _chemical_label(precursor_name)
        + rf" $\times$ {partner_plot_label} "
        + (r"[(molecules/cm$^3$)$^2$]" if concentration_unit == "molecules/cm3" else "[pptv$^2$]")
    )

    with plt.ioff():
        kinetic_fig, kinetic_ax, regression_results = plot_kinetic_relationship(
            grouped_means=combined_means,
            x_column=kinetic_x_column,
            y_compound=product_name,
            x_label=x_label,
            show="kinetics" in shown_plots,
            savepath=kinetic_savepath,
            label_intervals=label_intervals,
            color_by_ap=color_by_ap,
            concentration_unit=concentration_unit,
        )

    regression_results["temperature_c"] = (
        regression_results["temperature"].map(plot_temperatures)
        if "temperature" in regression_results else pd.Series(dtype=float)
    )
    with plt.ioff():
        k_temperature_fig, k_temperature_ax = plot_k_temperature(
            regression_results, concentration_unit,
            temperature_order=list(plot_temperatures),
            show="k(T)" in shown_plots, savepath=k_temperature_savepath,
            title=_chemical_label(product_name) + " versus " + _chemical_label(precursor_name)
                  + " × " + partner_plot_label,
        )

    return (
        trace_figures,
        combined_means,
        regression_results,
        kinetic_fig,
        kinetic_ax,
        k_temperature_fig,
        k_temperature_ax,
    )
