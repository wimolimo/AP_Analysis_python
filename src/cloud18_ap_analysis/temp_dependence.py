"""Temperature-dependence traces, interval means, and kinetic relationship plots.

All datetimes in this module are UTC. CLOUD18 timestamps are already UTC and
are therefore kept as timezone-naive ``datetime`` values throughout.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
import re
from typing import Optional, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.odr import Model, ODR, RealData
from scipy.stats import t as student_t

from . import config as cfg
from .loading import load_data, read_trace_timeseries
from .models import Channel
from .plotting import (_draw_stages, _draw_fan_speed_shading,
                       _shade_intervals, _interval_label, _smooth_channel, _prepare_savepath,
                       _maximize_window, _append_to_filename)
from .stages import load_stages
from .units import (convert_channel_concentration, concentration_label,
                    concentration_yscale, validate_concentration_unit, STANDARD_PRESSURE_PA)

TimeRange = tuple[datetime, datetime]


def _signal_name(selection, adduct):
    if isinstance(selection, tuple):
        if len(selection) != 2 or not all(isinstance(s, str) and s.strip() for s in selection):
            raise ValueError("Use (instrument, exact channel name) for a signal")
        return f"{selection[0]}:{selection[1]}"
    if not isinstance(selection, str):
        raise TypeError("A signal must be a name or an (instrument, channel) tuple")
    return _prepare_trace_compounds([selection], adduct)[0]


def _read_selected_signals(run, names, time_range, target_unit, temperature_c, pressure_pa):
    """Read exact prepared channels, retaining each instrument's sample clock."""
    from .prepared import prepared_channel_info
    sources = run.get("sources", {})
    catalogs, series, units, sample_times = {}, [], set(), {}
    for name in names:
        instrument, separator, channel = name.partition(":")
        if not separator:
            instrument, channel = "ptr", name
        if instrument not in sources:
            raise KeyError(f"{run.get('run_id', 'run')}: instrument {instrument!r} is unavailable; use prepared runs")
        if instrument not in catalogs:
            catalogs[instrument] = prepared_channel_info(sources[instrument])
        if channel not in catalogs[instrument]:
            raise KeyError(f"{instrument}: exact channel {channel!r} not found; use the name stored in H5")
        unit = catalogs[instrument][channel]["unit"]
        trace = load_data(sources[instrument], time_range, channels=[channel])[channel]
        if not len(trace.time) or not np.any(np.isfinite(trace.values)):
            raise ValueError(f"{run.get('run_id', 'run')}: {name} has no valid data in this run window")
        if target_unit != "native" and unit != target_unit:
            trace = convert_channel_concentration(trace, unit, target_unit,
                        temperature_c=temperature_c, pressure_pa=pressure_pa)
        units.add(unit if target_unit == "native" else target_unit)
        series.append(pd.Series(trace.values, index=trace.time, name=name))
        sample_times[name] = pd.to_datetime(trace.time, unit="ms")
    if len(units) > 1:
        raise ValueError("Mixed instrument units require concentration_unit='molecules/cm3'")
    frame = pd.concat(series, axis=1).sort_index()
    frame.insert(0, "Time", pd.to_datetime(frame.index, unit="ms"))
    frame = frame.reset_index(drop=True)
    frame.attrs["sample_times"] = sample_times
    frame.attrs["native_unit"] = next(iter(units))
    return frame


def _available_temperature_runs(runs, trace_names, partner_name, color_by_ap):
    """Filter prepared runs by selected measurements, without masking typos.

    A channel must exist in at least one candidate run. Individual runs may
    then be skipped for missing channels or finite samples in their window.
    Legacy raw-file runs retain their existing validation in temp_dep.
    """
    from .prepared import prepared_channel_info
    requests = []
    for name in trace_names:
        instrument, separator, channel = name.partition(":")
        requests.append((instrument, channel) if separator else ("ptr", name))
    horus_channels = {"HO2": "HO2_ppt", "RO2": "RO2_crosstalk_ppt"}
    if partner_name.strip().upper() in horus_channels:
        requests.append(("horus", horus_channels[partner_name.strip().upper()]))
    if color_by_ap:
        requests.append(("vocus", "AP"))
    requests = list(dict.fromkeys(requests))
    prepared = [run for run in runs if "sources" in run]
    if not prepared:
        return runs
    catalogs = {}
    for index, run in enumerate(prepared):
        _validate_time_ranges([run["time_range"]])
        for instrument, _ in requests:
            if instrument in run["sources"] and (index, instrument) not in catalogs:
                catalogs[index, instrument] = prepared_channel_info(run["sources"][instrument])
    # Validate all names before skipping anything, so a typo is never silent.
    for instrument, channel in requests:
        if not any(channel in catalogs.get((i, instrument), {}) for i in range(len(prepared))):
            raise KeyError(f"Exact channel {instrument}:{channel} not found in any candidate prepared run; check the instrument and channel name")
    excluded = set()
    for index, run in enumerate(prepared):
        reasons = []
        for instrument, channel in requests:
            if channel not in catalogs.get((index, instrument), {}):
                reasons.append(f"missing {instrument}:{channel}")
                continue
            trace = load_data(run["sources"][instrument], run["time_range"], channels=[channel])[channel]
            if not np.any(np.isfinite(trace.values)):
                reasons.append(f"no valid {instrument}:{channel} samples in this window")
        if reasons:
            excluded.add(run["run_id"])
            print(f"Skipping {run['run_id']}: " + "; ".join(reasons))
    available = [run for run in runs if run["run_id"] not in excluded]
    if not available:
        raise ValueError("No temperature runs have data for all selected signals")
    for group in dict.fromkeys(run["fit_group"] for run in runs if run["pooled"]):
        members = [run for run in runs if run["fit_group"] == group]
        if any(run["run_id"] in excluded for run in members):
            retained = [run["run_id"] for run in members if run["run_id"] not in excluded]
            print(f"Combined fit {group}: available members " + (", ".join(retained) or "none"))
    return available


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

    Reaction partners ``"HO2"`` and ``"RO2"`` are implemented. The selection is kept
    explicit so additional partners such as O3 can later be added without
    changing the public ``temp_dep`` interface.

    Parameters
    ----------
    means:
        AP interval-mean table.
    horus_path:
        Path to the HORUS HOx/ROx file used for either partner.
    time_range:
        Overall UTC loading interval.
    important_time_ranges:
        UTC averaging intervals T1, T2, ... .
    precursor_compound:
        AP precursor column to multiply by the partner mean.
    reaction_partner:
        Chemical reaction partner: ``"HO2"`` or ``"RO2"``.

    Returns
    -------
    means_with_partner, kinetic_x_column
        Updated interval-mean table and the name of the derived kinetic
        x-column.
    """
    partner = reaction_partner.strip().upper()
    if reaction_partner in means.columns:
        out = means.copy()
        return _add_kinetic_product(out, precursor_compound, reaction_partner)

    partner_channels = {"HO2": "HO2_ppt", "RO2": "RO2_crosstalk_ppt"}
    if partner in partner_channels:
        partner_channel = partner_channels[partner]
        partner_label = partner
        partner_data = load_data(horus_path, time_range, channels=[partner_channel])
    else:
        raise NotImplementedError(
            f"Reaction partner {reaction_partner!r} is not implemented yet. "
            "Currently supported: 'HO2', 'RO2'."
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
    return _add_kinetic_product(out, precursor_compound, partner_channel)


def _add_kinetic_product(out, precursor_compound, partner_channel):
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
    fit_intercept: bool = True,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Plot interval means with errors and fit linear ODR models.

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

    fit_intercept : bool, default True
        True fits y = k*x + b; False fixes b=0. Both x/y errors are used.
        A free intercept needs at least three points with distinct x values;
        fixed zero needs at least two points and a nonzero x.

    Returns
    -------
    figure, axis, regression_results
        The combined kinetic plot and a table containing the linear
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
    ap_cmap = plt.get_cmap("turbo").copy()
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
        candidates = grouped_means.loc[grouped_means["temperature"] == temperature]
        valid = grouped_means.loc[
            grouped_means["temperature"] == temperature,
            [x_column, y_compound, x_error_column, y_error_column]
            + (["range_index"] if label_intervals else [])
            + (["AP_pptv"] if color_by_ap else [])
            + (["interval_label"] if "interval_label" in grouped_means.columns else []),
        ].replace([np.inf, -np.inf], np.nan).dropna(
            subset=[x_column, y_compound, x_error_column, y_error_column]
            + (["range_index"] if label_intervals else [])
        )

        # ODR requires strictly positive standard deviations. A zero SEM means
        # that the interval provides no usable estimate of measurement spread.
        valid = valid.loc[
            (valid[x_error_column] > 0) & (valid[y_error_column] > 0)
        ]

        finite_means = np.isfinite(candidates[[x_column, y_compound]]).all(axis=1)
        usable_errors = (np.isfinite(candidates[[x_error_column, y_error_column]]).all(axis=1)
                         & (candidates[x_error_column] > 0) & (candidates[y_error_column] > 0))
        missing_means = int((~finite_means).sum())
        missing_errors = int((finite_means & ~usable_errors).sum())
        excluded = len(candidates) - len(valid)
        reasons = []
        if missing_means:
            reasons.append(f"{missing_means} missing/nonfinite means")
        if missing_errors:
            reasons.append(f"{missing_errors} missing/nonpositive uncertainties")
        missing_labels = excluded - missing_means - missing_errors
        if missing_labels:
            reasons.append(f"{missing_labels} missing interval labels")
        print(f"Kinetics {temperature}: excluded {excluded} of {len(candidates)} interval points; "
              f"{len(valid)} retained" + (" (" + "; ".join(reasons) + ")" if reasons else "") + ".")
        if "run_id" in candidates and candidates["run_id"].nunique() > 1:
            usable = finite_means & usable_errors
            if label_intervals:
                usable &= np.isfinite(candidates["range_index"])
            for run_id, run_points in candidates.groupby("run_id", sort=False):
                run_mask = candidates["run_id"] == run_id
                retained = int(usable.loc[run_mask].sum())
                print(f"  {run_id}: excluded {len(run_points) - retained} of "
                      f"{len(run_points)} interval points; {retained} retained.")

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
            for xi, yi, interval_tag in zip(x, y, valid["interval_label"] if "interval_label" in valid else valid["range_index"].map(_interval_label)):
                ax.annotate(
                    interval_tag, (xi, yi),
                    xytext=(5, 5), textcoords="offset points",
                    color=color, fontsize=cfg.PAPER_LABEL_SIZE,
                    ha="left", va="bottom", zorder=5,
                )

        minimum_points = 3 if fit_intercept else 2
        identifiable = np.ptp(x) > 0 if fit_intercept else np.any(x != 0.0)
        if len(valid) < minimum_points or not identifiable:
            requirement = "distinct x values" if fit_intercept else "a non-zero x value"
            print(f"Skipping ODR fit at {temperature}: need at least {minimum_points} "
                  f"valid points and {requirement} (fit_intercept={fit_intercept}).")
            continue

        def linear_model(beta: np.ndarray, x_values: np.ndarray) -> np.ndarray:
            return beta[0] * x_values + (beta[1] if fit_intercept else 0.0)

        if fit_intercept:
            centered_x = x - x.mean()
            initial_slope = float(np.dot(centered_x, y - y.mean()) / np.dot(centered_x, centered_x))
            beta0 = [initial_slope, float(y.mean() - initial_slope * x.mean())]
        else:
            initial_slope = float(np.dot(x, y) / np.dot(x, x))
            beta0 = [initial_slope]
        odr_result = ODR(
            RealData(x, y, sx=x_error, sy=y_error),
            Model(linear_model),
            beta0=beta0,
        ).run()

        intercept = float(odr_result.beta[1]) if fit_intercept else 0.0
        intercept_stderr = float(odr_result.sd_beta[1]) if fit_intercept else 0.0
        slope = float(odr_result.beta[0])
        slope_stderr = float(odr_result.sd_beta[0])
        if concentration_unit == "molecules/cm3":
            positive_x = x[x > 0]
            x_fit = (np.geomspace(positive_x.min() / 1.05, positive_x.max() * 1.05, 200)
                     if positive_x.size else np.array([]))
        else:
            x_fit = np.linspace(0.0, x.max() * 1.05, 200)
        y_fit = slope * x_fit + intercept
        if concentration_unit == "molecules/cm3":
            y_fit = np.where(y_fit > 0, y_fit, np.nan)

        ax.plot(
            x_fit,
            y_fit,
            color=color,
            linewidth=cfg.PAPER_LINE_WIDTH,
            zorder=2,
        )

        fitted_y = slope * x + intercept
        ss_res = float(np.sum((y - fitted_y) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        centered_r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        uncentered_ss_tot = float(np.sum(y ** 2))
        uncentered_r_squared = (
            1.0 - ss_res / uncentered_ss_tot
            if uncentered_ss_tot > 0
            else np.nan
        )
        degrees_of_freedom = len(valid) - (2 if fit_intercept else 1)
        p_value = (
            float(2.0 * student_t.sf(abs(slope / slope_stderr), degrees_of_freedom))
            if degrees_of_freedom > 0 and slope_stderr > 0
            else np.nan
        )

        regression_rows.append({
            "temperature": temperature,
            "slope_unit": "cm3/molecule" if concentration_unit == "molecules/cm3" else "1/pptv",
            "slope": slope,
            "intercept": intercept,
            "intercept_stderr": intercept_stderr,
            "intercept_unit": "molecules/cm3" if concentration_unit == "molecules/cm3" else "pptv",
            "fit_intercept": fit_intercept,
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
        if fit_intercept:
            b, b_se = result["intercept"], result["intercept_stderr"]
            b_magnitude = max(abs(b), abs(b_se))
            b_exponent = int(np.floor(np.log10(b_magnitude))) if np.isfinite(b_magnitude) and b_magnitude > 0 else 0
            b_scale = 10.0 ** b_exponent
            b_unit = (r"\mathrm{molecules\,cm^{-3}}" if concentration_unit == "molecules/cm3"
                      else r"\mathrm{pptv}")
            fit_labels[str(result["temperature"])] += (
                "\n" + rf"$b=({b / b_scale:.3g} \pm {b_se / b_scale:.2g})"
                + rf"\times 10^{{{b_exponent}}}\ {b_unit}$"
            )
    handles, labels = ax.get_legend_handles_labels()
    labels = [fit_labels.get(label, label) for label in labels]
    if color_by_ap and missing_ap:
        handles.append(Line2D([], [], marker="o", linestyle="none", color="0.6"))
        labels.append("AP unavailable")
    ax.legend(handles, labels, title="Temperature; fit ±1σ")
    _style_axis(ax)
    fig.tight_layout()

    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
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
    fan_speed_data=None,
    smoothing_time_traces: timedelta | None = None,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Create the time-trace figure for one temperature condition."""
    plot_time = pd.to_datetime(data[time_col], errors="coerce")
    sample_times = data.attrs.get("sample_times", {})

    def samples(compound):
        mask = plot_time.isin(sample_times[compound]) if compound in sample_times else np.ones(len(data), dtype=bool)
        times, values = plot_time[mask], data.loc[mask, compound]
        if smoothing_time_traces is None:
            return times, values
        time_ms = times.to_numpy(dtype="datetime64[ms]").astype("int64").astype(float)
        smoothed = _smooth_channel(Channel(compound, time_ms, values.to_numpy(dtype=float)),
                                   smoothing_time_traces)
        return times, smoothed.values

    if overlay:
        fig, ax = plt.subplots(figsize=cfg.TRACE_FIGSIZE, layout="constrained")
        axes = [ax]
        _shade_intervals(ax, important_time_ranges, label_intervals=label_intervals, xlims=time_range)
        for compound in trace_compounds:
            ax.plot(
                *samples(compound),
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
            layout="constrained",
        )
        axes = list(axes_array[:, 0])
        for ax, compound in zip(axes, trace_compounds):
            _shade_intervals(ax, important_time_ranges, label_intervals=label_intervals, xlims=time_range)
            ax.plot(*samples(compound), linewidth=cfg.PAPER_LINE_WIDTH)
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

    for index, ax in enumerate(axes):
        ax.set_xlim(*time_range)
        fan_strips = []
        if fan_speed_data is not None:
            fan_strips = _draw_fan_speed_shading(
                ax, fan_speed_data, ymin=0.0, ymax=0.035, alpha=0.65, zorder=4,
            ) or []
        if fan_strips and index == 0:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles + [Patch(facecolor="gray", alpha=0.65)],
                      labels + ["Fan speed 12%"], loc="upper left",
                      handlelength=1.5, borderaxespad=0.2)
        if stages is not None:
            _draw_stages(ax, stages, xlims=time_range,
                         label_y=0.05 if fan_strips else 0.02)

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

    if smoothing_time_traces is not None:
        savepath, averaging_label = _prepare_savepath(savepath, smoothing_time_traces)
        axes[0].set_title(f"Moving average: {averaging_label.removeprefix('avg')}")
    if cfg.MAXIMIZE_TIME_TRACE_WINDOWS:
        fig.set_size_inches(*cfg.MAXIMIZED_FIGURE_SIZE, forward=True)
        if show:
            _maximize_window(fig)
    # Keep constrained layout active: GUI maximization can resize again after
    # show(), so a one-shot tight_layout would leave stale outer margins.
    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
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
    connected_points = []
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
        connected_points.append((row["temperature_c"], row["slope"]))
    if len(connected_points) > 1:
        connected_points.sort(key=lambda point: point[0])
        temperatures, slopes = zip(*connected_points)
        ax.plot(temperatures, slopes, color="0.65", linewidth=1.0,
                zorder=1, label="_nolegend_")
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
    ax.text(0.02, 0.98, "Error bars: ±1σ", transform=ax.transAxes,
            ha="left", va="top", fontsize=cfg.PAPER_LABEL_SIZE)
    _style_axis(ax)
    fig.tight_layout()
    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
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


def _prepare_run_groups(temperature_runs, combine_runs=None):
    """Copy run inputs and validate explicit pooling of same-temperature runs."""
    runs = [dict(run) for run in temperature_runs]
    by_id = {}
    for run in runs:
        run_id = str(run.get("run_id", run["temperature"]))
        if run_id in by_id:
            raise ValueError("Each measurement needs a unique run_id: " + run_id)
        run["run_id"] = run_id
        run["fit_group"] = str(run["temperature"])
        run["pooled"] = False
        by_id[run_id] = run
    if combine_runs is not None and not isinstance(combine_runs, dict):
        raise TypeError("combine_runs must be None or a dict mapping labels to run ID lists")
    assigned = set()
    for label, members in (combine_runs or {}).items():
        if not isinstance(label, str) or not label.strip():
            raise ValueError("Combined fit labels must be nonempty strings")
        if not isinstance(members, (list, tuple)) or len(members) < 2:
            raise ValueError("Each combined fit needs at least two run IDs")
        selected = []
        for run_id in members:
            if run_id not in by_id:
                raise ValueError(f"Unknown run ID in combine_runs: {run_id}")
            if run_id in assigned:
                raise ValueError(f"Run appears more than once in combine_runs: {run_id}")
            assigned.add(run_id)
            selected.append(by_id[run_id])
        temperatures = [_nominal_temperature_c(run) for run in selected]
        if not all(np.isfinite(t) and np.isclose(t, temperatures[0]) for t in temperatures):
            raise ValueError("Combined runs must have the same nominal temperature")
        for run_number, run in enumerate(selected, start=1):
            run["fit_group"] = label
            run["pooled"] = True
            run["combined_run_number"] = run_number
    # Identical display labels must never accidentally pool separate measurements.
    owners = {}
    for run in runs:
        label = run["fit_group"]
        if label in owners and not (run["pooled"] and owners[label]["pooled"]):
            raise ValueError("Fit labels must be unique unless explicitly combined: " + label)
        owners[label] = run
    return runs


def temp_dep(
    temperature_runs: Sequence[dict],
    compounds: Sequence[str | tuple[str, str]],
    horus_path: str,
    precursor_compound: str | tuple[str, str],
    product_compound: str | tuple[str, str],
    reaction_partner: str | tuple[str, str] = "HO2",
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
    fit_intercept: bool = True,
    fan_speed_path: Optional[str] = None,
    combine_runs: Optional[dict[str, Sequence[str]]] = None,
    high_res: bool | None = None,
    smoothing_time_traces: timedelta | None = None,
    trace_savepath: Optional[str] = None,
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

    Fit ``y = k*x + b`` (default) or ``y = k*x`` for each temperature using ODR
    with both x and y uncertainties. The kinetic legend displays each slope
    k and its standard error (1 SE). All datetimes must be timezone-naive UTC,
    e.g. ``datetime(2025, 9, 25, 10, 0)``; interval endpoints are inclusive.

    combine_runs optionally pools interval points for a single kinetic fit, e.g.
    ``{"-50 °C (combined)": ["m50_1", "m50_2"]}``. Each input run must have
    a unique ``run_id`` (defaults to its temperature label). None keeps runs
    separate. Only runs with the same nominal temperature may be pooled.
    Trace plots stay separate; each interval keeps its x/y errors and local t_i.
    Pooled kinetic labels use t^i_j: i is the 1-based run position in the
    combine_runs list and j is its original interval number. The combined fit produces one k(T)
    point; slopes are fitted from pooled points, not averaged. Returned means
    retain run_id/run_label; temperature identifies the fit group.

    Prepared runs are checked for finite samples of all selected trace and
    kinetic channels (plus HORUS partners and AP coloring when requested).
    Runs lacking coverage are skipped with a reason; names absent from every
    candidate prepared run raise an error. Combined fits retain available
    members, preserving their original interval labels. No available run is
    an error. Known campaign runs loaded by make_temperature_runs use the
    current intervals.py rather than the preparation-time H5 snapshot.
    high_res=True/False selects the prepared high/low-resolution PTR source
    for all runs, including tuple selections and interval uncertainties.
    None keeps the sources supplied by make_temperature_runs. Raw CSV callers
    choose resolution in make_temperature_runs(..., high_res=...). Missing
    prepared resolutions raise an error rather than substituting another one.

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
        Multiple trace files are combined chronologically; duplicate timestamps
        use the last listed file so overlapping data is not counted twice.

        * ``trace_path`` (required, str, Path, or list of paths): calibrated compound trace
          CSV readable by ``read_trace_timeseries``; values are assumed pptv.
        * ``time_range`` (required, tuple of datetime): ``(start, stop)`` for
          loading and displaying that run; stop must be >= start.
        * ``important_time_ranges`` (required, nonempty sequence of tuples):
          ``[(start1, stop1), (start2, stop2), ...]`` for averaging. Choose
          intervals inside ``time_range``. Their input order defines t_1,
          t_2, etc.; numbering restarts for each temperature.
        * ``temperature_path`` (optional, str or Path): measurement file
          readable by ``load_data`` with a temperature channel in Celsius.
          Used only for ``concentration_unit="molecules/cm3"``.
        * ``temperature_channel`` (optional, str): temperature column name
          in that file; default ``"TE_Calib_3"``.
        * ``temperature_c`` (optional, float): constant Celsius temperature,
          e.g. ``-50.0``, as an alternative to measured temperature.
          In molecules/cm3 mode, supply either ``temperature_path`` or
          ``temperature_c``. If both are supplied, the file takes precedence.

    smoothing_time_traces : timedelta or None, default None
        Centered moving-average window for displayed traces only, applied
        after unit conversion, independently on each instrument's sample times.
        Interval means, SEMs and kinetic fits always use unsmoothed data.
        Trace filenames receive an averaging suffix such as _avg10min.

    compounds : sequence of str
        Required, nonempty names of trace compounds to load and plot, e.g.
        ``["C10H15O4.NH3H+", "C10H16O4.NH3H+"]``. Bare formulas such as
        ``"C10H15O4"`` receive the ``adduct`` suffix. Both the precursor and
        product must occur in this list after adduct expansion when using strings.
        Tuples ``(instrument, exact H5 channel)`` are also accepted and do not
        receive an adduct suffix. Tuple kinetic selections are added automatically.
    horus_path : str or Path
        Required HORUS measurement filename, readable by ``load_data``.
        Must contain ``HO2_ppt`` for HO2 or ``RO2_crosstalk_ppt`` for RO2;
        those measurements are interpreted as pptv.
    precursor_compound : str or (instrument, channel) tuple
        Required precursor name, e.g. ``"C10H15O4.NH3H+"`` or ``"C10H15O4"``.
        Its interval mean is multiplied by the partner mean to form x.
    product_compound : str or (instrument, channel) tuple
        Required product name, e.g. ``"C10H16O4.NH3H+"`` or ``"C10H16O4"``.
        Its interval mean forms y. Uses the same adduct rules as the precursor.
    reaction_partner : str or (instrument, channel) tuple, default "HO2"
        Supports ``"HO2"`` and ``"RO2"`` (case-insensitive; surrounding
        whitespace is ignored), or an exact channel selected with a tuple.
        Tuple selections require prepared ``sources`` in each run, supplied by
        ``make_temperature_runs(..., prepared=True)``. Mixed pptv and molecules/cm3
        sources require ``concentration_unit="molecules/cm3"``; calibrated LToF
        values are kept unchanged. Missing channels or empty coverage raise errors.
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
    trace_savepath : str or Path or None, default None
        Template for trace figures; None means do not save them. Supports
        {run_id}, {precursor_compound}, {reaction_partner}, and {product_compound}. Compounds use
        formulas without ion suffixes. If multiple runs share a template without
        {run_id}, the run ID is appended to prevent overwriting. Smoothing adds
        its suffix afterward. Run dictionaries no longer set trace filenames.
    kinetic_savepath : str or Path or None, default None
        Filename for the combined kinetic figure; None means do not save it.
        The extension selects the format, e.g. .png or .pdf. All three save paths
        support precursor/partner/product placeholders in directories and filenames.
        Missing parent directories are created automatically when saving.
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
        ``"native"`` keeps compound and reaction-partner concentrations in pptv.
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

    fit_intercept : bool, default True
        True fits both slope k and intercept b (y = k*x + b); False fixes b=0.
        The legend reports fitted parameters ±1 SE. Results include intercept,
        intercept_stderr and intercept_unit. Fixed b has zero uncertainty.
        Free-intercept fits require three valid points with distinct x values;
        otherwise points remain visible and the fit is skipped with a message.
        The k(T) plot always uses the fitted slope and its standard error.

    fan_speed_path : str or Path or None, default None
        Optional fan-speed CSV containing Fan_BOT_Speed. Values below 50%
        identify the 12% operating mode, as in run_plotting. Adds a thin gray
        strip at the bottom of each trace axis and one legend entry. Stage
        labels sit above the strip; averaging-interval shading is unchanged.
        None disables the strip. This does not affect kinetic calculations.

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
        One row per fitted temperature: slope, slope_stderr (1 SE),
        intercept, intercept_stderr, intercept_unit, fit_intercept, slope_unit, centered/uncentered R-squared, p-value,
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
    and an x-array that is not entirely zero for a fixed-zero intercept.
    A free intercept requires three usable points with distinct x values.

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
    if smoothing_time_traces is not None:
        if not isinstance(smoothing_time_traces, timedelta):
            raise TypeError("smoothing_time_traces must be a timedelta or None")
        if smoothing_time_traces.total_seconds() <= 0:
            raise ValueError("smoothing_time_traces must be positive")
    validate_concentration_unit(concentration_unit)
    if color_by_ap and ap_path is None:
        raise ValueError("ap_path is required when color_by_ap=True.")
    if not temperature_runs:
        raise ValueError("temperature_runs must contain at least one condition.")
    if not compounds:
        raise ValueError("compounds must contain at least one compound.")

    plt.rcParams.update(cfg.PAPER_RC_PARAMS)
    trace_compounds = list(dict.fromkeys(_signal_name(c, adduct) for c in compounds))
    explicit_signals = any(isinstance(c, tuple) for c in [*compounds, precursor_compound, product_compound, reaction_partner])
    # Expand filename placeholders using the formula without its ion adduct.
    precursor_formula = (precursor_compound[1] if isinstance(precursor_compound, tuple)
                         else precursor_compound).split(".", 1)[0]
    partner_formula = (reaction_partner[1] if isinstance(reaction_partner, tuple)
                       else reaction_partner.strip()).split(".", 1)[0]
    product_formula = (product_compound[1] if isinstance(product_compound, tuple)
                       else product_compound.strip()).split(".", 1)[0]
    if kinetic_savepath is not None:
        kinetic_savepath = (str(kinetic_savepath)
                           .replace("{precursor_compound}", precursor_formula)
                           .replace("{reaction_partner}", partner_formula)
                           .replace("{product_compound}", product_formula))
    if k_temperature_savepath is not None:
        k_temperature_savepath = (str(k_temperature_savepath)
                                 .replace("{precursor_compound}", precursor_formula)
                                 .replace("{reaction_partner}", partner_formula)
                                 .replace("{product_compound}", product_formula))
    precursor_name = _signal_name(precursor_compound, adduct)
    product_name = _signal_name(product_compound, adduct)
    partner_name = (_signal_name(reaction_partner, adduct)
                    if isinstance(reaction_partner, tuple) else reaction_partner)
    for selection, name in ((precursor_compound, precursor_name), (product_compound, product_name)):
        if name not in trace_compounds:
            if isinstance(selection, tuple):
                trace_compounds.append(name)
            else:
                raise ValueError(f"Compound {name!r} is not in compounds.")
    if isinstance(reaction_partner, tuple) and partner_name not in trace_compounds:
        trace_compounds.append(partner_name)

    if stages_path is not None and stages is not None:
        load_stages(stages_path)

    if high_res is not None:
        from .prepared import select_prepared_resolution
        selected_runs = []
        for original in temperature_runs:
            run = dict(original)
            if "sources" in run:
                run["sources"] = dict(run["sources"])
                run["sources"]["ptr"] = select_prepared_resolution(run["sources"]["ptr"], high_res)
                run["trace_path"] = run["sources"]["ptr"]
            else:
                run["trace_path"] = select_prepared_resolution(run["trace_path"], high_res)
            selected_runs.append(run)
        temperature_runs = selected_runs
    temperature_runs = _prepare_run_groups(temperature_runs, combine_runs)
    temperature_runs = _available_temperature_runs(
        temperature_runs, trace_compounds, partner_name, color_by_ap)
    plot_temperatures = {}
    for run in temperature_runs:
        label = run["fit_group"]
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
        if "stages_path" in run:
            load_stages(run["stages_path"])
        time_range = _validate_time_ranges([run["time_range"]])[0]
        intervals = _validate_time_ranges(run["important_time_ranges"])
        if not intervals:
            raise ValueError(
                f"important_time_ranges must not be empty for {temperature}."
            )

        temperature_c = None
        if concentration_unit == "molecules/cm3":
            if "temperature_path" in run:
                temperature_name = run.get("temperature_channel", "TE_Calib_3")
                temperature_c = load_data(run["temperature_path"], channels=[temperature_name])[temperature_name]
            else:
                temperature_c = run.get("temperature_c")
        if explicit_signals:
            data = _read_selected_signals(run, trace_compounds, time_range,
                                          concentration_unit, temperature_c, pressure_pa)
            if concentration_unit == "native" and data.attrs["native_unit"] != "pptv":
                raise ValueError("Instrument selections in native mode must use pptv; use concentration_unit='molecules/cm3' for LToF")
        else:
            data = read_trace_timeseries(run["trace_path"], trace_compounds, time_range)
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

        if concentration_unit == "molecules/cm3" and not explicit_signals:
            if temperature_c is None:
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
            horus_path=run.get("horus_path", horus_path),
            time_range=time_range,
            important_time_ranges=intervals,
            precursor_compound=precursor_name,
            reaction_partner=partner_name,
            concentration_unit=concentration_unit,
            temperature_c=temperature_c,
            pressure_pa=pressure_pa,
        )
        if color_by_ap:
            ap_source = run.get("ap_path", ap_path)
            ap_data = load_data(ap_source, time_range, channels=["AP"])
            if "AP" not in ap_data.channels:
                raise KeyError("AP channel not found in " + str(ap_path))
            ap_means = calculate_channel_interval_means(
                convert_channel_concentration(ap_data["AP"], "ppbv", "pptv"),
                intervals, column_name="AP_pptv",
            )
            means = means.merge(ap_means, on="range_index", how="left",
                                validate="one_to_one")
        means.insert(0, "temperature", run["fit_group"])
        means["run_id"] = run["run_id"]
        means["run_label"] = temperature
        means["interval_label"] = [
            (rf"$t^{{{run['combined_run_number']}}}_{{{int(i)}}}$"
             if run["pooled"] else _interval_label(i))
            for i in means["range_index"]
        ]
        all_means.append(means)

        if kinetic_x_column is None:
            kinetic_x_column = current_x_column
        elif kinetic_x_column != current_x_column:
            raise RuntimeError("Inconsistent kinetic x-column between runs.")

        fan_source = run.get("fan_speed_path", fan_speed_path) if fan_speed_path is not None else None
        fan_speed_data = (load_data(fan_source, time_range, channels=["Fan_BOT_Speed"])
                          if fan_source is not None else None)
        current_trace_savepath = None
        if trace_savepath is not None:
            current_trace_savepath = (str(trace_savepath)
                                     .replace("{run_id}", str(run["run_id"]))
                                     .replace("{precursor_compound}", precursor_formula)
                                     .replace("{reaction_partner}", partner_formula)
                                     .replace("{product_compound}", product_formula))
            if "{run_id}" not in str(trace_savepath) and len(temperature_runs) > 1:
                current_trace_savepath = _append_to_filename(current_trace_savepath, run["run_id"])
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
                    savepath=current_trace_savepath,
                    stages=stages,
                    label_intervals=label_intervals,
                    concentration_unit=concentration_unit,
                    yscale=yscale,
                    fan_speed_data=fan_speed_data,
                    smoothing_time_traces=smoothing_time_traces,
                )
            )

    combined_means = pd.concat(all_means, ignore_index=True)
    combined_means.attrs["concentration_unit"] = concentration_unit
    combined_means.attrs["AP_pptv_unit"] = "pptv"

    partner_label = partner_name.strip().upper()
    partner_plot_label = {"HO2": r"HO$_2$", "RO2": r"RO$_2$"}.get(partner_label, _chemical_label(partner_name))
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
            fit_intercept=fit_intercept,
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
