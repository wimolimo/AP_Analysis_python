"""Steady-state radical loss-budget calculation and plotting.

Bounds assume direct AP ozonolysis only and steady state; central kinetic
coefficients are used without calibration or kinetic uncertainty propagation.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as cfg
from cloud18_ap_analysis.interpolation import resample
from cloud18_ap_analysis.loading import load_data
from cloud18_ap_analysis.prepared import prepared_channel_info
from cloud18_ap_analysis.temp_dependence import (calculate_interval_means, _style_axis,
    _plot_temperature_trace, _prepare_run_groups)
from cloud18_ap_analysis.stages import load_stages
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from cloud18_ap_analysis.units import (convert_channel_concentration,
                                      AVOGADRO_CONSTANT, GAS_CONSTANT)

def _align(channel, times, max_interpolation_gap_s):
    """Reuse resample, with an additional configurable finite-data gap limit."""
    out = resample(channel, times)
    finite = np.isfinite(channel.time) & np.isfinite(channel.values)
    measured = np.sort(np.asarray(channel.time)[finite])
    if len(measured) < 2:
        return out
    right = np.searchsorted(measured, times, side="left").clip(0, len(measured) - 1)
    left = (right - 1).clip(0, len(measured) - 1)
    gap = ((measured[right] - measured[left] > max_interpolation_gap_s * 1000)
           & (times > measured[left]) & (times < measured[right]))
    out.values[gap] = np.nan
    return out


def _read_channel(run, instrument, name):
    source = run["sources"][instrument]
    catalog = prepared_channel_info(source)
    if name not in catalog:
        raise KeyError(f"{run['run_id']}: missing {instrument}:{name}")
    return load_data(source, run["time_range"], channels=[name])[name], catalog[name]["unit"]


def _native_trace_frame(channels, *, pressure_pa, max_interpolation_gap_s):
    """Build full-run traces on each channel's own finite sample timestamps.

    Concentrations are never interpolated. Only temperature is aligned to the
    native timestamps for mixing-ratio conversion. Long gaps between usable
    measurements receive a NaN line break. sample_times excludes timestamps
    introduced solely by joining other instruments.
    """
    temperature = channels["temperature_c"][0]
    traces, sample_times = [], {}
    for name, (channel, unit) in channels.items():
        if name == "temperature_c" or unit == "molecules/cm3":
            values = np.asarray(channel.values, dtype=float)
        else:
            temp = _align(temperature, channel.time, max_interpolation_gap_s)
            fraction = {"ppt": 1e-12, "pptv": 1e-12, "ppb": 1e-9, "ppbv": 1e-9}[unit]
            values = channel.values * fraction * pressure_pa / (
                GAS_CONSTANT * (temp.values + 273.15)) * 1e-6 * AVOGADRO_CONSTANT
        finite = np.isfinite(channel.time) & np.isfinite(values)
        trace = pd.Series(values[finite], index=pd.to_datetime(
            channel.time[finite], unit="ms"), name=name).groupby(level=0).mean().sort_index()
        # Connect actual samples, but never draw across a long missing period.
        gaps = np.flatnonzero((trace.index[1:] - trace.index[:-1])
                              > pd.Timedelta(seconds=max_interpolation_gap_s))
        if gaps.size:
            breaks = trace.index[gaps] + (trace.index[gaps + 1] - trace.index[gaps]) / 2
            trace = trace.reindex(trace.index.union(breaks))
        sample_times[name] = trace.index
        traces.append(trace)
    frame = pd.concat(traces, axis=1).sort_index().rename_axis("Time").reset_index()
    frame.attrs["sample_times"] = sample_times
    return frame


def _separate_interval_budget(run, channels, *, pressure_pa, wall_lifetime_s,
                              dilution_rate_s, production_yield_max,
                              max_interpolation_gap_s):
    """Average native samples independently, then convert and calculate rates.

    All finite concentration samples, including zeros and negatives, contribute
    to the means. Missing timestamps are not reconstructed.
    Derived-rate SEMs are NaN because uncertainties are not propagated.
    Returned traces span the full run at native concentration timestamps.
    """
    tables = {}
    ranges = run["important_time_ranges"]
    for name, (channel, unit) in channels.items():
        values = np.asarray(channel.values, dtype=float).copy()
        values[~np.isfinite(values)] = np.nan
        if name == "temperature_c" and np.any(values + 273.15 <= 0):
            raise ValueError("Nonpositive absolute temperature")
        frame = pd.DataFrame({"Time": pd.to_datetime(channel.time, unit="ms"), name: values})
        table = calculate_interval_means(frame, [name], "Time", ranges)
        tables[name] = table
    means = tables["R"].copy()
    for name, table in tables.items():
        if name != "R":
            for column in (name, f"{name}_n_points", f"{name}_sem"):
                means[column] = table[column]
    temp_k = means["temperature_c"] + 273.15
    fractions = {"ppt": 1e-12, "pptv": 1e-12, "ppb": 1e-9, "ppbv": 1e-9}
    air = pressure_pa / (8.31446261815324 * temp_k) * 1e-6 * 6.02214076e23
    for name in ("R", "AP", "O3", "HO2"):
        unit = channels[name][1]
        if unit == "molecules/cm3":
            factor = 1.0
        else:
            if unit not in fractions:
                raise ValueError(f"Unsupported concentration unit: {unit}")
            factor = fractions[unit] * air
        if name == "AP":
            means["AP_pptv"] = means[name] * factor / air / 1e-12
            means["AP_pptv_sem"] = means[f"{name}_sem"] * factor / air / 1e-12
            means["AP_pptv_n_points"] = means[f"{name}_n_points"]
        means[name] *= factor
        means[f"{name}_sem"] *= factor
    means["P_max"] = production_yield_max * 1.4e-15 * np.exp(-833.0 / temp_k) * means["AP"] * means["O3"]
    means["wall_loss"] = means["R"] / wall_lifetime_s
    means["dilution_loss"] = dilution_rate_s * means["R"]
    means["HO2_loss"] = 2.91e-13 * np.exp(1300.0 / temp_k) * (1 - np.exp(-0.245 * 10)) * means["HO2"] * means["R"]
    means["known_loss"] = means["wall_loss"] + means["dilution_loss"] + means["HO2_loss"]
    means["residual"] = means["P_max"] - means["known_loss"]
    counts = pd.concat([t[f"{name}_n_points"] for name, t in tables.items()], axis=1).min(axis=1)
    for name in ("P_max", "wall_loss", "dilution_loss", "HO2_loss", "known_loss", "residual"):
        means[f"{name}_sem"] = np.nan
        means[f"{name}_n_points"] = counts
    means["input_n_points_min"] = counts
    frame = _native_trace_frame(channels, pressure_pa=pressure_pa,
                                max_interpolation_gap_s=max_interpolation_gap_s)
    return means, frame


def calculate_loss_budget(
    run, *, pressure_pa=97000.0, wall_lifetime_s=480.0,
    dilution_rate_s=2e-4, production_yield_max=1.0,
    radical_channel="C10H15O4.NH3H+", temperature_channel="TE_Calib_3",
    ozone_channel="o3", max_interpolation_gap_s=1800.0,
    return_timeseries=False, average_separately=False,
):
    """Return interval means, component losses, and conditional bounds.

    Uses shared UTC intervals in the prepared run. All finite concentration
    samples, including zeros and negatives, are retained for every channel.
    Missing/nonfinite values are skipped. Intervals with nonpositive mean R
    are still averaged and reported, but excluded from fits.
    average_separately=False (default) aligns samples to radical timestamps,
    calculates instantaneous rates, then averages. Interpolation cannot bridge
    gaps beyond max_interpolation_gap_s.
    True averages each instrument's native samples over the same UTC intervals,
    converts concentrations using mean temperature, then calculates rates from
    those means. No interpolation is used for interval calculations. Each input
    needs at least two valid samples; derived-rate SEMs are not propagated.
    Returned traces in this mode contain full-run input concentrations on each
    channel's native timestamps rather than derived rates. Only temperature is
    interpolated for trace unit conversion; interval calculations use mean T.
    """
    if not all(np.isfinite(v) for v in (pressure_pa, wall_lifetime_s,
            dilution_rate_s, production_yield_max, max_interpolation_gap_s)):
        raise ValueError("Settings must be finite")
    if wall_lifetime_s <= 0 or pressure_pa <= 0 or dilution_rate_s < 0 or not 0 <= production_yield_max <= 1:
        raise ValueError("Invalid physical settings")
    if not 0 < max_interpolation_gap_s <= 1800:
        raise ValueError("Gap must be in (0, 1800] seconds")
    temperature, unit = _read_channel(run, "temperature", temperature_channel)
    if unit != "degC":
        raise ValueError(f"Expected degC for {temperature_channel}, got {unit}")
    radical, radical_unit = _read_channel(run, "ptr", radical_channel)
    if average_separately:
        channels = dict(R=(radical, radical_unit),
                        temperature_c=(temperature, unit),
                        AP=_read_channel(run, "vocus", "AP"),
                        O3=_read_channel(run, "ozone", ozone_channel),
                        HO2=_read_channel(run, "horus", "HO2_ppt"))
        means, frame = _separate_interval_budget(
            run, channels, pressure_pa=pressure_pa, wall_lifetime_s=wall_lifetime_s,
            dilution_rate_s=dilution_rate_s, production_yield_max=production_yield_max,
            max_interpolation_gap_s=max_interpolation_gap_s)
    else:
        times = radical.time
        temp = _align(temperature, times, max_interpolation_gap_s)
        temp_k = temp.values + 273.15
        if np.any(np.isfinite(temp_k) & (temp_k <= 0)):
            raise ValueError("Nonpositive absolute temperature")

        def concentration(channel, source_unit):
            # Align first so every conversion and rate uses the same measured T.
            channel = _align(channel, times, max_interpolation_gap_s)
            if source_unit == "molecules/cm3":
                return channel.values
            return convert_channel_concentration(
                channel, source_unit, "molecules/cm3", temperature_c=temp,
                pressure_pa=pressure_pa).values

        r = concentration(radical, radical_unit)
        ap_channel, ap_unit = _read_channel(run, "vocus", "AP")
        ap = concentration(ap_channel, ap_unit)
        ap_pptv = _align(convert_channel_concentration(
            ap_channel, ap_unit, "pptv"), times, max_interpolation_gap_s).values
        ozone = concentration(*_read_channel(run, "ozone", ozone_channel))
        ho2 = concentration(*_read_channel(run, "horus", "HO2_ppt"))
        # Retain finite values of either sign for every concentration channel.
        valid = (np.isfinite(r) & np.isfinite(ap)
                 & np.isfinite(ozone) & np.isfinite(ho2)
                 & np.isfinite(temp_k))
        k_ap_o3 = 1.4e-15 * np.exp(-833.0 / temp_k)
        k_ro2_ho2 = 2.91e-13 * np.exp(1300.0 / temp_k) * (1 - np.exp(-0.245 * 10))
        production = production_yield_max * k_ap_o3 * ap * ozone
        wall = r / wall_lifetime_s
        dilution = dilution_rate_s * r
        ho2_loss = k_ro2_ho2 * ho2 * r
        known = wall + dilution + ho2_loss
        quantities = dict(R=r, AP_pptv=ap_pptv, AP=ap, O3=ozone, HO2=ho2, temperature_c=temp.values,
                          P_max=production, wall_loss=wall, dilution_loss=dilution,
                          HO2_loss=ho2_loss, known_loss=known, residual=production-known)
        frame = pd.DataFrame({name: np.where(valid, values, np.nan)
                              for name, values in quantities.items()})
        frame.insert(0, "Time", pd.to_datetime(times, unit="ms"))
        means = calculate_interval_means(frame, list(quantities), "Time",
                                         run["important_time_ranges"])
    means["status"] = np.select(
        [means.get("input_n_points_min", means["R_n_points"]) < 2, means["R"] <= 0,
         means["residual"] < 0],
        ["insufficient_data", "nonpositive_mean_R", "negative_budget"], default="ok")
    means.insert(0, "run_id", run["run_id"])
    means.insert(1, "temperature", run["temperature"])
    means["nominal_temperature_c"] = run["plot_temperature_c"]
    means["production_yield_max"] = production_yield_max
    means["wall_lifetime_s"] = wall_lifetime_s
    means["dilution_rate_s-1"] = dilution_rate_s
    means["pressure_pa"] = pressure_pa
    means["average_separately"] = average_separately
    means["max_interpolation_gap_s"] = max_interpolation_gap_s
    return (means, frame) if return_timeseries else means


def fit_loss_budget(means, *, fit_intercept=False):
    """Fit y=lambda*x by unweighted least squares for each nominal temperature.

    Explicit fit groups are pooled. Negative budgets are retained. The reported SE
    describes interval scatter with x treated as fixed; it is not a propagated
    measurement uncertainty. A fitted slope is not a strict pointwise bound.
    """
    rows = []
    for label, candidates in means.groupby("temperature", sort=False):
        temperature_c = float(candidates["nominal_temperature_c"].iloc[0])
        valid = candidates[candidates["status"].isin(["ok", "negative_budget"])]
        valid = valid[np.isfinite(valid[["R", "residual"]]).all(axis=1) & (valid["R"] > 0)]
        x = valid["R"].to_numpy(dtype=float)
        y = valid["residual"].to_numpy(dtype=float)
        slope = se = intercept = intercept_se = np.nan
        parameters = 2 if fit_intercept else 1
        identifiable = np.ptp(x) > 0 if len(x) and fit_intercept else np.any(x != 0)
        if len(x) > parameters and identifiable:
            # Scale x to keep the two-column regression well conditioned.
            scale = float(np.max(np.abs(x)))
            design = np.column_stack([x/scale, np.ones(len(x))]) if fit_intercept else (x/scale)[:, None]
            beta = np.linalg.lstsq(design, y, rcond=None)[0]
            variance = float(np.sum((y-design@beta)**2) / (len(x)-parameters))
            covariance = variance*np.linalg.inv(design.T@design)
            slope, se = float(beta[0]/scale), float(np.sqrt(covariance[0, 0])/scale)
            intercept = float(beta[1]) if fit_intercept else 0.0
            intercept_se = float(np.sqrt(covariance[1, 1])) if fit_intercept else 0.0
        rows.append(dict(temperature=label, temperature_c=temperature_c,
                         slope=slope, slope_stderr=se, slope_unit="s^-1", intercept=intercept,
                         intercept_stderr=intercept_se, intercept_unit="molecules cm^-3 s^-1",
                         fit_intercept=fit_intercept,
                         n_points=len(valid), n_excluded=len(candidates)-len(valid),
                         run_ids=",".join(candidates["run_id"].unique()),
                         fit_method="unweighted least squares"))
    return pd.DataFrame(rows)


def plot_loss_budget(means, *, production_yield_max=1.0, label_intervals=False,
                     regression_results=None, color_by_ap=False, fit_intercept=False):
    """Plot temperature fits using the shared kinetic figure layout and style."""
    plt.rcParams.update(cfg.PAPER_RC_PARAMS)
    fig, ax = plt.subplots(figsize=cfg.KINETIC_FIGSIZE)
    results = fit_loss_budget(means, fit_intercept=fit_intercept) if regression_results is None else regression_results
    drawable = means[means["status"].isin(["ok", "negative_budget"])]
    ap_values = drawable["AP_pptv"].to_numpy(dtype=float) if color_by_ap else np.array([])
    ap_values = ap_values[np.isfinite(ap_values)]
    norm = None
    if ap_values.size:
        lo, hi = ap_values.min(), ap_values.max()
        margin = max(abs(lo)*0.05, 0.01) if lo == hi else 0.0
        norm = Normalize(lo-margin, hi+margin)
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_bad("0.6")
    for index, result in results.iterrows():
        group = means[(means["temperature"] == result["temperature"])
                      & means["status"].isin(["ok", "negative_budget"])]
        group = group[np.isfinite(group[["R", "residual"]]).all(axis=1) & (group["R"] > 0)]
        if group.empty:
            continue
        color = cfg.PAPER_LINE_COLORS[index % len(cfg.PAPER_LINE_COLORS)]
        slope, se = result["slope"], result["slope_stderr"]
        label = result["temperature"]
        if np.isfinite(slope):
            magnitude = max(abs(slope), abs(se))
            exponent = int(np.floor(np.log10(magnitude))) if magnitude > 0 else 0
            scale = 10.0**exponent
            label += ("\n" + rf"$\lambda=({slope/scale:.3g} \pm {se/scale:.2g})"
                      + rf"\times 10^{{{exponent}}}\ \mathrm{{s^{{-1}}}}$")
            x_fit = np.linspace(0.0, group["R"].max()*1.05, 200)
            if result["fit_intercept"]:
                intercept = result["intercept"]
                intercept_se = result["intercept_stderr"]
                intercept_magnitude = max(abs(intercept), abs(intercept_se))
                intercept_exponent = int(np.floor(np.log10(intercept_magnitude))) if intercept_magnitude > 0 else 0
                intercept_scale = 10.0**intercept_exponent
                label += ("\n" + rf"$b=({intercept/intercept_scale:.3g} \pm {intercept_se/intercept_scale:.2g})"
                          + rf"\times 10^{{{intercept_exponent}}}\ \mathrm{{molecules\,cm^{{-3}}\,s^{{-1}}}}$")
            ax.plot(x_fit, slope*x_fit + result["intercept"], color=color, linewidth=cfg.PAPER_LINE_WIDTH, zorder=2)
        marker = ("o", "s", "^", "D", "v", "P", "X")[index % 7] if color_by_ap else "o"
        if color_by_ap:
            ax.scatter(group["R"], group["residual"],
                       c=np.ma.masked_invalid(group["AP_pptv"].to_numpy(dtype=float)),
                       cmap=cmap, norm=norm, plotnonfinite=True,
                       marker=marker, s=cfg.KINETIC_MARKER_SIZE,
                       edgecolors=color, linewidths=0.8, zorder=3)
            ax.plot([], [], marker=marker, color=color, markerfacecolor="white",
                    label=label, linewidth=cfg.PAPER_LINE_WIDTH)
        else:
            ax.plot(group["R"], group["residual"], linestyle="none", marker=marker,
                    markersize=np.sqrt(cfg.KINETIC_MARKER_SIZE), markeredgewidth=0,
                    color=color, label=label, zorder=3)
        if label_intervals:
            for _, row in group.iterrows():
                ax.annotate(row.get("interval_label", f"$t_{{{row['range_index']}}}$"),
                            (row["R"], row["residual"]), xytext=(5, 5),
                            textcoords="offset points", color=color,
                            fontsize=cfg.PAPER_LABEL_SIZE, ha="left", va="bottom")
    if color_by_ap and norm is not None:
        fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Mean chamber AP [pptv]")
    ax.axhline(0, color="0.75", linewidth=cfg.PAPER_LINE_WIDTH, zorder=1)
    # Linear axes retain zero and negative loss budgets.
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$\mathrm{C}_{10}\mathrm{H}_{15}\mathrm{O}_4$ [molecules cm$^{-3}$]")
    ax.set_ylabel(r"$P_{\mathrm{R}}-L_{\mathrm{known}}$ [molecules cm$^{-3}$ s$^{-1}$]")
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    _style_axis(ax)
    fig.tight_layout()
    return fig, ax



def loss_budget(
    temperature_runs, *, pressure_pa=97000.0, wall_lifetime_s=480.0,
    dilution_rate_s=2e-4, production_yield_max=1.0,
    radical_channel="C10H15O4.NH3H+", temperature_channel="TE_Calib_3",
    ozone_channel="o3", max_interpolation_gap_s=1800.0,
    label_intervals=False, color_by_ap=False, overlay=False,
    fit_intercept=False, combine_runs=None, average_separately=False,
    show=True, interactive=True, smoothing_time_traces=None,
    stages="type", fan_speed=True, trace_yscale="linear",
    savepath=None, trace_savepath=None, lambda_temperature_savepath=None,
    table_savepath=None, fit_savepath=None, print_results=True,
):
    """Return traces, means, fits, loss figure/axis, and lambda(T) figure/axis.

    show accepts True, False, 'traces', 'loss_kinetics', 'k_loss(T)', or a list.
    Saving is independent of show. combine_runs uses the same explicit mapping
    as temp_dep; None keeps runs separate. fit_intercept=True fits an offset.
    interactive=False suppresses all windows. Each nominal temperature gets
    an unweighted fit; explicit combine_runs groups pool measurements. Only fitted
    temperature slopes are calculated and printed, never interval y/x ratios.
    Trace tags and loss-point tags share the same interval numbering.
    average_separately=True averages native measurements independently over each
    interval, then converts using mean temperature and calculates rates. False
    (default) averages instantaneous rates after alignment to radical timestamps.
    Both modes retain finite zero and negative samples for every concentration
    channel when averaging; only missing/nonfinite samples are skipped.
    Time traces plot each concentration at its own timestamps over the full run,
    independently of average_separately. Only temperature is interpolated for
    trace unit conversion. See calculate_loss_budget for sample counts.
    """
    available = {"traces", "loss_kinetics", "k_loss(T)"}
    if isinstance(show, bool):
        selected = available if show else set()
    elif isinstance(show, str):
        selected = {show}
    elif isinstance(show, (list, tuple)) and all(isinstance(name, str) for name in show):
        selected = set(show)
    else:
        raise TypeError("show must be a bool, plot name, or list of plot names")
    if selected - available:
        raise ValueError("show must select from 'traces', 'loss_kinetics', 'k_loss(T)'")
    runs = _prepare_run_groups(list(temperature_runs), combine_runs)
    if not runs:
        raise ValueError("temperature_runs must not be empty")
    traces, tables = [], []
    for run in runs:
        means, frame = calculate_loss_budget(
            run, pressure_pa=pressure_pa, wall_lifetime_s=wall_lifetime_s,
            dilution_rate_s=dilution_rate_s, production_yield_max=production_yield_max,
            radical_channel=radical_channel, temperature_channel=temperature_channel,
            ozone_channel=ozone_channel, max_interpolation_gap_s=max_interpolation_gap_s,
            return_timeseries=True,
            average_separately=average_separately)
        means["run_label"] = run["temperature"]
        means["temperature"] = run["fit_group"]
        means["interval_label"] = [
            (rf"$t^{{{run['combined_run_number']}}}_{{{i}}}$" if run["pooled"]
             else rf"$t_{{{i}}}$") for i in means["range_index"]]
        tables.append(means)
        trace_path = (str(trace_savepath).format(run_id=run["run_id"])
                      if trace_savepath is not None else None)
        if "traces" in selected or trace_path is not None:
            if stages is not None:
                load_stages(run["stages_path"])
            fans = (load_data(run["fan_speed_path"], run["time_range"],
                              channels=["Fan_BOT_Speed"]) if fan_speed else None)
            if "sample_times" not in frame.attrs:
                channels = dict(
                    R=_read_channel(run, "ptr", radical_channel),
                    temperature_c=_read_channel(run, "temperature", temperature_channel),
                    AP=_read_channel(run, "vocus", "AP"),
                    O3=_read_channel(run, "ozone", ozone_channel),
                    HO2=_read_channel(run, "horus", "HO2_ppt"))
                frame = _native_trace_frame(channels, pressure_pa=pressure_pa,
                                           max_interpolation_gap_s=max_interpolation_gap_s)
            trace_frame = frame.rename(columns={"R": radical_channel})
            trace_frame.attrs["sample_times"] = {
                radical_channel if name == "R" else name: times
                for name, times in frame.attrs["sample_times"].items()}
            trace_fig, trace_axes = _plot_temperature_trace(
                trace_frame, [radical_channel, "HO2", "AP", "O3"], "Time",
                run["time_range"], run["important_time_ranges"], run["temperature"],
                overlay=overlay, show=interactive and "traces" in selected,
                savepath=trace_path, stages=stages, label_intervals=label_intervals,
                concentration_unit="molecules/cm3", yscale=trace_yscale,
                fan_speed_data=fans, smoothing_time_traces=smoothing_time_traces)
            traces.append((run["run_id"], trace_fig, trace_axes))
    means = pd.concat(tables, ignore_index=True)
    regression_results = fit_loss_budget(means, fit_intercept=fit_intercept)
    fig, ax = plot_loss_budget(means, production_yield_max=production_yield_max,
                              label_intervals=label_intervals, color_by_ap=color_by_ap,
                              regression_results=regression_results, fit_intercept=fit_intercept)
    if print_results:
        print(regression_results.to_string(index=False))
        print("Temperature fits: " + ("y = lambda*x + b" if fit_intercept else "y = lambda*x") + " (unweighted).")
        print("SE describes interval scatter; kinetic/calibration uncertainties are not propagated.")
        print("Intervals with insufficient data or nonpositive mean R excluded; negative budgets retained. Fits are not strict pointwise bounds.")
    for output, table in [(fit_savepath, regression_results), (table_savepath, means)]:
        if output is not None:
            path = Path(output)
            path.parent.mkdir(parents=True, exist_ok=True)
            table.to_csv(path, index=False)
    if savepath is not None:
        path = Path(savepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
    if interactive and "loss_kinetics" in selected:
        plt.show()
    else:
        plt.close(fig)
    lambda_fig = lambda_ax = None
    if "k_loss(T)" in selected or lambda_temperature_savepath is not None:
        lambda_fig, lambda_ax = plt.subplots(figsize=cfg.KINETIC_FIGSIZE)
        for index, row in regression_results.iterrows():
            if not np.isfinite([row["slope"], row["slope_stderr"]]).all():
                continue
            lambda_ax.errorbar(row["temperature_c"], row["slope"], yerr=row["slope_stderr"],
                fmt="o", markersize=np.sqrt(cfg.KINETIC_MARKER_SIZE), capsize=3,
                color=cfg.PAPER_LINE_COLORS[index % len(cfg.PAPER_LINE_COLORS)],
                label=row["temperature"])
        ordered = regression_results.sort_values("temperature_c")
        lambda_ax.plot(ordered["temperature_c"], ordered["slope"],
                       color="0.65", linewidth=cfg.PAPER_LINE_WIDTH, zorder=1)
        lambda_ax.set_xlabel(r"Temperature [$^\circ$C]")
        lambda_ax.set_ylabel(r"Fitted loss rate $\lambda$ [s$^{-1}$]")
        lambda_ax.set_title("Additional loss rate versus temperature")
        lambda_ax.legend()
        _style_axis(lambda_ax)
        lambda_fig.tight_layout()
        if lambda_temperature_savepath is not None:
            path = Path(lambda_temperature_savepath)
            path.parent.mkdir(parents=True, exist_ok=True)
            lambda_fig.savefig(path, dpi=cfg.PAPER_SAVE_DPI, bbox_inches="tight")
        if interactive and "k_loss(T)" in selected:
            plt.show()
        else:
            plt.close(lambda_fig)
    return traces, means, regression_results, fig, ax, lambda_fig, lambda_ax
