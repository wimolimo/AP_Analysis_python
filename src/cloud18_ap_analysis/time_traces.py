"""Plot selected compound time traces with CLOUD stage information."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

from .loading import load_data
from .stages import load_stages
from .plotting import plot_data, _plot_window_kind
from .models import Channel, Dataset
from .trace_selection import available_channels, resolve_channel
from .prepared import PreparedSource, prepared_channel_info
from .units import (convert_dataset_concentrations, concentration_label,
                    concentration_yscale, validate_concentration_unit, STANDARD_PRESSURE_PA)


def _plot_trace_data(*args, **kwargs):
    with _plot_window_kind("time_traces"):
        return plot_data(*args, **kwargs)


def plot_time_traces(
    trace_path: str | Path | list[str | Path] | dict,
    compounds: list[str] | dict[str, list[str]],
    time_range: tuple[datetime, datetime],
    stages_path: str | Path | None = None,
    savepath: str | Path | None = None,
    smoothing: timedelta | float | None = None,
    interactive: bool = True,
    stages: str | bool | None = "onlytype",
    stacked: bool = False,
    ylabel: str | None = None,
    ylims=None,
    yscale: str | None = None,
    fan_speed_path: str | Path | None = None,
    important_time_ranges=None,
    label_intervals: bool = False,
    concentration_unit: str = "native",
    source_unit: str = "pptv",
    temperature_path: str | Path | None = None,
    temperature_channel: str = "TE_Calib_3",
    temperature_c=None,
    pressure_pa=STANDARD_PRESSURE_PA,
    adduct_preferences=None,
    source_units=None,
    high_res: bool | None = None,
):
    """
    Plot time traces of selected compounds.

    Parameters
    ----------
    trace_path
        Path to one trace CSV, or a nonempty list/tuple of files with the same
        channels. Files are combined chronologically before conversion/smoothing;
        duplicate timestamps use the last listed file (no double counting).

    compounds
        Channel names, or {instrument: [neutral formula or full channel name]}.
        For grouped selection, trace_path is a mapping of instrument to source.
        Two different units share a panel with separate y axes (cps on the right).
        stacked=True gives one panel per signal.

    adduct_preferences
        Optional {instrument: [adducts in priority order]} for grouped selection.
        Ambiguous matches without a preferred ion raise an actionable error.

    source_units
        Optional {instrument: unit} for raw grouped sources. Prepared channels
        use their stored units. Counts cannot be converted to concentration.

    high_res
        True/False selects prepared high/low-resolution PTR measurements.
        None keeps the supplied source. Other instruments are unchanged.
        Raw CSV callers select paths through runs.py instead.

    time_range
        Tuple ``(start, stop)`` defining the displayed time range.

    stages_path
        Path to the CLOUD ``stages.txt`` file. If supplied, the stages
        are loaded before plotting.

    savepath
        Optional output path for the figure.

    smoothing
        Optional smoothing window. Can be a ``timedelta`` or seconds.

    interactive
        If True, show the Matplotlib figure interactively.

    stages
        Stage-label mode passed to ``plot_data``.
        ``"onlytype"`` shows only the stage type.

    stacked
        If False, all selected compounds are overlaid in one panel.
        If True, each compound gets its own panel.

    ylabel
        Optional axis-label override. None builds a label from the concentration unit.

    ylims
        Optional y-axis limits.

    yscale
        None chooses log for molecules/cm3 and linear for native mixing ratios.
        Explicit "log" or "linear" overrides this. Log displays positive values only.

    fan_speed_path
        Optional fan-speed CSV file (same format as run_plotting). Load
        Fan_BOT_Speed for time_range and shade periods below 50 percent gray.
        Applies to overlay and stacked traces. None disables shading.
        Uses a bottom strip when important_time_ranges is supplied.

    important_time_ranges
        Optional list of UTC (start, stop) datetime pairs to shade. None hides
        intervals. Visual references only; no interval averaging is performed.

    label_intervals
        True adds t_i labels in original input order on the first subplot only
        (or the overlay axes); False shows shading only.
    concentration_unit
        "native" (default) keeps the input mixing ratio; "molecules/cm3"
        converts selected compounds before smoothing using the shared temp_dep conversion.

    source_unit
        Input concentration unit: "pptv" (default) or "ppbv"; "ppt"/"ppb"
        aliases are accepted. Input must be calibrated concentration, not counts.

    temperature_path, temperature_channel
        Optional temperature CSV and channel (default "TE_Calib_3", Celsius).
        Measurements are interpolated to trace timestamps without extrapolation.
        Takes precedence over temperature_c; only loaded for conversion.

    temperature_c
        Constant Celsius temperature or temperature Channel, required for conversion
        if temperature_path is absent. Missing temperature coverage yields NaN.

    pressure_pa
        Pressure in Pa (constant or Channel); default 101325 Pa, normal pressure.
    """

    validate_concentration_unit(concentration_unit)
    if high_res is not None:
        from .prepared import select_prepared_resolution
        if isinstance(trace_path, dict):
            trace_path = dict(trace_path)
            if "ptr" in trace_path:
                trace_path["ptr"] = select_prepared_resolution(trace_path["ptr"], high_res)
        else:
            trace_path = select_prepared_resolution(trace_path, high_res)

    if isinstance(compounds, dict):
        return _plot_instruments(
            trace_path, compounds, time_range, adduct_preferences or {}, source_units or {},
            concentration_unit, temperature_path, temperature_channel, temperature_c, pressure_pa,
            stacked, ylabel, stages_path, fan_speed_path,
            dict(xlims=time_range, ylims=ylims, yscale=concentration_yscale(concentration_unit, yscale),
                 smoothing=smoothing, stages=stages, important_time_ranges=important_time_ranges,
                 label_intervals=label_intervals, savepath=savepath, interactive=interactive),
        )

    # Load only the requested time period
    data = load_data(trace_path, time_range, channels=compounds)

    # Load CLOUD stages so plotting.py can draw them
    if stages_path is not None:
        load_stages(stages_path)

    # Make sure requested compounds actually exist
    missing = [
        compound
        for compound in compounds
        if compound not in data.channels
    ]

    if missing:
        available = "\n".join(sorted(data.channels))

        raise KeyError(
            "The following requested compounds were not found:\n"
            + "\n".join(f"  - {compound}" for compound in missing)
            + "\n\nAvailable channels:\n"
            + available
        )

    if concentration_unit == "molecules/cm3":
        if temperature_path is not None:
            temperature_c = load_data(temperature_path, channels=[temperature_channel])[temperature_channel]
        data = convert_dataset_concentrations(
            data, {name: source_unit for name in compounds}, concentration_unit,
            temperature_c=temperature_c, pressure_pa=pressure_pa,
        )
    if ylabel is None:
        native_label = {"ppt": "pptv", "ppb": "ppbv"}.get(source_unit, source_unit)
        ylabel = "Concentration [" + concentration_label(concentration_unit, native_label) + "]"
    yscale = concentration_yscale(concentration_unit, yscale)

    fan_speed_data = (
        load_data(fan_speed_path, time_range, channels=["Fan_BOT_Speed"])
        if fan_speed_path is not None else None
    )

    _plot_trace_data(
        data,
        channels=compounds,
        stacked=stacked,
        xlims=time_range,
        ylims=ylims,
        yscale=yscale,
        ylabels=ylabel,
        smoothing=smoothing,
        stages=stages,
        fan_speed_data=fan_speed_data,
        important_time_ranges=important_time_ranges,
        label_intervals=label_intervals,
        savepath=savepath,
        interactive=interactive,
    )


def _plot_instruments(sources, requests, time_range, preferences, source_units,
                      target_unit, temperature_path, temperature_channel, temperature_c,
                      pressure_pa, stacked, ylabel, stages_path, fan_speed_path, options):
    if not isinstance(sources, dict):
        raise TypeError("Grouped compounds require trace_path=run['sources']")
    selected = []
    for instrument, compounds in requests.items():
        if not compounds:
            continue
        if isinstance(compounds, str):
            raise TypeError(f"{instrument}: use a list of compound names")
        if instrument not in sources:
            raise KeyError(f"Instrument {instrument!r} is not available in this prepared run. "
                           "Check the run ID and rebuild with prepare_runs.py if needed.")
        source = sources[instrument]
        available = available_channels(source)
        names = []
        for request in compounds:
            try:
                name = resolve_channel(request, available, preferences.get(instrument, ()))
            except (KeyError, ValueError) as exc:
                raise type(exc)(f"{instrument}: {exc}") from exc
            if name not in names:
                names.append(name)
        info = prepared_channel_info(source) if isinstance(source, PreparedSource) else {}
        units = {name: info[name]["unit"] if info else source_units.get(instrument, "unknown")
                 for name in names}
        units = {name: {"ppt": "pptv", "ppb": "ppbv"}.get(unit, unit)
                 for name, unit in units.items()}
        if target_unit != "native" and any(unit not in {"pptv", "ppbv", target_unit} for unit in units.values()):
            raise ValueError(f"{instrument}: cannot convert {set(units.values())} to {target_unit}; "
                             "use concentration_unit='native' until the signals are calibrated.")
        selected.append((instrument, source, names, units))
    if not selected:
        raise ValueError("Select at least one compound")
    if target_unit != "native" and temperature_path is not None:
        temperature_c = load_data(temperature_path, channels=[temperature_channel])[temperature_channel]
    groups = {}
    for instrument, source, names, units in selected:
        data = load_data(source, time_range, channels=names)
        if target_unit != "native":
            data = convert_dataset_concentrations(data,
                                                 {name: unit for name, unit in units.items()
                                                  if unit != target_unit}, target_unit,
                                                 temperature_c=temperature_c, pressure_pa=pressure_pa)
        for name in names:
            label = f"{instrument.upper()} · {name}"
            print(f"Selected {label}")
            unit = units[name] if target_unit == "native" else target_unit
            group = (unit, label) if stacked else (unit,)
            ch = data[name]
            groups.setdefault(group, {})[label] = Channel(label, ch.time.copy(), ch.values.copy())
    if stages_path is not None:
        load_stages(stages_path)
    options["fan_speed_data"] = (load_data(fan_speed_path, time_range, channels=["Fan_BOT_Speed"])
                                  if fan_speed_path is not None else None)
    if not stacked and len(groups) == 2:
        # Keep concentration left and count rate right regardless of input order.
        keys = sorted(groups, key=lambda key: key[0] == "cps")
        left, right = [list(groups[key]) for key in keys]
        combined = Dataset({name: channel for key in keys
                            for name, channel in groups[key].items()})
        axis_labels = [ylabel or ("Signal" if key[0] in {"cps", "unknown"} else "Concentration")
                       + " [" + concentration_label(target_unit, key[0]) + "]" for key in keys]
        # The existing multi-dataset plotter supports a dual-axis channel tuple.
        return _plot_trace_data([combined], channels=[(left, right)],
                         ylabels=[tuple(axis_labels)], **options)
    panels = [Dataset(channels) for channels in groups.values()]
    labels = [ylabel or ("Signal" if key[0] in {"cps", "unknown"} else "Concentration")
              + " [" + concentration_label(target_unit, key[0]) + "]" for key in groups]
    if len(panels) == 1:
        return _plot_trace_data(panels[0], channels=list(panels[0].channels), ylabels=labels[0], **options)
    # plot_data's list interface keeps independent time axes and all stage/fan annotations.
    if options["ylims"] is not None:
        options["ylims"] = [options["ylims"]] * len(panels)
    return _plot_trace_data(panels, channels=[list(panel.channels) for panel in panels],
                     ylabels=labels, **options)
