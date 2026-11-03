"""Plot selected compound time traces with CLOUD stage information."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

from .loading import load_data
from .stages import load_stages
from .plotting import plot_data
from .units import (convert_dataset_concentrations, concentration_label,
                    concentration_yscale, validate_concentration_unit, STANDARD_PRESSURE_PA)


def plot_time_traces(
    trace_path: str | Path | list[str | Path],
    compounds: list[str],
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
        Compound/channel names that should be shown.

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

    # Load only the requested time period
    data = load_data(trace_path, time_range)

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
            temperature_c = load_data(temperature_path)[temperature_channel]
        data = convert_dataset_concentrations(
            data, {name: source_unit for name in compounds}, concentration_unit,
            temperature_c=temperature_c, pressure_pa=pressure_pa,
        )
    if ylabel is None:
        native_label = {"ppt": "pptv", "ppb": "ppbv"}.get(source_unit, source_unit)
        ylabel = "Concentration [" + concentration_label(concentration_unit, native_label) + "]"
    yscale = concentration_yscale(concentration_unit, yscale)

    fan_speed_data = (
        load_data(fan_speed_path, time_range)
        if fan_speed_path is not None else None
    )

    plot_data(
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