"""Plot selected compound time traces with CLOUD stage information."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

from .loading import load_data
from .stages import load_stages
from .plotting import plot_data


def plot_time_traces(
    trace_path: str | Path,
    compounds: list[str],
    time_range: tuple[datetime, datetime],
    stages_path: str | Path | None = None,
    savepath: str | Path | None = None,
    smoothing: timedelta | float | None = None,
    interactive: bool = True,
    stages: str | bool | None = "onlytype",
    stacked: bool = False,
    ylabel: str = "Signal [cps]",
    ylims=None,
    yscale: str = "log",
    fan_speed_path: str | Path | None = None,
):
    """
    Plot time traces of selected compounds.

    Parameters
    ----------
    trace_path
        Path to the trace CSV file.

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
        Label for the y-axis.

    ylims
        Optional y-axis limits.

    yscale
        Y-axis scale: "log" (default) or "linear". Log scales display
        only positive values.

    fan_speed_path
        Optional fan-speed CSV file (same format as run_plotting). Load
        Fan_BOT_Speed for time_range and shade periods below 50 percent gray.
        Applies to overlay and stacked traces. None disables shading.
    """

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
        savepath=savepath,
        interactive=interactive,
    )