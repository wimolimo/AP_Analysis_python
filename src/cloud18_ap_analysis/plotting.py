from __future__ import annotations

import os
import sys
from pathlib import Path
from datetime import datetime, timedelta
import calendar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap, LogNorm#
import matplotlib.patches as mpatches
import itertools

from .models import Channel, Dataset, SMPSData
from .loading import load_data
from .stages import get_stages
from .config import DPI, PAPER_LINE_COLORS, PAPER_UNSTACKED_SIZE, STAGE_SETTINGS

import pandas as pd
import warnings


def set_plot_backend(backend: str):
    """Select a Matplotlib backend before plotting.

    ``interactive`` plotting normally works with the GUI backend selected by
    Matplotlib. Call this only if you need to force a particular backend.
    """
    backend = backend.lower()
    if backend in {"interactive", "gui"}:
        plt.ion()
    elif backend in {"noninteractive", "static"}:
        plt.ioff()
    else:
        plt.switch_backend(backend)


def _channel_datetimes(ch):
    return [
        datetime.utcfromtimestamp(t / 1000.0)
        if np.isfinite(t)
        else datetime(1970, 1, 1)
        for t in ch.time
    ]


def _utcify_datetime(dt):
    return dt


def _utcify_range(xlims):
    if xlims is None:
        return None
    return tuple(_utcify_datetime(dt) for dt in xlims)


def _datetime_to_unix_ms(dt):
    return (calendar.timegm(dt.timetuple()) + dt.microsecond / 1_000_000.0) * 1000.0


def _present(fig, interactive=False, savepath=None):
    if interactive:
        plt.show(block=True)
    else:
        filepath = Path(savepath or "display.png")
        filepath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(filepath, dpi=DPI, bbox_inches="tight")
        print(f"Figure saved to {filepath}")
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(filepath))
            elif sys.platform == "darwin":
                os.system(f'open "{filepath}"')
            else:
                os.system(f'xdg-open "{filepath}" >/dev/null 2>&1 &')
        except Exception:
            pass
    return None


def _get_next_colors(n):
    """Gets the next `n` colors from the global color generator."""
    if n <= 0:
        return []
    # Take the next n colors from the generator
    return [next(_COLOR_GENERATOR) for _ in range(n)]

def _reset_color_cycle():
    """Resets the global color generator to the beginning of the list."""
    global _COLOR_GENERATOR
    _COLOR_GENERATOR = itertools.cycle(PAPER_LINE_COLORS)


def _append_to_filename(path, tag):
    if not tag:
        return path
    p = Path(path)
    return str(p.with_name(f"{p.stem}_{tag}{p.suffix}"))


def _smoothing_label(window):
    if window is None:
        return ""
    seconds = _window_seconds(window)
    ms = round(seconds * 1000)
    if ms >= 3600000 and ms % 3600000 == 0:
        return f"avg{ms // 3600000}h"
    if ms >= 60000 and ms % 60000 == 0:
        return f"avg{ms // 60000}min"
    if ms >= 1000 and ms % 1000 == 0:
        return f"avg{ms // 1000}s"
    return f"avg{ms}ms"


def _window_seconds(window):
    if isinstance(window, timedelta):
        return window.total_seconds()
    if isinstance(window, (int, float)):
        # Numeric smoothing values are interpreted as seconds.
        return float(window)
    raise TypeError("smoothing must be a datetime.timedelta, seconds, or None")


def _smooth_channel(ch, window):
    half_w = _window_seconds(window) * 1000.0 / 2.0
    out = np.full(len(ch.time), np.nan)

    valid = np.isfinite(ch.values) & np.isfinite(ch.time)
    vi = np.flatnonzero(valid)
    if len(vi) == 0:
        return Channel(ch.name, ch.time.copy(), out)

    all_idx = np.flatnonzero(np.isfinite(ch.time))
    all_idx = all_idx[np.argsort(ch.time[all_idx])]
    valid_idx = vi[np.argsort(ch.time[vi])]
    times_valid = ch.time[valid_idx]

    left = 0
    right = -1
    running_sum = 0.0
    running_count = 0

    for i in all_idx:
        t = ch.time[i]
        hi = t + half_w
        lo = t - half_w

        while right + 1 < len(valid_idx) and times_valid[right + 1] <= hi:
            right += 1
            running_sum += ch.values[valid_idx[right]]
            running_count += 1

        while left <= right and times_valid[left] < lo:
            running_sum -= ch.values[valid_idx[left]]
            running_count -= 1
            left += 1

        if running_count:
            out[i] = running_sum / running_count

    return Channel(ch.name, ch.time.copy(), out)


def _maybe_smooth(channels, smoothing):
    if smoothing is None:
        return channels
    return [_smooth_channel(ch, smoothing) for ch in channels]


def _resolve_channels(data, requested):
    if requested is None or requested == []:
        return [data.channels[k] for k in sorted(data.channels)]
    return [x if isinstance(x, Channel) else data.channels[str(x)] for x in requested]


def _padded_ylims(ymin, ymax):
    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return None
    if ymin == ymax:
        if ymin == 0:
            return (-1, 1)
        margin = abs(ymin) * 0.05
    else:
        margin = (ymax - ymin) * 0.05
    return ymin - margin, ymax + margin


def _compute_ylims(channels, xlims=None):
    xlims = _utcify_range(xlims)
    ys = []
    for ch in channels:
        mask = np.isfinite(ch.values)
        if xlims is not None:
            lo = _datetime_to_unix_ms(xlims[0])
            hi = _datetime_to_unix_ms(xlims[1])
            mask &= (ch.time >= lo) & (ch.time <= hi)
        ys.extend(ch.values[mask])
    if not ys:
        return None
    return _padded_ylims(float(np.min(ys)), float(np.max(ys)))


def _apply_ylims(ax, yl, channels, xlims):
    if yl is not None:
        ax.set_ylim(*yl)
    else:
        computed = _compute_ylims(channels, xlims)
        if computed:
            ax.set_ylim(*computed)


def _get_stage_label(stage, key):
    if key is None or key == "stage":
        return stage["stage"]
    if key == "type":
        return f'{stage["stage"]}: {stage["type"]}'
    if key == "onlytype":
        return stage["type"]
    if key == "description":
        return stage["description"] or stage["stage"]
    if key == "comments":
        return stage["comments"] or stage["stage"]
    if key == "details":
        d = "\n".join(stage["details"])
        return d or stage["stage"]
    if key == "full":
        return " - ".join(x for x in [stage["stage"], stage["type"], stage["description"]] if x)
    return stage["stage"]


def _print_stage_details(stage):
    lines = [
        f"Stage {stage['stage']}",
        f"Time: {stage['time'].isoformat(sep=' ') if stage['time'] else 'n/a'}",
        f"Type: {stage['type'] or 'n/a'}",
        f"Description: {stage['description'] or 'n/a'}",
        f"Comments: {stage['comments'] or 'n/a'}",
    ]
    if stage["details"]:
        lines.append("Details:")
        lines.extend(f"  {detail}" for detail in stage["details"])
    print("\n" + "\n".join(lines))


def _draw_stages(ax, stages, xlims=None, print_details=False):
    data = get_stages()
    if not data:
        print("Warning: no stages loaded. Call load_stages(filepath) first.")
        return
    xlims = _utcify_range(xlims)
    label_key = None if stages is True else stages
    x_span = (xlims[1] - xlims[0]).total_seconds() if xlims else 60
    dt_offset = timedelta(seconds=x_span * STAGE_SETTINGS["text_offset"])

    ymin, ymax = ax.get_ylim()
    for stage in data:
        t = stage["time"]
        if t is None:
            continue
        if xlims and not (xlims[0] <= t <= xlims[1]):
            continue
        if print_details:
            _print_stage_details(stage)
        label = _get_stage_label(stage, label_key)
        if len(label) > STAGE_SETTINGS["max_length"]:
            label = label[:STAGE_SETTINGS["max_length"]] + "…"
        ax.axvline(
            t, color=STAGE_SETTINGS["line_color"],
            linewidth=STAGE_SETTINGS["line_width"],
            linestyle=STAGE_SETTINGS["line_style"],
        )
        ax.text(
            t + dt_offset,
            ymin + 0.02 * (ymax - ymin),
            label,
            rotation=90,
            fontsize=STAGE_SETTINGS["font_size"],
            color=STAGE_SETTINGS["font_color"],
            ha="left",
            va="bottom",
        )

def _draw_fan_speed_shading(ax, data, fan_speed_col="Fan_BOT_Speed", threshold=50):
    """
    Draws shaded gray regions on an axis where fan speed is above a threshold.
    This version uses pure NumPy to find the contiguous blocks.

    Args:
        ax: The matplotlib Axes object to draw on.
        data: The Dataset object containing the channels.
        fan_speed_col: The name of the fan speed Channel in the Dataset.
        threshold: The fan speed percentage to use as a threshold.
    """
    # 1. Safely get the fan speed Channel
    if fan_speed_col not in data.channels:
        warnings.warn(f"Fan speed channel '{fan_speed_col}' not found. Cannot draw shading.")
        return

    fan_channel = data.channels[fan_speed_col]
    if fan_channel.time.size < 2:  # Need at least 2 points to find blocks
        return

    # 2. Convert time array to datetimes for plotting with axvspan
    # This is the only part where Pandas is useful here.
    time_as_datetime = pd.to_datetime(fan_channel.time, unit='ms')
    n_times = len(time_as_datetime)
    
    # 3. Find contiguous blocks using pure NumPy
    is_low = fan_channel.values < threshold

    # Find the indices where the state changes
    # Pad with False at both ends to catch changes at the very start/end
    padded_is_low = np.concatenate(([False], is_low, [False]))
    # Use diff() to find changes: +1 means start of a low block, -1 means end
    diffs = np.diff(padded_is_low.astype(np.int8))
    
    start_indices = np.where(diffs == 1)[0]
    end_indices = np.where(diffs == -1)[0]

    # 4. Draw a shaded region for each start/end pair
    for start_idx, end_idx in zip(start_indices, end_indices):
        # Get the corresponding start and end times
        start_time = time_as_datetime[start_idx]
        # The end_idx from diff corresponds to the start of the "low" period,
        # so it's the correct exclusive end time for our high period span.
        if end_idx < n_times:
            end_time = time_as_datetime[end_idx]
        else:
            end_time = time_as_datetime[-1]

        ax.axvspan(
            start_time,
            end_time,
            color='gray',
            alpha=0.2,
            zorder=0,  # Draw behind the plot lines
            linewidth=0
        )


def _style_axis(ax):
    ax.grid(True, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.75, length=3)


def _plot_on_axis(ax, channels, colors, offset=0):
    entries = []
    for j, ch in enumerate(channels):
        line, = ax.plot(_channel_datetimes(ch), ch.values,
                        label=ch.name, color=colors[offset + j])
        entries.append((line, ch.name))
    return entries


def _parse_panel_ylims(yl, is_dual):
    if yl is None:
        return None, None
    if not is_dual:
        if isinstance(yl, tuple) and len(yl) == 2 and not all(isinstance(x, (int, float)) for x in yl):
            return yl[0], None
        return yl, None
    if isinstance(yl, tuple) and len(yl) == 2:
        if all(isinstance(x, (int, float)) for x in yl):
            return yl, None
        return yl[0], yl[1]
    return yl, None


def _prepare_savepath(savepath, smoothing):
    tag = _smoothing_label(smoothing)
    return (_append_to_filename(savepath, tag) if savepath and tag else savepath), tag


def _configure_figure(fig):
    fig.tight_layout()

def _is_smps_data(data):
    return isinstance(data, SMPSData) or all(
        hasattr(data, attr) for attr in ("time", "diameters", "dNdlogDp")
    )

def _smps_time_to_mpl(time_values):
    """Convert SMPS datetime/object time array to Matplotlib date numbers."""
    out = []
    valid = []

    for t in time_values:
        if isinstance(t, datetime):
            out.append(mdates.date2num(t))
            valid.append(True)
        elif isinstance(t, np.datetime64):
            if np.isnat(t):
                out.append(np.nan)
                valid.append(False)
            else:
                py_dt = t.astype("datetime64[ms]").astype(datetime)
                out.append(mdates.date2num(py_dt))
                valid.append(True)
        else:
            out.append(np.nan)
            valid.append(False)

    return np.asarray(out, dtype=float), np.asarray(valid, dtype=bool)

def _prepare_smps_for_plot(smps):
    """Prepare SMPS data for pcolormesh.

    SMPSData.dNdlogDp is expected as time × diameter bins.
    Matplotlib wants C as diameter × time, so this transposes it.
    """
    time_num, valid_time = _smps_time_to_mpl(smps.time)

    diameters = np.asarray(smps.diameters, dtype=float)
    conc = np.asarray(smps.dNdlogDp, dtype=float)

    if conc.ndim != 2:
        raise ValueError("SMPS dNdlogDp must be a 2D array")

    n_time = len(time_num)
    n_diam = len(diameters)

    if conc.shape == (n_time, n_diam):
        conc = conc.T
    elif conc.shape == (n_diam, n_time):
        pass
    else:
        raise ValueError(
            "SMPS dNdlogDp shape does not match time/diameters. "
            f"Got {conc.shape}, expected {(n_time, n_diam)} or {(n_diam, n_time)}."
        )

    valid_diam = np.isfinite(diameters) & (diameters > 0)

    time_num = time_num[valid_time]
    diameters = diameters[valid_diam]
    conc = conc[np.ix_(valid_diam, valid_time)]

    conc_masked = np.ma.masked_where(
        (~np.isfinite(conc)) | (conc < 0),
        conc,
    )

    positive = conc[np.isfinite(conc) & (conc > 0)]

    if positive.size:
        q01, q99 = np.nanquantile(positive, [0.01, 0.99])
        vmin = max(1.0, float(q01))
        vmax = float(q99)
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin * 10.0
    else:
        vmin, vmax = 1.0, 1e4

    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)

    return time_num, diameters, conc_masked, norm

def _plot_smps_on_axis(fig, ax, smps, ylabel=None, cbar_label=None,
                       xlims=None, ylims=None):
    time_num, diameters, conc, norm = _prepare_smps_for_plot(smps)

    ax.set_ylabel(ylabel or "Diameter [nm]")
    ax.set_yscale("log")

    original_cmap = plt.get_cmap("inferno")
    new_colors = original_cmap(np.linspace(0, 1, 256))

    # Set the first color in the list (for the bottom of the scale) to black
    new_colors[0] = (0, 0, 0, 1) # RGBA for black
    cmap = ListedColormap(new_colors)
    cmap.set_bad("white")

    if len(time_num) == 0 or len(diameters) == 0 or conc.size == 0:
        ax.text(
            0.5, 0.5, "No SMPS data",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
        return None
    mesh = ax.pcolormesh(
        time_num,
        diameters,
        conc,
        shading="auto",
        cmap=cmap,
        norm=norm,
    )

    cax = ax.inset_axes([1.012, 0.0, 0.025, 1.0], transform=ax.transAxes)
    cax.set_in_layout(False)

    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(cbar_label or "dN/dlogDp [cm⁻³]")

    ax.xaxis_date()

    if xlims:
        ax.set_xlim(*xlims)

    if ylims is not None:
        ax.set_ylim(*ylims)

    return mesh

def _plot_channels(channels, title="Overview", savepath=None, stacked=False,
                   interactive=False, xlims=None, ylims=None, smoothing=None,
                   stages=None):
    xlims = _utcify_range(xlims)
    channels = _maybe_smooth(channels, smoothing)
    savepath, tag = _prepare_savepath(savepath, smoothing)
    if tag:
        title = f"{title} ({tag})"

    colors = _get_next_colors(len(channels))

    if stacked:
        fig, axes = plt.subplots(
            len(channels), 1,
            figsize=(PAPER_UNSTACKED_SIZE[0], 2 * max(1, len(channels))),
            sharex=True,
            squeeze=False,
        )
        axes = axes[:, 0]
        for i, (ax, ch) in enumerate(zip(axes, channels)):
            ax.plot(_channel_datetimes(ch), ch.values, color=colors[i])
            if color_by_fan_speed and not _is_smps_data(data):
                draw_fan_speed_shading(ax, ch.values, xlims=xlims)
            ax.set_ylabel(ch.name)
            ax.set_title(title if i == 0 else "")
            _style_axis(ax)
            if i == len(channels) - 1:
                ax.set_xlabel("Time (UTC)")
            if xlims:
                ax.set_xlim(*xlims)
            yl = ylims[i] if isinstance(ylims, (list, tuple)) and not (len(ylims) == 2 and all(isinstance(x, (int, float)) for x in ylims)) else ylims
            _apply_ylims(ax, yl, [ch], xlims)
            if stages is not None:
                _draw_stages(ax, stages, xlims)
    else:
        fig, ax = plt.subplots(figsize=PAPER_UNSTACKED_SIZE)
        _plot_on_axis(ax, channels, colors)
        ax.set_title(title)
        ax.set_xlabel("Time (UTC)")
        _style_axis(ax)
        ax.legend(loc="upper right")
        if xlims:
            ax.set_xlim(*xlims)
        _apply_ylims(ax, ylims, channels, xlims)
        if stages is not None:
            _draw_stages(ax, stages, xlims)

    _configure_save_and_present(fig, interactive, savepath)
    return None


def _configure_save_and_present(fig, interactive, savepath=None):
    if savepath is not None:
        path = Path(savepath)
        path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"Figure saved to {path}")

    if interactive:
        plt.show(block=True)
    else:
        plt.close(fig)


def plot_channel(ch, ylabel=None, title=None, savepath=None, interactive=False,
                 xlims=None, ylims=None, smoothing=None):
    xlims = _utcify_range(xlims)
    ch = _maybe_smooth([ch], smoothing)[0]
    savepath, tag = _prepare_savepath(savepath, smoothing)
    title = title or ch.name
    ylabel = ylabel or ch.name
    if tag:
        title = f"{title} ({tag})"

    fig, ax = plt.subplots(figsize=PAPER_UNSTACKED_SIZE)
    ax.plot(_channel_datetimes(ch), ch.values)
    ax.set_title(title)
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel(ylabel)
    _style_axis(ax)
    if xlims:
        ax.set_xlim(*xlims)
    _apply_ylims(ax, ylims, [ch], xlims)
    _configure_save_and_present(fig, interactive, savepath)


def plot_data(data_or_path, channels=None, stacked=False, savepath=None,
              interactive=False, xlims=None, ylims=None, ylabels=None,
              smoothing=None, stages=None, fan_speed_data=None):
    _reset_color_cycle()
    if isinstance(data_or_path, (list, tuple)):
        return plot_datasets(
            data_or_path, channels=channels, savepath=savepath,
            interactive=interactive, xlims=xlims, ylims=ylims,
            ylabels=ylabels, smoothing=smoothing, stages=stages,
            fan_speed_data=fan_speed_data,
        )
    data = load_data(data_or_path) if isinstance(data_or_path, (str, Path)) else data_or_path

    if _is_smps_data(data):
        return plot_datasets(
            [data],
            channels=[None],
            savepath=savepath,
            interactive=interactive,
            xlims=xlims,
            ylims=[ylims] if ylims is not None else None,
            ylabels=[ylabels] if ylabels is not None else None,
            smoothing=smoothing,
            stages=stages,
        )

    selected = _resolve_channels(data, channels)
    title = Path(data_or_path).name if isinstance(data_or_path, (str, Path)) else "Data"
    return _plot_channels(
        selected, title=title, stacked=stacked, savepath=savepath,
        interactive=interactive, xlims=xlims, ylims=ylims,
        smoothing=smoothing, stages=stages, fan_speed_data=fan_speed_data,
    )


def plot_channels(channels, title="Overview", savepath=None, stacked=False,
                  interactive=False, xlims=None, ylims=None, fan_speed_data=None):
    return _plot_channels(
        channels, title=title, savepath=savepath, stacked=stacked,
        interactive=interactive, xlims=xlims, ylims=ylims, fan_speed_data=fan_speed_data
    )


def plot_all(filepath, stacked=False, savepath=None, interactive=False,
             xlims=None, ylims=None):
    return plot_data(filepath, stacked=stacked, savepath=savepath,
                     interactive=interactive, xlims=xlims, ylims=ylims)


def plot_datasets(datasets, channels=None, savepath=None, interactive=False,
                  xlims=None, ylims=None, ylabels=None, smoothing=None,
                  stages=None, fan_speed_data=None):
    """Plot multiple Dataset and/or SMPSData objects in stacked panels.

    Normal Dataset panels become line plots.

    SMPSData panels become particle-size heatmaps using:
        SMPSData.time
        SMPSData.diameters
        SMPSData.dNdlogDp

    For SMPS labels, use either:

        ylabels[i] = "Diameter [nm]"

    or:

        ylabels[i] = ("Diameter [nm]", "dN/dlogDp [cm⁻³]")
    """
    xlims = _utcify_range(xlims)
    resolved = [load_data(d) if isinstance(d, (str, Path)) else d for d in datasets]

    n = len(resolved)
    if n == 0:
        return None

    if channels is None:
        specs = [None] * n
    elif isinstance(channels, (list, tuple)) and len(channels) == n and all(
        c is None or isinstance(c, (list, tuple)) for c in channels
    ):
        specs = list(channels)
    else:
        specs = [channels] * n

    fig, axes = plt.subplots(
        n, 1,
        figsize=(PAPER_UNSTACKED_SIZE[0], 2 * n),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes[:, 0]

    for i, (data, spec) in enumerate(zip(resolved, specs)):
        ax = axes[i]
        _style_axis(ax)

        if i == n - 1:
            ax.set_xlabel("Time (UTC)")

        if fan_speed_data is not None and not _is_smps_data(data):
            _draw_fan_speed_shading(ax, fan_speed_data)

        panel_ylims = ylims[i] if isinstance(ylims, list) and i < len(ylims) else ylims

        ylabel = ""
        right_ylabel = ""

        if ylabels is not None and i < len(ylabels):
            if isinstance(ylabels[i], tuple):
                ylabel = ylabels[i][0] if len(ylabels[i]) >= 1 else ""
                right_ylabel = ylabels[i][1] if len(ylabels[i]) >= 2 else ""
            else:
                ylabel = ylabels[i]

        # ── SMPS heatmap panel ─────────────────────────────────────
        if _is_smps_data(data):
            _plot_smps_on_axis(
                fig,
                ax,
                data,
                ylabel=ylabel or "Diameter [nm]",
                cbar_label=right_ylabel or "dN/dlogDp [cm⁻³]",
                xlims=xlims,
                ylims=panel_ylims,
            )

            if stages is not None:
                _draw_stages(ax, stages, xlims, print_details=(i == n - 1))

            continue

        # ── Normal Dataset line panel ──────────────────────────────
        if not isinstance(data, Dataset):
            raise NotImplementedError(
                f"Unsupported dataset type: {type(data).__name__}. "
                "Expected Dataset or SMPSData."
            )

        ax.set_ylabel(ylabel or "")

        dual = isinstance(spec, tuple) and len(spec) == 2
        ax2 = None

        if dual:
            left = _maybe_smooth(_resolve_channels(data, spec[0]), smoothing)
            right = _maybe_smooth(_resolve_channels(data, spec[1]), smoothing)

            if not right:
                colors = _get_next_colors(len(left))
                _plot_on_axis(ax, left, colors)
                if left:
                    ax.legend(loc="upper left")
            else:
                total = len(left) + len(right)
                colors = _get_next_colors(total)

                _plot_on_axis(ax, left, colors, 0)

                ax2 = ax.twinx()
                _style_axis(ax2)
                ax2.set_ylabel(right_ylabel or "")

                _plot_on_axis(ax2, right, colors, len(left))

                if left:
                    ax.legend(loc="upper left")
                if right:
                    ax2.legend(loc="upper right")

            if xlims:
                ax.set_xlim(*xlims)
                if ax2 is not None:
                    ax2.set_xlim(*xlims)

            left_yl, right_yl = _parse_panel_ylims(panel_ylims, is_dual=True)

            _apply_ylims(ax, left_yl, left, xlims)

            if ax2 is not None:
                _apply_ylims(ax2, right_yl, right, xlims)

        else:
            selected = _maybe_smooth(_resolve_channels(data, spec), smoothing)

            colors = _get_next_colors(len(selected))
            _plot_on_axis(ax, selected, colors)

            if selected:
                ax.legend(loc="upper left")

            if xlims:
                ax.set_xlim(*xlims)

            _apply_ylims(ax, panel_ylims, selected, xlims)

        if stages is not None:
            _draw_stages(ax, stages, xlims, print_details=(i == n - 1))

    if fan_speed_data is not None:
        ax0 = axes[0] # Target the first subplot's axis

        # Check if a legend already exists on the first plot (it might be SMPS)
        if ax0.get_legend() is not None:
            # Get the handles and labels that are already in the legend
            current_handles, current_labels = ax0.get_legend_handles_labels()

            # Create our proxy artist for the gray area.
            # Use the same alpha as in _draw_fan_speed_shading for consistency.
            fan_patch = mpatches.Patch(color='gray', alpha=0.2)
            fan_label = 'Fan Speed 12%'

            # Prepend the new entry to the existing lists
            new_handles = [fan_patch] + current_handles
            new_labels = [fan_label] + current_labels

            # Redraw the legend on the first subplot with the combined list
            ax0.legend(new_handles, new_labels, loc="upper left")

    if smoothing is not None:
        tag = _smoothing_label(smoothing)
        if savepath:
            savepath = _append_to_filename(savepath, tag)
        axes[0].set_title(tag)

    _configure_save_and_present(fig, interactive, savepath)
    return None