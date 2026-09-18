"""Time-trace and interval-mean analysis for temperature-dependence work."""

from datetime import datetime, timezone
from typing import List, Optional, Tuple, Sequence

import matplotlib.pyplot as plt
import pandas as pd
import re

from .loading import read_trace_timeseries
from . import config as cfg


def _to_plot_datetime(value):
    """Convert a datetime/Timestamp to timezone-naive UTC for matplotlib."""
    if isinstance(value, pd.Timestamp):
        if value.tzinfo is not None:
            return value.tz_convert("UTC").tz_localize(None)
        return value

    if isinstance(value, datetime):
        if value.tzinfo is not None:
            return value.astimezone(timezone.utc).replace(tzinfo=None)
        return value

    return value


def _normalize_time_ranges(
    important_time_ranges: Sequence[Tuple[datetime, datetime]],
) -> list[Tuple[datetime, datetime]]:
    """Validate and convert interval boundaries to timezone-naive UTC."""
    normalized = []

    for start, stop in important_time_ranges:
        if not isinstance(start, datetime) or not isinstance(stop, datetime):
            raise TypeError(
                "important_time_ranges must contain (datetime, datetime) pairs"
            )

        start = _to_plot_datetime(start)
        stop = _to_plot_datetime(stop)

        if stop < start:
            raise ValueError(
                "An important_time_ranges interval has end < start"
            )

        normalized.append((start, stop))

    return normalized


def calculate_interval_means(
    data: pd.DataFrame,
    compounds: List[str],
    time_col: str,
    important_time_ranges: Sequence[Tuple[datetime, datetime]],
) -> pd.DataFrame:
    """
    Calculate the mean intensity of every requested compound in every
    important time interval.

    Returns a DataFrame with one row per interval and columns:
        range_index, start, end, n_points, <compound 1>, <compound 2>, ...
    """
    if not important_time_ranges:
        raise ValueError("important_time_ranges must not be empty")

    # Use a datetime copy for interval selection. The original data is not changed.
    time_values = pd.to_datetime(data[time_col], errors="coerce", utc=True)

    rows = []

    for i, (start, stop) in enumerate(important_time_ranges, start=1):
        start_utc = pd.Timestamp(start)
        stop_utc = pd.Timestamp(stop)

        if start_utc.tzinfo is None:
            start_utc = start_utc.tz_localize("UTC")
        else:
            start_utc = start_utc.tz_convert("UTC")

        if stop_utc.tzinfo is None:
            stop_utc = stop_utc.tz_localize("UTC")
        else:
            stop_utc = stop_utc.tz_convert("UTC")

        mask = time_values.between(start_utc, stop_utc)
        selected = data.loc[mask, compounds]

        row = {
            "range_index": i,
            "start": start,
            "end": stop,
            "n_points": len(selected),
        }

        for compound in compounds:
            values = pd.to_numeric(selected[compound], errors="coerce")
            row[compound] = values.mean()

        if len(selected) == 0:
            print(
                f"Warning: important time range {i} contains no trace points: "
                f"{start} – {stop}"
            )

        rows.append(row)

    return pd.DataFrame(rows)



def _chemical_label(formula: str) -> str:
    """Convert a trace-column formula to a compact matplotlib chemical label."""
    formula = str(formula).strip()
    parts = formula.split(".", 1)

    def format_part(part: str) -> str:
        # Element counts -> subscripts, terminal charge -> superscript.
        charge = ""
        if part.endswith("+") or part.endswith("-"):
            charge = part[-1]
            part = part[:-1]

        formatted = re.sub(r"([A-Z][a-z]?)(\d+)", r"\1$_{\2}$", part)
        if charge:
            formatted += rf"$^{{{charge}}}$"
        return formatted

    label = format_part(parts[0])
    if len(parts) == 2:
        label += r"$\cdot$" + format_part(parts[1])
    return label


def _apply_nature_axes(ax: plt.Axes) -> None:
    """Apply restrained publication-style axis formatting."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(direction="out")


def plot_kinetic_relationship(
    means: pd.DataFrame,
    x_compound: str,
    y_compound: str,
    show: bool = True,
    savepath: Optional[str] = None,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Plot interval-mean intensity of ``y_compound`` versus ``x_compound``.

    Each important time range contributes one point.
    """
    if x_compound not in means.columns:
        raise ValueError(f"{x_compound!r} is not present in the interval means")

    if y_compound not in means.columns:
        raise ValueError(f"{y_compound!r} is not present in the interval means")

    valid = means[[x_compound, y_compound]].dropna().copy()

    if valid.empty:
        raise RuntimeError(
            "No valid interval means are available for the requested "
            "compound pair."
        )

    fig, ax = plt.subplots(figsize=getattr(cfg, "KINETIC_FIGSIZE", (3.5, 3.0)))

    ax.scatter(
        valid[x_compound],
        valid[y_compound],
        s=getattr(cfg, "KINETIC_MARKER_SIZE", 24),
        color=getattr(cfg, "KINETIC_MARKER_COLOR", "black"),
        edgecolors="none",
        zorder=3,
    )

    # Label each point with the corresponding important time range.
    for _, row in valid.join(
        means[["range_index", "start", "end"]]
    ).iterrows():
        ax.annotate(
            f"T{int(row['range_index'])}",
            (row[x_compound], row[y_compound]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=getattr(cfg, "PAPER_FONT_SIZE", 7),
        )

    ax.set_xlabel(_chemical_label(x_compound) + " (ppt)")
    ax.set_ylabel(_chemical_label(y_compound) + " (ppt)")
    _apply_nature_axes(ax)

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=getattr(cfg, "PAPER_SAVE_DPI", 600), bbox_inches="tight")

    if show:
        plt.show()

    return fig, ax


def temp_dep(
    trace_path: str,
    compounds: List[str],
    time_range: Tuple[datetime, datetime],
    adduct: str = ".NH3H+",
    show: bool = True,
    savepath: Optional[str] = None,
    overlay: bool = False,
    important_time_ranges: Optional[List[Tuple[datetime, datetime]]] = None,
    x_compound: Optional[str] = None,
    y_compound: Optional[str] = None,
    kinetic_savepath: Optional[str] = None,
) -> tuple[
    plt.Figure,
    list[plt.Axes],
    pd.DataFrame,
    pd.DataFrame,
    Optional[plt.Figure],
    Optional[plt.Axes],
]:
    """
    Plot selected ion time traces and calculate interval means.

    ``important_time_ranges`` define the intervals used for the mean
    calculation. Every requested compound gets one mean per interval.

    If ``x_compound`` and ``y_compound`` are both given, a second figure is
    created showing the interval-mean values of y_compound versus x_compound.

    Returns
    -------
    trace_fig, trace_axes, data, means, kinetic_fig, kinetic_ax
    """
    plt.rcParams.update(getattr(cfg, "PAPER_RC_PARAMS", {}))

    if not compounds:
        raise ValueError("compounds must contain at least one compound")

    if important_time_ranges is None:
        important_time_ranges = []

    # The trace columns contain the ionic/adducted species.
    if adduct and not adduct.startswith("."):
        adduct = "." + adduct

    trace_compounds = [
        c if "." in c else f"{c}{adduct}"
        for c in (str(compound).strip() for compound in compounds)
    ]

    data = read_trace_timeseries(
        path=trace_path,
        compounds=trace_compounds,
        time_range=time_range,
    )

    time_col = "Time"
    if time_col not in data.columns:
        time_candidates = [
            c for c in data.columns if c.lower() == "time"
        ]
        if not time_candidates:
            raise RuntimeError("Loaded trace data contains no Time column.")
        time_col = time_candidates[0]

    compounds = trace_compounds

    if data.empty:
        raise RuntimeError(
            "No trace data found inside the requested time range: "
            f"{time_range[0].isoformat()} to {time_range[1].isoformat()}"
        )

    # Matplotlib gets a timezone-naive UTC copy; loaded data stays untouched.
    plot_time = data[time_col]
    if isinstance(plot_time.dtype, pd.DatetimeTZDtype):
        plot_time = plot_time.dt.tz_convert("UTC").dt.tz_localize(None)

    plot_important_ranges = _normalize_time_ranges(important_time_ranges)

    def _shade_important_ranges(ax):
        for i, (start, stop) in enumerate(plot_important_ranges, start=1):
            ax.axvspan(
                start,
                stop,
                color=getattr(cfg, "IMPORTANT_RANGE_COLOR", "#8E7CC3"),
                alpha=getattr(cfg, "IMPORTANT_RANGE_ALPHA", 0.12),
                zorder=0,
            )

            midpoint = start + (stop - start) / 2
            ax.text(
                midpoint,
                0.97,
                f"T{i}",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=getattr(cfg, "PAPER_FONT_SIZE", 7),
                fontweight="bold",
                color="black",
            )

    if overlay:
        trace_fig, ax = plt.subplots(figsize=getattr(cfg, "TRACE_FIGSIZE", (7.2, 3.0)))
        trace_axes = [ax]

        _shade_important_ranges(ax)

        for compound in compounds:
            ax.plot(
                plot_time,
                data[compound],
                label=_chemical_label(compound),
                linewidth=getattr(cfg, "PAPER_LINE_WIDTH", 1.0),
            )

        ax.set_ylabel("Signal (ppt)")
        ax.set_xlabel("Time (UTC)")
        ax.legend(
            frameon=False,
            handlelength=1.5,
            borderaxespad=0.2,
        )
        _apply_nature_axes(ax)

    else:
        trace_fig, axes_array = plt.subplots(
            nrows=len(compounds),
            ncols=1,
            sharex=True,
            squeeze=False,
            figsize=(getattr(cfg, "TRACE_FIGSIZE", (7.2, 3.0))[0], max(1.8 * len(compounds), 2.5)),
        )
        trace_axes = list(axes_array[:, 0])

        for ax, compound in zip(trace_axes, compounds):
            _shade_important_ranges(ax)

            ax.plot(
                plot_time,
                data[compound],
                linewidth=getattr(cfg, "PAPER_LINE_WIDTH", 1.0),
            )
            ax.set_ylabel("Signal (ppt)")
            ax.text(
                0.01, 0.95,
                _chemical_label(compound),
                transform=ax.transAxes,
                ha="left", va="top",
            )
            _apply_nature_axes(ax)

        trace_axes[-1].set_xlabel("Time (UTC)")

    trace_fig.tight_layout()

    if savepath:
        trace_fig.savefig(savepath, dpi=getattr(cfg, "PAPER_SAVE_DPI", 600), bbox_inches="tight")

    # Calculate interval means independently of the trace plot.
    if plot_important_ranges:
        means = calculate_interval_means(
            data=data,
            compounds=compounds,
            time_col=time_col,
            important_time_ranges=plot_important_ranges,
        )
    else:
        means = pd.DataFrame(
            columns=["range_index", "start", "end", "n_points"] + compounds
        )

    kinetic_fig = None
    kinetic_ax = None

    if (x_compound is None) != (y_compound is None):
        raise ValueError(
            "Provide both x_compound and y_compound, or neither."
        )

    if x_compound is not None and y_compound is not None:
        # Allow users to refer to the base formula when the default adduct
        # is used.
        if "." not in x_compound:
            x_compound = f"{x_compound}{adduct}"
        if "." not in y_compound:
            y_compound = f"{y_compound}{adduct}"

        kinetic_fig, kinetic_ax = plot_kinetic_relationship(
            means=means,
            x_compound=x_compound,
            y_compound=y_compound,
            show=show,
            savepath=kinetic_savepath,
        )

    if show and kinetic_fig is None:
        plt.show()

    return (
        trace_fig,
        trace_axes,
        data,
        means,
        kinetic_fig,
        kinetic_ax,
    )
