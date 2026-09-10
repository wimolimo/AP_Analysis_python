# src/mdplot.py
import re
from dataclasses import dataclass
from typing import Optional, Tuple, Iterable, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timezone

# Pull plotting settings from config.py (same folder); provide safe defaults
try:
    import config as cfg  # same folder as this file
except Exception:  # pragma: no cover
    class cfg:  # fallback defaults
        FIGSIZE = (9, 5.5)
        CMAP = "viridis"
        POINT_SIZE = 12.0
        ALPHA = 0.8
        LOG_COLOR = True
        MIN_INTENSITY = 0.0
        DPI = 150
        BACKGROUND = "white"


# -----------------------------
# Public API
# -----------------------------

@dataclass
class MDPlotResult:
    data: pd.DataFrame   # Data used for plotting: SumFormula, Mass, AggIntensity, NominalMass, MassDefect
    fig: plt.Figure
    ax: plt.Axes


def md_plot(
    composition_path: str,
    trace_path: str,
    time_range: Tuple[datetime, datetime],  # datetime only
    agg: str = "sum",                       # "sum" or "mean"
    clip_neg: bool = False,                 # clip negative intensities to 0 before aggregating
    # plotting knobs
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show: bool = True,                      # call plt.show() automatically
    savepath: Optional[str] = None,         # if set, save the figure here
    stream: bool = True,                    # streaming read for huge CSVs
    chunksize: int = 200_000,
    assume_sorted_unix: bool = True,
) -> MDPlotResult:
    """
    Build a classical mass defect (MD = m - round(m)) plot using datetime-only time ranges.
    - composition_path: composition txt (with 'number of header rows', columns include SumFormula, Mass)
    - trace_path: huge trace csv (with 'number of header rows', columns Time, unixTime, and SumFormula columns)
    - time_range: (datetime_start, datetime_end); naive datetimes assumed UTC
    - agg: 'sum' or 'mean' aggregation within the time window
    - clip_neg: clip negative intensities to 0 before aggregating
    - xlim, ylim: axis limits (tuples), e.g. xlim=(50, 1000), ylim=(-0.2, 0.2)
    - show: if True, call plt.show()
    - savepath: if provided, save the figure to this path (PNG, PDF, etc.)
    - stream: stream the CSV by time window and only needed columns
    """
    # Validate/normalize time range
    t0_dt, t1_dt = _ensure_datetimes(time_range)
    if t1_dt < t0_dt:
        raise ValueError("time_range end < start")

    # Load small composition table
    comp = read_composition_table(composition_path)  # ['SumFormula','Mass']
    comp["SumFormula"] = comp["SumFormula"].astype(str).str.strip()

    # Peek at trace header without loading data
    header = read_trace_header(trace_path)
    time_col = header["time_col"]
    unix_col = header["unix_col"]
    all_trace_formulas = header["intensity_cols"]

    # Intersect formulas present in both files
    common = sorted(set(comp["SumFormula"]).intersection(all_trace_formulas))
    if not common:
        raise RuntimeError(
            "No overlap between composition SumFormula and trace columns. "
            "Ensure names match exactly (e.g., 'CH2O.H+' in both files)."
        )

    # Aggregate over requested window
    if stream:
        agg_series = aggregate_window_streaming(
            trace_path=trace_path,
            t0_dt=t0_dt,
            t1_dt=t1_dt,
            formulas=common,
            how=agg,
            clip_neg=clip_neg,
            assume_sorted_unix=assume_sorted_unix,
            chunksize=chunksize,
            skiprows=_detect_and_skip_header_rows(trace_path),
            sep=_detect_delimiter(trace_path),
        )
    else:
        traces = read_trace_table(trace_path)
        agg_series = aggregate_window_unix_dt(traces, t0_dt, t1_dt, how=agg, clip_neg=clip_neg)
        agg_series = agg_series.reindex(common)

    agg_series = agg_series.dropna()

    # Build plotting table
    df_plot = comp[comp["SumFormula"].isin(common)].merge(
        agg_series.rename("AggIntensity"), left_on="SumFormula", right_index=True, how="inner"
    )
    df_plot = df_plot[df_plot["AggIntensity"] > getattr(cfg, "MIN_INTENSITY", 0.0)]
    if df_plot.empty:
        raise RuntimeError("No data to plot after time-window aggregation and filtering.")

    # Classical mass defect: MD = m - round(m)
    nominal, defect = classical_mass_defect(df_plot["Mass"].values)
    df_plot["NominalMass"] = nominal
    df_plot["MassDefect"] = defect

    # Colors
    if getattr(cfg, "LOG_COLOR", True):
        color = np.log10(df_plot["AggIntensity"].values + 1.0)
        cbar_label = "log10(intensity in window + 1)"
    else:
        color = df_plot["AggIntensity"].values
        cbar_label = "Intensity in window"

    # Sizes
    if getattr(cfg, "SIZE_BY_INTENSITY", True):
        sizes = _compute_sizes_from_intensity(df_plot["AggIntensity"].values)
    else:
        sizes = float(getattr(cfg, "POINT_SIZE", 12.0))

    # Plot
    fig, ax = plt.subplots(figsize=getattr(cfg, "FIGSIZE", (9, 5.5)))
    sc = ax.scatter(
        df_plot["NominalMass"].values,
        df_plot["MassDefect"].values,
        c=color,
        s=sizes,
        cmap=getattr(cfg, "CMAP", "viridis"),
        alpha=getattr(cfg, "ALPHA", 0.8),
        edgecolors="none",
    )
    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Mass defect")
    ax.set_title(f"Mass Defect Plot, UTC ∈ [{t0_dt.isoformat()} — {t1_dt.isoformat()}]")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(cbar_label)

    # Optional: size legend so readers can interpret marker sizes
    if getattr(cfg, "SIZE_BY_INTENSITY", True) and getattr(cfg, "SIZE_LEGEND", True):
        intens = df_plot["AggIntensity"].values
        # Determine legend levels
        levels = getattr(cfg, "SIZE_LEGEND_LEVELS", None)
        if levels is None:
            # Use representative percentiles
            qs = [10, 50, 90]
            levels = [np.nanpercentile(intens, q) for q in qs]
        # Build legend handles
        handles = []
        labels = []
        tmp_sizes = _compute_sizes_from_intensity(np.array(levels, dtype=float))
        for lvl, s in zip(levels, tmp_sizes):
            h = ax.scatter([], [], s=s, color="gray", alpha=0.6, edgecolors="none")
            handles.append(h)
            # Human-friendly label; adapt if using LOG_COLOR etc.
            labels.append(f"{lvl:.2g}")
        leg_title = getattr(cfg, "SIZE_LEGEND_TITLE", "Intensity (relative)")
        size_leg = ax.legend(handles, labels, title=leg_title, frameon=True, loc="upper right")
        ax.add_artist(size_leg)

    # Limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        # default signed MD limits
        ax.set_ylim(-0.55, 0.55)

    plt.tight_layout()

    # Save and/or show
    if savepath:
        fig.savefig(
            savepath,
            dpi=getattr(cfg, "DPI", 150),
            facecolor=getattr(cfg, "BACKGROUND", "white"),
            bbox_inches="tight",
        )
    if show:
        plt.show()

    return MDPlotResult(data=df_plot, fig=fig, ax=ax)


# -----------------------------
# Internals (file reading, streaming, utilities)
# -----------------------------

def _ensure_datetimes(time_range: Tuple[datetime, datetime]) -> Tuple[datetime, datetime]:
    if (not isinstance(time_range, tuple)) or len(time_range) != 2:
        raise TypeError("time_range must be a tuple of two datetime objects")
    t0, t1 = time_range
    if not isinstance(t0, datetime) or not isinstance(t1, datetime):
        raise TypeError("time_range must contain datetime objects only")
    if t0.tzinfo is None:
        t0 = t0.replace(tzinfo=timezone.utc)
    if t1.tzinfo is None:
        t1 = t1.replace(tzinfo=timezone.utc)
    return t0, t1


def _as_unix_seconds(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp()

def _compute_sizes_from_intensity(intensity: np.ndarray) -> np.ndarray:
    """
    Map intensity to marker sizes using config settings.
    - Supports linear or log scaling.
    - Uses robust percentile clipping to avoid outlier domination.
    Returns an array of sizes suitable for matplotlib scatter s=...
    """
    # Guard
    x = np.asarray(intensity, dtype=float)
    # Select transform
    scale = getattr(cfg, "SIZE_SCALE", "log").lower()
    if scale == "log":
        # log10(1 + I)
        x_t = np.log10(np.maximum(x, 0.0) + 1.0)
    elif scale == "linear":
        x_t = np.maximum(x, 0.0)
    else:
        raise ValueError("config.SIZE_SCALE must be 'log' or 'linear'")

    # Robust clipping
    plow = float(getattr(cfg, "SIZE_PLOW", 5.0))
    phigh = float(getattr(cfg, "SIZE_PHIGH", 95.0))
    lo = np.nanpercentile(x_t, plow) if np.isfinite(x_t).any() else 0.0
    hi = np.nanpercentile(x_t, phigh) if np.isfinite(x_t).any() else 1.0
    if not np.isfinite(lo): lo = 0.0
    if not np.isfinite(hi) or hi <= lo: hi = lo + 1.0  # avoid zero span

    # Normalize to 0..1
    norm = (np.clip(x_t, lo, hi) - lo) / (hi - lo)

    # Map to sizes
    smin = float(getattr(cfg, "SIZE_MIN", 8.0))
    smax = float(getattr(cfg, "SIZE_MAX", 80.0))
    sizes = smin + norm * (smax - smin)
    return sizes


def _detect_and_skip_header_rows(path: str) -> int:
    n_header_rows = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        line = f.readline()
        if not line:
            return 0
        low = line.strip().lower()
        if "number of header rows" in low:
            tokens = re.split(r"[:\t, ]+", low)
            for token in reversed(tokens):
                try:
                    n_header_rows = int(float(token))
                    break
                except Exception:
                    continue
    return n_header_rows


def _detect_delimiter(path: str) -> str:
    skiprows = _detect_and_skip_header_rows(path)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for _ in range(skiprows):
            f.readline()
        header = f.readline()
    tabs = header.count("\t")
    commas = header.count(",")
    return "\t" if tabs > commas else ","


def _read_table_with_header_skip(path: str) -> pd.DataFrame:
    skiprows = _detect_and_skip_header_rows(path)
    sep = _detect_delimiter(path)
    return pd.read_csv(path, skiprows=skiprows, sep=sep, engine="c")


def read_composition_table(path: str) -> pd.DataFrame:
    df = _read_table_with_header_skip(path)
    df.columns = [c.strip() for c in df.columns]
    colmap = {c.lower(): c for c in df.columns}
    if "sumformula" not in colmap or "mass" not in colmap:
        raise ValueError(f"Composition file must contain 'SumFormula' and 'Mass'. Found: {df.columns.tolist()}")
    out = df[[colmap["sumformula"], colmap["mass"]]].copy()
    out.columns = ["SumFormula", "Mass"]
    out["SumFormula"] = out["SumFormula"].astype(str).str.strip()
    out["Mass"] = pd.to_numeric(out["Mass"], errors="coerce")
    out = out.dropna(subset=["SumFormula", "Mass"])
    return out


def read_trace_header(path: str) -> dict:
    skiprows = _detect_and_skip_header_rows(path)
    sep = _detect_delimiter(path)
    df0 = pd.read_csv(path, skiprows=skiprows, sep=sep, engine="c", nrows=0)
    cols = [c.strip() for c in df0.columns]
    colmap = {c.lower(): c for c in cols}
    if "time" not in colmap or "unixtime" not in colmap:
        raise ValueError("Trace file must contain 'Time' and 'unixTime' columns.")
    time_col = colmap["time"]
    unix_col = colmap["unixtime"]
    intensity_cols = [c for c in cols if c not in {time_col, unix_col}]
    return {
        "time_col": time_col,
        "unix_col": unix_col,
        "intensity_cols": intensity_cols,
        "sep": sep,
        "skiprows": skiprows,
    }


def read_trace_table(path: str) -> pd.DataFrame:
    header = read_trace_header(path)
    df = pd.read_csv(
        path,
        skiprows=header["skiprows"],
        sep=header["sep"],
        engine="c",
    )
    df[header["unix_col"]] = pd.to_numeric(df[header["unix_col"]], errors="coerce")
    for c in header["intensity_cols"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=[header["unix_col"]])
    return df


def classical_mass_defect(m: Iterable[float]):
    """
    Classic MD only: MD = m - round(m)
    Returns (nominal_mass_array, defect_array)
    """
    m = np.asarray(m, dtype=float)
    nominal = np.rint(m).astype(int)
    defect = m - nominal
    return nominal, defect


def aggregate_window_unix_dt(
    df_tr: pd.DataFrame,
    t0_dt: datetime,
    t1_dt: datetime,
    how: str = "sum",
    clip_neg: bool = False,
) -> pd.Series:
    colmap = {c.lower(): c for c in df_tr.columns}
    time_col = colmap["time"]
    unix_col = colmap["unixtime"]
    intensity_cols = [c for c in df_tr.columns if c not in {time_col, unix_col}]
    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)
    sel = df_tr[(df_tr[unix_col] >= t0s) & (df_tr[unix_col] <= t1s)]
    if sel.empty:
        return pd.Series(dtype=float)
    if clip_neg:
        sel = sel.copy()
        sel[intensity_cols] = sel[intensity_cols].clip(lower=0)
    if how == "sum":
        agg = sel[intensity_cols].sum(axis=0)
    elif how == "mean":
        agg = sel[intensity_cols].mean(axis=0)
    else:
        raise ValueError("how must be 'sum' or 'mean'")
    return agg


def aggregate_window_streaming(
    trace_path: str,
    t0_dt: datetime,
    t1_dt: datetime,
    formulas: List[str],
    how: str = "sum",
    clip_neg: bool = False,
    assume_sorted_unix: bool = True,
    chunksize: int = 200_000,
    skiprows: Optional[int] = None,
    sep: Optional[str] = None,
) -> pd.Series:
    if skiprows is None:
        skiprows = _detect_and_skip_header_rows(trace_path)
    if sep is None:
        sep = _detect_delimiter(trace_path)

    header = read_trace_header(trace_path)
    unix_col = header["unix_col"]

    usecols = [unix_col] + formulas
    dtype_map = {unix_col: "float64"}
    dtype_map.update({f: "float64" for f in formulas})

    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)

    if how not in {"sum", "mean"}:
        raise ValueError("how must be 'sum' or 'mean'")

    acc = None
    rows_in_window = 0

    # For early-stop we need unfiltered mins/maxes per chunk; read unixTime twice:
    unix_usecols = [unix_col]  # small
    for chunk_unix, chunk in zip(
        pd.read_csv(
            trace_path, skiprows=skiprows, sep=sep, engine="c",
            usecols=unix_usecols, chunksize=chunksize, dtype={unix_col: "float64"}
        ),
        pd.read_csv(
            trace_path, skiprows=skiprows, sep=sep, engine="c",
            usecols=usecols, chunksize=chunksize, dtype=dtype_map
        )
    ):
        # Early-stop check with unfiltered bounds
        min_u = float(chunk_unix[unix_col].min())
        max_u = float(chunk_unix[unix_col].max())

        # Filter the second chunk (with intensities)
        chunk = chunk[pd.to_numeric(chunk[unix_col], errors="coerce").between(t0s, t1s)]
        if not chunk.empty:
            if clip_neg:
                chunk[formulas] = chunk[formulas].clip(lower=0)
            if acc is None:
                acc = chunk[formulas].sum(axis=0)
            else:
                acc = acc.add(chunk[formulas].sum(axis=0), fill_value=0.0)
            rows_in_window += len(chunk)

        # Early stop if sorted and we are past t1
        if assume_sorted_unix and max_u >= t1s:
            break
        # Also skip ahead quickly if this entire chunk is before t0
        if assume_sorted_unix and max_u < t0s:
            continue
        # If unsorted, we just keep scanning all chunks

    if acc is None:
        return pd.Series(index=formulas, dtype=float)

    if how == "mean":
        acc = acc / max(rows_in_window, 1)

    return acc.reindex(formulas)