# src/mdplot.py
import re
from dataclasses import dataclass
from typing import Optional, Tuple, Iterable, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timezone

from matplotlib.colors import (
    BoundaryNorm,
    LinearSegmentedColormap,
    ListedColormap,
)
import matplotlib.patches as mpatches
from .loading import _ensure_datetimes

# Pull plotting settings from the package config module; provide safe defaults
try:
    from . import config as cfg
except Exception:  # pragma: no cover
    class cfg:  # fallback defaults
        FIGSIZE = (9, 5.5)
        CMAP = "plasma"
        POINT_SIZE = 12.0
        ALPHA = 0.8
        LOG_COLOR = True
        MIN_INTENSITY = 0.0
        DPI = 150
        BACKGROUND = "white"




from .loading import (
    read_trace_header, read_composition_table, read_trace_table,
    aggregate_window_unix_dt, aggregate_window_streaming,
    _detect_and_skip_header_rows, _detect_delimiter,
    classical_mass_defect, element_counts,
    select_preferred_adducts,
)

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
    clip_neg: bool = False,                 # clip negative intensities to 0 before aggregating
    mass_range: Optional[Tuple[float, float]] = None, # e.g., (100, 350)
    color_by: Optional[str] = None,        # "OC" or "On"; None -> config.COLOR_BY
    reference_families: Optional[List[str]] = None,   # e.g. ["C10H16Ox", "C10H17Ox"]
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
    - clip_neg: clip negative intensities to 0 before aggregating
    - mass_range: range of masses to include in the plot, e.g. (100, 350)
    - reference_families: list of homologous families, e.g. ["C10H16Ox", "C9H16Ox"];
      each is drawn as a dashed line through the actual plotted compounds
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
    header = read_trace_header(trace_path)
    all_trace_formulas = header["intensity_cols"]

    comp = read_composition_table(composition_path, available=all_trace_formulas)
    comp["SumFormula"] = comp["SumFormula"].astype(str).str.strip()

    # MD plots use one representative ion per neutral compound, preserving
    # the historical preference order from config.ADDUCT_PRIORITY.
    comp = select_preferred_adducts(comp)

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
            clip_neg=clip_neg,
            assume_sorted_unix=assume_sorted_unix,
            chunksize=chunksize,
            skiprows=_detect_and_skip_header_rows(trace_path),
            sep=_detect_delimiter(trace_path),
        )
    else:
        traces = read_trace_table(trace_path)
        agg_series = aggregate_window_unix_dt(traces, t0_dt, t1_dt, clip_neg=clip_neg)
        agg_series = agg_series.reindex(common)

    agg_series = agg_series.dropna()

    # Build plotting table
    df_plot = comp[comp["SumFormula"].isin(common)].merge(
        agg_series.rename("AggIntensity"), left_on="SumFormula", right_index=True, how="inner"
    )

    # Filter by mass range, if provided
    if mass_range is not None:
        if not (isinstance(mass_range, (list, tuple)) and len(mass_range) == 2 and mass_range[0] < mass_range[1]):
            raise ValueError("mass_range must be a tuple or list of two numbers (min_mass, max_mass).")
        
        n_before_mass_filter = len(df_plot)
        m_min, m_max = mass_range
        
        # Filter based on the NeutralMass column
        df_plot = df_plot[
            (df_plot["NeutralMass"] >= m_min) & (df_plot["NeutralMass"] <= m_max)
        ].copy()
        
        n_dropped_mass = n_before_mass_filter - len(df_plot)
        if n_dropped_mass:
            print(f"[md_plot] dropped {n_dropped_mass}/{n_before_mass_filter} compounds "
                  f"outside mass range [{m_min:.1f}, {m_max:.1f}]")

    i_cut = float(getattr(cfg, "MIN_INTENSITY", 0.005))

    n_before = len(df_plot)
    df_plot = df_plot[df_plot["AggIntensity"] >= i_cut].copy()
    n_dropped = n_before - len(df_plot)
    if n_dropped:
        print(f"[md_plot] dropped {n_dropped}/{n_before} compounds with "
              f"intensity < {i_cut:.3g} ppt")

    if df_plot.empty:
        raise RuntimeError("No data to plot after time-window aggregation and filtering.")

    # Classical mass defect: MD = m - round(m)
    nominal, defect = classical_mass_defect(df_plot["NeutralMass"].values)
    df_plot["NominalMass"] = nominal
    df_plot["MassDefect"] = defect

    # Color mapping (continuous O:C ratio or discrete number of O atoms)
    mode = color_by if color_by is not None else getattr(cfg, "COLOR_BY", "OC")
    cspec = _make_color_spec(df_plot, mode)

    # Sizes
    if getattr(cfg, "SIZE_BY_INTENSITY", True):
        intensities = df_plot["AggIntensity"].values

        sizes = _compute_sizes_from_intensity(
            intensities,
            reference_intensity=intensities,
        )
    else:
        sizes = float(getattr(cfg, "POINT_SIZE", 12.0))

    # Plot
    fig, ax = plt.subplots(figsize=getattr(cfg, "MD_FIG_SIZE", (9, 5.5)))

    scatter_kwargs = dict(
        c=cspec.values,
        s=sizes,
        cmap=cspec.cmap,
        alpha=getattr(cfg, "ALPHA", 0.7),
        edgecolors="black",
        linewidth=1.0,
    )

    # BoundaryNorm (discrete) and vmin/vmax (continuous) are mutually exclusive
    if cspec.norm is not None:
        scatter_kwargs["norm"] = cspec.norm
    else:
        scatter_kwargs["vmin"] = cspec.vmin
        scatter_kwargs["vmax"] = cspec.vmax

    sc = ax.scatter(
        df_plot["NominalMass"].values,
        df_plot["MassDefect"].values,
        **scatter_kwargs,
    )
    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.set_xlabel(getattr(cfg, "XLABEL", "Nominal neutral mass"))
    ax.set_ylabel(getattr(cfg, "YLABEL", "Mass defect"))
    ax.set_title(f"Mass Defect Plot, UTC ∈ [{t0_dt.isoformat()} — {t1_dt.isoformat()}]")

    if cspec.ticks is not None:
        cbar = plt.colorbar(sc, ax=ax, ticks=cspec.ticks, spacing="proportional")
        if cspec.tick_labels is not None:
            cbar.ax.set_yticklabels(cspec.tick_labels)
    else:
        cbar = plt.colorbar(sc, ax=ax)

    cbar.set_label(cspec.label)
    # Optional: size legend so readers can interpret marker sizes
    if getattr(cfg, "SIZE_BY_INTENSITY", True) and getattr(cfg, "SIZE_LEGEND", True):
        intens = df_plot["AggIntensity"].values

        levels = [v for v in _size_legend_levels(intens) if v >= i_cut]
        level_sizes = _compute_sizes_from_intensity(
            np.array(levels, dtype=float),
            reference_intensity=intens,
        )

        handles = [
            ax.scatter([], [], s=s, color="0.35", alpha=0.6, edgecolors="black")
            for s in level_sizes
        ]
        labels = [_format_size_label(v) for v in levels]

        dummy_handle = mpatches.Patch(color='none')
        dummy_label = ""
        handles.insert(1, dummy_handle)
        labels.insert(1, dummy_label)
        handles.insert(0, dummy_handle)
        labels.insert(0, dummy_label)
        size_leg = ax.legend(
            handles,
            labels,
            title=getattr(cfg, "SIZE_LEGEND_TITLE", "Mean intensity [ppt]"),
            loc=getattr(cfg, "SIZE_LEGEND_LOC", "upper right"),
            frameon=True,
            framealpha=0.9,
            labelspacing=1.5,        # room for the big circles
            handletextpad=3.5,
            borderpad=1.8,
            scatterpoints=1,
            fontsize=10,
            title_fontsize=9,
        )
        ax.add_artist(size_leg)

    # Limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        # default signed MD limits
        ax.set_ylim(-0.55, 0.55)

    families = reference_families
    if families is None:
        families = getattr(cfg, "REFERENCE_FAMILIES", [])

    x_end = getattr(cfg, "REF_LINE_X_END", 320.0)
    x_start = getattr(cfg, "REF_LINE_X_START", None)

    for family in families:
        fam = (_family_members(df_plot, family)
               .dropna(subset=["On", "MassDefect"])
               .sort_values("NominalMass"))
        if len(fam) < 2:
            continue

        x = fam["NominalMass"].to_numpy(dtype=float)
        y = fam["MassDefect"].to_numpy(dtype=float)

        xs, ys = list(x), list(y)

        # Extend right, along the slope of the last segment
        if x_end is not None and x_end > x[-1]:
            s = (y[-1] - y[-2]) / (x[-1] - x[-2])
            xs.append(float(x_end))
            ys.append(y[-1] + s * (x_end - x[-1]))

        # Optional: extend left, along the slope of the first segment
        if x_start is not None and x_start < x[0]:
            s = (y[1] - y[0]) / (x[1] - x[0])
            xs.insert(0, float(x_start))
            ys.insert(0, y[0] + s * (x_start - x[0]))

        ax.plot(xs, ys,
                ls="--",
                lw=getattr(cfg, "REF_LINE_WIDTH", 1.1),
                color=getattr(cfg, "REF_LINE_COLOR", "0.5"),
                zorder=0.5)

        ax.annotate(
            family, (xs[-1], ys[-1]),
            textcoords="offset points", xytext=(5, 0),
            fontsize=12, color="0.35", ha="left", va="center",
            annotation_clip=False,
        )

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

@dataclass
class _ColorSpec:
    values: np.ndarray
    cmap: object
    norm: Optional[object]
    vmin: Optional[float]
    vmax: Optional[float]
    label: str
    ticks: Optional[np.ndarray] = None
    tick_labels: Optional[List[str]] = None

def _get_oc_cmap():
    name = getattr(cfg, "OC_CMAP", "viridis")
    if name == "custom_oc":
        cmap = LinearSegmentedColormap.from_list(
            "custom_oc",
            getattr(cfg, "OC_COLORS",
                    ["#7f3b08", "#fdb863", "#ffffbf", "#91bfdb", "#313695"]),
            N=256,
        )
    else:
        cmap = plt.get_cmap(name).copy()
    cmap.set_bad(getattr(cfg, "OC_NAN_COLOR", "lightgray"))
    return cmap


def _make_color_spec(df_plot: pd.DataFrame, color_by: str) -> _ColorSpec:
    """
    Build the color mapping for the MD plot.

    color_by:
      "OC" / "oc_ratio"  -> continuous O:C ratio
      "On" / "o_number"  -> discrete number of oxygen atoms
    """
    key = str(color_by).strip().lower().replace(":", "").replace("_", "")

    # ---------- continuous O:C ----------
    if key in {"oc", "ocratio", "oc-ratio"}:
        values = np.ma.masked_invalid(df_plot["OC"].values.astype(float))
        return _ColorSpec(
            values=values,
            cmap=_get_oc_cmap(),
            norm=None,
            vmin=getattr(cfg, "OC_VMIN", None),
            vmax=getattr(cfg, "OC_VMAX", None),
            label="O:C ratio",
        )

    # ---------- discrete number of oxygens ----------
    if key in {"on", "o", "onumber", "numo", "nox"}:
        raw = df_plot["On"].values.astype(float)

        omin = int(getattr(cfg, "ON_MIN", 0))
        omax = getattr(cfg, "ON_MAX", None)
        if omax is None:
            finite = raw[np.isfinite(raw)]
            omax = int(np.nanmax(finite)) if finite.size else omin
        omax = int(omax)
        if omax < omin:
            omax = omin

        clip_high = bool(getattr(cfg, "ON_CLIP_HIGH", True))
        vals = raw.copy()
        if clip_high:
            vals = np.where(np.isfinite(vals), np.clip(vals, omin, omax), vals)

        values = np.ma.masked_invalid(vals)

        n_bins = omax - omin + 1
        base = plt.get_cmap(getattr(cfg, "ON_CMAP", "viridis"))
        colors = base(np.linspace(0.0, 1.0, n_bins))
        cmap = ListedColormap(colors)
        cmap.set_bad(getattr(cfg, "ON_NAN_COLOR", "lightgray"))

        boundaries = np.arange(omin - 0.5, omax + 1.5, 1.0)
        norm = BoundaryNorm(boundaries, cmap.N)

        ticks = np.arange(omin, omax + 1, 1.0)
        tick_labels = [str(int(t)) for t in ticks]
        finite_raw = raw[np.isfinite(raw)]
        if clip_high and finite_raw.size and finite_raw.max() > omax:
            tick_labels[-1] = f"\u2265{omax}"
        return _ColorSpec(
            values=values,
            cmap=cmap,
            norm=norm,
            vmin=None,
            vmax=None,
            label=getattr(cfg, "ON_LABEL", "Number of O atoms"),
            ticks=ticks,
            tick_labels=tick_labels,
        )

    raise ValueError("color_by must be 'OC' (O:C ratio) or 'On' (number of O atoms)")

_formula_token = re.compile(r"([A-Z][a-z]?)(\d*)")

_FAMILY_RE = re.compile(r"^([A-Za-z0-9()]*?)O?x$")


def _family_members(comp: pd.DataFrame, family: str) -> pd.DataFrame:
    """
    'C10H16Ox' -> all rows whose CompoundFormula is C10H16 + n oxygens
    (n >= 0), sorted by oxygen count.
    """
    m = _FAMILY_RE.match(family.strip())
    if not m:
        raise ValueError(f"Bad family spec {family!r}; expected e.g. 'C10H16Ox'")
    base = element_counts(m.group(1))
    base.pop("O", None)

    def matches(f):
        try:
            c = element_counts(f)
        except ValueError:
            return False
        c.pop("O", None)
        return c == base

    sel = comp[comp["CompoundFormula"].apply(matches)]
    return sel.sort_values("On")




def _compute_sizes_from_intensity(x, reference_intensity=None) -> np.ndarray:
    """
    Absolute power-law size mapping, no clipping:

        s / SIZE_ANCHOR_S = (I / SIZE_ANCHOR_I) ** SIZE_EXP

    Non-positive intensities return NaN — they must be filtered out
    upstream, not silently floored.
    """
    i_ref = float(getattr(cfg, "SIZE_ANCHOR_I", 100.0))
    s_ref = float(getattr(cfg, "SIZE_ANCHOR_S", 45.0))
    p     = float(getattr(cfg, "SIZE_EXP", 1.0))

    x = np.asarray(x, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(x > 0.0, s_ref * np.power(x / i_ref, p), np.nan)
    return s

def _size_legend_levels(reference: np.ndarray) -> List[float]:
    """
    Pick nice intensity values for the size legend.

    "decade"   -> whole powers of ten inside [vmin, vmax]  (recommended for log)
    "geom"     -> geometrically spaced values (exponential distribution)
    "explicit" -> whatever config.SIZE_LEGEND_LEVELS says
    """
    explicit = getattr(cfg, "SIZE_LEGEND_LEVELS", None)
    if explicit:
        return [float(v) for v in explicit]

    vmin, vmax = _size_limits(reference)
    scale = getattr(cfg, "SIZE_SCALE", "log").lower()
    style = str(getattr(cfg, "SIZE_LEGEND_STYLE", "decade")).lower()
    nmax = max(int(getattr(cfg, "SIZE_LEGEND_NMAX", 4)), 2)

    if scale != "log":
        return list(np.linspace(vmin, vmax, nmax))

    if style == "decade":
        e0 = int(np.ceil(np.log10(vmin) - 1e-9))
        e1 = int(np.floor(np.log10(vmax) + 1e-9))
        exps = list(range(e0, e1 + 1))

        if len(exps) < 2:                      # range narrower than one decade
            return list(np.geomspace(vmin, vmax, nmax))

        if len(exps) > nmax:                   # thin out, keeping both endpoints
            idx = np.unique(np.round(np.linspace(0, len(exps) - 1, nmax)).astype(int))
            exps = [exps[i] for i in idx]

        return [10.0 ** e for e in exps]

    # "geom"
    return list(np.geomspace(vmin, vmax, nmax))


def _format_size_label(v: float) -> str:
    """Compact label: plain numbers in the readable range, 10^n for decades."""
    if not np.isfinite(v) or v <= 0:
        return f"{v:g}"

    e = np.log10(v)
    is_decade = abs(e - round(e)) < 1e-9

    if 1e-2 <= v < 1e4:
        if v >= 100:
            return f"{v:.0f}"
        if v >= 1:
            return f"{v:.3g}"
        return f"{v:.3g}".rstrip("0").rstrip(".")

    if is_decade:
        return rf"$10^{{{int(round(e))}}}$"

    mant = v / 10.0 ** np.floor(e)
    return rf"${mant:.1f}\times10^{{{int(np.floor(e))}}}$"

    """
    Classic MD only: MD = m - round(m)
    Returns (nominal_mass_array, defect_array)
    """
    m = np.asarray(m, dtype=float)
    nominal = np.rint(m).astype(int)
    defect = m - nominal
    return nominal, defect

