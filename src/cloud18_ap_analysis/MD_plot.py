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
    df_plot = df_plot[df_plot["AggIntensity"] > i_cut].copy()
    n_dropped = n_before - len(df_plot)
    if n_dropped:
        print(f"[md_plot] dropped {n_dropped}/{n_before} compounds with "
              f"intensity <= {i_cut:.3g} ppt")

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

        levels = [v for v in _size_legend_levels(intens) if v > i_cut]
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

# Ion/adduct masses to subtract from the exported Mass column
# Your examples suggest Mass includes the ion/adduct.
ELEMENT_MASSES = {
    "H": 1.00782503223,   "H(2)":  2.01410177812,
    "C": 12.0,            "C(13)": 13.00335483507,
    "N": 14.00307400443,  "N(15)": 15.00010889889,
    "O": 15.99491461957,  "O(18)": 17.99915961286,
    "S": 31.9720711744,   "S(34)": 33.967867004,
    "P": 30.97376199842,  "F": 18.99840316273,
    "Cl": 34.968852682,   "Cl(37)": 36.965902602,
    "Br": 78.9183376,     "I": 126.9044719,
}

_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\(\d+\))?(\d*)")
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


def split_sumformula(sumformula: str) -> Tuple[str, str]:
    """
    Split something like:
        CH2O.H+       -> ("CH2O", "H+")
        C2H2O.NH3H+   -> ("C2H2O", "NH3H+")
        C10H16O2      -> ("C10H16O2", "")
    """
    s = str(sumformula).strip()
    if "." not in s:
        return s, ""
    compound, adduct = s.split(".", 1)
    return compound.strip(), adduct.strip()

def element_counts(formula: str) -> dict:
    """'C(13)C9H16O3' -> {'C(13)': 1, 'C': 9, 'H': 16, 'O': 3}"""
    s = str(formula).strip()
    counts, pos = {}, 0
    for m in _TOKEN_RE.finditer(s):
        if m.start() != pos:
            raise ValueError(f"Cannot parse {formula!r} at {pos}")
        pos = m.end()
        sym, iso, num = m.groups()
        token = f"{sym}{iso}" if iso else sym
        counts[token] = counts.get(token, 0) + (int(num) if num else 1)
    if pos != len(s):
        raise ValueError(f"Cannot parse {formula!r} at {pos}")
    return counts


def exact_mass(formula: str) -> float:
    total = 0.0
    for token, n in element_counts(formula).items():
        if token not in ELEMENT_MASSES:
            return np.nan
        total += ELEMENT_MASSES[token] * n
    return total


def _element_total(counts: dict, symbol: str) -> int:
    """Sum an element across isotopes: 'C' matches 'C' and 'C(13)'."""
    return sum(n for t, n in counts.items()
               if t == symbol or t.startswith(f"{symbol}("))@dataclass
def oxygen_to_carbon_ratio(formula: str) -> float:
    """
    Calculate O:C ratio from the neutral compound formula.
    Returns NaN if C=0.
    """
    counts = element_counts(formula)
    c = counts.get("C", 0)
    o = counts.get("O", 0)

    if c == 0:
        return np.nan

    return o / c

def oxygen_count(formula: str) -> float:
    """
    Number of oxygen atoms in the neutral compound formula.
    Returns NaN if the formula cannot be parsed at all.
    """
    counts = element_counts(formula)
    if not counts:
        return np.nan
    return float(counts.get("O", 0))

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


def read_composition_table(path: str, available: Optional[Iterable[str]] = None) -> pd.DataFrame:
    df = _read_table_with_header_skip(path)
    df.columns = [c.strip() for c in df.columns]
    colmap = {c.lower(): c for c in df.columns}

    if "sumformula" not in colmap or "mass" not in colmap:
        raise ValueError(
            f"Composition file must contain 'SumFormula' and 'Mass'. "
            f"Found: {df.columns.tolist()}"
        )

    out = df[[colmap["sumformula"], colmap["mass"]]].copy()
    out.columns = ["SumFormula", "Mass"]

    out["SumFormula"] = out["SumFormula"].astype(str).str.strip()
    out["Mass"] = pd.to_numeric(out["Mass"], errors="coerce")
    out = out.dropna(subset=["SumFormula", "Mass"])

    # Keep only ions the trace file can actually provide, so dedup never
    # discards a compound whose preferred adduct has no trace column.
    if available is not None:
        out = out[out["SumFormula"].isin(set(available))]
        if out.empty:
            raise RuntimeError("No composition SumFormula matches any trace column.")

    # Split e.g. CH2O.H+ -> CH2O + H+
    split = out["SumFormula"].apply(split_sumformula)
    out["CompoundFormula"] = split.apply(lambda x: x[0])
    out["Adduct"] = split.apply(lambda x: x[1])

    # Subtract adduct mass to get neutral compound mass
    out["NeutralMass"] = out["CompoundFormula"].apply(exact_mass)

    # O:C ratio from neutral compound formula
    out["OC"] = out["CompoundFormula"].apply(oxygen_to_carbon_ratio)

    # Number of oxygen atoms (for discrete coloring)
    out["On"] = out["CompoundFormula"].apply(oxygen_count)

    ADDUCT_PRIORITY = getattr(cfg, "ADDUCT_PRIORITY", {"NH3H+": 0, "H+": 1, "": 2})

    out["_prio"] = out["Adduct"].map(ADDUCT_PRIORITY).fillna(99)
    out = (out.sort_values(["CompoundFormula", "_prio"])
            .drop_duplicates("CompoundFormula", keep="first")
            .drop(columns="_prio"))

    print(df.head())

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
    print(df.head())

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
    return sel[intensity_cols].mean(axis=0)


def aggregate_window_streaming(
    trace_path: str,
    t0_dt: datetime,
    t1_dt: datetime,
    formulas: List[str],
    clip_neg: bool = False,
    assume_sorted_unix: bool = True,
    chunksize: int = 200_000,
    skiprows: Optional[int] = None,
    sep: Optional[str] = None,
) -> pd.Series:
    """
    Aggregate trace intensities inside a time window without loading the full CSV.

    The trace file can be very large, so this function reads it in chunks and
    only keeps rows whose ``unixTime`` lies inside ``[t0_dt, t1_dt]``. The
    returned series is indexed by ``formulas`` and contains the aggregated
    intensity for each compound column.

    Workflow:
    - detect the trace file header layout and delimiter
    - read the ``unixTime`` column plus the requested compound columns in chunks
    - filter each chunk to the requested time window
    - optionally clip negative intensities to zero
    - accumulate the mean intensity across all selected rows

    Parameters
    ----------
    trace_path:
        Path to the trace CSV/TXT file.
    t0_dt, t1_dt:
        Start and end of the aggregation window.
    formulas:
        Compound columns to aggregate.
    clip_neg:
        If True, negative values are clipped to zero before aggregation.
    assume_sorted_unix:
        If True, stop scanning once the file is past the end of the window.
    chunksize:
        Number of rows to read per chunk.
    skiprows, sep:
        Optional header settings; if omitted they are detected automatically.
    """
    if skiprows is None:
        skiprows = _detect_and_skip_header_rows(trace_path)
    if sep is None:
        sep = _detect_delimiter(trace_path)

    header = read_trace_header(trace_path)
    unix_col = header["unix_col"]

    # Read only the timestamp column and the compounds we care about.
    usecols = [unix_col] + formulas
    dtype_map = {unix_col: "float64"}
    dtype_map.update({f: "float64" for f in formulas})

    # Convert the datetime window to unix seconds for numeric comparison.
    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)

    acc = None
    rows_in_window = 0

    # To support early stopping, we read the unixTime column separately so we
    # can inspect each chunk's bounds before processing the intensity columns.
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
        # Chunk bounds tell us whether this chunk is entirely before or after
        # the requested window when the file is sorted by time.
        min_u = float(chunk_unix[unix_col].min())
        max_u = float(chunk_unix[unix_col].max())

        # Keep only rows inside the requested time window.
        chunk = chunk[pd.to_numeric(chunk[unix_col], errors="coerce").between(t0s, t1s)]
        if not chunk.empty:
            # Optionally suppress negative intensities before aggregating.
            if clip_neg:
                chunk[formulas] = chunk[formulas].clip(lower=0)

            # Accumulate the column-wise sum for this chunk. The mean is
            # computed once at the end so chunk boundaries do not bias it.
            chunk_sum = chunk[formulas].sum(axis=0)
            if acc is None:
                acc = chunk_sum
            else:
                acc = acc.add(chunk_sum, fill_value=0.0)
            rows_in_window += len(chunk)

        # Stop early once the file has moved past the end of the time window.
        if assume_sorted_unix and max_u >= t1s:
            break
        # Skip chunks that are entirely before the window.
        if assume_sorted_unix and max_u < t0s:
            continue
        # If unsorted, we keep scanning all chunks.

    if acc is None:
        return pd.Series(index=formulas, dtype=float)

    acc = acc / max(rows_in_window, 1)

    print(acc.head(15))

    return acc.reindex(formulas)