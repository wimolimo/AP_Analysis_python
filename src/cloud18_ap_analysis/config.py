"""Central configuration for CLOUD18 AP analysis and plotting."""

# =============================================================================
# General publication style
# =============================================================================

DPI = 300

# Existing overview/diagnostic plots in plotting.py
PAPER_UNSTACKED_SIZE = (10.0, 4.0)

# Shared line palette used by plotting.py
PAPER_LINE_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#E69F00",
    "#56B4E9",
    "#CC79A7",
    "#F0E442",
    "#000000",
]

# Common styling for publication figures
PAPER_FONT_SIZE = 7
PAPER_LABEL_SIZE = 8
PAPER_LINE_WIDTH = 1.0
PAPER_SAVE_DPI = DPI

PAPER_RC_PARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": PAPER_FONT_SIZE,
    "axes.labelsize": PAPER_LABEL_SIZE,
    "axes.titlesize": PAPER_LABEL_SIZE,
    "axes.linewidth": 0.8,
    "xtick.labelsize": PAPER_FONT_SIZE,
    "ytick.labelsize": PAPER_FONT_SIZE,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.fontsize": PAPER_FONT_SIZE,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
}


# =============================================================================
# Stage annotations
# =============================================================================

STAGE_SETTINGS = {
    "line_color": "0.6",
    "line_width": 1.0,
    "line_style": "--",
    "font_size": 10,
    "font_color": "0.4",
    "text_offset": 0.003,
    "max_length": 50,
}


# =============================================================================
# Temperature-dependence / kinetic plots
# =============================================================================

# Time traces are usually wider; kinetic relationship plots are compact.
TRACE_FIGSIZE = (10.0, 4.0)
KINETIC_FIGSIZE = (4.0, 4.0)

# Important averaging intervals T1, T2, ...
IMPORTANT_RANGE_COLOR = "#8E7CC3"
IMPORTANT_RANGE_ALPHA = 0.12

# Interval-mean kinetic relationship plot
KINETIC_MARKER_SIZE = 24
KINETIC_MARKER_COLOR = "black"


# =============================================================================
# Mass-defect analysis
# =============================================================================

MD_FIG_SIZE = (14.0, 8.0)

# Preferred ion/adduct when several traces represent the same neutral compound.
# Lower number = higher priority.
ADDUCT_PRIORITY = {
    "NH3H+": 0,
    "H+": 1,
    "": 2,
}

# Axis labels
XLABEL = "Compound m/z [Da]"
YLABEL = "Mass defect [Da]"

# Default color variable: "OC" or "On"
COLOR_BY = "On"

# --- continuous O:C color scale ---
OC_CMAP = "plasma"  # use "custom_oc" to use OC_COLORS below
OC_COLORS = [
    "#7f3b08",
    "#d95f02",
    "#fdb863",
    "#ffffbf",
    "#91bfdb",
    "#4575b4",
    "#313695",
]
OC_VMIN = 0.0
OC_VMAX = 1.2
OC_NAN_COLOR = "lightgray"

# --- discrete number-of-oxygen-atoms color scale ---
ON_CMAP = "turbo"
ON_MIN = 0
ON_MAX = None          # None -> determine maximum automatically from data
ON_CLIP_HIGH = False
ON_NAN_COLOR = "lightgray"
ON_LABEL = "Number of O atoms"

# --- marker size from mean intensity ---
SIZE_BY_INTENSITY = True
SIZE_SCALE = "log"     # "log" or "linear"
ALPHA = 0.7
MIN_INTENSITY = 0.01   # ppt

SIZE_ANCHOR_I = 1.0    # intensity [ppt] at the size anchor
SIZE_ANCHOR_S = 400.0  # marker area [pt^2] at SIZE_ANCHOR_I
SIZE_EXP = 0.5         # 1 = linear, 0.5 = square-root scaling

# --- marker-size legend ---
SIZE_LEGEND = True
SIZE_LEGEND_TITLE = "Mean intensity [ppt]"
SIZE_LEGEND_STYLE = "geom"  # "decade", "geom", or "explicit"
SIZE_LEGEND_NMAX = 3
SIZE_LEGEND_LEVELS = [100, 1, 0.01]
SIZE_LEGEND_LOC = "upper right"

# --- homologous-family reference lines ---
REFERENCE_FAMILIES = [
    "C10H15Ox",
    "C10H16Ox",
    "C10H17Ox",
    "C10H18Ox",
]

REF_LINE_X_START = None
REF_LINE_X_END = 320.0
REF_LINE_COLOR = "0.5"
REF_LINE_WIDTH = 1.1
