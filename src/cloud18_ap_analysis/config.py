from pathlib import Path

PAPER_PLOT_WIDTH = 780
PAPER_STACKED_PANEL_HEIGHT = 200
PAPER_UNSTACKED_SIZE = (PAPER_PLOT_WIDTH / 100, 400 / 100)
DPI = 300

PAPER_LINE_COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#E69F00",
    "#56B4E9", "#CC79A7", "#F0E442", "#000000",
]

STAGE_SETTINGS = {
    "line_color": "0.6",
    "line_width": 1,
    "line_style": "--",
    "font_size": 10,
    "font_color": "0.4",
    "text_offset": 0.003,
    "max_length": 50,
}

OUTPUT_DIR = Path("output")

# mass defect plot settings
MD_FIG_SIZE = (14, 8)

# Default coloring mode: "OC" or "On"
COLOR_BY = "On"

# --- continuous O:C scale ---
OC_CMAP = "plasma"          # or "custom_oc" if you use the custom list
OC_COLORS = [
    "#7f3b08", "#d95f02", "#fdb863", "#ffffbf",
    "#91bfdb", "#4575b4", "#313695",
]
OC_VMIN = 0.0
OC_VMAX = 1.2
OC_NAN_COLOR = "lightgray"

# --- discrete number-of-oxygens scale ---
ON_CMAP = "turbo"         # base cmap that gets sampled into discrete bins
ON_MIN = 0                  # lowest O bin
ON_MAX = None                 # highest O bin, None -> auto from data
ON_CLIP_HIGH = False         # values > ON_MAX go into the top bin, labelled ">=12"
ON_NAN_COLOR = "lightgray"
ON_LABEL = "Number of O atoms"
ADDUCT_PRIORITY = {"NH3H+": 0, "H+": 1, "": 2}   # lower = preferred
REFERENCE_FAMILIES = ["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox"]  # e.g. ["C10H16Ox", "C10H17Ox"]


# Axis / label styling optional
XLABEL = "Compound m/z [Da]"
YLABEL = "Mass defect [Da]"

# --- size mapping ---
SIZE_BY_INTENSITY = True
SIZE_SCALE = "log"          # "log" or "linear"
ALPHA = 0.7
MIN_INTENSITY = 0.005      # ppt

SIZE_ANCHOR_I = 1.0     # ppt
SIZE_ANCHOR_S = 400.0      # pt² at the anchor intensity
SIZE_EXP = 0.5            # power for mapping intensity to marker size (1.0 = linear, 0.5 = sqrt, etc.)

# --- size legend ---
SIZE_LEGEND = True
SIZE_LEGEND_TITLE = "Mean intensity [ppt]"
SIZE_LEGEND_STYLE = "geom"   # "decade" | "geom" | "explicit"
SIZE_LEGEND_NMAX = 3           # max entries for "decade", exact count for "geom"
SIZE_LEGEND_LEVELS = [0.01, 1, 100]      # e.g. [0.01, 0.1, 1, 100] -> forces "explicit"
SIZE_LEGEND_LOC = "upper right"

# reference lines
REF_LINE_X_END = 320.0     # extend every family line to this nominal mass
REF_LINE_X_START = None    # optional: fixed left end; None = start at first member
REF_LINE_COLOR = "0.5"
REF_LINE_WIDTH = 1.1