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

FIGSIZE = (9, 5.5)
CMAP = "viridis"
POINT_SIZE = 12.0
ALPHA = 0.8
LOG_COLOR = True       # use log10(intensity + 1) for color
MIN_INTENSITY = 0.0    # filter out weak peaks
DPI = 150              # when saving
BACKGROUND = "white"   # facecolor for saved figs (optional)

# Size settings for marker size based on intensity
SIZE_BY_INTENSITY = True     # enable size ~ intensity
SIZE_SCALE = "linear"           # "log" or "linear"
SIZE_MIN = 8.0               # smallest marker size (points^2 in matplotlib)
SIZE_MAX = 80.0              # largest marker size
SIZE_PLOW = 5.0              # lower percentile for robust clipping
SIZE_PHIGH = 95.0            # upper percentile for robust clipping
SIZE_LEGEND = True
SIZE_LEGEND_LEVELS = None    # None = use percentiles [10, 50, 90] of intensity in window
SIZE_LEGEND_TITLE = "Intensity (relative)"