# CLOUD18 AP Analysis — Python

Python conversion of the original Julia `CLOUD18_AP_Analysis` project.

## Current scope

This first conversion includes:

- normal paired-column CSV loading
- custom `NumberOfHeaderRows` loading
- AP/VOCUS CSV loading
- `Channel` and `Dataset` data structures
- linear interpolation / resampling
- time-window filtering
- moving time-window smoothing
- stage-log parsing and annotations
- single-channel and multi-channel plots
- stacked plots
- multiple datasets
- dual y-axes
- interactive **or** non-interactive plotting via `interactive=True/False`

**SMPS loading and SMPS heatmap plotting are intentionally omitted for now.**
They will be implemented separately later.

## Installation

On Windows:

```powershell
python -m venv .venv
.venv\Scriptsctivate
pip install -e .
```

## Basic usage

```python
from datetime import datetime, timedelta
from cloud18_ap_analysis import load_data, plot_data, load_stages

data = load_data(
    r"C:\path	o\data.csv",
    time_range=(
        datetime(2025, 10, 2, 4),
        datetime(2025, 10, 2, 13),
    ),
)

load_stages(r"C:\path	o\stages.txt")

plot_data(
    data,
    channels=["TE_Calib_3", "DewPoint_C"],
    smoothing=timedelta(minutes=10),
    stages="onlytype",
    interactive=True,
)
```

## Interactive vs. non-interactive

The public plotting functions all expose the same switch:

```python
plot_data(data, interactive=True)
```

With `interactive=True`, Matplotlib opens a normal GUI figure with zoom, pan,
save, etc.

With `interactive=False`, the figure is saved to `savepath` (or `display.png`
when no path is supplied) and the saved image is opened with the operating
system's default image viewer.

For a script that produces several figures, use `interactive=False` to avoid
blocking on every `plt.show()`.

## Notes on smoothing

A smoothing value can be a `datetime.timedelta`:

```python
smoothing=timedelta(minutes=10)
```

or a number of seconds:

```python
smoothing=600
```

The implementation uses the same centered, time-based averaging concept as the
Julia version; it does not simply average every N samples.

## Project structure

```text
CLOUD18_AP_Analysis/
├── src/cloud18_ap_analysis/
├── scripts/
├── tests/
├── config/
├── pyproject.toml
├── README.md
└── .gitignore
```
