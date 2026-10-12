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


## Concentration units

Both `scripts/run_plotting.py` and `scripts/run_temp_dep.py` have a
`concentration_unit` setting near the top:

```python
concentration_unit = "molecules/cm3"  # or "native" to restore mixing-ratio units
pressure_pa = 101325.0          # standard atmospheric pressure
```

The shared `cloud18_ap_analysis.units` module converts copies of the loaded
channels. `run_plotting` converts O3 (source ppbv) and HO2/OH/RO2 (source pptv);
`temp_dep` converts the selected compound traces and HO2 before calculating
interval means, SEMs, and the precursor-times-HO2 product. AP and IP always
use pptv (1000 times the VOCUS file's ppbv values), including the AP colorbar.

Conversion uses `mixing_fraction * pressure_pa / (R * temperature_K) * 1e-6 * N_A`
with Avogadro's constant `N_A = 6.02214076e23 molecules/mol`.
Measured `TE_Calib_3` temperature is linearly interpolated onto each channel's
UTC timestamps before smoothing or averaging. No extrapolation is performed:
points outside condition coverage become NaN. Temperature and pressure are
treated as known; their uncertainties are not propagated. The source files
and loaded datasets are not changed.

For `temp_dep`, each temperature run provides `temperature_path` and optionally
`temperature_channel` (default `TE_Calib_3`). Alternatively, provide an explicit
constant `temperature_c`. These inputs are only used in molecules/cm3 mode.

Time traces default to log for molecules/cm3 and linear for native pptv/ppbv.
Override with `concentration_scale` in `run_plotting`, or `trace_yscale` in
`run_temp_dep` (passed as `yscale` to `temp_dep`). Temperature and AP/IP panels
remain linear; SMPS diameter retains its logarithmic axis. Nonpositive
concentrations cannot be displayed on logarithmic axes.

Kinetic axes are logarithmic for molecules/cm3: y is molecules/cm3 and x is (molecules/cm3)^2.
Native mode uses linear axes in pptv and pptv^2. Axis scaling does not change
the original-space ODR fit; nonpositive values cannot appear on log axes. The fit uses an exact all-zero check so valid nonzero
values are accepted; regression results include `slope_unit`. Tests verify
that equivalent data expressed in either unit give equivalent slopes and
uncertainties without manual fit scaling.

Returned mean-table column names remain source identifiers (for example
`HO2_ppt`), even in molecules/cm3 mode; `means.attrs["concentration_unit"]` records
the chosen mode. AP color data have explicit `AP_pptv`, `AP_pptv_sem`, and
`AP_pptv_n_points` columns. Column names do not infer conversion automatically.
