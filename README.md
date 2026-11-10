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

## Prepared HDF5 run files

The four run scripts now read one compressed file per run from
`output/prepared/`. Original exports stay unchanged. Install/update the project
dependencies with `python -m pip install -e .`, then prepare the files:

```powershell
python scripts/prepare_runs.py
# Or prepare just one run:
python scripts/prepare_runs.py m25
# Rebuild after changing input files or stages (not display ranges):
python scripts/prepare_runs.py m25 --overwrite
```

Then run `run_plotting.py`, `run_time_traces.py`, `run_temp_dep.py`, or
`run_md_plot.py` as before. Select `run_id` and plot options in those scripts.
Preparation keeps all samples from every explicitly configured input file,
independent of `RUNS[run_id]["time_range"]`. That setting now controls the
default display/analysis range, read fresh by `load_prepared_run` at plot time.
Change it or override `time_range` in a script without rebuilding. Actual stored
coverage is available as `run["prepared_time_range"]`; instruments may have
different coverage or gaps. Files are not automatically discovered. Rebuild
older cropped H5 files once to recover full configured-source coverage.

Each file contains PTR, temperature, ozone, HORUS, VOCUS, fan-speed and SMPS
data, matching PTR composition tables, stages, and averaging intervals.
Instrument groups preserve independent UTC timestamps and native float64 values.
PTR timestamps consistently use the numeric `unixTime` clock, including when
the export's text timestamps mix different fractional-second formats.
There is no preparation-time interpolation, smoothing, concentration conversion,
or downsampling. Shared time axes within each instrument are stored once.
All measured temperature samples are retained for concentration conversion.
The 30-minute interpolation gap limit still applies when converting conditions.

All source stage records are retained, including those without timestamps,
with original stage numbers, types, descriptions, comments and details. Existing
plot label options are preserved. H5 retains the preparation-time averaging
intervals as provenance. `load_prepared_run` uses the current `intervals.py`
for known campaign runs, so edits immediately affect `run_plotting.py`,
`run_time_traces.py`, and `run_temp_dep.py` without rebuilding measurements.
Restart/rerun the script after editing. Stage-file edits still require rebuilding.

For calculations or a new plot, read just the channels and time window needed:

```python
from cloud18_ap_analysis.prepared import load_prepared_run, prepared_channel_info
from cloud18_ap_analysis import load_data, load_stages

run = load_prepared_run("m25")
data = load_data(run["sources"]["horus"], run["time_range"], channels=["HO2_ppt"])
load_stages(run["stages_path"])
info = prepared_channel_info(run["sources"]["horus"])  # units, coverage, missing counts
```

Channel selection reads only those channel arrays; window selection uses an
on-disk binary search and reads only the requested slices. It does not load the
whole file. Channels are qualified by instrument, so matching names from different
instruments do not collide. SMPS stays a separate time-by-diameter matrix.

To add future data, edit `source_specs(run_id)` in `runs.py`. Each `SourceSpec`
declares files, explicit format and native units. Legacy boundary-padding
settings are unnecessary because preparation now retains full source coverage.
Supported formats are `trace`, `paired`, `custom`, `ap`, and `smps`. New export
formats need one reader in the preparation layer. Units not configured are
marked `unknown`; they must be specified before using new channels in calculations.
Current concentration analyses expect PTR/HORUS in pptv and VOCUS/ozone in ppbv.
Matching composition paths are configured separately in `composition_paths()`.

Files are merged chronologically per channel, with the last listed file/row
winning at an identical timestamp (including an explicitly missing value).
Different channel sets are supported; absent channels do not erase data from
other files. Coverage and duplicate counts are stored in each instrument's
report. File paths, sizes, modification times, SHA-256 fingerprints, preparation
settings and schema version are saved for reproducibility. Rebuild explicitly
when sources change; analysis does not silently refresh prepared files.
Rebuilds use a temporary file and replace the previous file only on success.

Preparation currently loads one instrument's configured files into memory at
a time. Reading prepared files is selective, but preparing much larger future
exports may warrant a chunked importer. Existing direct-file library calls remain
supported for comparisons or workflows that do not use prepared files.


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
points outside condition coverage become NaN. Interpolation also leaves NaN
strictly inside gaps longer than 30 minutes between valid condition measurements;
measured endpoints are retained and gaps of exactly 30 minutes are allowed.
This limit applies to measured temperature and pressure channels. Temperature and pressure are
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

## Time traces from multiple instruments

Overview and time-trace figures use an enlarged canvas and open maximized. Control them in
`config.py` with `MAXIMIZE_OVERVIEW_WINDOWS` and `MAXIMIZE_TIME_TRACE_WINDOWS`
(both True by default). Maximization happens before layout and saving, so saved
figures use the enlarged canvas too. Batch/unsupported backends use
`MAXIMIZED_FIGURE_SIZE` (19.2 × 10.8 inches by default). GUI backends use their
actual maximized canvas after processing the resize. Window controls remain
visible. The same time-trace flag also controls the trace figures produced by
`temp_dep`; its kinetics and k(T) figures and MD figures are unaffected. Set
either flag to False to retain ordinary dimensions for that plot type.

### PTR time resolution

Each run in `runs.py` has `low_res_trace_paths` and `high_res_trace_paths`.
Matching composition filenames are derived for each resolution. Preparation
stores low-resolution PTR under `ptr`/`composition`, and high-resolution PTR
under `ptr_high`/`composition_high` in the same run H5. Large high-resolution
CSV files are prepared in chunks, without smoothing or downsampling.

```python
high_res = True  # False selects low resolution
run = load_prepared_run(run_id, high_res=high_res)
temperature_runs = make_temperature_runs(run_ids, prepared=True, high_res=high_res)
```

All four plot scripts expose `high_res`; set it to True for high resolution
or False for low resolution. `run_plotting.py` currently plots
non-PTR instruments, so its panels are unaffected. The functions `temp_dep`,
`plot_time_traces`, and `md_plot` also accept `high_res=True/False` to switch
prepared inputs; omitting it keeps the supplied resolution. For direct CSV
workflows, select paths in `runs.py`/`make_temperature_runs` instead.
`md_plot` switches the composition table together with the trace source.
No instrument other than PTR is changed. Unconfigured/unprepared resolutions
raise an error; there is no silent fallback.

Rebuild existing H5 files once to add both resolutions:

```powershell
python scripts/prepare_runs.py m50_1 m25 m10 m50_2 --overwrite
```

For `m50_1`, the available high-resolution export is **V2**, while the existing
low-resolution export is **V4**. These are explicitly configured separately;
switching also changes the export version. Other high-resolution inputs use V1.

For `temp_dep`, select a measurement using `(instrument, exact H5 channel)`:

```python
precursor_compound=("ptr", "C10H15O4.NH3H+"),
product_compound=("ltof", "C20H30O6"),
reaction_partner=("ltof", "C10H15O4"),
```

Smooth only the temperature-dependence trace display with
`smoothing_time_traces=timedelta(minutes=10)` (import `timedelta` from
`datetime`), or use None for raw traces. This is a centered moving average
after unit conversion, applied independently to each signal on its original
sample times. Interval means, SEMs, kinetics and k(T) results always use raw
measurements. Smoothed trace filenames receive a suffix such as `_avg10min`;
kinetic/k(T) filenames are unchanged.

Set trace output directly in the `temp_dep(...)` call, alongside the kinetic
save paths:

```python
trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}.png",
```

None disables trace saving. Templates support the run ID and the precursor
and partner/product formulas without ion suffixes; smoothing appends `_avg10min`, etc.
For multiple runs, a missing `{run_id}` is appended automatically to prevent
overwriting. Trace filenames are no longer supplied by `runs.py`.
All three save paths expand `{precursor_compound}`, `{reaction_partner}`, and
`{product_compound}` in both folder and file names. Missing folders are created
automatically, so reactions can be grouped under a folder such as
`output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/`.

Use `make_temperature_runs(run_ids, prepared=True)` and
`concentration_unit="molecules/cm3"` to compare PTR and LToF. Only requested
channels are read; PTR mixing ratios are converted using run temperature and
pressure, while LToF concentrations remain unchanged. Exact names retain any
ion text already present in the export; no suffix is inferred or removed.
Tuple kinetic selections are automatically included in the interval means
and trace plots. `compounds` also accepts these tuples. Returned mean columns
include the instrument, e.g. `ltof:C20H30O6`. Original string selections and
the HORUS partners `"HO2"` and `"RO2"` remain supported. `run_temp_dep.py`
offers all four runs to each plot. Prepared runs missing a required signal or
finite samples in their window are skipped with the reason printed. PTR-only
plots can therefore include all four runs, while the current LToF export
limits plots requiring LToF to `m50_1` and `m25`. A channel absent from every
candidate run raises an error rather than silently skipping a typo. Combined
fits retain available members and report exclusions; if no run remains, the
analysis raises an error. Required HORUS partners and AP coloring are included
in the availability check.

In `scripts/run_time_traces.py`, select compounds separately for each instrument:

```python
run_id = "m50_1"  # Processed LTOF data is available for all four runs
compounds = {
    "ptr": ["C20H30O6", "C20H30O8"],
    "ltof": ["C20H30O6", "C20H30O8"],
}
```

Pass `trace_path=run["sources"]` to `plot_time_traces`. Neutral formulas match
by elemental counts against the part before the channel's dot; `CH3COOH` and
`C2H4O2` therefore match. Use standard capitalization. An exact full channel
such as `C20H30O6.NO3-` overrides selection. Channels without a neutral formula
cannot be resolved automatically.

`DEFAULT_ADDUCT_PREFERENCES` in `runs.py` sets the PTR preference. Override
per run with e.g. `"adduct_preferences": {"ltof": ["NO3-", "H+"]}`. A unique
match needs no preference; ambiguous matches raise an error listing candidates.
Preferences select existing channels, without correcting export labels. They are
read at plot time and do not require rebuilding prepared files.

Prepared channels supply their own units. Signals with the same unit share a
single axis. PTR pptv and LTOF molecules/cm3 overlay in one panel, with pptv on the left
and molecules/cm3 on the right. Selecting only
one unit gives a single y axis. More than two units use separate panels.
`stacked=True` gives each signal its own panel. Legends show the instrument and
actual channel. Original timestamps are retained. Use `concentration_unit="native"`
for counts, which require calibration before conversion to molecules/cm3.
Raw sources can supply `source_units={"ptr": "pptv", "ltof": "molecules/cm3"}`.
The LTOF source is the processed `06112026.csv` campaign export, read directly
with its shared UTC `time` column and metadata comments. Values are already
calibrated molecules/cm3 and are stored unchanged. Rebuild prepared runs after
changing this source: `python scripts/prepare_runs.py m50_1 m25 m10 m50_2 --overwrite`.
The original single-source list of full channel names remains supported.
