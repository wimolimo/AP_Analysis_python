"""Plot selected CLOUD18 compound time traces."""

# Add source files in runs.py, then rebuild with scripts/prepare_runs.py.
from datetime import datetime, timedelta
from cloud18_ap_analysis.prepared import load_prepared_run
from cloud18_ap_analysis.runs import RUNS, DEFAULT_ADDUCT_PREFERENCES

from cloud18_ap_analysis import plot_time_traces


# ---------------------------------------------------------------------------
# Temperature / Compounds
# ---------------------------------------------------------------------------

run_id = "m25"  # m50_1, m25, m10, m50_2

compounds = {
    "ptr": [
        # first generation RO2 oxidation products
        # "C10H14O3",       # RO2 + RO2 -> R-CO-R' (ketone)
        # "C10H15O3",       # RO2 + RO2 -> RO' (alkoxyl radical), but probably RO2 from RO2-RO2 reaction (we just do not know)
        # "C10H15O5",
        # "C10H16O5",
        "C10H15O4",       # AP + O3 -> RO2
        # "C10H16O3",       # RO2 + RO2 -> ROH (alcohol)
        "C10H16O4",       # RO2 + HO2 -> ROOH (hydroperoxide)
        # "C10H15O6",       # RO2 + O2 -> RO2 (autoxidation product)
        "C10H17O3",       # AP + OH -> RO2
        # "C10H17O4",
        # "C10H17O6",       
        "C10H18O3",       # RO2 + HO2 -> ROOH (hydroperoxide)
        # "C10H18O4",       # ROOOH
        # "C10H18O8",
        # "C10H16O8",
        # "C10H16O10",
    ],

    # Processed LTOF concentrations (molecules/cm3) are available for all runs.
    "ltof": [
        # "C20H30O6",
        # "C20H30O8",
        # "C19H28O5",
        # "C19H28O7",
        # "C20H34O6",
        # "C20H34O8",
        # "C20H34O10",
        # "C20H32O7",

    ],
}

# ROOOH?



high_res = False  # True selects high-resolution PTR traces and matching composition.
run = load_prepared_run(run_id, high_res=high_res)
trace_path = run["sources"]
time_range = (datetime(2025,10,1,21,0,0), datetime(2025,10,2,13,0,0))  # Override here for a narrower plot window. run["time_range"]
fan_speed_path = run["sources"]["fan"]
stages_path = run["stages_path"]

important_time_ranges_md = [
    (datetime(2025, 10, 2, 0, 0, 0), datetime(2025, 10, 2, 1, 0, 0)),
    (datetime(2025, 10, 2, 2, 0, 0), datetime(2025, 10, 2, 4, 0, 0)),
    (datetime(2025, 10, 2, 7, 0, 0), datetime(2025, 10, 2, 9, 0, 0)),
]

plot_time_traces(
    trace_path=trace_path,
    compounds=compounds,
    adduct_preferences={**DEFAULT_ADDUCT_PREFERENCES,
                        **RUNS[run_id].get("adduct_preferences", {})},
    time_range=time_range,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # Set None to disable low-fan-speed shading
    important_time_ranges=important_time_ranges_md,  # run["important_time_ranges"],  # None disables intervals
    label_intervals=True,
    smoothing=timedelta(minutes=10),  # None disables smoothing

    # Plot appearance
    stacked=False,
    concentration_unit="native", #"molecules/cm3",  # "native" restores pptv
    temperature_path=run["temperature_path"],
    pressure_pa=97000.0,
    yscale="log",  # Automatic: log for molecules/cm3, linear for native

    # CLOUD stage labels
    stages="type",

    # Output
    savepath=f"output/{run_id}/time_traces.png",
    interactive=True,
)
