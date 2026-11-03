"""Plot selected CLOUD18 compound time traces."""

# Trace paths accept a single filename or a list; append additional parts to the list.
from datetime import datetime, timedelta
from cloud18_ap_analysis.runs import RUNS, SHARED_PATHS

from cloud18_ap_analysis import plot_time_traces


# ---------------------------------------------------------------------------
# Temperature / Compounds
# ---------------------------------------------------------------------------

run_id = "m10"  # m50_1, m25, m10, m50_2

compounds = [
    # first generation RO2 oxidation products
    "C10H14O3.NH3H+",       # RO2 + RO2 -> R-CO-R' (ketone)
    "C10H15O3.NH3H+",       # RO2 + RO2 -> RO' (alkoxyl radical), but probably RO2 from RO2-RO2 reaction (we just do not know)
    "C10H15O5.NH3H+",
    "C10H16O5.NH3H+",
    "C10H15O4.NH3H+",       # AP + O3 -> RO2
    "C10H16O3.NH3H+",       # RO2 + RO2 -> ROH (alcohol)
    "C10H16O4.NH3H+",       # RO2 + HO2 -> ROOH (hydroperoxide)
    # "C10H15O6.NH3H+",       # RO2 + O2 -> RO2 (autoxidation product)
    "C10H17O3.NH3H+",       # AP + OH -> RO2
    "C10H17O4.NH3H+",
    "C10H17O6.NH3H+",       
    "C10H18O3.NH3H+",       # RO2 + HO2 -> ROOH (hydroperoxide)
    "C10H18O4.NH3H+",       # ROOOH
    "C10H18O8.NH3H+",
    # "C20H30O6.NH3H+",
    # "C20H30O8.NH3H+",
    # "C19H28O5.NH3H+",
    # "C19H28O7.NH3H+",
    # "C20H34O4.NH3H+",
    # "C20H34O6.NH3H+",
]

# ROOOH?



run = RUNS[run_id]
trace_path = run["trace_paths"]
time_range = run["time_range"] #(datetime(2025,10,1,21,0,0), datetime(2025,10,2,13,0,0) # Override here for a narrower plot window.
fan_speed_path = SHARED_PATHS["fan_speed_path"]
stages_path = SHARED_PATHS["stages_path"]

plot_time_traces(
    trace_path=trace_path,
    compounds=compounds,
    time_range=time_range,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # Set None to disable low-fan-speed shading
    important_time_ranges=run["important_time_ranges"],  # None disables intervals
    label_intervals=True,
    smoothing=timedelta(minutes=30),  # None disables smoothing

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
