"""Plot selected CLOUD18 compound time traces."""

from datetime import datetime, timedelta

from cloud18_ap_analysis import plot_time_traces


# ---------------------------------------------------------------------------
# Temperature / Compounds
# ---------------------------------------------------------------------------

T = -10

compounds = [
    "C10H15O4.NH3H+",
    "C10H16O4.NH3H+",
    "C10H17O3.NH3H+",
    "C10H18O3.NH3H+",
]


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# -50°C data
if T == -50:
    trace_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-50C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1.csv"

    time_range = (
        datetime(2025, 9, 24, 0, 0, 0),
        datetime(2025, 9, 28, 0, 0, 0),
    )

# -25°C data
if T == -25:
    trace_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

    time_range = (
        datetime(2025, 9, 28, 0, 0, 0),
        datetime(2025, 10, 3, 0, 0, 0),
    )

# -10°C data
if T == -10:
    trace_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-10C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-10C_V1_low_res_part1.csv"

    time_range = (
        datetime(2025, 10, 8, 0, 0, 0),
        datetime(2025, 10, 12, 0, 0, 0),
    )

# -50_2°C data
if T == -50.2:
    trace_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-10C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-10C_V1_low_res_part1.csv"

    time_range = (
        datetime(2025, 10, 11, 0, 0, 0),
        datetime(2025, 10, 14, 12, 0, 0),
    )

fan_speed_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\fan_speed.csv"

stages_path = (
    r"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"
)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

plot_time_traces(
    trace_path=trace_path,
    compounds=compounds,
    time_range=time_range,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # Set None to disable low-fan-speed shading

    # Plot appearance
    stacked=False,
    ylabel="Intensity [ppt]",
    smoothing=None,
    yscale="log",

    # CLOUD stage labels
    stages="onlytype",

    # Output
    savepath=f"output/{T}C/time_traces.png",
    interactive=True,
)