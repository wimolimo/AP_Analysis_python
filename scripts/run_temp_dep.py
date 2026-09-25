"""Run the multi-temperature kinetic comparison."""

from datetime import datetime

from cloud18_ap_analysis import temp_dep


# ---------------------------------------------------------------------------
# -50 °C
# ---------------------------------------------------------------------------

time_range_m50C = (
    datetime(2025, 9, 24, 0, 0, 0),
    datetime(2025, 9, 28, 0, 0, 0),
)

important_time_ranges_m50C = [
    (datetime(2025, 9, 25, 10, 0, 0), datetime(2025, 9, 25, 20, 0, 0)),
    (datetime(2025, 9, 26, 1, 30, 0), datetime(2025, 9, 26, 6, 0, 0)),
    (datetime(2025, 9, 26, 12, 0, 0), datetime(2025, 9, 26, 16, 0, 0)),
    (datetime(2025, 9, 26, 19, 0, 0), datetime(2025, 9, 27, 0, 0, 0)),
    (datetime(2025, 9, 27, 4, 0, 0), datetime(2025, 9, 27, 4, 30, 0)),
]

trace_path_m50C = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-50C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part3.csv"


# ---------------------------------------------------------------------------
# -25 °C
# ---------------------------------------------------------------------------
# Fill these in once you have identified the -25 °C AP period and intervals.

time_range_m25C = (
    datetime(2025, 9, 30, 0, 0, 0),
    datetime(2025, 10, 3, 0, 0, 0),
)

important_time_ranges_m25C = [
    (datetime(2025, 9, 30, 2, 40, 0), datetime(2025, 9, 30, 3, 40, 0)),
    (datetime(2025, 9, 30, 8, 10, 0), datetime(2025, 9, 30, 9, 10, 0)),
    (datetime(2025, 9, 30, 11, 40, 0), datetime(2025, 9, 30, 12, 30, 0)),
    (datetime(2025, 9, 30, 14, 30, 0), datetime(2025, 9, 30, 15, 20, 0)),
    (datetime(2025, 9, 30, 16, 10, 0), datetime(2025, 9, 30, 17, 0, 0)),
    (datetime(2025, 9, 30, 18, 30, 0), datetime(2025, 9, 30, 21, 0, 0)),
    (datetime(2025, 9, 30, 22, 0, 0), datetime(2025, 9, 30, 23, 30, 0)),
    (datetime(2025, 10, 1, 2, 10, 0), datetime(2025, 10, 1, 3, 15, 0)),
    (datetime(2025, 10, 1, 6, 0, 0), datetime(2025, 10, 1, 7, 30, 0)),
    (datetime(2025, 10, 1, 15, 20, 0), datetime(2025, 10, 1, 16, 20, 0)),
]

trace_path_m25C = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

# ---------------------------------------------------------------------------
# -10 °C
# ---------------------------------------------------------------------------
# Fill these in once you have identified the -10 °C AP period and intervals.

# time_range_m10C = (
#     datetime(2025, 10, 12, 0, 0, 0),
#     datetime(2025, 10, 15, 0, 0, 0),
# )

# important_time_ranges_m25C = [
#     (datetime(2025, 10, 1, 0, 30, 0), datetime(2025, 10, 1, 3, 0, 0)),
#     (datetime(2025, 10, 1, 6, 0, 0), datetime(2025, 10, 1, 7, 30, 0)),
# ]

# trace_path_m25C = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"


# ---------------------------------------------------------------------------
# Shared inputs
# ---------------------------------------------------------------------------

horus_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"

stages_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"

compounds = [
    "C10H15O4.NH3H+",
    "C10H16O4.NH3H+",
    "C10H17O3.NH3H+",
    "C10H18O3.NH3H+",
]


temperature_runs = [
    {
        "temperature": "-50 °C",
        "trace_path": trace_path_m50C,
        "time_range": time_range_m50C,
        "important_time_ranges": important_time_ranges_m50C,
        "trace_savepath": "output/temp_dep/temp_dep_m50C.png",
    },

    {
        "temperature": "-25 °C",
        "trace_path": trace_path_m25C,
        "time_range": time_range_m25C,
        "important_time_ranges": important_time_ranges_m25C,
        "trace_savepath": "output/temp_dep/temp_dep_m25C.png",
    },
]


trace_figures, means, regression_results, kinetic_fig, kinetic_ax = temp_dep(
    temperature_runs=temperature_runs,
    compounds=compounds,
    horus_path=horus_path,
    precursor_compound="C10H15O4.NH3H+",
    product_compound="C10H16O4.NH3H+",
    reaction_partner="HO2",
    overlay=True,
    stages_path=stages_path,
    stages="onlytype",
    label_intervals=True,  # Matching t_i tags on shaded intervals and kinetic points
    kinetic_savepath="output/temp_dep/temp_dep_kinetic_all_temperatures.png",
)

print(regression_results)
