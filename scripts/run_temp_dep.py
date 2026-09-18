from datetime import datetime
from cloud18_ap_analysis import temp_dep



time_range = (
        datetime(2025, 9, 24, 0, 0, 0),
        datetime(2025, 9, 28, 0, 0, 0),
    )
important_time_ranges = [
    (   datetime(2025, 9, 24, 18, 0, 0), datetime(2025, 9, 24, 21, 0, 0),),
    (   datetime(2025, 9, 25, 10, 0, 0), datetime(2025, 9, 25, 20, 0, 0),),
    (   datetime(2025, 9, 26, 1, 30, 0), datetime(2025, 9, 26, 6, 0, 0),),
    (   datetime(2025, 9, 26, 12, 0, 0), datetime(2025, 9, 26, 16, 0, 0),),
    (   datetime(2025, 9, 26, 19, 0, 0), datetime(2025, 9, 27, 0, 0, 0),),
    (   datetime(2025, 9, 27, 3, 30, 0), datetime(2025, 9, 27, 4, 30, 0),),
]

trace_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-50C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part3.csv"
HORUS_path = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"
compounds = [
        "C10H15O4.NH3H+",
        "C10H16O4.NH3H+",
        "C10H17O3.NH3H+",
        "C10H18O3.NH3H+",
    ]

temp_dep(
    trace_path=trace_path,
    time_range=time_range,
    compounds=compounds,
    overlay=True,
    important_time_ranges=important_time_ranges,
    savepath="output/temp_dep/temp_dep_plot.png",
    x_compound="C10H16O4.NH3H+",
    y_compound="C10H15O4.NH3H+",
    correlation_savepath="output/temp_dep/temp_dep_correlation.png",
)