"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime, timedelta
from cloud18_ap_analysis import load_data, load_stages, plot_data, md_plot

time_range = (
    datetime(2025, 10, 1, 22, 0, 0),
    datetime(2025, 10, 2, 13, 0, 0),
)

# file_name_temp = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\temp_dew_frost_-25C.csv"
# file_name_o3 = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\O3_-25C.csv"
# file_name_OH = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"
# file_name_AP = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\VOCUS_data\AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv"
# stages = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"

# data_temp = load_data(file_name_temp, time_range)
# data_o3 = load_data(file_name_o3, time_range)
# data_OH = load_data(file_name_OH, time_range)
# data_AP = load_data(file_name_AP, time_range)
# data_SMPS = load_data(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\SMPS_data\CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt", time_range)
# load_stages(stages)

# plot_data(
#     [data_temp, data_o3, data_OH, data_AP, data_SMPS],
#     channels=[
#         (["TE_Calib_3"], ["DewPoint_C"]),
#         [],
#         (["OH_ppt"], ["HO2_ppt"]),
#         (["AP"], ["IP"]),
#         [],
#     ],
#     ylabels=[
#         ("Temperature [°C]", "Dew Point [°C]"),
#         "O₃ [ppb]",
#         ("OH [ppt]", "HO₂ [ppt]"),
#         ("AP [ppb]", "IP [ppb]"),
#         "Diameter [nm]",
#     ],
#     xlims=time_range,
#     savepath="output/-25C/O3_temp_OH_AP.png",
#     smoothing=timedelta(minutes=1),
#     interactive=True,
#     stages="onlytype",
# )

time_range = (
    datetime(2025, 10, 2, 2, 0, 0),
    datetime(2025, 10, 2, 4, 0, 0),
)

composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.txt"
trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv"

md_plot(composition_file, trace_file, time_range=time_range, ylim=(-0.05, 0.45), savepath="md_plot.png")