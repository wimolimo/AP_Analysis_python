"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime, timedelta
from cloud18_ap_analysis import load_data, load_stages, plot_data

time_range = (
    datetime(2025, 11, 24, 7, 0, 0),
    datetime(2025, 12, 2, 9, 0, 0),
)

file_name_temp = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-50C\results\temp_dew_frost_-50C.csv"
file_name_o3 = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-50C\results\O3_-50C.csv"
file_name_OH = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"
file_name_AP = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\VOCUS_data\AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv"
stages = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"

data_temp = load_data(file_name_temp, time_range)
data_o3 = load_data(file_name_o3, time_range)
data_OH = load_data(file_name_OH, time_range)
data_AP = load_data(file_name_AP, time_range)
data_SMPS = load_data(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\SMPS_data\CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt", time_range)
load_stages(stages)

plot_data(
    [data_temp, data_o3, data_OH, data_AP, data_SMPS],
    channels=[
        (["TE_Calib_3"], ["DewPoint_C"]),
        [],
        (["RO2_crosstalk_ppt", "HO2_ppt"], ["OH_ppt"]),
        (["AP"], ["IP"]),
        [],
    ],
    ylabels=[
        ("Temperature [°C]", "Dew Point [°C]"),
        "O₃ [ppb]",
        ("RO₂ and HO₂ [ppt]", "OH [ppt]"),
        ("AP [ppb]", "IP [ppb]"),
        "Diameter [nm]",
    ],
    xlims=time_range,
    savepath="output/-50C/O3_temp_OH_AP.png",
    smoothing=timedelta(minutes=10),
    interactive=True,
    stages="onlytype",
)

# time_range = (
#     datetime(2025, 11, 25, 0, 0, 0),
#     datetime(2025, 12, 1, 40, 0, 0),
# )

# # high resolution data
# # composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.txt"
# # trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv"

# # low resolution data
# composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.txt"
# trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

# md_plot(composition_file, trace_file, time_range=time_range, xlim=(100, 350), ylim=(0.02, 0.205), savepath="output/-25C/md_plot_2-4.png",
#         stream=True, color_by="On", reference_families=["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox", "C9H16Ox", "C9H14Ox", "C10H20Ox"])

# time_range = (
#     datetime(2025, 10, 2, 0, 0, 0),
#     datetime(2025, 10, 2, 1, 0, 0),
# )

# # high resolution data
# # composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.txt"
# # trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv"

# # low resolution data
# composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.txt"
# trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

# md_plot(composition_file, trace_file, time_range=time_range, xlim=(100, 350), ylim=(0.02, 0.205), savepath="output/-25C/md_plot_0-1.png",
#         stream=True, color_by="On", reference_families=["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox", "C9H16Ox", "C9H14Ox", "C10H20Ox"])

# time_range = (
#     datetime(2025, 10, 2, 7, 0, 0),
#     datetime(2025, 10, 2, 9, 0, 0),
# )

# # high resolution data
# # composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.txt"
# # trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv"

# # low resolution data
# composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.txt"
# trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

# md_plot(composition_file, trace_file, time_range=time_range, xlim=(100, 350), ylim=(0.02, 0.205), savepath="output/-25C/md_plot_7-9.png",
#         stream=True, color_by="On", reference_families=["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox", "C9H16Ox", "C9H14Ox", "C10H20Ox"])