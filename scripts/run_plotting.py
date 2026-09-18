"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime, timedelta
from cloud18_ap_analysis import load_data, load_stages, plot_data

# -50C data
T = -50
time_range = (
    datetime(2025, 9, 24, 0, 0, 0),
    datetime(2025, 9, 28, 0, 0, 0),
)

file_name_temp = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\{T}C\results\temp_dew_frost_{T}C.csv"
file_name_o3 = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\{T}C\results\O3_{T}C.csv"
file_name_OH = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"
file_name_AP = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\VOCUS_data\AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv"
file_name_fan = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\fan_speed.csv"
stages = rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"

data_temp = load_data(file_name_temp, time_range)
data_o3 = load_data(file_name_o3, time_range)
data_OH = load_data(file_name_OH, time_range)
data_AP = load_data(file_name_AP, time_range)
data_SMPS = load_data(rf"C:\Users\c7441399\Documents\CLOUD18\AP_run\SMPS_data\CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt", time_range)
data_fan = load_data(file_name_fan, time_range)
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
    savepath=f"output/{T}C/O3_temp_OH_AP.png",
    smoothing=timedelta(minutes=10),
    interactive=True,
    stages="onlytype",
    fan_speed_data=data_fan,
)

# -25C data
# time_range = (
#     datetime(2025, 9, 28, 0, 0, 0),
#     datetime(2025, 10, 3, 0, 0, 0),
# )

# file_name_temp = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\temp_dew_frost_-25C.csv"
# file_name_o3 = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\O3_-25C.csv"
# file_name_OH = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt"
# file_name_AP = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\VOCUS_data\AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv"
# file_name_fan = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\fan_speed.csv"
# stages = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\CLOUD_log_stages\stages.txt"

# data_temp = load_data(file_name_temp, time_range)
# data_o3 = load_data(file_name_o3, time_range)
# data_OH = load_data(file_name_OH, time_range)
# data_AP = load_data(file_name_AP, time_range)
# data_SMPS = load_data(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\SMPS_data\CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt", time_range)
# data_fan = load_data(file_name_fan, time_range)
# load_stages(stages)

# plot_data(
#     [data_temp, data_o3, data_OH, data_AP, data_SMPS],
#     channels=[
#         (["TE_Calib_3"], ["DewPoint_C"]),
#         [],
#         (["RO2_crosstalk_ppt", "HO2_ppt"], ["OH_ppt"]),
#         (["AP"], ["IP"]),
#         [],
#     ],
#     ylabels=[
#         ("Temperature [°C]", "Dew Point [°C]"),
#         "O₃ [ppb]",
#         ("RO₂ and HO₂ [ppt]", "OH [ppt]"),
#         ("AP [ppb]", "IP [ppb]"),
#         "Diameter [nm]",
#     ],
#     xlims=time_range,
#     savepath="output/-25C/O3_temp_OH_AP.png",
#     smoothing=timedelta(minutes=10),
#     interactive=True,
#     stages="onlytype",
#     fan_speed_data=data_fan,
# )
