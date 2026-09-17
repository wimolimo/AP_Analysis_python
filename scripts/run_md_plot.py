"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime
from cloud18_ap_analysis import md_plot

time_range = (
    datetime(2025, 10, 2, 2, 0, 0),
    datetime(2025, 10, 2, 4, 0, 0),
)

# high resolution data
# composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.txt"
# trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv"

# low resolution data
composition_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3compositions_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.txt"
trace_file = r"C:\Users\c7441399\Documents\CLOUD18\AP_run\-25C\results\ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv"

md_plot(composition_file, trace_file, time_range=time_range, xlim=(80, 370), ylim=(0.02, 0.205), savepath="output/-25C/md_plot_2-4.png",
        stream=True, color_by="On", reference_families=["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox", "C9H16Ox", "C9H14Ox", "C10H20Ox"],
        mass_range=(100, 350))

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