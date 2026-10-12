"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime, timedelta
from cloud18_ap_analysis import load_data, load_stages, plot_data
from cloud18_ap_analysis.units import (convert_dataset_concentrations,
    concentration_label, concentration_yscale)

# Concentration display: "native" keeps O3 in ppbv and radicals in pptv.
concentration_unit = "molecules/cm3"
pressure_pa = 101325.0  # Standard atmospheric pressure
concentration_scale = None  # Automatic; or override with "linear" / "log"

# -50C data
T = -10
# time_range = (
#     datetime(2025, 9, 24, 0, 0, 0),
#     datetime(2025, 9, 28, 0, 0, 0),
# )
# time_range = (
#     datetime(2025, 9, 28, 0, 0, 0),
#     datetime(2025, 10, 3, 0, 0, 0),
# )
time_range = (
    datetime(2025, 10, 8, 0, 0, 0),
    datetime(2025, 10, 12, 0, 0, 0),
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

# Keep the raw datasets intact. Convert before plot smoothing/averaging.
# Full temperature coverage avoids losing interpolation endpoints at time_range.
temperature_c = (load_data(file_name_temp)["TE_Calib_3"]
                 if concentration_unit == "molecules/cm3" else None)
plot_o3 = convert_dataset_concentrations(
    data_o3, {"o3": "ppbv", "o3_psi": "ppbv"}, concentration_unit,
    temperature_c=temperature_c, pressure_pa=pressure_pa,
)
plot_radicals = convert_dataset_concentrations(
    data_OH, {"RO2_crosstalk_ppt": "pptv", "HO2_ppt": "pptv", "OH_ppt": "pptv"},
    concentration_unit, temperature_c=temperature_c, pressure_pa=pressure_pa,
)
plot_ap = convert_dataset_concentrations(data_AP, {"AP": "ppbv", "IP": "ppbv"}, "pptv")
radical_unit = concentration_label(concentration_unit)
ozone_unit = concentration_label(concentration_unit, native="ppbv")
trace_scale = concentration_yscale(concentration_unit, concentration_scale)

plot_data(
    [data_temp, plot_o3, plot_radicals, plot_ap, data_SMPS],
    channels=[
        (["TE_Calib_3"], ["DewPoint_C"]),
        [],
        (["RO2_crosstalk_ppt", "HO2_ppt"], ["OH_ppt"]),
        (["AP"], ["IP"]),
        [],
    ],
    ylabels=[
        ("Temperature [°C]", "Dew Point [°C]"),
        rf"O$_3$ [{ozone_unit}]",
        (rf"RO$_2$ and HO$_2$ [{radical_unit}]", f"OH [{radical_unit}]"),
        ("AP [pptv]", "IP [pptv]"),
        "Diameter [nm]",
    ],
    xlims=time_range,
    savepath=f"output/{T}C/O3_temp_OH_AP.png",
    smoothing=timedelta(minutes=10),
    interactive=True,
    stages="onlytype",
    fan_speed_data=data_fan,
    yscale=["linear", trace_scale, trace_scale, "linear", "log"],
)

