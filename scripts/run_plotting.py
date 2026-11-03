"""Example CLOUD18 plotting script.

Edit the paths below for your local data.
"""

from datetime import datetime, timedelta
from cloud18_ap_analysis.runs import RUNS, SHARED_PATHS
from cloud18_ap_analysis import load_data, load_stages, plot_data
from cloud18_ap_analysis.units import (convert_dataset_concentrations,
    concentration_label, concentration_yscale)

# Concentration display: "native" keeps O3 in ppbv and radicals in pptv.
concentration_unit = "native" #"molecules/cm3"
pressure_pa = 97000.0  # Standard atmospheric pressure
concentration_scale = None  # Automatic; or override with "linear" / "log"

run_id = "m25"  # m50_1, m25, m10, m50_2
run = RUNS[run_id]
time_range =(datetime(2025,10,1, 21,0,0), datetime(2025,10,2, 13,0,0)) # run["time_range"] # Override here for a narrower plot window.
file_name_temp = run["temperature_path"]
file_name_o3 = run["ozone_path"]
file_name_OH = SHARED_PATHS["horus_path"]
file_name_AP = SHARED_PATHS["ap_path"]
file_name_fan = SHARED_PATHS["fan_speed_path"]
stages = SHARED_PATHS["stages_path"]

data_temp = load_data(file_name_temp, time_range)
data_o3 = load_data(file_name_o3, time_range)
data_OH = load_data(file_name_OH, time_range)
data_AP = load_data(file_name_AP, time_range)
data_SMPS = load_data(SHARED_PATHS["smps_path"], time_range)
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
    savepath=f"output/{run_id}/O3_temp_OH_AP.png",
    smoothing=timedelta(minutes=10),
    interactive=True,
    stages="type",
    fan_speed_data=data_fan,
    important_time_ranges=None,#run["important_time_ranges"],  # None disables intervals
    label_intervals=True,
    yscale=["linear", trace_scale, trace_scale, "linear", "log"],
)

