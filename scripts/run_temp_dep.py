"""Run the multi-temperature kinetic comparison."""
from cloud18_ap_analysis import temp_dep
from cloud18_ap_analysis.runs import make_temperature_runs, SHARED_PATHS

# tau_H2SO4 = 480 s at T = 291 K, p = 970 mbar (standard CLOUD pressure)

run_ids = ["m50_1", "m25", "m10", "m50_2"]
temperature_runs = make_temperature_runs(run_ids)
# Pool interval points into one fit; set None for a separate fit per measurement.
combine_runs = {"-50 °C (combined)": ["m50_1", "m50_2"]}
concentration_unit = "molecules/cm3"
pressure_pa = 97000.0
trace_yscale = None
horus_path = SHARED_PATHS["horus_path"]
fan_speed_path = SHARED_PATHS["fan_speed_path"]
stages_path = SHARED_PATHS["stages_path"]
ap_path = SHARED_PATHS["ap_path"]
compounds = ['C10H14O4.NH3H+', 'C10H15O4.NH3H+', 'C10H16O4.NH3H+', 'C10H17O3.NH3H+', 'C10H18O3.NH3H+']

# O3 oxidation products
trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
    temperature_runs=temperature_runs,
    combine_runs=combine_runs,              # None for seperate fits
    compounds=compounds,
    horus_path=horus_path,
    precursor_compound="C10H15O4.NH3H+",    # RO2
    product_compound="C10H16O4.NH3H+",
    reaction_partner="HO2",
    overlay=True,
    show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
    concentration_unit=concentration_unit,
    pressure_pa=pressure_pa,
    yscale=trace_yscale,
    color_by_ap=True,  # Set False to color points by temperature instead
    ap_path=ap_path,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
    stages="onlytype",
    label_intervals=True,  # Matching t_i tags on shaded intervals and kinetic points
    fit_intercept=True,  # False fixes the kinetic fit intercept at zero
    kinetic_savepath="output/temp_dep/temp_dep_kinetic_{precursor_compound}_HO2.png",
    k_temperature_savepath="output/temp_dep/temp_dep_kT_{precursor_compound}_HO2.png",
)

print(regression_results)

# OH oxidation products
# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,
#     compounds=compounds,
#     horus_path=horus_path,
#     precursor_compound="C10H17O3.NH3H+",    # RO2
#     product_compound="C10H18O3.NH3H+",      # RO2 + HO2 -> ROOH
#     reaction_partner="HO2",
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="onlytype",
#     label_intervals=True,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     kinetic_savepath="output/temp_dep/temp_dep_kinetic_C10H18O3_HO2.png",
#     k_temperature_savepath="output/temp_dep/temp_dep_kT_C10H18O3_HO2.png",
# )

# print(regression_results)