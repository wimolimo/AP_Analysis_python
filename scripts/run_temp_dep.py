"""Run the multi-temperature kinetic comparison."""
from datetime import timedelta
from cloud18_ap_analysis import temp_dep
from cloud18_ap_analysis.runs import make_temperature_runs

# tau_H2SO4 = 480 s at T = 291 K, p = 970 mbar (standard CLOUD pressure)

# Each plot uses all runs with data for its selected signals.
run_ids = ["m50_1", "m25", "m10", "m50_2"]
high_res = True  # Use dense PTR samples for interval means and uncertainties.
temperature_runs = make_temperature_runs(run_ids, prepared=True, high_res=high_res)
# Pool interval points into one fit; set None for a separate fit per measurement.
combine_runs = {"-50 °C (combined)": ["m50_1", "m50_2"]}
concentration_unit = "molecules/cm3"
pressure_pa = 97000.0
trace_yscale = None
smoothing_time_traces = timedelta(minutes=10)  # None keeps raw traces; fits always use raw data.
# Each temperature run supplies its own instrument/stage references.
horus_path = temperature_runs[0]["horus_path"]
fan_speed_path = temperature_runs[0]["fan_speed_path"]
stages_path = None
ap_path = temperature_runs[0]["ap_path"]



############################## C10H15O4 + HO2 -> C10H16O4 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("horus", "HO2_ppt")]

# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,              # None for seperate fits
#     compounds=compounds,
#     horus_path=horus_path,
#     # (instrument, exact channel name in H5); no ion suffix is inferred.
#     precursor_compound=("ptr", "C10H15O4.NH3H+"),    # RO2
#     product_compound=("ptr", "C10H16O4.NH3H+"),
#     reaction_partner="HO2",     # ("ptr", "C10H15O4.NH3H+"),
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     smoothing_time_traces=smoothing_time_traces,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)

################################### C10H17O3 + HO2 -> C10H18O3 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("horus", "HO2_ppt")]


# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,
#     compounds=compounds,
#     horus_path=horus_path,
#     precursor_compound=("ptr", "C10H17O3.NH3H+"),    # RO2
#     product_compound=("ptr", "C10H18O3.NH3H+"),      # RO2 + HO2 -> ROOH
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
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)

############################## C10H15O4 + RO2 -> C20H30O6 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("horus", "RO2_crosstalk_ppt")] + [("ltof", "C20H30O6")]

trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
    temperature_runs=temperature_runs,
    combine_runs=combine_runs,              # None for seperate fits
    compounds=compounds,
    horus_path=horus_path,
    # (instrument, exact channel name in H5); no ion suffix is inferred.
    precursor_compound=("ptr", "C10H15O4.NH3H+"),    # RO2
    product_compound=("ltof", "C20H30O6"),
    reaction_partner="RO2",     # ("ptr", "C10H15O4.NH3H+"),
    overlay=True,
    show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
    concentration_unit=concentration_unit,
    pressure_pa=pressure_pa,
    yscale=trace_yscale,
    smoothing_time_traces=smoothing_time_traces,
    color_by_ap=True,  # Set False to color points by temperature instead
    ap_path=ap_path,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
    stages="type",
    label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
    fit_intercept=True,  # False fixes the kinetic fit intercept at zero
    trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
    kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
    k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
)

print(regression_results)

################################ C10H15O4 + C10H15O4 -> C20H30O6 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("ltof", "C20H30O6")]

trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
    temperature_runs=temperature_runs,
    combine_runs=combine_runs,              # None for seperate fits
    compounds=compounds,
    horus_path=horus_path,
    # (instrument, exact channel name in H5); no ion suffix is inferred.
    precursor_compound=("ptr", "C10H15O4.NH3H+"),    # RO2
    product_compound=("ltof", "C20H30O6"),
    reaction_partner=("ptr", "C10H15O4.NH3H+"),
    overlay=True,
    show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
    concentration_unit=concentration_unit,
    pressure_pa=pressure_pa,
    yscale=trace_yscale,
    smoothing_time_traces=smoothing_time_traces,
    color_by_ap=True,  # Set False to color points by temperature instead
    ap_path=ap_path,
    stages_path=stages_path,
    fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
    stages="type",
    label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
    fit_intercept=True,  # False fixes the kinetic fit intercept at zero
    trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
    kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
    k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
)

print(regression_results)

################################ C10H15O4 + C10H15O6 -> C20H30O8 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("ltof", "C20H30O6"), ("ltof", "C10H15O6")]

# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,              # None for seperate fits
#     compounds=compounds,
#     horus_path=horus_path,
#     # (instrument, exact channel name in H5); no ion suffix is inferred.
#     precursor_compound=("ptr", "C10H15O4.NH3H+"),    # RO2
#     product_compound=("ltof", "C20H30O8"),
#     reaction_partner=("ltof", "C10H15O6"),
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     smoothing_time_traces=smoothing_time_traces,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)

################################ C10H15O6 + RO2 -> C20H30O8 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("ltof", "C20H30O6"), ("ltof", "C10H15O6")]

# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,              # None for seperate fits
#     compounds=compounds,
#     horus_path=horus_path,
#     # (instrument, exact channel name in H5); no ion suffix is inferred.
#     precursor_compound=("ltof", "C10H15O6"),    # RO2
#     product_compound=("ltof", "C20H30O8"),
#     reaction_partner="RO2",
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     smoothing_time_traces=smoothing_time_traces,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)

################################ C10H15O4 + C10H15O8 -> C20H30O10 ##################################

compounds = [("ptr", name) for name in ['C10H15O4.NH3H+', 'C10H16O4.NH3H+',
    'C10H17O3.NH3H+', 'C10H18O3.NH3H+',]] + [("ltof", "C20H30O10"), ("ltof", "C10H15O8")]

# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,              # None for seperate fits
#     compounds=compounds,
#     horus_path=horus_path,
#     # (instrument, exact channel name in H5); no ion suffix is inferred.
#     precursor_compound=("ptr", "C10H15O4.NH3H+"),    # RO2
#     product_compound=("ltof", "C20H30O10"),
#     reaction_partner=("ltof", "C10H15O8"),
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     smoothing_time_traces=smoothing_time_traces,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)

################################### C10H15O6 + C10H15O6 -> C20H30O10 ##################################

compounds = [("ltof", "C20H30O6"), ("ltof", "C10H15O6")]

# trace_figures, means, regression_results, kinetic_fig, kinetic_ax, k_temperature_fig, k_temperature_ax = temp_dep(
#     temperature_runs=temperature_runs,
#     combine_runs=combine_runs,              # None for seperate fits
#     compounds=compounds,
#     horus_path=horus_path,
#     # (instrument, exact channel name in H5); no ion suffix is inferred.
#     precursor_compound=("ltof", "C10H15O6"),    # RO2
#     product_compound=("ltof", "C20H30O10"),
#     reaction_partner=("ltof", "C10H15O6"),
#     overlay=True,
#     show=["kinetics", "k(T)"],  # "traces", "kinetics", "k(T)", a list of these, or True for all
#     concentration_unit=concentration_unit,
#     pressure_pa=pressure_pa,
#     yscale=trace_yscale,
#     smoothing_time_traces=smoothing_time_traces,
#     color_by_ap=True,  # Set False to color points by temperature instead
#     ap_path=ap_path,
#     stages_path=stages_path,
#     fan_speed_path=fan_speed_path,  # None disables the low-fan-speed strip
#     stages="type",
#     label_intervals=False,  # Matching t_i tags on shaded intervals and kinetic points
#     fit_intercept=True,  # False fixes the kinetic fit intercept at zero
#     trace_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/trace_{run_id}_{precursor_compound}_{reaction_partner}_no_tags.png",
#     kinetic_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kinetic_{precursor_compound}_{reaction_partner}_no_tags.png",
#     k_temperature_savepath="output/temp_dep/{product_compound}_{precursor_compound}_{reaction_partner}/kT_{precursor_compound}_{reaction_partner}.png",
# )

# print(regression_results)


