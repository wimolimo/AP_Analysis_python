"""Example CLOUD18 plotting script.

Choose a prepared run and plotting options below; configure inputs in runs.py.
"""

from datetime import datetime
from cloud18_ap_analysis import md_plot
from cloud18_ap_analysis.prepared import load_prepared_run
import cloud18_ap_analysis.config as cfg


run_id = "m25"
instrument = "ltof"  # "ptr" or "ltof": selects data and all settings below.

reference_families_ptr = ["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox"]
reference_families_ltof = ["C19H28Ox", "C20H30Ox", "C20H32Ox", "C20H34Ox"]
mass_range_ptr = (100, 350)
mass_range_ltof = (250, 650)
xlim_ptr = (80, 370)
xlim_ltof = (230, 670)
ylim_ptr = (0.02, 0.205)
ylim_ltof = (-0.02, 0.32)  # Set None for automatic limits.
# Threshold, marker scale and legend levels are configured separately in config.py:
# MIN_INTENSITY_PTR/LTOF, SIZE_ANCHOR_I_PTR/LTOF, SIZE_ANCHOR_S_PTR/LTOF,
# SIZE_EXP_PTR/LTOF and SIZE_LEGEND_LEVELS_PTR/LTOF.

#LTOF formulas without PI adduct, PTR formulas with NH3H+ adduct.

settings = {
    "ptr": dict(reference_families=reference_families_ptr, mass_range=mass_range_ptr,
                xlim=xlim_ptr, ylim=ylim_ptr),
    "ltof": dict(reference_families=reference_families_ltof, mass_range=mass_range_ltof,
                 xlim=xlim_ltof, ylim=ylim_ltof),
}
if instrument not in settings:
    raise ValueError("instrument must be 'ptr' or 'ltof'")
plot_settings = settings[instrument]
high_res = False  # True selects high-resolution PTR traces and matching composition.
run = load_prepared_run(run_id, high_res=high_res if instrument == "ptr" else False)

time_range = (
    datetime(2025, 10, 2, 0, 0, 0),
    datetime(2025, 10, 2, 1, 0, 0),
)

time_ranges = [(datetime(2025, 10, 2, 0, 0, 0), datetime(2025, 10, 2, 1, 0, 0)), 
               (datetime(2025, 10, 2, 2, 0, 0), datetime(2025, 10, 2, 4, 0, 0)),
               (datetime(2025, 10, 2, 7, 0, 0), datetime(2025, 10, 2, 9, 0, 0))]

# PTR needs its matching composition table; LTOF formulas are in its channel names.
composition_file = run["composition_path"] if instrument == "ptr" else None
trace_file = run["sources"][instrument]
threshold = getattr(cfg, f"MIN_INTENSITY_{instrument.upper()}")
threshold_tag = f"{threshold:g}{'ppt' if instrument == 'ptr' else 'molec_cm3'}"

for time_range in time_ranges:

    window_tag = f"{time_range[0]:%Y%m%dT%H%M}-{time_range[1]:%Y%m%dT%H%M}"
    md_plot(
        composition_file, trace_file, time_range=time_range,
        instrument=instrument, **plot_settings,
        savepath=f"output/{run_id}/md_plot_{instrument}_{threshold_tag}_{window_tag}.png",
        stream=True, color_by="On",
    )
