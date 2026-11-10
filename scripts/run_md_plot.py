"""Example CLOUD18 plotting script.

Choose a prepared run and plotting options below; configure inputs in runs.py.
"""

from datetime import datetime
from cloud18_ap_analysis import md_plot
from cloud18_ap_analysis.prepared import load_prepared_run

run_id = "m25"
high_res = False  # True selects high-resolution PTR traces and matching composition.
run = load_prepared_run(run_id, high_res=high_res)

time_range = (
    datetime(2025, 10, 2, 0, 0, 0),
    datetime(2025, 10, 2, 1, 0, 0),
)

time_ranges = [(datetime(2025, 10, 2, 0, 0, 0), datetime(2025, 10, 2, 1, 0, 0)), 
               (datetime(2025, 10, 2, 2, 0, 0), datetime(2025, 10, 2, 4, 0, 0)),
               (datetime(2025, 10, 2, 7, 0, 0), datetime(2025, 10, 2, 9, 0, 0))]

# Composition and trace data come from the same prepared run.
composition_file = run["composition_path"]
trace_file = run["trace_paths"]

for time_range in time_ranges:

    md_plot(composition_file, trace_file, time_range=time_range, xlim=(80, 370), ylim=(0.02, 0.205), savepath=f"output/{run_id}/md_plot_10ppq_{time_range[0].hour}h-{time_range[1].hour}h.png",
            stream=True, color_by="On", reference_families=["C10H15Ox", "C10H16Ox", "C10H17Ox", "C10H18Ox"],
            mass_range=(100, 350))
