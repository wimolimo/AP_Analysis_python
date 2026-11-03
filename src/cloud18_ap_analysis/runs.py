"""Shared measurement settings. Select m50_1, m25, m10, or m50_2.

Edit paths/default windows here and averaging intervals in intervals.py.
Plot appearance and optional time-window overrides stay in the run scripts.
"""
from datetime import datetime
from pathlib import Path
from copy import deepcopy
from .intervals import IMPORTANT_TIME_RANGES

DATA_ROOT = Path(r"C:\Users\c7441399\Documents\CLOUD18\AP_run")

RUNS = {
    "m50_1": {
        "label": '-50 °C (run 1)',
        "temperature_c": -50.0,
        "folder": '-50C',
        "trace_paths": [
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part3.csv",
        ],
        "temperature_path": DATA_ROOT / "-50C/results/temp_dew_frost_-50C.csv",
        "ozone_path": DATA_ROOT / "-50C/results/O3_-50C.csv",
        "time_range": (datetime(2025, 9, 24, 0, 0), datetime(2025, 9, 28, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m50_1"],
    },
    "m25": {
        "label": '-25 °C',
        "temperature_c": -25.0,
        "folder": '-25C',
        "trace_paths": [
            DATA_ROOT / "-25C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv",
        ],
        "temperature_path": DATA_ROOT / "-25C/results/temp_dew_frost_-25C.csv",
        "ozone_path": DATA_ROOT / "-25C/results/O3_-25C.csv",
        "time_range": (datetime(2025, 9, 30, 0, 0), datetime(2025, 10, 3, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m25"],
    },
    "m10": {
        "label": '-10 °C',
        "temperature_c": -10.0,
        "folder": '-10C',
        "trace_paths": [
            DATA_ROOT / "-10C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-10C_V1_low_res_part1.csv",
        ],
        "temperature_path": DATA_ROOT / "-10C/results/temp_dew_frost_-10C.csv",
        "ozone_path": DATA_ROOT / "-10C/results/O3_-10C.csv",
        "time_range": (datetime(2025, 10, 8, 0, 0), datetime(2025, 10, 12, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m10"],
    },
    "m50_2": {
        "label": '-50 °C (run 2)',
        "temperature_c": -50,
        "folder": '-50C_2',
        "trace_paths": [
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_low_res_part1.csv",
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_low_res_part2.csv",
        ],
        "temperature_path": DATA_ROOT / "-50C_2/results/temp_dew_frost_-50C_2.csv",
        "ozone_path": DATA_ROOT / "-50C_2/results/O3_-50C_2.csv",
        "time_range": (datetime(2025, 10, 11, 10, 0), datetime(2025, 10, 14, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m50_2"],
    },
}

SHARED_PATHS = {
    "horus_path": DATA_ROOT / "HORUS/HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt",
    "fan_speed_path": DATA_ROOT / "fan_speed.csv",
    "stages_path": DATA_ROOT / "CLOUD_log_stages/stages.txt",
    "ap_path": DATA_ROOT / "VOCUS_data/AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv",
    "smps_path": DATA_ROOT / "SMPS_data/CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt",
}

def make_temperature_runs(run_ids):
    """Build independent temp_dep inputs from selected RUNS entries."""
    result = []
    for run_id in run_ids:
        run = deepcopy(RUNS[run_id])
        result.append(dict(run_id=run_id, temperature=run["label"],
                           plot_temperature_c=run["temperature_c"],
                           trace_path=run["trace_paths"], temperature_path=run["temperature_path"],
                           time_range=run["time_range"], important_time_ranges=run["important_time_ranges"],
                           trace_savepath=f"output/temp_dep/temp_dep_{run_id}.png"))
    return result
