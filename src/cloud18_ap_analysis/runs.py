"""Shared measurement settings. Select m50_1, m25, m10, or m50_2.

Edit paths/default windows here and averaging intervals in intervals.py.
Plot appearance and optional time-window overrides stay in the run scripts.
"""
from datetime import datetime
from pathlib import Path
from copy import deepcopy
from .intervals import IMPORTANT_TIME_RANGES

DATA_ROOT = Path(r"C:\Users\c7441399\Documents\CLOUD18\AP_run")
LTOF_PROCESSED_PATH = DATA_ROOT / "LTOF_UFRA/LToF_UFRA_CLOUD18_ParameterizationRun_220925to131025_06112026.csv"

# Selection preferences apply only when a neutral compound has multiple channels.
# Override per run below when the instrument's reagent ion changes.
DEFAULT_ADDUCT_PREFERENCES = {"ptr": ["NH3H+", "H+", ""]}

RUNS = {
    "m50_1": {
        "label": '-50 °C (run 1)',
        "temperature_c": -50.0,
        "folder": '-50C',
        # Available high-resolution export is V2; low-resolution remains V4.
        "high_res_trace_paths": [
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V2_high_res_part1.csv",
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V2_high_res_part2.csv",
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V2_high_res_part3.csv",
        ],
        "low_res_trace_paths": [
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part1.csv",
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part2.csv",
            DATA_ROOT / "-50C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V4_low_res_part3.csv",
        ],
        "temperature_path": DATA_ROOT / "-50C/results/temp_dew_frost_-50C.csv",
        # Empty means select a unique match; otherwise require an explicit channel.
        # Set e.g. ["NO3-", "H+"] when that ion preference is appropriate.
        "adduct_preferences": {"ltof": []},
        "ozone_path": DATA_ROOT / "-50C/results/O3_-50C.csv",
        "time_range": (datetime(2025, 9, 24, 0, 0), datetime(2025, 9, 28, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m50_1"],
    },
    "m25": {
        "label": '-25 °C',
        "temperature_c": -25.0,
        "folder": '-25C',
        "high_res_trace_paths": [
            DATA_ROOT / "-25C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1_high_res.csv",
        ],
        "low_res_trace_paths": [
            DATA_ROOT / "-25C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-25C_V1.csv",
        ],
        "temperature_path": DATA_ROOT / "-25C/results/temp_dew_frost_-25C.csv",
        "ozone_path": DATA_ROOT / "-25C/results/O3_-25C.csv",
        "time_range": (datetime(2025, 9, 29, 19, 0), datetime(2025, 10, 2, 13, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m25"],
    },
    "m10": {
        "label": '-10 °C',
        "temperature_c": -10.0,
        "folder": '-10C',
        "high_res_trace_paths": [
            DATA_ROOT / "-10C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-10C_V1_high_res.csv",
        ],
        "low_res_trace_paths": [
            DATA_ROOT / "-10C/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-10C_V1_low_res_part1.csv",
        ],
        "temperature_path": DATA_ROOT / "-10C/results/temp_dew_frost_-10C.csv",
        "ozone_path": DATA_ROOT / "-10C/results/O3_-10C.csv",
        "time_range": (datetime(2025, 10, 8, 8, 0), datetime(2025, 10, 11, 12, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m10"],
    },
    "m50_2": {
        "label": '-50 °C (run 2)',
        "temperature_c": -50,
        "folder": '-50C_2',
        "high_res_trace_paths": [
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_high_res_part1.csv",
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_high_res_part2.csv",
        ],
        "low_res_trace_paths": [
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_low_res_part1.csv",
            DATA_ROOT / "-50C_2/results/ptr3traces_UIBK_oVOCs_AP_run_CLOUD18_inletLossCorr_T-50C_V1_low_res_part2.csv",
        ],
        "temperature_path": DATA_ROOT / "-50C_2/results/temp_dew_frost_-50C_2.csv",
        "ozone_path": DATA_ROOT / "-50C_2/results/O3_-50C_2.csv",
        "time_range": (datetime(2025, 10, 11, 10, 0), datetime(2025, 10, 14, 0, 0)),
        "important_time_ranges": IMPORTANT_TIME_RANGES["m50_2"],
    },
}

# Backward-compatible default for direct CSV callers: low resolution.
for _run in RUNS.values():
    _run["trace_paths"] = _run["low_res_trace_paths"]

SHARED_PATHS = {
    "horus_path": DATA_ROOT / "HORUS/HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V2.txt",
    "fan_speed_path": DATA_ROOT / "fan_speed.csv",
    "stages_path": DATA_ROOT / "CLOUD_log_stages/stages.txt",
    "ap_path": DATA_ROOT / "VOCUS_data/AP_IP_tracer_calibrated_VOCUS_BGCorrect.csv",
    "smps_path": DATA_ROOT / "SMPS_data/CLOUD18_SMPS_merged_20nm_final_remove_clogging.txt",
}

def source_specs(run_id):
    """All inputs and native units used to prepare a run; add instruments here."""
    from .prepared import SourceSpec
    run = RUNS[run_id]
    sources = {
        "ptr": SourceSpec(run["low_res_trace_paths"], "trace", default_unit="pptv"),
        "temperature": SourceSpec([run["temperature_path"]], "paired",
                                  default_unit="degC", padding_minutes=30),
        "ozone": SourceSpec([run["ozone_path"]], "paired", default_unit="ppbv"),
        "horus": SourceSpec([SHARED_PATHS["horus_path"]], "custom",
                            units={"HO2_ppt": "pptv", "OH_ppt": "pptv", "RO2_crosstalk_ppt": "pptv",
                                   "OH_det_lim_ppt": "pptv", "OH_close_under_LOD": "flag"}),
        "vocus": SourceSpec([SHARED_PATHS["ap_path"]], "ap", default_unit="ppbv"),
        "fan": SourceSpec([SHARED_PATHS["fan_speed_path"]], "paired",
                          default_unit="%", padding_minutes=30),
        "smps": SourceSpec([SHARED_PATHS["smps_path"]], "smps", default_unit="cm^-3"),
    }
    sources["ltof"] = SourceSpec([LTOF_PROCESSED_PATH], "ltof_processed",
                                 default_unit="molecules/cm3")
    if run.get("high_res_trace_paths"):
        sources["ptr_high"] = SourceSpec(run["high_res_trace_paths"], "trace_stream", default_unit="pptv")
    return sources


def composition_paths(run_id, high_res=False):
    """Matching PTR composition exports; override here if naming differs."""
    return [p.with_name(p.name.replace("ptr3traces", "ptr3compositions")).with_suffix(".txt")
            for p in RUNS[run_id]["high_res_trace_paths" if high_res else "low_res_trace_paths"]]


def make_temperature_runs(run_ids, prepared=False, high_res=False):
    """Build independent temp_dep inputs from selected RUNS entries."""
    result = []
    for run_id in run_ids:
        if prepared:
            from .prepared import load_prepared_run
            run = load_prepared_run(run_id, high_res=high_res)
        else:
            run = deepcopy(RUNS[run_id])
            run["trace_paths"] = run["high_res_trace_paths" if high_res else "low_res_trace_paths"]
            if not run["trace_paths"]:
                raise ValueError(f"No requested PTR resolution configured for {run_id}")
        result.append(dict(run_id=run_id, temperature=run["label"],
                           plot_temperature_c=run["temperature_c"],
                           trace_path=run["trace_paths"], temperature_path=run["temperature_path"],
                           time_range=run["time_range"], important_time_ranges=run["important_time_ranges"]))
        if prepared:
            result[-1].update(horus_path=run["sources"]["horus"], ap_path=run["sources"]["vocus"],
                              fan_speed_path=run["sources"]["fan"], stages_path=run["stages_path"],
                              sources=run["sources"])
    return result
