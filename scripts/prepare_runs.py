"""Build run files: python scripts/prepare_runs.py [m25 ...] [--overwrite]."""
import argparse

from cloud18_ap_analysis.prepared import prepare_run
from cloud18_ap_analysis.runs import RUNS, SHARED_PATHS, source_specs, composition_paths


# Play-button settings. Set RUN_IDS = list(RUNS) to prepare all runs.
RUN_IDS = ["m25", "m10", "m50_2"]
OVERWRITE = True  # Rebuild the prepared run to include newly added instruments.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_ids", nargs="*", help="Run IDs; omit to use RUN_IDS above")
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE,
                        help="Rebuild existing files")
    args = parser.parse_args()
    unknown = set(args.run_ids).difference(RUNS)
    if unknown:
        parser.error(f"Unknown run IDs: {sorted(unknown)}")
    for run_id in args.run_ids or RUN_IDS:
        prepare_run(run_id, RUNS[run_id], source_specs(run_id), SHARED_PATHS["stages_path"],
                    {"composition": composition_paths(run_id),
                     "composition_high": composition_paths(run_id, high_res=True)},
                    overwrite=args.overwrite)


if __name__ == "__main__":
    main()
