"""Configure and run the steady-state C10H15O4 loss-budget plot."""
from datetime import timedelta
from pathlib import Path
import argparse

from cloud18_ap_analysis.loss_budget import loss_budget
from cloud18_ap_analysis.runs import make_temperature_runs

# Select prepared runs; their averaging intervals come from intervals.py (UTC).
run_ids = ["m50_1", "m25", "m10", "m50_2"]
high_res = True
output_dir = Path(__file__).resolve().parents[1] / "output" / "loss_budget"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--no-save", action="store_true")
    args = parser.parse_args()

    temperature_runs = make_temperature_runs(run_ids, prepared=True, high_res=high_res)
    trace_figures, means, regression_results, fig, ax, lambda_fig, lambda_ax = loss_budget(
        temperature_runs=temperature_runs,
        combine_runs={"-50 °C (combined)": ["m50_1", "m50_2"]},  # None keeps runs separate.
        fit_intercept=True,  # True fits y = lambda*x + b.
        pressure_pa=97000.0,
        wall_lifetime_s=480.0,
        dilution_rate_s=2e-4,
        production_yield_max=1.0,  # Direct ozonolysis only; conditional upper bound.
        radical_channel="C10H15O4.NH3H+",
        temperature_channel="TE_Calib_3",
        ozone_channel="o3",
        max_interpolation_gap_s=1800.0,
        average_separately=True,  # Average native samples first, then calculate rates; False aligns first using interpolation.
        label_intervals=False,  # Toggle matching t_i tags on traces and loss points.
        color_by_ap=True,       # AP colorbar; marker shapes identify temperature.
        overlay=True,
        show=["loss_kinetics", "k_loss(T)"],  # True, False, or any subset of ["traces", "loss_kinetics", "k_loss(T)"].
        smoothing_time_traces=timedelta(minutes=10),  # Display only; fits use raw data.
        stages="type",
        fan_speed=True,
        trace_yscale="log",
        interactive=not args.no_show,
        savepath=None if args.no_save else output_dir / "loss_budget_intercept.png",
        trace_savepath=None if args.no_save else str(output_dir / "trace_{run_id}.png"),
        lambda_temperature_savepath=None if args.no_save else output_dir / "lambda_vs_T_intercept.png",
        fit_savepath=None,  # Print fit results without saving a CSV.
        table_savepath=None,  # Print interval results without saving a CSV.
    )
