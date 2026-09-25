"""Public API for the CLOUD18 atmospheric-physics analysis package."""

from .models import Channel, Dataset, SMPSData
from .loading import load_data, read_trace_timeseries
from .processing import load_data_old, merge_log_files
from .interpolation import resample
from .stages import parse_cloud_log, load_stages
from .plotting import (
    plot_data,
    plot_all,
    plot_channels,
    plot_channel,
    plot_datasets,
    set_plot_backend,
)
from .MD_plot import md_plot
from .temp_dependence import calculate_interval_means, plot_kinetic_relationship, temp_dep
from .time_traces import plot_time_traces

__all__ = [
    "Channel", "Dataset", "SMPSData",
    "load_data", "read_trace_timeseries", "load_data_old", "merge_log_files",
    "resample", "parse_cloud_log", "load_stages",
    "plot_data", "plot_all", "plot_channels", "plot_channel",
    "plot_datasets", "set_plot_backend", "md_plot",
    "calculate_interval_means", "plot_kinetic_relationship", "temp_dep",
]
