from pathlib import Path
import pandas as pd
import numpy as np

from .models import Channel, Dataset


def load_data_old(filepath: str | Path):
    """Compatibility implementation of the old pair-column loader."""
    df = pd.read_csv(filepath, header=0)
    channels = {}
    for i in range(0, len(df.columns) - 1, 2):
        time_col = pd.to_numeric(df.iloc[:, i], errors="coerce").to_numpy(float)
        values = pd.to_numeric(df.iloc[:, i + 1], errors="coerce").to_numpy(float)
        mask = np.isfinite(time_col) & np.isfinite(values)
        name = str(df.columns[i + 1])
        channels[name] = Channel(name, time_col[mask], values[mask])
    return Dataset(channels)


def merge_log_files(filepaths, output_path):
    with open(output_path, "w", encoding="utf-8") as out:
        for path in filepaths:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                out.write(f.read())
            out.write("\n")
