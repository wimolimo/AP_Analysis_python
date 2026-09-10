from __future__ import annotations

from pathlib import Path
from datetime import datetime
import calendar
import numpy as np
import pandas as pd

from .models import Channel, Dataset, SMPSData


def _datetime_to_unix_ms(dt: datetime) -> float:
    return (calendar.timegm(dt.timetuple()) + dt.microsecond / 1_000_000.0) * 1000.0


def _range_mask(time_ms, time_range):
    if time_range is None:
        return np.ones(len(time_ms), dtype=bool)
    lo = _datetime_to_unix_ms(time_range[0])
    hi = _datetime_to_unix_ms(time_range[1])
    return (time_ms >= lo) & (time_ms <= hi)

def _unix_seconds_to_datetime(value):
    """Convert unix seconds to datetime.

    Missing/NaN values are returned as np.nan so downstream code can filter them.
    """
    try:
        value = float(value)
    except Exception:
        return np.nan

    if not np.isfinite(value):
        return np.nan

    return datetime.utcfromtimestamp(round(value))


def _utcify_datetime(value):
    if isinstance(value, datetime):
        return value
    return value

def _filter_time_range_dt(timestamps, time_range):
    """Return boolean mask for optional datetime range."""
    if time_range is None:
        return np.ones(len(timestamps), dtype=bool)

    start, stop = (_utcify_datetime(time_range[0]), _utcify_datetime(time_range[1]))
    mask = []

    for t in timestamps:
        if isinstance(t, datetime):
            t = _utcify_datetime(t)
            mask.append(start <= t <= stop)
        else:
            mask.append(False)

    return np.asarray(mask, dtype=bool)

def _numeric_time_to_ms(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    if not np.any(finite):
        return values

    out = values.astype(float, copy=True)
    magnitude = float(np.nanmedian(np.abs(out[finite])))

    # Infer epoch unit by magnitude.
    if magnitude >= 1e17:      # nanoseconds
        out[finite] = out[finite] / 1_000_000.0
    elif magnitude >= 1e14:    # microseconds
        out[finite] = out[finite] / 1_000.0
    elif magnitude >= 1e11:    # milliseconds
        pass
    elif magnitude >= 1e8:     # seconds
        out[finite] = out[finite] * 1000.0

    return out


def _timestamps_to_unix_ms(series: pd.Series, preferred_format: str | None = None) -> np.ndarray:
    numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    finite_ratio = float(np.isfinite(numeric).mean()) if len(numeric) else 0.0

    # Some exports carry epoch timestamps instead of formatted datetime strings.
    if finite_ratio > 0.95:
        return _numeric_time_to_ms(numeric)

    if preferred_format is not None:
        timestamps = pd.to_datetime(series, format=preferred_format, errors="coerce")
    else:
        timestamps = pd.to_datetime(series, errors="coerce")

    # Fallback in case the strict format misses a variant.
    if timestamps.isna().all():
        timestamps = pd.to_datetime(series, errors="coerce")

    # Normalize to millisecond precision first so conversion is unit-stable
    # across pandas datetime resolutions (s/us/ns).
    time_ms = timestamps.to_numpy(dtype="datetime64[ms]").astype("int64").astype(float)
    invalid = timestamps.isna().to_numpy()
    time_ms[invalid] = np.nan
    return time_ms


def load_data(filepath: str | Path, time_range=None):
    """Load one of the data formats used by the CLOUD18 analysis.

    SMPS data are intentionally not supported in this first Python version.
    """
    filepath = Path(filepath)
    with filepath.open("r", encoding="utf-8", errors="replace") as f:
        first_line = f.readline().strip()

    if first_line.lower().startswith("numberofheaderrows"):
        return _load_custom_format(filepath, first_line, time_range)
    if first_line.lower().startswith("number of header rows"):
        return load_smps_format(filepath, first_line, time_range)
    if filepath.name.startswith("AP"):
        return _load_ap_format(filepath, time_range)
    return _load_normal_format(filepath, time_range)


def _load_normal_format(filepath, time_range=None):
    df = pd.read_csv(filepath)
    channels = {}
    print("loading normal format data from", filepath)

    # pairing columns as (time, value).
    for i in range(0, len(df.columns) - 1, 2):
        time_raw = pd.to_numeric(df.iloc[:, i], errors="coerce").to_numpy(float)
        values = pd.to_numeric(df.iloc[:, i + 1], errors="coerce").to_numpy(float)

        valid_time = np.isfinite(time_raw)
        time_vals = time_raw[valid_time]
        data_vals = values[valid_time]

        mask = _range_mask(time_vals, time_range)
        name = str(df.columns[i + 1])
        channels[name] = Channel(name, time_vals[mask], data_vals[mask])

    return Dataset(channels)


def _read_header_names(lines, header_line_index):
    return [x.strip() for x in lines[header_line_index].split(",")]


def _load_custom_format(filepath, first_line, time_range=None):
    n_header = int(first_line.split(":", 1)[1].strip())
    lines = filepath.read_text(encoding="utf-8", errors="replace").splitlines()
    col_names = _read_header_names(lines, n_header - 1)  # Julia is 1-based.

    print(f"Loading custom format data from {filepath}")

    df = pd.read_csv(
        filepath,
        header=None,
        skiprows=n_header,
        names=col_names,
        na_values=["NaN"],
    )

    time_ms = _timestamps_to_unix_ms(df.iloc[:, 0], preferred_format="%Y-%m-%d %H:%M:%S")

    mask = _range_mask(time_ms, time_range) & np.isfinite(time_ms)
    channels = {}
    for name in col_names[1:]:
        values = pd.to_numeric(df[name], errors="coerce").to_numpy(float)
        channels[name] = Channel(name, time_ms[mask], values[mask])
    return Dataset(channels)


def _load_ap_format(filepath, time_range=None):
    df = pd.read_csv(filepath)
    time_ms = _timestamps_to_unix_ms(df.iloc[:, 0], preferred_format="%Y-%m-%d %H:%M:%S.%f")
    print(f"Loading AP format data from {filepath}")

    mask = _range_mask(time_ms, time_range) & np.isfinite(time_ms)
    channels = {}
    for name in df.columns[1:]:
        values = pd.to_numeric(df[name], errors="coerce").to_numpy(float)
        channels[str(name)] = Channel(str(name), time_ms[mask], values[mask])
    return Dataset(channels)

def load_smps_format(filepath, first_line=None, time_range=None):
    """Load SMPS CSV format.

    Expected file structure:
    - first line contains something like "...: N"
      where N is the number of header lines before the column-name line.
    - column-name line is line N+1 in Julia 1-based indexing.
    - first column is unix time in seconds.
    - remaining columns are diameter-bin names and must parse as floats.
    """
    filepath = Path(filepath)
    print(f"Loading SMPS data from {filepath}")

    if first_line is None:
        with filepath.open("r", encoding="utf-8", errors="replace") as f:
            first_line = f.readline()

    try:
        n_header = int(first_line.split(":", 1)[1].strip())
    except Exception as exc:
        raise ValueError(
            f"Could not parse SMPS header count from first line: {first_line!r}"
        ) from exc

    lines_all = filepath.read_text(encoding="utf-8", errors="replace").splitlines()

    try:
        header_line = lines_all[n_header]
    except IndexError as exc:
        raise ValueError(
            f"SMPS file does not contain expected header line at index {n_header}"
        ) from exc

    col_names = [x.strip() for x in header_line.split(",")]
    n_cols = len(col_names)

    if n_cols < 2:
        raise ValueError("SMPS file must contain time column plus at least one diameter bin")

    df = pd.read_csv(
        filepath,
        header=None,
        skiprows=n_header + 1,
        usecols=range(n_cols),
        na_values=["NaN", "nan", ""],
    )

    df.columns = col_names

    time_col = col_names[0]
    raw_time = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=float)

    timestamps = np.asarray(
        [_unix_seconds_to_datetime(t) for t in raw_time],
        dtype=object,
    )

    mask = _filter_time_range_dt(timestamps, time_range)
    timestamps = timestamps[mask]

    bin_names = col_names[1:]

    try:
        diameters = np.asarray([float(name.strip()) for name in bin_names], dtype=float)
    except Exception as exc:
        raise ValueError(
            "Could not parse SMPS diameter bin names as floats. "
            f"Got bin names: {bin_names[:10]}..."
        ) from exc

    mat = df[bin_names].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    # Apply time mask.
    mat = mat[mask, :]

    return SMPSData("SMPS", timestamps, diameters, mat)