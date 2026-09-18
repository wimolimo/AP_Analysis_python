from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Tuple, Iterable, List
from datetime import datetime, timezone
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


##########################################
# data loading for composition and trace files
##########################################
try:
    from . import config as cfg
except Exception:
    class cfg:
        ADDUCT_PRIORITY = {"NH3H+": 0, "H+": 1, "": 2}

# Ion/adduct masses to subtract from the exported Mass column
# Your examples suggest Mass includes the ion/adduct.
ELEMENT_MASSES = {
    "H": 1.00782503223,   "H(2)":  2.01410177812,
    "C": 12.0,            "C(13)": 13.00335483507,
    "N": 14.00307400443,  "N(15)": 15.00010889889,
    "O": 15.99491461957,  "O(18)": 17.99915961286,
    "S": 31.9720711744,   "S(34)": 33.967867004,
    "P": 30.97376199842,  "F": 18.99840316273,
    "Cl": 34.968852682,   "Cl(37)": 36.965902602,
    "Br": 78.9183376,     "I": 126.9044719,
}

_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\(\d+\))?(\d*)")

def split_sumformula(sumformula: str) -> Tuple[str, str]:
    """
    Split something like:
        CH2O.H+       -> ("CH2O", "H+")
        C2H2O.NH3H+   -> ("C2H2O", "NH3H+")
        C10H16O2      -> ("C10H16O2", "")
    """
    s = str(sumformula).strip()
    if "." not in s:
        return s, ""
    compound, adduct = s.split(".", 1)
    return compound.strip(), adduct.strip()

def element_counts(formula: str) -> dict:
    """'C(13)C9H16O3' -> {'C(13)': 1, 'C': 9, 'H': 16, 'O': 3}"""
    s = str(formula).strip()
    counts, pos = {}, 0
    for m in _TOKEN_RE.finditer(s):
        if m.start() != pos:
            raise ValueError(f"Cannot parse {formula!r} at {pos}")
        pos = m.end()
        sym, iso, num = m.groups()
        token = f"{sym}{iso}" if iso else sym
        counts[token] = counts.get(token, 0) + (int(num) if num else 1)
    if pos != len(s):
        raise ValueError(f"Cannot parse {formula!r} at {pos}")
    return counts


def exact_mass(formula: str) -> float:
    total = 0.0
    for token, n in element_counts(formula).items():
        if token not in ELEMENT_MASSES:
            return np.nan
        total += ELEMENT_MASSES[token] * n
    return total


def _element_total(counts: dict, symbol: str) -> int:
    """Sum an element across isotopes: 'C' matches 'C' and 'C(13)'."""
    return sum(n for t, n in counts.items()
               if t == symbol or t.startswith(f"{symbol}("))@dataclass

def oxygen_to_carbon_ratio(formula: str) -> float:
    """
    Calculate O:C ratio from the neutral compound formula.
    Returns NaN if C=0.
    """
    counts = element_counts(formula)
    c = counts.get("C", 0)
    o = counts.get("O", 0)

    if c == 0:
        return np.nan

    return o / c

def oxygen_count(formula: str) -> float:
    """
    Number of oxygen atoms in the neutral compound formula.
    Returns NaN if the formula cannot be parsed at all.
    """
    counts = element_counts(formula)
    if not counts:
        return np.nan
    return float(counts.get("O", 0))

def _ensure_datetimes(time_range: Tuple[datetime, datetime]) -> Tuple[datetime, datetime]:
    if (not isinstance(time_range, tuple)) or len(time_range) != 2:
        raise TypeError("time_range must be a tuple of two datetime objects")
    t0, t1 = time_range
    if not isinstance(t0, datetime) or not isinstance(t1, datetime):
        raise TypeError("time_range must contain datetime objects only")
    if t0.tzinfo is None:
        t0 = t0.replace(tzinfo=timezone.utc)
    if t1.tzinfo is None:
        t1 = t1.replace(tzinfo=timezone.utc)
    return t0, t1


def _as_unix_seconds(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp()

def _detect_and_skip_header_rows(path: str) -> int:
    n_header_rows = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        line = f.readline()
        if not line:
            return 0
        low = line.strip().lower()
        if "number of header rows" in low:
            tokens = re.split(r"[:\t, ]+", low)
            for token in reversed(tokens):
                try:
                    n_header_rows = int(float(token))
                    break
                except Exception:
                    continue
    return n_header_rows


def _detect_delimiter(path: str) -> str:
    skiprows = _detect_and_skip_header_rows(path)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for _ in range(skiprows):
            f.readline()
        header = f.readline()
    tabs = header.count("\t")
    commas = header.count(",")
    return "\t" if tabs > commas else ","


def _read_table_with_header_skip(path: str) -> pd.DataFrame:
    skiprows = _detect_and_skip_header_rows(path)
    sep = _detect_delimiter(path)
    return pd.read_csv(path, skiprows=skiprows, sep=sep, engine="c")


def read_composition_table(path: str, available: Optional[Iterable[str]] = None) -> pd.DataFrame:
    df = _read_table_with_header_skip(path)
    df.columns = [c.strip() for c in df.columns]
    colmap = {c.lower(): c for c in df.columns}

    if "sumformula" not in colmap or "mass" not in colmap:
        raise ValueError(
            f"Composition file must contain 'SumFormula' and 'Mass'. "
            f"Found: {df.columns.tolist()}"
        )

    out = df[[colmap["sumformula"], colmap["mass"]]].copy()
    out.columns = ["SumFormula", "Mass"]

    out["SumFormula"] = out["SumFormula"].astype(str).str.strip()
    out["Mass"] = pd.to_numeric(out["Mass"], errors="coerce")
    out = out.dropna(subset=["SumFormula", "Mass"])

    # Keep only ions the trace file can actually provide, so dedup never
    # discards a compound whose preferred adduct has no trace column.
    if available is not None:
        out = out[out["SumFormula"].isin(set(available))]
        if out.empty:
            raise RuntimeError("No composition SumFormula matches any trace column.")

    # Split e.g. CH2O.H+ -> CH2O + H+
    split = out["SumFormula"].apply(split_sumformula)
    out["CompoundFormula"] = split.apply(lambda x: x[0])
    out["Adduct"] = split.apply(lambda x: x[1])

    # Subtract adduct mass to get neutral compound mass
    out["NeutralMass"] = out["CompoundFormula"].apply(exact_mass)

    # O:C ratio from neutral compound formula
    out["OC"] = out["CompoundFormula"].apply(oxygen_to_carbon_ratio)

    # Number of oxygen atoms (for discrete coloring)
    out["On"] = out["CompoundFormula"].apply(oxygen_count)

    # Keep all adducted ions that are present in the trace file.
    # Selection/deduplication is deliberately left to the analysis module
    # that consumes this table (e.g. MD plot vs. time-trace analysis).

    print(df.head())

    return out


def select_preferred_adducts(
    df: pd.DataFrame,
    priority: Optional[dict] = None,
) -> pd.DataFrame:
    """
    Select one adduct per neutral compound using an explicit priority order.

    This reproduces the previous MD-plot behavior: among the adducts that
    are present in the already-filtered composition table, prefer NH3H+,
    then H+, then the neutral entry, unless a different priority mapping
    is supplied.

    Parameters
    ----------
    df:
        Composition table returned by ``read_composition_table``.
    priority:
        Mapping ``adduct -> priority`` where lower numbers are preferred.
        If omitted, ``config.ADDUCT_PRIORITY`` is used.

    Returns
    -------
    pd.DataFrame
        A copy containing one row per ``CompoundFormula``.
    """
    required = {"CompoundFormula", "Adduct"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Composition table is missing required columns: "
            + ", ".join(sorted(missing))
        )

    if priority is None:
        priority = getattr(
            cfg,
            "ADDUCT_PRIORITY",
            {"NH3H+": 0, "H+": 1, "": 2},
        )

    out = df.copy()
    out["_prio"] = out["Adduct"].map(priority).fillna(99)

    return (
        out.sort_values(["CompoundFormula", "_prio"])
        .drop_duplicates("CompoundFormula", keep="first")
        .drop(columns="_prio")
    )


def filter_adducts(
    df: pd.DataFrame,
    adducts: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    Keep only explicitly requested adducts.

    ``Adduct`` values are stored without the leading dot, e.g. ``NH3H+``
    or ``H+``. Passing ``None`` returns an unchanged copy.
    """
    if adducts is None:
        return df.copy()

    requested = {str(a).lstrip(".").strip() for a in adducts}
    if "Adduct" not in df.columns:
        raise ValueError("Composition table does not contain an 'Adduct' column.")

    return df[df["Adduct"].isin(requested)].copy()


def read_trace_header(path: str) -> dict:
    skiprows = _detect_and_skip_header_rows(path)
    sep = _detect_delimiter(path)
    df0 = pd.read_csv(path, skiprows=skiprows, sep=sep, engine="c", nrows=0)
    cols = [c.strip() for c in df0.columns]
    colmap = {c.lower(): c for c in cols}
    if "time" not in colmap or "unixtime" not in colmap:
        raise ValueError("Trace file must contain 'Time' and 'unixTime' columns.")
    time_col = colmap["time"]
    unix_col = colmap["unixtime"]
    intensity_cols = [c for c in cols if c not in {time_col, unix_col}]
    return {
        "time_col": time_col,
        "unix_col": unix_col,
        "intensity_cols": intensity_cols,
        "sep": sep,
        "skiprows": skiprows,
    }


def read_trace_table(path: str) -> pd.DataFrame:
    header = read_trace_header(path)
    df = pd.read_csv(
        path,
        skiprows=header["skiprows"],
        sep=header["sep"],
        engine="c",
    )
    print(df.head())

    df[header["unix_col"]] = pd.to_numeric(df[header["unix_col"]], errors="coerce")
    for c in header["intensity_cols"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=[header["unix_col"]])
    return df

def read_trace_timeseries(
    path: str,
    compounds: List[str],
    time_range: Tuple[datetime, datetime],
    chunksize: int = 200_000,
) -> pd.DataFrame:
    """
    Load selected compound time traces inside a requested time window.

    The trace file is read in chunks and filtered using the ``unixTime``
    column, so the complete trace file does not need to be loaded into memory.
    """
    if not compounds:
        raise ValueError("compounds must contain at least one compound")

    t0_dt, t1_dt = _ensure_datetimes(time_range)
    if t1_dt < t0_dt:
        raise ValueError("time_range end < start")

    compounds = [str(c).strip() for c in compounds]

    header = read_trace_header(path)
    time_col = header["time_col"]
    unix_col = header["unix_col"]
    available = set(header["intensity_cols"])

    missing = [c for c in compounds if c not in available]
    if missing:
        raise ValueError(
            "The following compounds were not found in the trace file: "
            + ", ".join(missing)
        )

    usecols = [time_col, unix_col] + compounds
    dtype_map = {unix_col: "float64"}
    dtype_map.update({c: "float64" for c in compounds})

    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)

    chunks = []
    for chunk in pd.read_csv(
        path,
        skiprows=header["skiprows"],
        sep=header["sep"],
        engine="c",
        usecols=usecols,
        chunksize=chunksize,
        dtype=dtype_map,
    ):
        unix = pd.to_numeric(chunk[unix_col], errors="coerce")
        selected = chunk.loc[unix.between(t0s, t1s)].copy()

        if not selected.empty:
            chunks.append(selected)

        if not unix.dropna().empty and unix.max() >= t1s:
            break

    if not chunks:
        return pd.DataFrame(columns=[time_col, unix_col] + compounds)

    out = pd.concat(chunks, ignore_index=True)

    parsed_time = pd.to_datetime(out[time_col], errors="coerce", utc=True)
    if parsed_time.notna().any():
        out[time_col] = parsed_time

    return out.sort_values(unix_col).reset_index(drop=True)

def classical_mass_defect(m: Iterable[float]):
    """
    Classic MD only: MD = m - round(m)
    Returns (nominal_mass_array, defect_array)
    """
    m = np.asarray(m, dtype=float)
    nominal = np.rint(m).astype(int)
    defect = m - nominal
    return nominal, defect
def aggregate_window_unix_dt(
        df_tr: pd.DataFrame,
        t0_dt: datetime,
        t1_dt: datetime,
        clip_neg: bool = False,
    ) -> pd.Series:
    colmap = {c.lower(): c for c in df_tr.columns}
    time_col = colmap["time"]
    unix_col = colmap["unixtime"]
    intensity_cols = [c for c in df_tr.columns if c not in {time_col, unix_col}]
    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)
    sel = df_tr[(df_tr[unix_col] >= t0s) & (df_tr[unix_col] <= t1s)]
    if sel.empty:
        return pd.Series(dtype=float)
    if clip_neg:
        sel = sel.copy()
        sel[intensity_cols] = sel[intensity_cols].clip(lower=0)
    return sel[intensity_cols].mean(axis=0)


def aggregate_window_streaming(
    trace_path: str,
    t0_dt: datetime,
    t1_dt: datetime,
    formulas: List[str],
    clip_neg: bool = False,
    assume_sorted_unix: bool = True,
    chunksize: int = 200_000,
    skiprows: Optional[int] = None,
    sep: Optional[str] = None,
) -> pd.Series:
    """
    Aggregate trace intensities inside a time window without loading the full CSV.

    The trace file can be very large, so this function reads it in chunks and
    only keeps rows whose ``unixTime`` lies inside ``[t0_dt, t1_dt]``. The
    returned series is indexed by ``formulas`` and contains the aggregated
    intensity for each compound column.

    Workflow:
    - detect the trace file header layout and delimiter
    - read the ``unixTime`` column plus the requested compound columns in chunks
    - filter each chunk to the requested time window
    - optionally clip negative intensities to zero
    - accumulate the mean intensity across all selected rows

    Parameters
    ----------
    trace_path:
        Path to the trace CSV/TXT file.
    t0_dt, t1_dt:
        Start and end of the aggregation window.
    formulas:
        Compound columns to aggregate.
    clip_neg:
        If True, negative values are clipped to zero before aggregation.
    assume_sorted_unix:
        If True, stop scanning once the file is past the end of the window.
    chunksize:
        Number of rows to read per chunk.
    skiprows, sep:
        Optional header settings; if omitted they are detected automatically.
    """
    if skiprows is None:
        skiprows = _detect_and_skip_header_rows(trace_path)
    if sep is None:
        sep = _detect_delimiter(trace_path)

    header = read_trace_header(trace_path)
    unix_col = header["unix_col"]

    # Read only the timestamp column and the compounds we care about.
    usecols = [unix_col] + formulas
    dtype_map = {unix_col: "float64"}
    dtype_map.update({f: "float64" for f in formulas})

    # Convert the datetime window to unix seconds for numeric comparison.
    t0s = _as_unix_seconds(t0_dt)
    t1s = _as_unix_seconds(t1_dt)

    acc = None
    rows_in_window = 0

    # To support early stopping, we read the unixTime column separately so we
    # can inspect each chunk's bounds before processing the intensity columns.
    unix_usecols = [unix_col]  # small
    for chunk_unix, chunk in zip(
        pd.read_csv(
            trace_path, skiprows=skiprows, sep=sep, engine="c",
            usecols=unix_usecols, chunksize=chunksize, dtype={unix_col: "float64"}
        ),
        pd.read_csv(
            trace_path, skiprows=skiprows, sep=sep, engine="c",
            usecols=usecols, chunksize=chunksize, dtype=dtype_map
        )
    ):
        # Chunk bounds tell us whether this chunk is entirely before or after
        # the requested window when the file is sorted by time.
        min_u = float(chunk_unix[unix_col].min())
        max_u = float(chunk_unix[unix_col].max())

        # Keep only rows inside the requested time window.
        chunk = chunk[pd.to_numeric(chunk[unix_col], errors="coerce").between(t0s, t1s)]
        if not chunk.empty:
            # Optionally suppress negative intensities before aggregating.
            if clip_neg:
                chunk[formulas] = chunk[formulas].clip(lower=0)

            # Accumulate the column-wise sum for this chunk. The mean is
            # computed once at the end so chunk boundaries do not bias it.
            chunk_sum = chunk[formulas].sum(axis=0)
            if acc is None:
                acc = chunk_sum
            else:
                acc = acc.add(chunk_sum, fill_value=0.0)
            rows_in_window += len(chunk)

        # Stop early once the file has moved past the end of the time window.
        if assume_sorted_unix and max_u >= t1s:
            break
        # Skip chunks that are entirely before the window.
        if assume_sorted_unix and max_u < t0s:
            continue
        # If unsorted, we keep scanning all chunks.

    if acc is None:
        return pd.Series(index=formulas, dtype=float)

    acc = acc / max(rows_in_window, 1)

    print(acc.head(15))

    return acc.reindex(formulas)