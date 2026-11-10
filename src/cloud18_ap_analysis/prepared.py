"""Compressed, versioned run files with independent instrument time axes.

No interpolation, smoothing, unit conversion or precision reduction is done
during preparation. A PreparedSource can be passed to the existing loaders.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .models import Channel, Dataset, SMPSData

SCHEMA_VERSION = 1
PREPARED_DIR = Path(__file__).resolve().parents[2] / "output" / "prepared"


def _h5py():
    # Keep direct CSV workflows usable without importing h5py.
    import h5py
    return h5py


def _json(value):
    def encode(item):
        if isinstance(item, (datetime, Path)):
            return item.isoformat() if isinstance(item, datetime) else str(item)
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(f"Cannot serialize {type(item).__name__}")
    return json.dumps(value, default=encode, ensure_ascii=False)


def _ms(value):
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.timestamp() * 1000


def _datetime(value):
    return datetime.fromtimestamp(float(value) / 1000, timezone.utc).replace(tzinfo=None)


@dataclass(frozen=True)
class SourceSpec:
    """Explicit input format and units for one instrument.

    Supported formats: trace, trace_stream, ltof_processed, paired, custom, ap,
    smps. trace_stream prepares large CSVs with bounded memory. Unknown units stay
    explicitly 'unknown'; add channel overrides before using them in calculations.
    Later files win at identical timestamps within each channel.
    """
    paths: list[Path]
    format: str
    units: dict[str, str] = field(default_factory=dict)
    default_unit: str = "unknown"
    padding_minutes: int = 0  # Legacy compatibility; full-source preparation needs no padding.


@dataclass(frozen=True)
class PreparedSource:
    path: Path
    instrument: str

    def load(self, channels=None, time_range=None):
        return read_prepared(self, channels=channels, time_range=time_range)


def _open(path):
    file = _h5py().File(path, "r")
    if file.attrs.get("schema_version") != SCHEMA_VERSION:
        file.close()
        raise ValueError(f"Unsupported prepared-file schema: {path}")
    return file


def _bounds(axis, time_range):
    """Binary search on disk: do not read the complete timestamp array."""
    if time_range is None:
        return slice(None)
    start, stop = map(_ms, time_range)
    if stop < start:
        raise ValueError("time_range end < start")

    def search(value, right=False):
        lo, hi = 0, len(axis)
        while lo < hi:
            mid = (lo + hi) // 2
            current = axis[mid]
            if current < value or (right and current == value):
                lo = mid + 1
            else:
                hi = mid
        return lo

    return slice(search(start), search(stop, right=True))


def _check_window(file, time_range):
    if time_range is None:
        return
    start, stop = map(_ms, time_range)
    if stop < start:
        raise ValueError("time_range end < start")
    if file.attrs.get("preparation_scope") == "all_source_samples":
        return
    window = file.attrs["time_range_ms"]
    if start < window[0] or stop > window[1]:
        raise ValueError("Requested window exceeds this older cropped prepared run. "
                         "Rebuild with scripts/prepare_runs.py to retain all source samples.")


def read_prepared(source, channels=None, time_range=None):
    """Read only selected channel datasets and the inclusive requested time slice."""
    with _open(source.path) as file:
        _check_window(file, time_range)
        group = file["instruments"][source.instrument]
        if group.attrs["kind"] == "smps":
            if channels is not None:
                raise ValueError("Select SMPS as a distribution, not scalar channels")
            rows = _bounds(group["time"], time_range)
            return SMPSData(source.instrument,
                            [_datetime(t) for t in group["time"][rows]],
                            group["diameters"][:], group["values"][rows, :])
        names = list(group["channels"])
        selected = names if channels is None else list(channels)
        missing = set(selected).difference(names)
        if missing:
            raise KeyError(f"Missing {source.instrument} channels: {sorted(missing)}")
        result = {}
        axes = {}
        for name in selected:
            item = group["channels"][name]
            axis_name = item.attrs["axis"]
            if axis_name not in axes:
                axis = group["axes"][axis_name]
                rows = _bounds(axis, time_range)
                axes[axis_name] = (rows, axis[rows])
            rows, times = axes[axis_name]
            result[name] = Channel(name, times, item[rows])
        return Dataset(result)


def prepared_channel_info(source):
    """Return names/units/coverage without loading measurement arrays."""
    with _open(source.path) as file:
        group = file["instruments"][source.instrument]
        return json.loads(group.attrs["report"])


def read_prepared_traces(source, compounds, time_range):
    if not compounds:
        raise ValueError("compounds must contain at least one compound")
    data = read_prepared(source, compounds, time_range)
    # PTR exports have a common time axis. Outer joins also support channel
    # sets that change between source files without inventing measurements.
    series = [pd.Series(ch.values, index=ch.time, name=name)
              for name, ch in data.channels.items()]
    frame = pd.concat(series, axis=1).sort_index()
    frame.insert(0, "unixTime", frame.index.to_numpy() / 1000)
    frame.insert(0, "Time", pd.to_datetime(frame.index.to_numpy(), unit="ms"))
    return frame.reset_index(drop=True)


def read_prepared_stages(source):
    with _open(source.path) as file:
        stages = json.loads(file["stages"].asstr()[()])
    for stage in stages:
        stage["time"] = datetime.fromisoformat(stage["time"]) if stage["time"] else None
    return stages


def read_prepared_composition(source, available=None):
    with _open(source.path) as file:
        key = source.instrument
        if key not in file:
            raise ValueError("No composition table stored; configure composition_paths and rebuild")
        frame = pd.DataFrame(json.loads(file[key].asstr()[()]))
    if available is not None:
        frame = frame[frame.SumFormula.isin(available)]
    return frame


def select_prepared_resolution(source, high_res):
    """Select a prepared PTR or composition source; never silently fall back."""
    if not isinstance(high_res, bool):
        raise TypeError("high_res must be True or False")
    if not isinstance(source, PreparedSource):
        raise TypeError("Resolution selection requires a PreparedSource; select raw paths in runs.py")
    composition = source.instrument.startswith("composition")
    key = ("composition_high" if high_res else "composition") if composition else ("ptr_high" if high_res else "ptr")
    with _open(source.path) as file:
        container = file if composition else file["instruments"]
        if key not in container:
            raise ValueError(f"{'High' if high_res else 'Low'}-resolution PTR {'composition' if composition else 'traces'} not prepared in {source.path.name}; configure runs.py and rebuild with scripts/prepare_runs.py")
    return PreparedSource(source.path, key)


def load_prepared_run(run_id, directory=PREPARED_DIR, *, high_res=False):
    """Return instrument references and current analysis intervals.

    H5 retains its preparation settings as provenance. Known campaign runs use
    the current runs.py time_range and intervals.py intervals at plot time.
    high_res selects ptr_high and composition_high; False selects the original
    ptr and composition datasets. Missing resolutions raise an error.
    """
    path = Path(directory) / f"{run_id}.h5"
    if not path.exists():
        raise FileNotFoundError(f"Prepare {run_id} first: python scripts/prepare_runs.py {run_id}")
    with _open(path) as file:
        run = json.loads(file["run"].asstr()[()])
        instruments = list(file["instruments"])
        window = file.attrs["time_range_ms"]
    run["prepared_time_range"] = tuple(_datetime(t) for t in window)
    run["time_range"] = tuple(datetime.fromisoformat(t) for t in run["time_range"])
    from .runs import RUNS
    if run_id in RUNS:
        run["time_range"] = tuple(RUNS[run_id]["time_range"])
    run["important_time_ranges"] = [tuple(datetime.fromisoformat(t) for t in pair)
                                    for pair in run["important_time_ranges"]]
    from .intervals import IMPORTANT_TIME_RANGES
    if run_id in IMPORTANT_TIME_RANGES:
        run["important_time_ranges"] = list(IMPORTANT_TIME_RANGES[run_id])
    ptr = select_prepared_resolution(PreparedSource(path, "ptr"), high_res)
    run["sources"] = {name: PreparedSource(path, name) for name in instruments if name != "ptr_high"}
    run["sources"]["ptr"] = ptr
    run["high_res"] = high_res
    run["stages_path"] = PreparedSource(path, "stages")
    run["composition_path"] = PreparedSource(path, "composition_high" if high_res else "composition")
    run["trace_paths"] = run["sources"]["ptr"]
    run["temperature_path"] = run["sources"]["temperature"]
    run["ozone_path"] = run["sources"]["ozone"]
    return run


def _fingerprint(path):
    path = Path(path)
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                sha256=digest.hexdigest())


def _load_source(spec, time_range):
    from .loading import (load_trace_format, _load_normal_format, _load_custom_format,
                          _load_ap_format, load_smps_format, _report_trace_overlaps,
                          load_processed_ltof)
    readers = dict(trace=load_trace_format, paired=_load_normal_format,
                   ap=_load_ap_format, smps=load_smps_format,
                   ltof_processed=load_processed_ltof)
    if spec.format not in {*readers, "custom"}:
        raise ValueError(f"Unknown source format {spec.format!r}")
    datasets = []
    for path in spec.paths:
        if spec.format == "custom":
            with Path(path).open(encoding="utf-8", errors="replace") as stream:
                first_line = stream.readline()
            data = _load_custom_format(Path(path), first_line, time_range)
        else:
            data = readers[spec.format](path, time_range=time_range)
        datasets.append(data)
    if not datasets:
        raise ValueError("Source paths must not be empty")
    if isinstance(datasets[0], SMPSData):
        if len(datasets) != 1:
            raise ValueError("Configure one merged SMPS file per source")
        return datasets[0], {}
    if len(datasets) > 1:
        _report_trace_overlaps(spec.paths, [np.concatenate([ch.time for ch in d.channels.values()])
                                          if d.channels else np.array([]) for d in datasets])
    channels, duplicates = {}, {}
    for name in dict.fromkeys(name for data in datasets for name in data.channels):
        parts = [data[name] for data in datasets if name in data.channels]
        frame = pd.DataFrame({"time": np.concatenate([p.time for p in parts]),
                              "value": np.concatenate([p.values for p in parts])})
        frame = frame[np.isfinite(frame.time)]
        original_count = len(frame)
        frame = frame.drop_duplicates("time", keep="last").sort_values("time")
        duplicates[name] = original_count - len(frame)
        channels[name] = Channel(name, frame.time.to_numpy(), frame.value.to_numpy())
    return Dataset(channels), duplicates


def _write_array(group, name, values):
    values = np.asarray(values, dtype="float64")
    return group.create_dataset(name, data=values, compression="gzip", compression_opts=4,
                                shuffle=True, chunks=True)


def _write_instrument(group, data, spec, duplicates):
    report = {}
    if isinstance(data, SMPSData):
        group.attrs["kind"] = "smps"
        times = np.array([_ms(t) if isinstance(t, datetime) else np.nan for t in data.time])
        valid = np.isfinite(times)
        indices = np.flatnonzero(valid)
        # Stable order preserves the explicit last-row-wins policy.
        frame = pd.DataFrame({"time": times[indices], "index": indices})
        frame = frame.drop_duplicates("time", keep="last").sort_values("time")
        indices = frame["index"].to_numpy(dtype=int)
        _write_array(group, "time", times[indices]).attrs["unit"] = "unix_ms_UTC"
        _write_array(group, "diameters", data.diameters).attrs["unit"] = "nm"
        _write_array(group, "values", data.dNdlogDp[indices]).attrs["unit"] = spec.default_unit
        report = dict(samples=len(indices), bins=len(data.diameters),
                      duplicates_removed=int(valid.sum() - len(indices)))
    else:
        group.attrs["kind"] = "channels"
        axes, channels = group.create_group("axes"), group.create_group("channels", track_order=True)
        axis_names = {}
        for name, channel in data.channels.items():
            if "/" in name:
                raise ValueError(f"Channel names cannot contain '/' in schema v1: {name}")
            key = hashlib.sha256(channel.time.tobytes()).hexdigest()
            if key not in axis_names:
                axis_name = str(len(axis_names))
                _write_array(axes, axis_name, channel.time).attrs["unit"] = "unix_ms_UTC"
                axis_names[key] = axis_name
            item = _write_array(channels, name, channel.values)
            item.attrs["axis"] = axis_names[key]
            item.attrs["unit"] = spec.units.get(name, spec.default_unit)
            report[name] = dict(unit=item.attrs["unit"], samples=len(channel.time),
                                missing=int((~np.isfinite(channel.values)).sum()),
                                start_utc=_datetime(channel.time[0]).isoformat() if len(channel.time) else None,
                                stop_utc=_datetime(channel.time[-1]).isoformat() if len(channel.time) else None,
                                duplicates_removed=duplicates[name])
    group.attrs["report"] = _json(report)
    return report


def _write_streamed_traces(group, spec, time_range, chunksize=10000):
    """Prepare large trace exports with bounded memory and last-row-wins merging.

    Stage each file on disk, then merge one channel at a time. This preserves
    float64 samples and supports overlapping files without loading all channels.
    """
    from .loading import read_trace_header, _range_mask, _report_trace_overlaps
    fd, staging_path = tempfile.mkstemp(prefix=".ptr-staging-", suffix=".h5",
                                         dir=Path(group.file.filename).parent)
    os.close(fd)
    try:
        with _h5py().File(staging_path, "w") as staging:
            names = None
            times = []
            for index, path in enumerate(spec.paths):
                header = read_trace_header(path)
                current = header["intensity_cols"]
                if names is None:
                    names = current
                if set(current) != set(names):
                    raise ValueError("Streamed trace files must contain the same channels")
                part = staging.create_group(str(index))
                for name in ["_time", *names]:
                    if "/" in name:
                        raise ValueError(f"Unsupported channel name: {name}")
                    part.create_dataset(name, shape=(0,), maxshape=(None,), dtype="float64",
                                        chunks=(8192,))
                count = 0
                reader = pd.read_csv(path, skiprows=header["skiprows"], sep=header["sep"],
                                     usecols=[header["unix_col"], *names], chunksize=chunksize)
                for frame in reader:
                    time_ms = pd.to_numeric(frame[header["unix_col"]], errors="coerce").to_numpy(dtype=float) * 1000
                    mask = np.isfinite(time_ms) & _range_mask(time_ms, time_range)
                    size = int(mask.sum())
                    if not size:
                        continue
                    stop = count + size
                    part["_time"].resize((stop,))
                    part["_time"][count:stop] = time_ms[mask]
                    for name in names:
                        item = part[name]
                        item.resize((stop,))
                        item[count:stop] = pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float)[mask]
                    count = stop
                times.append(part["_time"][:])
                print(f"Staged {Path(path).name}: {count} high-resolution samples", flush=True)
            if names is None:
                raise ValueError("Source paths must not be empty")
            _report_trace_overlaps(spec.paths, times)
            all_times = np.concatenate(times)
            # Reverse lookup implements the same last-listed-file/row policy.
            merged_times, reversed_indices = np.unique(all_times[::-1], return_index=True)
            indices = len(all_times) - 1 - reversed_indices
            duplicates = len(all_times) - len(merged_times)
            group.attrs["kind"] = "channels"
            axes = group.create_group("axes")
            _write_array(axes, "0", merged_times).attrs["unit"] = "unix_ms_UTC"
            channels = group.create_group("channels", track_order=True)
            report = {}
            for name in names:
                values = np.concatenate([staging[str(i)][name][:] for i in range(len(times))])[indices]
                item = _write_array(channels, name, values)
                item.attrs["axis"] = "0"
                item.attrs["unit"] = spec.units.get(name, spec.default_unit)
                report[name] = dict(unit=item.attrs["unit"], samples=len(merged_times),
                                    missing=int((~np.isfinite(values)).sum()),
                                    start_utc=_datetime(merged_times[0]).isoformat() if len(merged_times) else None,
                                    stop_utc=_datetime(merged_times[-1]).isoformat() if len(merged_times) else None,
                                    duplicates_removed=duplicates)
            group.attrs["report"] = _json(report)
            return report
    finally:
        Path(staging_path).unlink(missing_ok=True)


def prepare_run(run_id, run, sources, stages_path, composition_paths=(),
                directory=PREPARED_DIR, overwrite=False):
    """Store all configured source samples, independently of the plot range.

    Build atomically; an unsuccessful rebuild leaves the previous file intact.
    """
    from .stages import parse_cloud_log
    from .loading import read_composition_table

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"{run_id}.h5"
    if target.exists() and not overwrite:
        raise FileExistsError(f"{target} exists; use --overwrite to rebuild")
    start, stop = run["time_range"]
    if _ms(stop) <= _ms(start):
        raise ValueError("Run stop must be after start")
    compositions = (composition_paths if isinstance(composition_paths, dict)
                    else {"composition": composition_paths})
    all_paths = list(dict.fromkeys(Path(p) for spec in sources.values() for p in spec.paths))
    all_paths += [Path(stages_path), *(Path(p) for paths in compositions.values() for p in paths)]
    fingerprints = [_fingerprint(p) for p in dict.fromkeys(all_paths)]
    fd, temporary = tempfile.mkstemp(prefix=f".{run_id}-", suffix=".h5", dir=directory)
    os.close(fd)
    try:
        with _h5py().File(temporary, "w", track_order=True) as file:
            file.attrs["schema_version"] = SCHEMA_VERSION
            file.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
            file.attrs["time_range_ms"] = [_ms(start), _ms(stop)]
            file.attrs["preparation_scope"] = "all_source_samples"
            file.attrs["duplicate_policy"] = "last listed file/row wins per channel"
            string_type = _h5py().string_dtype("utf-8")
            def store_json(name, value):
                file.create_dataset(name, data=_json(value), dtype=string_type)
            settings = {k: run[k] for k in ("label", "temperature_c", "time_range", "important_time_ranges")}
            settings["run_id"] = run_id
            store_json("run", settings)
            store_json("provenance", fingerprints)
            instruments = file.create_group("instruments")
            reports = {}
            for name, spec in sources.items():
                group = instruments.create_group(name)
                group.attrs["source_spec"] = _json(vars(spec))
                if spec.format == "trace_stream":
                    reports[name] = _write_streamed_traces(group, spec, None)
                else:
                    data, duplicates = _load_source(spec, None)
                    reports[name] = _write_instrument(group, data, spec, duplicates)
                print(f"Prepared {run_id}/{name}", flush=True)
            # Store actual source coverage separately from the chosen plot range.
            bounds = []
            for group in instruments.values():
                axes = [group["time"]] if group.attrs["kind"] == "smps" else group["axes"].values()
                for axis in axes:
                    if len(axis):
                        bounds.append((float(axis[0]), float(axis[-1])))
            if bounds:
                file.attrs["time_range_ms"] = [min(a for a, b in bounds), max(b for a, b in bounds)]
            store_json("stages", parse_cloud_log(stages_path))
            for key, paths in compositions.items():
                if paths:
                    tables = [read_composition_table(p) for p in paths]
                    composition = pd.concat(tables, ignore_index=True).drop_duplicates("SumFormula", keep="last")
                    store_json(key, composition.to_dict(orient="records"))
            store_json("report", reports)
        for original in fingerprints:
            stat = Path(original["path"]).stat()
            if stat.st_size != original["size"] or stat.st_mtime_ns != original["mtime_ns"]:
                raise RuntimeError(f"Source changed during preparation: {original['path']}")
        os.replace(temporary, target)
    finally:
        if Path(temporary).exists():
            Path(temporary).unlink()
    print(f"Saved {target} ({target.stat().st_size / 1e6:.2f} MB)")
    return target
