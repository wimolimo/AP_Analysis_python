"""Resolve neutral formulas against instrument-specific trace channels."""
import re
from .loading import element_counts, read_trace_header
from .prepared import PreparedSource, prepared_channel_info


def available_channels(source):
    if isinstance(source, PreparedSource):
        return list(prepared_channel_info(source))
    paths = source if isinstance(source, (list, tuple)) else [source]
    return list(dict.fromkeys(name for path in paths
                             for name in read_trace_header(path)["intensity_cols"]))


def formula_key(formula):
    try:
        counts = element_counts(formula)
        return tuple(sorted(counts.items())) if counts else None
    except ValueError:
        return None


def resolve_channel(request, available, preferences=()):
    # A full channel name is an explicit override, including pandas duplicate suffixes.
    if "." in request:
        if request in available:
            return request
        raise KeyError(f"Channel {request!r} not found")
    key = formula_key(request)
    candidates = [name for name in available
                  if key is not None and formula_key(name.split(".", 1)[0]) == key]
    if len(candidates) == 1:
        return candidates[0]
    for adduct in preferences:
        preferred = [name for name in candidates
                     if re.sub(r"\.\d+$", "", name.partition(".")[2]) == adduct]
        if len(preferred) == 1:
            return preferred[0]
    if candidates:
        raise ValueError(f"Ambiguous compound {request!r}: {', '.join(candidates)}. "
                         "Set an adduct preference or enter the full channel name.")
    raise KeyError(f"No channel matches compound {request!r}. "
                   "Use standard element capitalization or an exact full channel name.")
