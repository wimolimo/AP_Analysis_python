import numpy as np
from .models import Channel


def resample(ch: Channel, new_times):
    """Linearly interpolate a channel onto ``new_times``.

    Values outside the original time range are NaN, matching Julia's
    ``extrapolation_bc=NaN`` behavior.
    """
    new_times = np.asarray(new_times, dtype=float)
    order = np.argsort(ch.time)
    x = ch.time[order]
    y = ch.values[order]

    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 2:
        values = np.full(new_times.shape, np.nan, dtype=float)
    else:
        values = np.interp(new_times, x[valid], y[valid])
        outside = (new_times < x[valid][0]) | (new_times > x[valid][-1])
        values[outside] = np.nan
    return Channel(ch.name, new_times, values)
