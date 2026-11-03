import numpy as np
from .models import Channel


def resample(ch: Channel, new_times):
    """Linearly interpolate a channel onto ``new_times``.

    Values outside the original time range or strictly inside gaps longer
    than 30 minutes between finite measurements are NaN. Measured endpoints
    are retained, and gaps of exactly 30 minutes can be interpolated.
    """
    new_times = np.asarray(new_times, dtype=float)
    order = np.argsort(ch.time)
    x = ch.time[order]
    y = ch.values[order]

    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 2:
        values = np.full(new_times.shape, np.nan, dtype=float)
    else:
        x = x[valid]
        y = y[valid]
        values = np.interp(new_times, x, y)
        outside = (new_times < x[0]) | (new_times > x[-1])
        values[outside] = np.nan
        # Channel timestamps are milliseconds. Check the surrounding finite
        # samples, so missing values cannot hide a long measurement gap.
        right = np.searchsorted(x, new_times, side="left")
        left = np.clip(right - 1, 0, len(x) - 1)
        right = np.clip(right, 0, len(x) - 1)
        long_gap = ((x[right] - x[left] > 30 * 60 * 1000)
                    & (new_times > x[left]) & (new_times < x[right]))
        values[long_gap] = np.nan
    return Channel(ch.name, new_times, values)
