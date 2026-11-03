from datetime import datetime, timedelta
import numpy as np

from cloud18_ap_analysis import Channel, Dataset, resample
from cloud18_ap_analysis.plotting import _smooth_channel


def test_dataset_indexing():
    ch = Channel("x", [0, 1, 2], [1, 2, 3])
    ds = Dataset({"x": ch})
    assert ds["x"].name == "x"


def test_resample():
    ch = Channel("x", [0, 10], [0, 10])
    out = resample(ch, [0, 5, 10, 20])
    assert np.allclose(out.values[:3], [0, 5, 10])
    assert np.isnan(out.values[3])


def test_centered_smoothing():
    ch = Channel("x", [0, 1, 2], [0, 10, 20])
    out = _smooth_channel(ch, timedelta(seconds=2))
    assert np.allclose(out.values, [5, 10, 15])


def test_resample_gap_limit():
    minute = 60_000
    # Unsorted input, a 30-minute gap, and a 31-minute gap with a missing sample.
    ch = Channel("x", np.array([61, 0, 45, 30]) * minute,
                 [61, 0, np.nan, 30])
    out = resample(ch, np.array([-1, 0, 15, 30, 45, 61, 62]) * minute)
    np.testing.assert_allclose(out.values, [np.nan, 0, 15, 30, np.nan, 61, np.nan],
                               equal_nan=True)
