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
