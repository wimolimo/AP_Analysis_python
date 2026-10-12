from datetime import datetime, timedelta
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import unittest

from cloud18_ap_analysis import Channel, Dataset
from cloud18_ap_analysis import plotting
from cloud18_ap_analysis.time_traces import plot_time_traces


def check_trace_scales(scale, mode):
    a = Channel('a', [0, 1000, 2000, 3000], [-1, 0, 1, 1000])
    b = Channel('b', a.time, [2, 3, 4, 5])
    data = Dataset({'a': a, 'b': b})
    figures = []
    with patch.object(plotting, '_configure_save_and_present', side_effect=lambda fig, *args: figures.append(fig)):
        if mode == 'single':
            plotting.plot_channel(a, yscale=scale)
        elif mode == 'dual':
            plotting.plot_data([data], channels=[(['a'], ['b'])], yscale=scale)
        elif mode == 'datasets':
            plotting.plot_data([data, data], yscale=scale)
        else:
            plotting.plot_data(data, stacked=mode == 'stacked', yscale=scale)
    try:
        for fig in figures:
            fig.canvas.draw()
            for ax in fig.axes:
                assert ax.get_yscale() == scale
                if scale == 'log':
                    assert 0 < ax.get_ylim()[0] < ax.get_ylim()[1]
            assert figures[0].axes[0].get_ylim()[0] < 1
            assert figures[0].axes[0].get_ylim()[1] > 1000
    finally:
        for fig in figures:
            plt.close(fig)


def check_time_trace_forwards_scale(scale):
    data = Dataset({'a': Channel('a', [0, 1000], [1, 10])})
    with patch('cloud18_ap_analysis.time_traces.load_data', return_value=data), patch('cloud18_ap_analysis.time_traces.plot_data') as plot:
        plot_time_traces('unused.csv', ['a'], (datetime(1970, 1, 1), datetime(1970, 1, 2)), yscale=scale)
    assert plot.call_args.kwargs['yscale'] == scale


def check_log_limits_use_visible_positive_values():
    ch = Channel('a', [0, 1000, 2000, 3000, 4000], [1e9, -1, 0, 1, 10])
    limits = plotting._compute_ylims([ch], (datetime(1970, 1, 1) + timedelta(seconds=1), datetime(1970, 1, 1) + timedelta(seconds=4)), yscale='log')
    assert 0 < limits[0] < 1
    assert 10 < limits[1] < 20


class TestPlottingScales(unittest.TestCase):
    def test_layouts(self):
        for scale in ('linear', 'log'):
            for mode in ('overlay', 'stacked', 'datasets', 'dual', 'single'):
                with self.subTest(scale=scale, mode=mode):
                    check_trace_scales(scale, mode)

    def test_wrapper(self):
        for scale in ('linear', 'log'):
            with self.subTest(scale=scale):
                check_time_trace_forwards_scale(scale)

    def test_visible_limits(self):
        check_log_limits_use_visible_positive_values()
