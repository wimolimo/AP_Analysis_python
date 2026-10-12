import calendar
import importlib
import unittest
from datetime import datetime, timedelta
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd

from cloud18_ap_analysis import Channel, Dataset

temp = importlib.import_module("cloud18_ap_analysis.temp_dependence")


class TestShowSelection(unittest.TestCase):
    def tearDown(self):
        plt.close("all")
        plt.ioff()

    def test_selections_preserve_results_and_saving(self):
        dates = [datetime(2025, 1, 1) + timedelta(minutes=i) for i in range(6)]
        times = [calendar.timegm(d.timetuple()) * 1000 for d in dates]
        data = pd.DataFrame({"Time": dates, "A": [1,2,3,4,5,6], "B": [3,5,7,9,11,13]})
        partner = Dataset({"HO2_ppt": Channel("HO2_ppt", times, [2,3,4,5,6,7])})
        runs = [dict(temperature=t, plot_temperature_c=-50 if t == "cold" else -25, trace_path="trace", time_range=(dates[0], dates[-1]),
                     important_time_ranges=[(dates[i], dates[i+1]) for i in (0,2,4)],
                     trace_savepath=t + ".png") for t in ("cold", "warm")]
        baseline = None
        cases = [(False, set()), ([], set()), (True, {"traces", "kinetics", "k(T)"}),
                 ("k(T)", {"k(T)"}), ("traces", {"traces"}), ("kinetics", {"kinetics"}),
                 (["traces", "kinetics"], {"traces", "kinetics"}),
                 (("kinetics",), {"kinetics"})]
        for selection, expected in cases:
            with self.subTest(show=selection):
                plt.close("all")
                plt.ion()  # Must not automatically expose unselected figures.
                snapshots = []
                def capture_show():
                    snapshots.append(set(plt.get_fignums()))
                    self.assertFalse(plt.isinteractive())
                with patch.object(temp, "read_trace_timeseries", return_value=data), patch.object(
                    temp, "load_data", return_value=partner
                ), patch.object(plt, "show", side_effect=capture_show), patch.object(Figure, "savefig") as save:
                    traces, means, fits, fig, ax, kt_fig, kt_ax = temp.temp_dep(
                        runs, ["A", "B"], "horus", "A", "B", adduct="",
                        show=selection, stages=None, kinetic_savepath="kinetics.png",
                    )
                self.assertTrue(plt.isinteractive())
                self.assertEqual(save.call_count, 3)
                self.assertEqual(len(traces), 2)
                self.assertEqual(len(fits), 2)
                if baseline is None:
                    baseline = (means, fits)
                else:
                    pd.testing.assert_frame_equal(means, baseline[0])
                    pd.testing.assert_frame_equal(fits, baseline[1])
                selected_figures = []
                if "traces" in expected:
                    selected_figures.extend(f for f, _ in traces)
                if "kinetics" in expected:
                    selected_figures.append(fig)
                if "k(T)" in expected:
                    selected_figures.append(kt_fig)
                self.assertEqual(len(snapshots), len(selected_figures))
                selected_ids = {f.number for f in selected_figures}
                self.assertTrue(all(snapshot <= selected_ids for snapshot in snapshots))
                self.assertEqual(set(plt.get_fignums()), selected_ids)
                # Unregistered figures are still usable after returning.
                fig.canvas.draw()
                kt_fig.canvas.draw()
                self.assertEqual(list(fits.temperature_c), [-50, -25])
                self.assertEqual(len(kt_ax.containers), 2)
                for container, row in zip(kt_ax.containers, fits.itertuples()):
                    self.assertEqual(container.lines[0].get_xdata()[0], row.temperature_c)
                    self.assertEqual(container.lines[0].get_ydata()[0], row.slope)
                for trace, _ in traces:
                    trace.canvas.draw()

    def test_invalid_selection_fails_before_loading(self):
        for selection, error in [("kinetic", ValueError), (["temperature"], ValueError),
                                 (None, TypeError), (1, TypeError), (["traces", False], TypeError)]:
            with self.subTest(show=selection), patch.object(temp, "read_trace_timeseries") as load:
                with self.assertRaises(error):
                    temp.temp_dep([], [], "unused", "A", "B", show=selection)
                load.assert_not_called()
