from datetime import datetime, timedelta
import unittest
from unittest.mock import patch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cloud18_ap_analysis import Channel, Dataset, plotting

class TestIntervalAnnotations(unittest.TestCase):
    def setUp(self):
        self.start = datetime(1970, 1, 1)
        self.window = (self.start, self.start + timedelta(seconds=3))
        self.ranges = [(self.start-timedelta(seconds=3), self.start-timedelta(seconds=2)),
                       (self.start+timedelta(seconds=1), self.start+timedelta(seconds=2))]
        self.data = Dataset({name: Channel(name, [0,1000,2000,3000], [1,2,3,4]) for name in ("a", "b")})
        self.fan = Dataset({"Fan_BOT_Speed": Channel("Fan_BOT_Speed", [0,1000,2000,3000], [12,12,100,100])})

    def tearDown(self):
        plt.close("all")

    def test_context_and_numbering(self):
        for ranges, expected_height in ((self.ranges, .035), (None, 1.0)):
            fig, ax = plt.subplots()
            plotting._draw_time_context(ax, ranges, True, self.window, self.fan)
            fan = [p for p in ax.patches if p.get_gid().startswith("fan-")][0]
            self.assertAlmostEqual(fan.get_height(), expected_height)
            plotting._add_fan_legend(ax)
            self.assertIn("Fan speed 12%", [t.get_text() for t in ax.get_legend().texts])
            if ranges:
                self.assertEqual([t.get_text() for t in ax.texts], [r"$t_{2}$"])
                self.assertEqual(len(ax.patches), 2)

    def test_heatmap_labels_without_fill(self):
        fig, ax = plt.subplots()
        plotting._draw_time_context(ax, self.ranges, True, self.window, self.fan, heatmap=True)
        self.assertEqual(len(ax.patches), 0)
        self.assertEqual([t.get_text() for t in ax.texts], [r"$t_{2}$"])

    def test_public_plot_paths(self):
        for stacked, multi in ((False, False), (True, False), (False, True)):
            figures = []
            with patch.object(plotting, "_configure_save_and_present", side_effect=lambda fig,*args: figures.append(fig)):
                plotting.plot_data([self.data,self.data] if multi else self.data,
                                   stacked=stacked, xlims=self.window, yscale="linear",
                                   fan_speed_data=self.fan, important_time_ranges=self.ranges,
                                   label_intervals=True)
            fig = figures[0]
            for size in ((12,8), (6,6)):
                fig.set_size_inches(*size)
                fig.canvas.draw()
                for index, ax in enumerate(fig.axes):
                    self.assertEqual([t.get_text() for t in ax.texts], [r"$t_{2}$"] if index == 0 else [])
                    self.assertTrue(any(p.get_gid()=="averaging-interval" for p in ax.patches))
                    self.assertTrue(any(p.get_gid()=="fan-strip" for p in ax.patches))
                    self.assertAlmostEqual(ax.get_position().x0, fig.axes[0].get_position().x0)
                    self.assertAlmostEqual(ax.get_position().x1, fig.axes[0].get_position().x1)
