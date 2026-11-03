import calendar
import importlib
import unittest
from datetime import datetime, timedelta
from unittest.mock import patch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from cloud18_ap_analysis import Channel, Dataset
from cloud18_ap_analysis.runs import RUNS, make_temperature_runs
mod = importlib.import_module("cloud18_ap_analysis.temp_dependence")

class TestRunGroups(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_registry_and_independent_inputs(self):
        self.assertEqual(set(RUNS), {"m50_1", "m25", "m10", "m50_2"})
        runs = make_temperature_runs(["m50_1", "m50_2"])
        self.assertEqual([r["plot_temperature_c"] for r in runs], [-50,-50])
        runs[0]["trace_path"].clear()
        self.assertTrue(RUNS["m50_1"]["trace_paths"])

    def test_validation(self):
        runs = make_temperature_runs(list(RUNS))
        for groups in ({"mixed": ["m50_1","m25"]}, {"bad": ["m50_1","missing"]},
                       {"bad": ["m50_1","m50_1"]}, {"bad": ["m50_1"]}):
            with self.assertRaises(ValueError):
                mod._prepare_run_groups(runs, groups)

    def test_pooled_fit_and_separate_traces(self):
        dates = [datetime(2025,1,1)+timedelta(minutes=i) for i in range(6)]
        times = [calendar.timegm(d.timetuple())*1000 for d in dates]
        data = pd.DataFrame({"Time": dates, "A": [1,2,3,4,5,6], "B": [3,5,7,9,11,13]})
        partner = Dataset({"HO2_ppt": Channel("HO2_ppt", times, [2,3,4,5,6,7])})
        runs = [dict(run_id=rid, temperature=rid, plot_temperature_c=-50, trace_path="unused",
                     time_range=(dates[0],dates[-1]), important_time_ranges=[(dates[i],dates[i+1]) for i in (0,2,4)])
                for rid in ("m50_1", "m50_2")]
        for groups, expected in ((None,2), ({"combined": ["m50_1","m50_2"]},1)):
            with patch.object(mod,"read_trace_timeseries",return_value=data), patch.object(mod,"load_data",return_value=partner):
                traces, means, fits, fig, ax, kt, kax = mod.temp_dep(
                    runs,["A","B"],"unused","A","B",adduct="",show=False,stages=None,
                    label_intervals=True,combine_runs=groups)
            self.assertEqual(len(traces),2)
            self.assertEqual(len(means),6)
            self.assertEqual(len(fits),expected)
            self.assertEqual(len(kax.containers),expected)
            self.assertTrue((fits.temperature_c == -50).all())
            self.assertEqual(set(means.run_id), {"m50_1","m50_2"})
            if groups:
                self.assertEqual([t.get_text() for t in ax.texts],
                                 [rf"$t^{{{run}}}_{{{interval}}}$"
                                  for run in (1, 2) for interval in (1, 2, 3)])
                self.assertEqual(fits.iloc[0]["n_points"], 6)
            else:
                self.assertEqual(list(fits.n_points), [3, 3])
            self.assertNotIn("fit_group",runs[0])
