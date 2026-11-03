from datetime import datetime, timedelta
import unittest
from unittest.mock import patch
import numpy as np
from cloud18_ap_analysis import Channel, Dataset
from cloud18_ap_analysis.time_traces import plot_time_traces

class TestTimeTraceUnits(unittest.TestCase):
    def test_conversion_and_native(self):
        data = Dataset({"a": Channel("a", [0,1000], [1,2])})
        temp = Dataset({"TE_Calib_3": Channel("TE_Calib_3", [0,1000], [0,0])})
        window = (datetime(1970,1,1), datetime(1970,1,1)+timedelta(seconds=1))
        for unit, scale in (("native", "linear"), ("molecules/cm3", "log")):
            with patch("cloud18_ap_analysis.time_traces.load_data", side_effect=lambda path,*args: temp if path == "temperature" else data) as load, patch("cloud18_ap_analysis.time_traces.plot_data") as plot:
                plot_time_traces("trace", ["a"], window, concentration_unit=unit, temperature_path="temperature")
                result = plot.call_args.args[0]["a"].values
                np.testing.assert_allclose(result, np.array([1,2]) * (2.686780111e7 if unit != "native" else 1), rtol=1e-8)
                self.assertEqual(plot.call_args.kwargs["yscale"], scale)
                self.assertIn("molecules" if unit != "native" else "pptv", plot.call_args.kwargs["ylabels"])
                self.assertEqual(load.call_count, 2 if unit != "native" else 1)
                np.testing.assert_array_equal(data["a"].values, [1,2])

    def test_missing_temperature(self):
        data = Dataset({"a": Channel("a", [0,1000], [1,2])})
        with patch("cloud18_ap_analysis.time_traces.load_data", return_value=data):
            with self.assertRaisesRegex(ValueError, "temperature_c is required"):
                plot_time_traces("trace", ["a"], (datetime(1970,1,1),datetime(1970,1,2)), concentration_unit="molecules/cm3")
