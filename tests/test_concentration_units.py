import calendar
import importlib
import unittest
from datetime import datetime, timedelta
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cloud18_ap_analysis import Channel, Dataset
from cloud18_ap_analysis import plotting
from cloud18_ap_analysis.units import (
    convert_channel_concentration, convert_dataset_concentrations,
    concentration_yscale, GAS_CONSTANT, AVOGADRO_CONSTANT,
)

temp = importlib.import_module("cloud18_ap_analysis.temp_dependence")


class TestConcentrationUnits(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_conversion_and_copy(self):
        ch = Channel("O3", [0, 1000], [1, 2])
        data = Dataset({"O3": ch, "AP": Channel("AP", ch.time, [3, 4])})
        converted = convert_dataset_concentrations(
            data, {"O3": "ppbv"}, "molecules/cm3", temperature_c=0,
        )
        expected = ch.values * 1e-9 * 101325 / (GAS_CONSTANT * 273.15) * 1e-6 * AVOGADRO_CONSTANT
        np.testing.assert_allclose(converted["O3"].values, expected, atol=0)
        np.testing.assert_array_equal(data["O3"].values, [1, 2])
        np.testing.assert_array_equal(converted["AP"].values, [3, 4])
        ap = convert_channel_concentration(data["AP"], "ppbv", "pptv")
        np.testing.assert_allclose(ap.values, [3000, 4000])
        with self.assertRaises(ValueError):
            convert_channel_concentration(ch, "ppbv", "molecules/cm3", temperature_c=-274)

    def test_known_number_density(self):
        ch = Channel("HO2", [0], [1])
        out = convert_channel_concentration(ch, "pptv", "molecules/cm3", temperature_c=-25)
        # Independent ideal-gas check: 1 pptv at -25 C and 1 atm.
        np.testing.assert_allclose(out.values, [2.9574611627553694e7], rtol=1e-10)

    def test_temperature_alignment(self):
        ch = Channel("x", [-1, 0, 500, 1000, 1001], [1] * 5)
        temperature = Channel("T", [0, 1000], [0, 20])
        out = convert_channel_concentration(ch, "pptv", "molecules/cm3", temperature_c=temperature)
        self.assertTrue(np.isnan(out.values[[0, -1]]).all())
        np.testing.assert_allclose(out.values[1:4],
            1e-12 * 101325 / (GAS_CONSTANT * np.array([273.15, 283.15, 293.15])) * 1e-6 * AVOGADRO_CONSTANT, atol=0)

    def test_panel_scales_and_override(self):
        self.assertEqual(concentration_yscale("native"), "linear")
        self.assertEqual(concentration_yscale("molecules/cm3"), "log")
        self.assertEqual(concentration_yscale("molecules/cm3", "linear"), "linear")
        ds = Dataset({"x": Channel("x", [0,1000], [1,10])})
        with patch.object(plotting, "_configure_save_and_present") as present:
            plotting.plot_data([ds, ds], channels=[None, (["x"], ["x"])],
                               yscale=["linear", "log"])
        fig = present.call_args.args[0]
        self.assertEqual([ax.get_yscale() for ax in fig.axes], ["linear", "log", "log"])
        fig.canvas.draw()

    def test_pipeline_converts_before_averaging(self):
        dates = [datetime(2025, 1, 1) + timedelta(minutes=i) for i in range(6)]
        ms = [calendar.timegm(d.timetuple()) * 1000 for d in dates]
        ranges = [(dates[i], dates[i+1]) for i in (0,2,4)]
        data = pd.DataFrame({"Time": dates, "A": [1,2,3,4,5,6], "B": [3,5,7,9,11,13]})
        temperature = Channel("TE_Calib_3", ms, [-50,-40,-30,-20,-10,0])
        partner = Dataset({"HO2_ppt": Channel("HO2_ppt", ms, [2,3,4,5,6,7])})
        ap = Dataset({"AP": Channel("AP", ms, [1,3,5,7,9,11])})
        runs = [dict(temperature="test", trace_path="trace", temperature_path="temp",
                     time_range=(dates[0], dates[-1]), important_time_ranges=ranges)]
        loads = {"temp": Dataset({"TE_Calib_3": temperature}), "horus": partner, "ap": ap}
        factor = 1e-12 * 101325 / (GAS_CONSTANT * (temperature.values + 273.15)) * 1e-6 * AVOGADRO_CONSTANT
        for mode in ("native", "molecules/cm3"):
            with patch.object(temp, "read_trace_timeseries", return_value=data), patch.object(
                temp, "load_data", side_effect=lambda path, *args: loads[path]
            ):
                traces, means, fits, fig, ax, kt_fig, kt_ax = temp.temp_dep(
                    runs, ["A", "B"], "horus", "A", "B", adduct="", show=False,
                    stages=None, label_intervals=True, color_by_ap=True, ap_path="ap",
                    concentration_unit=mode,
                )
            f = factor if mode == "molecules/cm3" else np.ones(6)
            expected = (data.A.to_numpy() * f).reshape(3,2)
            np.testing.assert_allclose(means.A, expected.mean(axis=1), atol=0)
            np.testing.assert_allclose(means.A_sem, expected.std(axis=1, ddof=1) / np.sqrt(2), atol=0)
            expected_ho2 = (partner["HO2_ppt"].values * f).reshape(3,2).mean(axis=1)
            np.testing.assert_allclose(means["A * HO2_ppt"], expected.mean(axis=1) * expected_ho2, atol=0)
            np.testing.assert_allclose(means.AP_pptv, [2000,6000,10000])
            self.assertEqual(traces[0][1][0].get_yscale(), concentration_yscale(mode))
            self.assertEqual(ax.get_yscale(), concentration_yscale(mode))
            self.assertEqual(ax.get_xscale(), concentration_yscale(mode))
            self.assertEqual(len(fits), 1)
            legend_text = " ".join(t.get_text() for t in ax.get_legend().get_texts())
            self.assertIn("k=(", legend_text)
            self.assertIn(r"\pm", legend_text)
            self.assertIn("molecule^{-1}" if mode == "molecules/cm3" else "pptv^{-1}", legend_text)
            fit = fits.iloc[0]
            scale = 10.0 ** int(np.floor(np.log10(max(abs(fit.slope), abs(fit.slope_stderr)))))
            self.assertIn(f"{fit.slope / scale:.3g}", legend_text)
            self.assertIn(f"{fit.slope_stderr / scale:.2g}", legend_text)
            self.assertEqual(len(ax.texts), 3)
            self.assertIn("pptv", fig.axes[1].get_ylabel())
            fig.canvas.draw()
        np.testing.assert_array_equal(data.A, [1,2,3,4,5,6])

    def test_odr_unit_equivalence(self):
        x = np.array([1.,3.,7.,12.,18.])
        y = np.array([2.1,5.7,14.5,23.2,37.])
        reference = None
        for celsius in (-50, -25, 0):
            factor = 1e-12 * 101325 / (GAS_CONSTANT * (273.15+celsius)) * 1e-6 * AVOGADRO_CONSTANT
            for f in (1., factor):
                means = pd.DataFrame(dict(temperature=["test"]*5, x=x*f*f, y=y*f,
                                          x_sem=np.full(5,.1)*f*f, y_sem=np.full(5,.2)*f))
                _, ax, fits = temp.plot_kinetic_relationship(
                    means, "x", "y", show=False,
                    concentration_unit="native" if f == 1 else "molecules/cm3",
                )
                self.assertEqual(ax.get_xscale(), "linear" if f == 1 else "log")
                self.assertEqual(fits.slope_unit.iloc[0], "1/pptv" if f == 1 else "cm3/molecule")
                result = np.array([fits.slope.iloc[0], fits.slope_stderr.iloc[0]]) * f
                if reference is None:
                    reference = result
                np.testing.assert_allclose(result, reference, rtol=1e-4, atol=0)


if __name__ == "__main__":
    unittest.main()
