import contextlib
import io
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
import numpy as np
from cloud18_ap_analysis.loading import load_data, read_trace_timeseries

class TestMultipleTraceFiles(unittest.TestCase):
    def test_both_loaders(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = [Path(folder)/name for name in ("first.csv", "second.csv")]
            for path, rows in zip(paths, ([ (1,10), (2,20)], [(2,99), (3,30)])):
                path.write_text("number of header rows: 1\nTime,unixTime,a\n" + "".join(
                    f"1970-01-01 00:00:0{t},{t},{v}\n" for t,v in rows))
            for loader in (lambda: load_data(paths), lambda: read_trace_timeseries(paths, ["a"], (datetime(1970,1,1),datetime(1970,1,2)))):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    result = loader()
                values = result["a"].values if hasattr(result, "channels") else result["a"].to_numpy()
                np.testing.assert_array_equal(values, [10,99,30])
                self.assertIn("WARNING", output.getvalue())
                self.assertIn("1 duplicate cross-file timestamp entries removed", output.getvalue())
                self.assertIn("second.csv", output.getvalue())

    def test_overlap_without_identical_timestamps(self):
        from cloud18_ap_analysis.loading import _report_trace_overlaps
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            _report_trace_overlaps(["a", "b"], [np.array([0,2000]), np.array([1000,3000])])
        self.assertIn("WARNING", output.getvalue())
        self.assertIn("0 identical timestamps", output.getvalue())

    def test_no_overlap(self):
        from cloud18_ap_analysis.loading import _report_trace_overlaps
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            _report_trace_overlaps(["a", "b"], [np.array([0,1000]), np.array([2000,3000])])
        self.assertIn("no overlapping", output.getvalue())

    def test_empty_lists(self):
        with self.assertRaises(ValueError):
            load_data([])
        with self.assertRaises(ValueError):
            read_trace_timeseries([], ["a"], (datetime(1970,1,1),datetime(1970,1,2)))
