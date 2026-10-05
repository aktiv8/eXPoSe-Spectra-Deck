"""A saved figure's zoom: encode / sanitise / the zoomed test (Tk-free).

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import viewzoom as vz  # noqa: E402


class TestEncode(unittest.TestCase):
    def test_round_trip_keeps_the_drawn_orientation(self):
        z = vz.encode({"C 1s": ((292.0, 280.0), (0.0, 900.0))}, "Binding")
        self.assertEqual(z, {"scale": "Binding", "panels": {
            "C 1s": {"x": [292.0, 280.0], "y": [0.0, 900.0]}}})
        self.assertEqual(vz.sanitise(z, "Binding"),
                         {"C 1s": {"x": [292.0, 280.0], "y": [0.0, 900.0]}})

    def test_nothing_zoomed_is_an_empty_zoom_not_a_missing_one(self):
        z = vz.encode({}, "Binding")
        self.assertEqual(vz.sanitise(z, "Binding"), {})          # full range
        self.assertIsNone(vz.sanitise(None, "Binding"))          # no key at all

    def test_bad_rows_are_dropped(self):
        z = {"scale": "Binding", "panels": {
            "ok": {"x": [1, 2], "y": [3, 4]},
            "nan": {"x": [float("nan"), 2], "y": [3, 4]},
            "inf": {"x": [1, float("inf")], "y": [3, 4]},
            "flat": {"x": [5, 5], "y": [3, 4]},
            "short": {"x": [1], "y": [3, 4]},
            "str": {"x": ["a", "b"], "y": [3, 4]},
            "bool": {"x": [True, False], "y": [3, 4]},
            "notdict": 7, "nokey": {"y": [3, 4]}}}
        self.assertEqual(list(vz.sanitise(z, "Binding")), ["ok"])

    def test_malformed_zoom_means_none(self):
        for bad in ("x", 3, [], {"scale": "Binding"},
                    {"scale": "Binding", "panels": []}):
            self.assertIsNone(vz.sanitise(bad, "Binding"), bad)

    def test_a_zoom_taken_under_another_scale_is_not_used(self):
        z = vz.encode({"C 1s": ((292.0, 280.0), (0.0, 900.0))}, "Binding")
        self.assertEqual(vz.sanitise(z, "Kinetic"), {})


class TestIsZoomed(unittest.TestCase):
    def test_same_limits_are_not_a_zoom(self):
        auto = ((1200.0, 0.0), (0.0, 1000.0))
        self.assertFalse(vz.is_zoomed(auto, ((1200.0, 0.0), (0.0, 1000.0))))
        self.assertFalse(vz.is_zoomed(auto, ((1200.0 + 1e-9, 0.0), (0.0, 1000.0))))

    def test_a_changed_limit_is_a_zoom(self):
        auto = ((1200.0, 0.0), (0.0, 1000.0))
        self.assertTrue(vz.is_zoomed(auto, ((300.0, 280.0), (0.0, 1000.0))))
        self.assertTrue(vz.is_zoomed(auto, ((1200.0, 0.0), (0.0, 400.0))))
        self.assertTrue(vz.is_zoomed(auto, ((0.0, 1200.0), (0.0, 1000.0))))


class TestFits(unittest.TestCase):
    def test_overlap_is_enough_whichever_way_round(self):
        self.assertTrue(vz.fits((1200.0, 0.0), [300.0, 280.0]))
        self.assertTrue(vz.fits((0.0, 1200.0), [280.0, 300.0]))
        self.assertTrue(vz.fits((300.0, 280.0), [310.0, 290.0]))   # partly

    def test_a_range_outside_the_data_is_not_used(self):
        self.assertFalse(vz.fits((300.0, 280.0), [900.0, 880.0]))
        self.assertFalse(vz.fits((300.0, 280.0), [280.0, 270.0]))   # touching only


if __name__ == "__main__":
    unittest.main()
