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


class TestFitY(unittest.TestCase):
    NAN = float("nan")

    def test_only_the_data_inside_the_window_count(self):
        xs = [0, 1, 2, 3, 4]
        ys = [100, 5, 7, 3, 200]
        lo, hi = vz.fit_y([(xs, ys)], (1, 3), pad=0.1)
        self.assertAlmostEqual(lo, 3 - 0.4)
        self.assertAlmostEqual(hi, 7 + 0.4)

    def test_window_given_high_to_low_and_several_lines(self):
        a = ([0, 1, 2], [1, 2, 3])
        b = ([0, 1, 2], [10, 20, 30])
        lo, hi = vz.fit_y([a, b], (2, 0), pad=0)
        self.assertEqual((lo, hi), (1, 30))

    def test_nan_outside_a_fit_region_is_ignored(self):
        lo, hi = vz.fit_y([([0, 1, 2], [self.NAN, 4, self.NAN])], (0, 2),
                          pad=0.5)
        self.assertLess(lo, 4)
        self.assertGreater(hi, 4)

    def test_flat_data_still_gets_a_range(self):
        lo, hi = vz.fit_y([([0, 1], [10, 10])], (0, 1))
        self.assertLess(lo, 10)
        self.assertGreater(hi, 10)

    def test_empty_window_gives_none(self):
        self.assertIsNone(vz.fit_y([([0, 1], [1, 2])], (5, 6)))
        self.assertIsNone(vz.fit_y([], (0, 1)))

    def test_log_axis_ignores_non_positive_and_pads_by_a_factor(self):
        lo, hi = vz.fit_y([([0, 1, 2], [0, 10, 1000])], (0, 2), pad=0.1,
                          log=True)
        self.assertLess(lo, 10)
        self.assertGreater(lo, 0)
        self.assertGreater(hi, 1000)


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
