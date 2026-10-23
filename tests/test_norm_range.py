"""The "Range" normalisation: the lower of two end-point intensities is 0 and
the higher 1; the data's own ends until the user chooses two energies.

Run:  python -m unittest tests.test_norm_range
"""

import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import panelview  # noqa: E402
import plots  # noqa: E402
from readers import Region  # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    HAVE_MPL = True
except Exception:                                   # pragma: no cover
    HAVE_MPL = False

XS = [300.0, 299.0, 298.0, 297.0, 296.0, 295.0, 294.0]     # descending (BE)
YS = [10.0, 11.0, 30.0, 80.0, 40.0, 14.0, 12.0]


def region(xs=XS, ys=YS):
    return Region(name="x", index=0, offset=0, energy=list(xs),
                  counts=list(ys))


class TestNormAffine(unittest.TestCase):
    def test_other_modes_only_divide(self):
        r = region()
        for mode in ("None", "Max = 1", "Area = 1"):
            off, div = plots.norm_affine(r, mode)
            self.assertEqual(off, 0.0)
            self.assertEqual(div, plots.norm_factor(r, mode))

    def test_automatic_ends_are_the_mean_of_the_outer_three(self):
        lo_end = (10 + 11 + 30) / 3          # 300, 299, 298
        hi_end = (40 + 14 + 12) / 3          # 296, 295, 294
        off, div = plots.norm_affine(region(), "Range")
        self.assertAlmostEqual(off, min(lo_end, hi_end))
        self.assertAlmostEqual(div, abs(lo_end - hi_end))

    def test_lower_point_is_zero_and_higher_is_one_on_either_side(self):
        for ys in (YS, YS[::-1]):                 # step up / step down
            r = region(XS, ys)
            off, div = plots.norm_affine(r, "Range", rng=(299.0, 295.0))
            ends = plots.range_ends(r, (299.0, 295.0))
            self.assertAlmostEqual((min(ends) - off) / div, 0.0)
            self.assertAlmostEqual((max(ends) - off) / div, 1.0)

    def test_chosen_energies_are_interpolated_and_clamped(self):
        r = region()
        self.assertEqual(plots.range_ends(r, (298.5, 295.5)), (20.5, 27.0))
        self.assertEqual(plots.range_ends(r, (500.0, 0.0)), (10.0, 12.0))

    def test_ascending_axis_gives_the_same_answer(self):
        a = plots.norm_affine(region(), "Range", rng=(298.5, 294.5))
        b = plots.norm_affine(region(XS[::-1], YS[::-1]), "Range",
                              rng=(298.5, 294.5))
        self.assertEqual(a, b)

    def test_equal_or_unreadable_ends_leave_the_spectrum_alone(self):
        self.assertEqual(plots.norm_affine(region(XS[:4], [5.0] * 4),
                                           "Range"), (0.0, 1.0))
        self.assertEqual(plots.norm_affine(region([], []), "Range"),
                         (0.0, 1.0))

    def test_panel_views_know_the_mode(self):
        self.assertIn("Range", panelview.NORMS)
        self.assertEqual(panelview.sanitise({"norm": "Range"}),
                         {"norm": "Range"})


@unittest.skipUnless(HAVE_MPL, "matplotlib not available")
class TestRangeDrawing(unittest.TestCase):
    def test_a_fit_follows_the_range_normalisation(self):
        import casafit
        from test_casafit import fitted_region
        r = fitted_region(n=61)
        cv = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                            *r.dwell_and_scans())
        off, div = plots.norm_affine(r, "Range")
        ax = Figure().add_subplot()
        fit = {"curves": cv, "colours": ["#111", "#c00", "#0a0", "#00c"],
               "show": {"envelope": True, "components": True,
                        "background": True}}
        plots.draw_stack(ax, [r], norm="Range", fit=fit)
        data = [(y - off) / div for y in r.counts]
        i = next(k for k, v in enumerate(cv[0].envelope) if v == v)
        env = (cv[0].envelope[i] - off) / div
        lines = [list(line.get_ydata()) for line in ax.lines]
        self.assertTrue(any(math.isclose(v[0], data[0]) for v in lines))
        self.assertTrue(any(math.isclose(v[i], env, abs_tol=1e-9)
                            for v in lines[1:]))
        # the two ends of the data are 0 and 1 (in some order)
        self.assertAlmostEqual(min(data[0], data[-1]), 0.0, delta=0.1)

    def test_chosen_end_points_are_drawn_and_other_views_run(self):
        ax = Figure().add_subplot()
        plots.draw_stack(ax, [region()], norm="Range", rng=(299.0, 295.0))
        self.assertEqual(len([ln for ln in ax.lines
                              if list(ln.get_xdata()) in ([299.0, 299.0],
                                                          [295.0, 295.0])]),
                         2)
        fig = Figure()
        zi = plots.viewdata.ZInfo([0.0, 1.0], "Trace order", "Trace order")
        plots.draw_heatmap(fig, fig.add_subplot(), [region(), region()], zi,
                           "Range")


if __name__ == "__main__":
    unittest.main()
