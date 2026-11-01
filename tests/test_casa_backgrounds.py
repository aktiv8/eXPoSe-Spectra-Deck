"""CasaXPS's own simple backgrounds (``lineshapes.casa_levels``,
``shirley_casa``, ``casa_step``, ``casa_poly_regression`` and the dispatcher's
``Zero`` / ``None`` / ``Min Limits`` / ``Step`` / ``Edge`` / ``Poly
Regression`` types), read from CasaXPS.exe.

``XPS_CASA_BG=<folder holding amy.vms and amy.csv>`` (a Thermo K-Alpha+ run
fitted once with every type, then "Export All to ASCII") compares every one
with the Background column CasaXPS exported."""
import csv
import os
import unittest

import casafit
import lineshapes as ls
import readers

try:
    import numpy as np
    HAVE_NP = True
except ImportError:                                        # pragma: no cover
    HAVE_NP = False

FOLDER = os.environ.get("XPS_CASA_BG")


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestEndLevels(unittest.TestCase):
    Y = np.arange(1.0, 21.0) if HAVE_NP else None          # 1..20

    def test_single_point_when_avg_is_one(self):
        self.assertEqual(ls.casa_levels(self.Y[5:15], 1), (6.0, 15.0))

    def test_symmetric_window_reaches_outside_the_region(self):
        # avg 3 = +-2 points: the window at the start of y[5:15] is 4..8 (mean 6)
        # and at the end 13..17 (mean 15) -- points outside the region count
        lo, hi = ls.casa_levels(self.Y[5:15], 3, self.Y[3:5], self.Y[15:17])
        self.assertAlmostEqual(lo, 6.0)
        self.assertAlmostEqual(hi, 15.0)
        # a lopsided spectrum shows the window is not one-sided
        y = np.array([100.0, 0, 0, 0, 0, 0, 0, 0, 0, 100.0])
        lo, hi = ls.casa_levels(y[2:8], 2, y[1:2], y[8:9])
        self.assertAlmostEqual(lo, 0.0)

    def test_all_or_nothing(self):
        # one end cannot be given its window: both ends are single points
        self.assertEqual(ls.casa_levels(self.Y[0:10], 3, None, self.Y[10:12]),
                         (1.0, 10.0))
        self.assertEqual(ls.casa_levels(self.Y[0:10], 3, self.Y[0:0], None),
                         (1.0, 10.0))

    def test_offsets_scale_the_levels(self):
        lo, hi = ls.casa_levels(self.Y[5:15], 1, st_offset=10.0, end_offset=50.0)
        self.assertAlmostEqual(lo, 6.0 * 0.9)
        self.assertAlmostEqual(hi, 15.0 * 0.5)


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestSimpleTypes(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(1000.0, 1010.0, 101)
        self.y = 100.0 + 20.0 * np.exp(-((self.x - 1005.0) / 1.0) ** 2)

    def test_shirley_ends_on_the_levels(self):
        b = ls.shirley_casa(self.y, 90.0, 120.0)
        self.assertAlmostEqual(b[-1], 120.0)
        self.assertLess(abs(b[0] - 90.0), 0.5)
        self.assertTrue(np.all(np.diff(b) >= -1e-9))

    def test_min_limits_is_flat_at_the_smaller_end(self):
        y = self.y.copy()
        y[0], y[-1] = 80.0, 150.0
        self.assertTrue(np.allclose(ls.background("Min Limits", y, 1), 80.0))
        y[0], y[-1] = 150.0, 80.0
        self.assertTrue(np.allclose(ls.background("Min Limits", y, 1), 80.0))

    def test_none_is_zero_for_other_sources_flat_for_casaxps(self):
        y = self.y.copy()
        y[0], y[-1] = 80.0, 150.0
        self.assertTrue(np.allclose(ls.background("None", y, 1), 0.0))
        self.assertTrue(np.allclose(
            ls.background("None", y, 1, ctx=(None, None)), 80.0))
        self.assertTrue(np.allclose(ls.background("Zero", y, 1), 0.0))

    def test_steps(self):
        p = (0.0, 0.0, 50.0, 1005.0, 0.5, 10.0)
        down = ls.background("Step Down", self.y, 1, x=self.x, params=p)
        up = ls.background("Step Up", self.y, 1, x=self.x, params=p)
        self.assertAlmostEqual(down[0], 60.0, 3)            # height + offset
        self.assertAlmostEqual(down[-1], 10.0, 3)
        self.assertAlmostEqual(up[0], 10.0, 3)
        self.assertAlmostEqual(up[-1], 60.0, 3)
        self.assertAlmostEqual(down[50], 35.0, 3)           # half height at p1
        # the Edge and Fixed forms are the same functions
        self.assertTrue(np.allclose(
            ls.background("Edge Down", self.y, 1, x=self.x, params=p), down))
        # a non-positive height or width is CasaXPS's 1e300 marker: no curve
        self.assertIsNone(ls.background(
            "Step Down", self.y, 1, x=self.x,
            params=(0.0, 0.0, 0.0, 1005.0, 0.5, 10.0)))
        self.assertIsNone(ls.background("Step Down", self.y, 1, params=p))

    def test_poly_regression_recovers_a_polynomial(self):
        t = (np.arange(60) - 20.0) / 59.0
        y = 1000.0 + 500.0 * t - 800.0 * t ** 2 + 300.0 * t ** 3
        b = ls.background("Poly Regression", y, 2, params=(0.0, 0.0, 4.0, 0.0))
        self.assertTrue(np.allclose(b[2:-2], y[2:-2], rtol=1e-9))
        self.assertTrue(np.array_equal(b[:2], y[:2]))        # not computed
        self.assertTrue(np.array_equal(b[-2:], y[-2:]))
        # too few points for the order, or the derivative option: no curve
        self.assertIsNone(ls.background("Poly Regression", y[:4], 1,
                                        params=(0.0, 0.0, 4.0, 0.0)))
        self.assertIsNone(ls.background("Poly Regression", y, 1,
                                        params=(0.0, 0.0, 4.0, 1.0)))

    def test_dispatcher_keeps_the_old_shirley_without_a_casa_context(self):
        a = ls.background("Shirley", self.y, 3)
        b = ls.shirley(self.y, 3)
        self.assertTrue(np.array_equal(a, b))


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestRegionsFromOtherReaders(unittest.TestCase):
    def test_only_a_casa_line_marks_a_region(self):
        reg = casafit._region(
            "CASA region (*C 1s*) (*None*) 1000 1010 1 3 0 0 0 0 0 0 "
            "(*NoTag*) 0 0 0")
        self.assertTrue(reg.casa)
        self.assertFalse(casafit.FitRegion(name="C 1s").casa)

    def test_none_background_of_a_casa_region_is_flat_not_zero(self):
        n = 40
        be = [290.0 - 0.1 * i for i in range(n)]
        counts = [1000.0 + 5.0 * i for i in range(n)]
        hv = 1486.6
        for casa, expect in ((True, 1000.0), (False, 0.0)):
            reg = casafit.FitRegion(name="C 1s", background="None",
                                    start_ke=hv - 290.0 + 0.5,
                                    end_ke=hv - 286.5, casa=casa)
            fit = casafit.Fit(regions=[reg], components=[casafit.FitComponent(
                name="p", shape="GL(30)", area=1.0, fwhm=1.0,
                pos_ke=hv - 288.0, region="C 1s")])
            cv = casafit.curves(fit, be, counts, hv, 1.0, 1)[0]
            bg = np.array([v for v in cv.background if v == v])
            self.assertTrue(np.allclose(bg, bg[0]))
            if casa:
                # the lower of the two end levels (low-KE end = high BE)
                self.assertLess(abs(bg[0] - min(counts[5], counts[34])), 1e-6)
            else:
                self.assertAlmostEqual(bg[0], 0.0)


def _blocks(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    h = rows[6]
    out = []
    for i, c in enumerate(h):
        if c == "Background":
            j = i
            while h[j] != "B.E.":
                j -= 1
            be, bg = [], []
            for r in rows[7:]:
                if len(r) > i and r[i]:
                    be.append(float(r[j]))
                    bg.append(float(r[i]))
            out.append((h[j + 1], np.array(be), np.array(bg)))
    return out


@unittest.skipUnless(HAVE_NP and FOLDER and os.path.isdir(FOLDER),
                     "XPS_CASA_BG not set")
class TestAgainstCasaXpsExport(unittest.TestCase):
    """amy.vms: 21 fitted regions (Step Down x9, Step Up, Edge Up, Edge Down,
    Min Limits x2, Poly Regression x3, Shirley x2, Zero, None), amy.csv the
    ASCII export of the same file."""

    @classmethod
    def setUpClass(cls):
        doc = readers.load_file(os.path.join(FOLDER, "amy.vms"))
        regs = [(r, f) for r in doc.regions if r.fit and r.fit.regions
                for f in r.fit.regions]
        blocks = _blocks(os.path.join(FOLDER, "amy.csv"))
        assert len(regs) == len(blocks)
        cls.cases = []
        for (r, f), (label, be, bg) in zip(regs, blocks):
            dwell, scans = r.dwell_and_scans()
            cvs = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                 dwell, scans)
            cv = next(c for c in cvs if c.fit_region is f)
            k = dwell * scans if dwell else 1.0
            ref = np.interp(np.array(r.energy), be[::-1], bg[::-1])
            mine = np.array(cv.background, dtype=float) / k \
                if cv.background is not None else None
            counts = np.array(r.counts, dtype=float) / k
            cls.cases.append((f.background, f.params, mine, ref,
                              counts.max() - counts.min()))

    def _diff(self, kind):
        out = []
        for name, params, mine, ref, rng in self.cases:
            if name != kind:
                continue
            self.assertIsNotNone(mine, name)
            ok = ~np.isnan(mine)
            out.append((float(np.abs(mine[ok] - ref[ok]).max() / rng), params))
        self.assertTrue(out, kind)
        return out

    def test_step_up_and_down_and_edges(self):
        # erfc is read from a 0.001-step table by CasaXPS: ~0.003 % of range
        for kind in ("Step Down", "Step Up", "Edge Down", "Edge Up"):
            for d, _p in self._diff(kind):
                self.assertLess(d, 1e-4, kind)

    def test_min_limits_zero_none_are_exact(self):
        for kind in ("Min Limits", "Zero", "None"):
            for d, _p in self._diff(kind):
                self.assertLess(d, 1e-8, kind)

    def test_shirley_is_exact_for_both_averaging_widths(self):
        diffs = self._diff("Shirley")
        self.assertEqual(len(diffs), 2)
        for d, _p in diffs:
            self.assertLess(d, 1e-7)

    def test_poly_regression(self):
        for d, p in self._diff("Poly Regression"):
            order = int(p[2])
            # order 14 reproduces CasaXPS to 0.05 counts; orders above ~30 are
            # numerically unstable in CasaXPS itself
            self.assertLess(d, 1e-6 if order <= 20 else 0.05, order)


if __name__ == "__main__":
    unittest.main()
