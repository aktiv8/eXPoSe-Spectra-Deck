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


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestShirleyFamily(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(1000.0, 1020.0, 201)
        peak = 400.0 * np.exp(-((self.x - 1010.0) / 1.5) ** 2)
        step = np.cumsum(peak[::-1])[::-1] * 0.02          # Shirley-like rise to low KE
        self.y = 100.0 + peak + step

    def test_extrap_with_no_window_is_a_flat_baseline_shirley(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        b = ls.casa_extrap_shirley(self.x, self.y, lo, hi, 0.0, 0.0, 0.0, 25.0)
        self.assertAlmostEqual(b[-1], hi, 6)               # flat baseline at the end level
        self.assertAlmostEqual(b[0], lo, delta=0.001 * lo)  # the step reaches the start level
        self.assertTrue(np.all(np.diff(b) <= 1e-9))        # falls with KE

    def test_extrap_fits_a_line_to_the_last_p0_ev(self):
        y = self.y.copy()
        slope = 5.0                                         # a sloping tail
        tail = self.x > 1015.5          # the window starts at the point below 1016
        y[tail] = 100.0 + slope * (self.x[tail] - 1016.0)
        lo, hi = float(y[0]), float(y[-1])
        b = ls.casa_extrap_shirley(self.x, y, lo, hi, 0.0, 0.0, 4.0, 25.0)
        # inside the window the background is the fitted line: slope 5
        k = np.nonzero(self.x > 1017.0)[0]
        self.assertAlmostEqual(float(np.polyfit(self.x[k], b[k], 1)[0]), slope, 6)
        # (k stops changing by 1e-4 after a few passes, so the step reaches the
        # start level only to ~1.5 %: the program stops there too)
        self.assertAlmostEqual(float(b[0]), lo, delta=0.02 * lo)

    def test_extrap_shift_moves_the_step_along_the_axis(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        b0 = ls.casa_extrap_shirley(self.x, self.y, lo, hi, 0.0, 0.0, 0.0, 25.0)
        b3 = ls.casa_extrap_shirley(self.x, self.y, lo, hi, 0.0, 0.0, 0.0, 25.0, 3.0)
        j = int(np.argmin(np.abs(b0 - 0.5 * (lo + hi))))
        k = int(np.argmin(np.abs(b3 - 0.5 * (lo + hi))))
        self.assertAlmostEqual(self.x[j] - self.x[k], 3.0, delta=0.2)

    def test_offset_shirley_shifts_and_blends(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        plain = ls.casa_offset_shirley(self.x, self.y, lo, hi, shift=0.0)
        moved = ls.casa_offset_shirley(self.x, self.y, lo, hi, shift=4.0)
        self.assertGreater(float(np.abs(plain - moved).max()), 0.05 * (lo - hi))
        j = int(np.argmin(np.abs(plain - 0.5 * (lo + hi))))
        k = int(np.argmin(np.abs(moved - 0.5 * (lo + hi))))
        # (the shifted values feed back into the loop, so the step moves by
        # about, not exactly, the shift)
        self.assertGreater(self.x[j] - self.x[k], 3.0)
        self.assertLess(self.x[j] - self.x[k], 6.0)
        # p2 = 1 is the straight line between the levels
        line = ls.casa_offset_shirley(self.x, self.y, lo, hi, p2=1.0, shift=4.0)
        self.assertTrue(np.allclose(line, np.linspace(lo, hi, len(self.y))))

    def test_compress_shirley_squeezes_into_the_start_of_the_region(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        c = ls.casa_offset_shirley(self.x, self.y, lo, hi, shift=6.0, compress=True)
        self.assertIsNotNone(c)
        self.assertAlmostEqual(c[-1], hi, 3)

    def test_spline_family_keeps_the_ends_and_adds_to_the_knots(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        base = np.linspace(lo, hi, len(self.y))
        same = ls.casa_spline_bg(base, [0.0, 0.0, 0.0, 0.0])
        self.assertTrue(np.allclose(same, base))             # a line through a line
        up = ls.casa_spline_bg(base, [0.0, 50.0, 0.0, 0.0])
        e = len(base) - 1
        self.assertEqual(up[0], base[0])
        self.assertEqual(up[-1], base[-1])
        k = int(e / 5.0 * 2)
        self.assertAlmostEqual(up[k], base[k] + 50.0, 6)     # the second knot
        self.assertIsNone(ls.casa_spline_bg(base[:5], [0.0] * 4))   # too short
        params = (0.0, 0.0, 0.0, 50.0, 0.0, 0.0)
        out = ls.background("Spline Linear", self.y, 1, x=self.x, params=params,
                            ctx=(None, None))
        self.assertAlmostEqual(out[int(e / 5.0 * 2)],
                               ls.background("Linear", self.y, 1,
                                             ctx=(None, None))[int(e / 5.0 * 2)]
                               + 50.0, 6)

    def test_l_shirley_is_a_line_when_the_start_offset_is_zero(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        b = ls.casa_l_shirley(self.x, self.y, lo, hi)
        self.assertTrue(np.allclose(b, np.linspace(lo, hi, len(self.y))))
        with_offset = ls.casa_l_shirley(self.x, self.y, lo, hi, st_offset=30.0)
        self.assertGreater(with_offset[0], b[0] * 0.69)
        self.assertAlmostEqual(with_offset[-1], hi, 3)

    def test_iterated_shirley_converges_to_the_shirley_loop(self):
        lo, hi = float(self.y[0]), float(self.y[-1])
        it = ls.casa_iterated_shirley(self.y, lo, hi, 0.0)
        sh = ls.shirley_casa(self.y, lo, hi)
        self.assertLess(float(np.abs(it - sh).max()), 0.02 * (lo - hi))
        one = ls.casa_iterated_shirley(self.y, lo, hi, -1.0)      # one pass
        self.assertGreater(float(np.abs(one - it).max()), 1.0)

    def test_the_dispatcher_reaches_every_member(self):
        params = (0.0, 0.0, 2.0, 25.0, 0.0, 0.0)
        for kind in ("Extrap Shirley", "Offset Shirley", "Compress Shirley",
                     "L Shirley", "Spline Shirley", "Spline Linear",
                     "Spline None", "Iterated Shirley"):
            out = ls.background(kind, self.y, 2, x=self.x, params=params,
                                ctx=(None, None))
            self.assertIsNotNone(out, kind)
            self.assertEqual(len(out), len(self.y))
        # without a CasaXPS context they are not drawn (other sources)
        self.assertIsNone(ls.background("Extrap Shirley", self.y, 2, x=self.x,
                                        params=params))


@unittest.skipUnless(HAVE_NP and FOLDER and os.path.isdir(FOLDER)
                     and os.path.exists(os.path.join(FOLDER, "extrap_shroley.vms")),
                     "XPS_CASA_BG with extrap_shroley.vms not set")
class TestExtrapShirleyAgainstExport(unittest.TestCase):
    """``extrap_shroley.vms`` / ``.csv``: 12 Ti 2p regions of Extrap Shirley
    (window 2 eV, shifts 0-10 eV, start offset 7.1 %) exported by CasaXPS."""

    def test_every_background_matches(self):
        doc = readers.load_file(os.path.join(FOLDER, "extrap_shroley.vms"))
        regs = [r for r in doc.regions if r.fit and r.fit.regions]
        blocks = _blocks(os.path.join(FOLDER, "extrap_shroley.csv"))
        self.assertEqual(len(regs), len(blocks))
        self.assertEqual(len(regs), 12)
        for r, (_label, be, bg) in zip(regs, blocks):
            dwell, scans = r.dwell_and_scans()
            cv = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                dwell, scans)[0]
            k = dwell * scans
            ref = np.interp(np.array(r.energy), be[::-1], bg[::-1])
            mine = np.array(cv.background, dtype=float) / k
            ok = ~np.isnan(mine)
            rng = (max(r.counts) - min(r.counts)) / k
            self.assertLess(float(np.abs(mine[ok] - ref[ok]).max() / rng), 2e-4)



@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestUTougaardTypes(unittest.TestCase):
    """U 2 / U 3 / U 4 Tougaard as the program draws them
    (``FUN_14043db50`` cases 0x39 / 0x44, ``FUN_1404472f0``)."""

    def setUp(self):
        self.x = np.linspace(1000.0, 1030.0, 301)
        peak = 5000.0 * np.exp(-((self.x - 1015.0) / 2.0) ** 2)
        self.y = 200.0 + peak + 300.0 * (self.x < 1015.0) * 0.5
        self.ctx = (None, None)

    def tail(self, y, base, kern, step):
        n = len(y)
        t = np.arange(n) * step
        yb = np.maximum(y - base, 0.0)
        return np.correlate(yb, kern(t), mode="full")[n - 1:] * step

    def test_u2_is_scaled_to_the_start_level_with_the_c_floor(self):
        step = float(self.x[1] - self.x[0])
        params = (0.0, 0.0, 999.0, -650.0, 0.0, 0.0)       # the stored B is not used
        b = ls.background("U 2 Tougaard", self.y, 1, x=self.x, params=params,
                          ctx=self.ctx)
        base = self.y[-1]
        tail = self.tail(self.y, base, lambda t: t / (650.0 + t * t) ** 2, step)
        scale = (self.y[0] - base) / tail[0]
        self.assertTrue(np.allclose(b, tail * scale + base))
        self.assertAlmostEqual(b[0], self.y[0], 6)
        # |C| under 0.5 is read as 0.5
        b2 = ls.background("U 2 Tougaard", self.y, 1, x=self.x,
                           params=(0.0, 0.0, 999.0, -0.1, 0.0, 0.0), ctx=self.ctx)
        b3 = ls.background("U 2 Tougaard", self.y, 1, x=self.x,
                           params=(0.0, 0.0, 999.0, -0.5, 0.0, 0.0), ctx=self.ctx)
        self.assertTrue(np.allclose(b2, b3))

    def test_u2_falls_back_to_681_when_the_start_is_below_the_end(self):
        y = self.y.copy()
        y[:20] = 50.0                                       # start under the end level
        step = float(self.x[1] - self.x[0])
        b = ls.background("U 2 Tougaard", y, 1, x=self.x,
                          params=(0.0, 0.0, 1.0, -650.0, 0.0, 0.0), ctx=self.ctx)
        tail = self.tail(y, y[-1], lambda t: t / (650.0 + t * t) ** 2, step)
        self.assertTrue(np.allclose(b, tail * 681.2 + y[-1]))

    def test_other_sources_keep_the_stored_b(self):
        params = (0.0, 0.0, 800.0, -650.0, 0.0, 0.0)
        old = ls.tougaard_u2(self.x, self.y, 800.0, 650.0, 1)
        self.assertTrue(np.allclose(
            ls.background("U 2 Tougaard", self.y, 1, x=self.x, params=params), old))

    def test_u3_uses_the_stored_c_d_t0_and_matches_the_start(self):
        params = (0.0, 0.0, 0.0, 551.0, 436.0, 3.0)         # "U 3 T[551,436,3]"
        b = ls.background("U 3 Tougaard", self.y, 1, x=self.x, params=params,
                          ctx=self.ctx)
        self.assertAlmostEqual(b[0], self.y[0], 6)
        self.assertTrue(np.all(b[-5:] >= self.y[-1] - 1e-9))
        poly = ls.background("U Poly Tougaard", self.y, 1, x=self.x,
                             params=(0.0,) * 6, ctx=self.ctx)
        self.assertTrue(np.allclose(b, poly))               # the same cross section

    def test_u4_uses_the_stored_b_c_d_and_the_t0_gate(self):
        step = float(self.x[1] - self.x[0])
        params = (0.0, 0.0, 300.0, 400.0, 150.0, 2.0)
        b = ls.background("U 4 Tougaard", self.y, 1, x=self.x, params=params,
                          ctx=self.ctx)
        base = self.y[-1]

        def kern(t):
            den = (400.0 - t * t) ** 2 + 150.0 * t * t
            return np.where(t >= 2.0, 300.0 * t / den, 0.0)
        self.assertTrue(np.allclose(b, self.tail(self.y, base, kern, step) + base))
        # no scale is derived from the data: B doubles the tail
        b2 = ls.background("U 4 Tougaard", self.y, 1, x=self.x,
                           params=(0.0, 0.0, 600.0, 400.0, 150.0, 2.0), ctx=self.ctx)
        self.assertTrue(np.allclose(b2 - base, 2.0 * (b - base)))
        # nothing below T0 contributes
        none = ls.background("U 4 Tougaard", self.y, 1, x=self.x,
                             params=(0.0, 0.0, 300.0, 400.0, 150.0, 1e6), ctx=self.ctx)
        self.assertTrue(np.allclose(none, base))

    def test_u4_follows_a_line_when_the_start_is_below_the_end(self):
        y = np.full(len(self.x), 100.0)
        y[-1] = 300.0
        y[0] = 100.0
        b = ls.background("U 4 Tougaard", y, 1, x=self.x,
                          params=(0.0, 0.0, 0.0, 400.0, 150.0, 0.0), ctx=self.ctx)
        self.assertTrue(np.allclose(b, np.linspace(100.0, 300.0, len(y))))


if __name__ == "__main__":
    unittest.main()
