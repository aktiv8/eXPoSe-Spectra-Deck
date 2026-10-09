"""The ``GL(m)T(k)`` / ``SGL(m)T(k)`` exponential tail and the invalid-shape
guards (``lineshapes.check_shape``), and how ``casafit.curves`` uses them.

The tail is Casa's own formula (Cookbook 2026, p.67): the base shape plus
``(1 - base) * exp(k x)`` on the low-KE (high binding energy) side, x in FWHM.
No real CasaXPS fit using it is known on disk, so these tests check the
formula's own properties, not a Casa curve (see the lineshapes docstring).

Run:  python -m unittest discover tests
"""

import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import casafit  # noqa: E402
import htmlbrowser  # noqa: E402
import lineshapes as ls  # noqa: E402

try:
    import numpy as np
    HAVE_NP = True
except ImportError:                                  # pragma: no cover
    HAVE_NP = False


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class TestTail(unittest.TestCase):
    X = np.linspace(-8, 8, 16001) if HAVE_NP else None

    def raw(self, shape, fwhm=1.0):
        return ls._raw_values(self.X, ls.parse_shape(shape), 0.0, fwhm)

    def test_k_is_read(self):
        self.assertEqual(ls.parse_shape("GL(30)T(1.5)")["tk"], 1.5)
        self.assertEqual(ls.parse_shape("SGL(20)T( 0.8 )")["tk"], 0.8)
        self.assertIsNone(ls.parse_shape("GL(30)")["tk"])
        self.assertIsNone(ls.parse_shape("GL(30)T(x)")["tk"])
        self.assertEqual(ls.parse_shape("GL(30)T(1.5)")["mix"], 30.0)

    def test_only_the_high_binding_energy_side_changes(self):
        for kind in ("GL(30)", "SGL(30)"):
            base, tail = self.raw(kind), self.raw(kind + "T(1.5)")
            mid = len(self.X) // 2
            self.assertEqual(float(tail[mid]), float(base[mid]))   # the peak
            above = self.X > 0                       # high KE = low binding
            self.assertTrue(np.array_equal(tail[above], base[above]))
            below = self.X < 0
            self.assertTrue(np.all(tail[below] >= base[below]))
            self.assertGreater(float((tail - base)[below].max()), 0.1)

    def test_it_is_continuous_at_the_position(self):
        tail = self.raw("GL(30)T(1.5)")
        mid = len(self.X) // 2
        self.assertAlmostEqual(float(tail[mid - 1]), float(tail[mid]), 3)

    def test_a_larger_k_fades_faster(self):
        far = np.searchsorted(self.X, -3.0)
        v = [float((self.raw(f"GL(30)T({k})") - self.raw("GL(30)"))[far])
             for k in (0.3, 1.5, 6.0)]
        self.assertTrue(v[0] > v[1] > v[2] >= 0)

    def test_the_extra_is_the_documented_exponential(self):
        base, tail = self.raw("GL(30)"), self.raw("GL(30)T(1.5)")
        want = (1.0 - base) * np.where(self.X < 0, np.exp(1.5 * self.X), 0.0)
        self.assertTrue(np.allclose(tail - base, want, atol=1e-12))

    def test_fwhm_scales_the_tail(self):
        a = ls._raw_values(self.X * 2, ls.parse_shape("GL(30)T(1.5)"), 0.0, 2.0)
        b = self.raw("GL(30)T(1.5)")
        self.assertTrue(np.allclose(a, b, atol=1e-12))

    def test_the_component_keeps_its_stored_area(self):
        ke = np.linspace(-40, 40, 8001)
        step = ke[1] - ke[0]
        for shape in ("GL(30)T(1.5)", "SGL(40)T(0.5)", "GL(30)T(0.1)"):
            area = float(ls.component_curve(ke, shape, 0.0, 1.0, 100.0).sum()
                         * step)
            # 0.1 has a tail longer than this window: the rest is outside it
            self.assertAlmostEqual(area, 100.0, delta=1.0 if "0.1)" not in shape
                                   else 12.0, msg=shape)

    def test_a_very_slow_tail_still_normalises_over_its_whole_extent(self):
        ke = np.linspace(-400, 400, 80001)
        c = ls.component_curve(ke, "GL(30)T(0.05)", 0.0, 1.0, 100.0)
        self.assertAlmostEqual(float(c.sum() * (ke[1] - ke[0])), 100.0,
                               delta=0.5)

    def test_it_is_still_flagged_approximate_and_only_for_gl_and_sgl(self):
        self.assertFalse(ls.is_exact("GL(30)T(1.5)"))
        la = ls._raw_values(self.X, ls.parse_shape("LA(1,1,0)T(1.5)"),
                            0.0, 1.0)
        plain = ls._raw_values(self.X, ls.parse_shape("LA(1,1,0)"), 0.0, 1.0)
        self.assertTrue(np.array_equal(la, plain))


class TestCheckShape(unittest.TestCase):
    def test_real_world_shapes_pass(self):
        for s in ("GL(30)", "SGL(50)", "GL(30)T(1.5)", "LA(1.2,5,8)",
                  "LA(1.53,243)", "LA(30)", "LF(0.5,0.6,45,180,1)",
                  "LF(1.1,1.2,75,200)", "DS(0.09,500)", "A(0.15,0.7,20)SGL(12)",
                  "TLA(1,20,0)", "QF(1,2)", "H(0.09,250)SGL(90)"):
            self.assertEqual(ls.check_shape(s), [], s)

    def test_the_testers_refusals(self):
        cases = {
            "LA(0,1,1)": "exponents",
            "LA(1,-1,1)": "exponents",
            "LA(1,1,-3)": "Gaussian index",
            "LF(1,1)": "four",
            "LF(1,1,10,-5)": "m cannot be negative",
            "LF(-1,1,10,5)": "exponents",
            "DS(1.0,5)": "below 1",
            "DS(-0.1,5)": "at least 0",
            "DS(0.5,-5)": "Gaussian index",
            "GL(30)T(0)": "positive decay",
            "GL(30)T(-2)": "positive decay",
            "GL(30)T(x)": "positive decay",
            "GL(130)": "between 0 and 100",
            "TLA(0,5,1)": "alpha and mu",
            "A(0.1,0.5)": "three parameters",
        }
        for shape, word in cases.items():
            msgs = ls.check_shape(shape)
            self.assertTrue(msgs and any(word in m for m in msgs),
                            f"{shape}: {msgs}")

    def test_t_on_another_shape_is_said_to_be_unreconstructed(self):
        self.assertIn("only reconstructed for GL and SGL",
                      " ".join(ls.check_shape("LA(1,1,0)T(2)")))

    def test_non_finite_numbers(self):
        self.assertEqual(ls.check_shape("LA(nan,1,1)"),
                         ["a parameter is not a finite number"])
        self.assertEqual(ls.check_shape("GL(inf)"),
                         ["a parameter is not a finite number"])

    def test_garbage_never_raises(self):
        for s in ("", None, "???", "(", "GL(", "LA()", "LF()", "DS()", 7):
            self.assertIsInstance(ls.check_shape(s), list)


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class TestCurvesUseTheGuards(unittest.TestCase):
    HV = 1486.71

    def fit(self, *shapes):
        reg = casafit.FitRegion(name="R", background="none",
                                start_ke=1020.0, end_ke=1034.0)
        comps = [casafit.FitComponent(name=f"c{i}", shape=s, area=500.0,
                                      fwhm=1.5, pos_ke=1024.0 + 4.0 * i,
                                      region="R")
                 for i, s in enumerate(shapes)]
        return casafit.Fit(regions=[reg], components=comps)

    def curves(self, fit):
        be = np.linspace(self.HV - 1034.0, self.HV - 1020.0, 281).tolist()
        counts = [10.0] * len(be)
        return casafit.curves(fit, be, counts, self.HV, None, 1)[0]

    def test_an_invalid_component_is_left_out_and_named(self):
        cv = self.curves(self.fit("GL(30)", "LA(0,1,1)"))
        self.assertEqual(len(cv.components), 2)           # keeps its place
        good, bad = cv.components
        self.assertGreater(float(np.nanmax(good[1])), 0.0)
        self.assertEqual(float(np.nanmax(bad[1])), 0.0)   # drawn as nothing
        self.assertTrue(cv.approximate)
        self.assertTrue(np.all(np.isfinite(np.array(cv.envelope)[
            ~np.isnan(np.array(cv.envelope))])))
        self.assertEqual(len(cv.notes), 1)
        self.assertIn("c1", cv.notes[0])
        self.assertIn("exponents must be positive", cv.notes[0])

    def test_valid_shapes_give_no_notes(self):
        cv = self.curves(self.fit("GL(30)", "GL(30)T(1.5)", "LA(1.2,5,8)"))
        self.assertEqual(cv.notes, [])
        self.assertTrue(cv.approximate)                  # T(k) and LA

    def test_the_page_says_so_too(self):
        row = {"approximate": True, "background": "none",
               "background_known": True, "scale_known": True,
               "shape_notes": ["c1 (LA(0,1,1)): LA's exponents must be "
                               "positive; left out of the fit"]}
        notes = htmlbrowser.fit_notes([row])
        self.assertTrue(any("left out of the fit" in n for n in notes))
        row.pop("shape_notes")
        self.assertFalse(any("left out" in n for n in htmlbrowser.fit_notes([row])))


if __name__ == "__main__":
    unittest.main()
