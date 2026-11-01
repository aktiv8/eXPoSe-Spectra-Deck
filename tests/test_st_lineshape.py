"""``ST(w, d)BASE``: a peak multiplied by its own Shirley step
(``lineshapes._st_table``, read from ``FUN_140409340`` of CasaXPS.exe).

``XPS_CASA_ST=<path of es_bmim-pyrrntf2mix_230911.vms>`` (the one real file
known to use it, an O 1s fit of ``ST(1)LA(1,1.04,393)`` over ``Min Limits``)
runs the real-fit test, which compares residuals (no CasaXPS export of an ST
component exists, so the shape is checked against the data, not a curve)."""
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

PATH = os.environ.get("XPS_CASA_ST")


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestStShape(unittest.TestCase):
    def curve(self, shape, fwhm=1.0, area=1000.0, half=30.0, n=60001):
        x = np.linspace(-half * fwhm, half * fwhm, n)
        return x, ls.component_curve(x + 1000.0, shape, 1000.0, fwhm, area)

    def test_position_is_the_maximum(self):
        for shape in ("ST(1)LA(50)", "ST(5)LA(1,1.04,393)", "ST(0.2)GL(30)"):
            x, y = self.curve(shape, fwhm=1.5, half=20.0)
            self.assertAlmostEqual(float(x[np.argmax(y)]), 0.0, delta=0.05, msg=shape)

    def test_the_tail_is_on_the_low_ke_side(self):
        x, y = self.curve("ST(1)LA(50)", fwhm=1.5, half=40.0)
        low, high = y[x < -1.5].sum(), y[x > 1.5].sum()
        self.assertGreater(low, 1.15 * high)
        x, b = self.curve("LA(50)", fwhm=1.5, half=40.0)
        self.assertAlmostEqual(b[x < -1.5].sum() / b[x > 1.5].sum(), 1.0, delta=0.05)

    def test_w_sets_the_asymmetry(self):
        def ratio(w):
            x, y = self.curve(f"ST({w})LA(50)", fwhm=1.5, half=40.0)
            return y[x < -1.5].sum() / y[x > 1.5].sum()
        # h = g * (I + T / w): a larger w weakens the constant, so the tail grows
        self.assertGreater(ratio(5.0), ratio(0.5))
        self.assertGreater(ratio(500.0), ratio(5.0))
        self.assertAlmostEqual(ratio(1e-3), 1.0, delta=0.05)  # w -> 0: symmetric

    def test_w_of_zero_leaves_the_base_alone(self):
        self.assertIs(ls._table_for(ls.parse_shape("ST(0)LA(50)")),
                      ls._table_for(ls.parse_shape("LA(50)")))

    def test_area_and_smoothing(self):
        t = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        x, y0 = self.curve("ST(1)LA(50)", fwhm=1.0, area=500.0, half=300.0, n=240001)
        self.assertAlmostEqual(float(t(y0, x)) / 500.0, 1.0, delta=0.02)
        _, y1 = self.curve("ST(1,60)LA(50)", fwhm=1.0, area=500.0, half=300.0, n=240001)
        self.assertAlmostEqual(float(t(y1, x)) / 500.0, 1.0, delta=0.02)
        self.assertLess(np.abs(np.diff(y1, 2)).max(), np.abs(np.diff(y0, 2)).max())

    def test_invalid_forms(self):
        self.assertIsNone(ls._table_for(ls.parse_shape("ST(-1)LA(50)")))
        self.assertIsNone(ls._table_for(ls.parse_shape("ST(1)")))
        self.assertIsNone(ls._table_for(ls.parse_shape("ST(1)SB(0)LA(50)")))
        self.assertIsNone(ls._table_for(ls.parse_shape("ST(1,2000)LA(50)")))
        self.assertFalse(ls.is_exact("ST(1)LA(50)"))        # checked by residual only


@unittest.skipUnless(HAVE_NP and PATH and os.path.isfile(PATH),
                     "XPS_CASA_ST not set")
class TestAgainstARealFit(unittest.TestCase):
    def test_the_o1s_fit_is_closer_with_the_tail_and_the_recentring(self):
        doc = readers.load_file(PATH)
        r = next(r for r in doc.regions if r.fit and any(
            c.shape.startswith("ST(") for c in r.fit.components))
        dwell, scans = r.dwell_and_scans()

        def resid():
            ls._TABLE_CACHE.clear()
            return casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                  dwell, scans)[0].residual_rms
        with_st = resid()
        self.assertLess(with_st, 0.012)                     # it was 0.078 undrawn
        shapes = [c.shape for c in r.fit.components]
        try:
            for c in r.fit.components:
                c.shape = c.shape[c.shape.index(")") + 1:]  # the plain LA
            plain = resid()
        finally:
            for c, s in zip(r.fit.components, shapes):
                c.shape = s
            ls._TABLE_CACHE.clear()
        self.assertLess(with_st, 0.8 * plain)


if __name__ == "__main__":
    unittest.main()
