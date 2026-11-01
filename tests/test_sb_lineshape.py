"""``SB(d)BASE``, CasaXPS's "Shirley background" lineshape: a peak shape
(``LA(50)`` ...) turned into the step a Shirley background would have under it
(``lineshapes._sb_table``, read from ``FUN_140409ce0`` of CasaXPS.exe).

``XPS_CASA_SB=<folder holding SB0.vms and sb0.csv>`` (a Thermo K-Alpha+ Ti 2p
scan fitted with ``Min Limits`` and ``SB(0)LA(50)`` / ``LA(70)`` / ``LA(90)``
components, then "Export All to ASCII") compares every component, background
and envelope with CasaXPS's own."""
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

FOLDER = os.environ.get("XPS_CASA_SB")
HV = 1486.68


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestSbShape(unittest.TestCase):
    def curve(self, shape, fwhm=1.0, area=1000.0, span=400.0, n=80001):
        x = np.linspace(-span * fwhm, span * fwhm, n)
        return x, ls.component_curve(x + 1000.0, shape, 1000.0, fwhm, area)

    def test_it_is_a_step_that_is_level_on_the_low_ke_side(self):
        x, y = self.curve("SB(0)LA(50)", span=200.0, n=40001)
        # (the Savitzky-Golay integration leaves a ripple of ~1e-5 of the
        # plateau; the table itself ends at +-300 units, beyond which it is 0)
        self.assertTrue(np.all(np.diff(y) <= 1e-4 * y.max()))   # never rises with KE
        self.assertGreater(y[0], 0.99 * y.max())                # level at low KE
        self.assertLess(y[-1], 0.02 * y.max())                  # gone at high KE

    def test_the_area_is_the_step_integrated_over_the_table(self):
        # the stated Area is the integral of the step over the table's +-300
        # units (a few per cent of it lies beyond a 400-fwhm window's tail)
        for shape in ("SB(0)LA(50)", "SB(0)LA(90)"):
            # (+-300 table units are +-300 / W fwhm: 485 for LA(50), 327 for LA(90))
            x, y = self.curve(shape, fwhm=1.5, area=2.5e6, span=700.0, n=280001)
            t = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
            self.assertAlmostEqual(float(t(y, x)) / 2.5e6, 1.0, delta=0.005)

    def test_the_peak_under_the_step_is_wider_than_the_stated_fwhm_by_1_over_w(self):
        # SB(0)LA(m) keeps the base table's measured width W, so the peak whose
        # integral the step is has a FWHM of fwhm / W
        for m in (50, 70, 90):
            base = ls.parse_shape(f"LA({m})")
            w = ls._table_for(base)[1]
            x = np.linspace(-20.0, 20.0, 40001)
            y = ls.component_curve(x + 1000.0, f"SB(0)LA({m})", 1000.0, 2.0, 1000.0)
            q = -np.gradient(y, x)
            half = q.max() / 2
            above = x[q >= half]
            # (the table is sampled every 0.025 units: ~2.5 % of this width)
            self.assertAlmostEqual((above[-1] - above[0]) / (2.0 / w), 1.0, delta=0.04)

    def test_the_tabulated_width_is_the_bases_not_remeasured(self):
        base = ls._table_for(ls.parse_shape("LA(50)"))
        sb = ls._table_for(ls.parse_shape("SB(0)LA(50)"))
        self.assertEqual(sb[1], base[1])

    def test_d_smooths_the_step_and_keeps_the_area(self):
        x, y0 = self.curve("SB(0)LA(50)", fwhm=2.0, span=600.0, n=120001)
        _, y1 = self.curve("SB(60)LA(50)", fwhm=2.0, span=600.0, n=120001)
        self.assertLess(np.abs(np.diff(y1, 2)).max(), np.abs(np.diff(y0, 2)).max())
        t = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        self.assertAlmostEqual(float(t(y1, x) / t(y0, x)), 1.0, delta=0.01)

    def test_exactness_and_invalid_forms(self):
        self.assertTrue(ls.is_exact("SB(0)LA(50)"))
        self.assertFalse(ls.is_exact("SB(2000)LA(50)"))     # d beyond the program's range
        self.assertFalse(ls.is_exact("SB(0)"))              # no base shape
        self.assertFalse(ls.is_exact("SB(0)SB(0)LA(50)"))   # not nested
        self.assertEqual(ls.check_shape("SB(0)LA(50)"), [])

    def test_a_peak_shape_is_unchanged_by_the_lookup_refactor(self):
        x = np.linspace(990.0, 1010.0, 2001)
        y = ls.component_curve(x, "LA(50)", 1000.0, 1.5, 100.0)
        self.assertAlmostEqual(float(np.trapezoid(y, x) if hasattr(np, "trapezoid")
                                     else np.trapz(y, x)), 100.0, delta=1.0)


def _blocks(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    head = rows[6]
    out, i = [], 0
    while i < len(head):
        if head[i] == "B.E.":
            j = i + 1
            while j < len(head) and head[j] != "":
                j += 1
            data = [[float(v) if v != "" else np.nan for v in r[i:j]]
                    for r in rows[7:] if len(r) > i and r[i] != ""]
            out.append((head[i:j], np.array(data)))
            i = j
        i += 1
    return out


@unittest.skipUnless(HAVE_NP and FOLDER and os.path.isdir(FOLDER),
                     "XPS_CASA_SB not set")
class TestAgainstCasaXpsExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        doc = readers.load_file(os.path.join(FOLDER, "SB0.vms"))
        cls.regs = [r for r in doc.regions if r.fit and r.fit.regions]
        cls.blocks = _blocks(os.path.join(FOLDER, "sb0.csv"))
        assert len(cls.regs) == len(cls.blocks) == 12

    def test_every_component_matches(self):
        seen = {}
        for (cols, a), r in zip(self.blocks, self.regs):
            ke = HV - a[:, 0]
            bg = a[:, -2]
            for k, c in enumerate(r.fit.components):
                lit = a[:, 2 + k] - bg
                mine = ls.component_curve(ke, c.shape, c.pos_ke, c.fwhm, c.area)
                err = float(np.abs(mine - lit).max() / lit.max())
                base = c.shape[:c.shape.index(")") + 1] if c.shape.startswith("SB") \
                    else "peak"
                self.assertLess(err, 5e-4, c.shape)
                seen[base] = seen.get(base, 0) + 1
        self.assertGreaterEqual(seen.get("SB(0)", 0), 20)

    def test_whole_fits_reproduce_the_background_and_envelope(self):
        for (cols, a), r in zip(self.blocks, self.regs):
            dwell, scans = r.dwell_and_scans()
            cv = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                dwell, scans)[0]
            k = dwell * scans
            be = a[:, 0]
            o = np.argsort(be)
            energy = np.array(r.energy)
            rng = (max(r.counts) - min(r.counts)) / k
            for mine, col in ((cv.background, -2), (cv.envelope, -1)):
                ref = np.interp(energy, be[o], a[:, col][o])
                m = np.array(mine) / k
                ok = ~np.isnan(m)
                self.assertLess(float(np.abs(m[ok] - ref[ok]).max() / rng), 2e-4)
            self.assertFalse(cv.approximate)


if __name__ == "__main__":
    unittest.main()
