"""``TB(C, D, T0, n[, scale])BASE``, CasaXPS's "Tougaard background" lineshape
(``lineshapes._tb_table``, the negative-``d`` branch of ``FUN_140409ce0``), and
the region-index rule it exposed (``casafit._casa_positions``).

``XPS_CASA_TB=<folder holding TB0.vms and tb0.csv>`` (a Thermo K-Alpha+ Ti 2p
scan fitted over ``Min Limits`` with nine different ``TB(...)LA(50)``
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

FOLDER = os.environ.get("XPS_CASA_TB")
HV = 1486.68


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestTbShape(unittest.TestCase):
    def curve(self, shape, fwhm=2.0, area=1000.0, half=60.0, n=24001):
        x = np.linspace(-half, half, n)
        return x, ls.component_curve(x + 1000.0, shape, 1000.0, fwhm, area)

    def test_the_base_is_built_on_the_coarser_grid(self):
        # the parser builds a TB base on a 0.1 step (coarse step 1.0), not 0.025
        t = ls._table_for(ls.parse_shape("TB(27,0.5,-650,1)LA(50)"))
        self.assertEqual(len(t[0]), 6001)
        self.assertEqual(t[4], (0.1, 10))
        # and the default grid is restored afterwards
        self.assertEqual(len(ls._table_for(ls.parse_shape("LA(50)"))[0]), 24001)
        self.assertEqual(ls._TABLE_STEP, 0.025)

    def test_no_passes_is_the_stretched_peak(self):
        # n of 0 (or less) applies no Tougaard pass: a symmetric peak at Position
        for n in (0, -1):
            x, y = self.curve(f"TB(700,251,-250,{n})LA(50)", half=40.0)
            self.assertAlmostEqual(float(x[np.argmax(y)]), 0.0, delta=0.05)
            k = len(x) // 2
            self.assertAlmostEqual(float(y[k - 2000]), float(y[k + 2000]), delta=0.01 * float(y.max()))

    def test_scale_stretches_the_peak(self):
        def width(shape):
            x, y = self.curve(shape, half=60.0)
            above = x[y >= y.max() / 2]
            return float(above[-1] - above[0])
        w1 = width("TB(700,251,-250,0)LA(50)")
        w2 = width("TB(700,251,-250,0,2)LA(50)")
        self.assertAlmostEqual(w2 / w1, 2.0, delta=0.05)

    def test_trailing_comma_is_the_four_number_form(self):
        a = ls._table_for(ls.parse_shape("TB(27,0.5,-650,1,)LA(50)"))
        b = ls._table_for(ls.parse_shape("TB(27,0.5,-650,1)LA(50)"))
        self.assertTrue(np.array_equal(a[0], b[0]))

    def test_t0_below_the_data_is_the_same_cross_section(self):
        # a T0 under every kernel distance gates nothing
        a = ls._table_for(ls.parse_shape("TB(27,0.5,-650,1)LA(50)"))
        b = ls._table_for(ls.parse_shape("TB(27,0.5,-250,1)LA(50)"))
        self.assertTrue(np.array_equal(a[0], b[0]))

    def test_each_pass_acts_again(self):
        one = ls._table_for(ls.parse_shape("TB(27,0.5,-650,1)LA(50)"))
        two = ls._table_for(ls.parse_shape("TB(27,0.5,-650,2)LA(50)"))
        self.assertGreater(float(np.abs(one[0] - two[0]).max()), 0.01 * float(one[0].max()))

    def test_refused_forms_leave_nothing(self):
        # more than six passes, no overlap with the cross section (T0 beyond the
        # data leaves the second pass empty), no base shape
        for shape in ("TB(1,1,1,7)LA(50)", "TB(729.745,551,436,3)LA(50)",
                      "TB(1,1,1,1)"):
            self.assertIsNone(ls._table_for(ls.parse_shape(shape)), shape)
        self.assertEqual(ls.check_shape("TB(1,1,1,7)LA(50)"),
                         ["TB's repeat count is at most 6"])
        self.assertTrue(ls.check_shape("TB(729.745,551,436,3)LA(50)"))
        self.assertTrue(ls.check_shape("TB(1,1,1,1)"))
        self.assertEqual(ls.check_shape("TB(27,0.5,-650,1,)LA(50)"), [])

    def test_exactness(self):
        self.assertTrue(ls.is_exact("TB(27,0.5,-650,1,)LA(50)"))
        self.assertFalse(ls.is_exact("TB(1,1,1,7)LA(50)"))
        self.assertFalse(ls.is_exact("TB(729.745,551,436,3)LA(50)"))


@unittest.skipUnless(HAVE_NP, "numpy needed")
class TestRegionIndexRule(unittest.TestCase):
    def test_a_limit_on_a_grid_point_follows_the_float_arithmetic(self):
        ks = 1008.6 + 0.1 * np.arange(602)
        pos = casafit._casa_positions(ks, 1010.1, 1037.1, np.arange(3))
        # (1037.1 - 1008.6) * 601 / 60.1 comes out just under 285: the program
        # ends the region at 1037.0, and starts it at 1010.1
        self.assertAlmostEqual(float(ks[pos[0]]), 1010.1, places=6)
        self.assertAlmostEqual(float(ks[pos[-1]]), 1037.0, places=6)

    def test_an_ordinary_limit_is_ceil_and_floor(self):
        ks = 1000.0 + 0.5 * np.arange(101)
        pos = casafit._casa_positions(ks, 1010.1, 1030.9, np.arange(3))
        self.assertEqual(float(ks[pos[0]]), 1010.5)
        self.assertEqual(float(ks[pos[-1]]), 1030.5)

    def test_an_irregular_axis_keeps_the_inclusive_window(self):
        ks = np.array([0.0, 1.0, 2.5, 3.0, 4.0, 6.0])
        fallback = np.array([1, 2, 3])
        self.assertTrue(np.array_equal(
            casafit._casa_positions(ks, 1.0, 3.0, fallback), fallback))


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
                     "XPS_CASA_TB not set")
class TestAgainstCasaXpsExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        doc = readers.load_file(os.path.join(FOLDER, "TB0.vms"))
        cls.regs = [r for r in doc.regions if r.fit and r.fit.regions]
        cls.blocks = _blocks(os.path.join(FOLDER, "tb0.csv"))
        assert len(cls.regs) == len(cls.blocks) == 12

    def test_every_tb_component_matches(self):
        seen, drawn_nothing = 0, 0
        for (cols, a), r in zip(self.blocks, self.regs):
            ke = HV - a[:, 0]
            bg = a[:, -2]
            for k, c in enumerate(r.fit.components):
                if not c.shape.startswith("TB("):
                    continue
                lit = a[:, 2 + k] - bg
                if float(np.abs(lit).max()) == 0.0:
                    # CasaXPS refuses the string and draws nothing
                    self.assertTrue(ls.check_shape(c.shape), c.shape)
                    drawn_nothing += 1
                    continue
                mine = ls.component_curve(ke, c.shape, c.pos_ke, c.fwhm, c.area)
                self.assertLess(float(np.abs(mine - lit).max() / lit.max()), 5e-4,
                                c.shape)
                seen += 1
        self.assertEqual(drawn_nothing, 1)
        self.assertEqual(seen, 11)

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


if __name__ == "__main__":
    unittest.main()
