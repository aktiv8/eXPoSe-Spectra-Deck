"""``LA`` / ``LF`` / ``DS`` against CasaXPS's own exported component curves.

Two CasaXPS "Export All to ASCII" files, one component per column (the
component column carries the background, so the peak is the column minus
``Background``), each with the Position / FWHM / Area / Lineshape of the
components in its header:

* ``LF shape.csv`` - Mo 3d, ``LF(0.7,1.5,75,100)`` and ``LF(0.5,1,45,200)``
* ``DS F and H lineshapes.csv`` - HOPG C 1s: ``DS(0.09,500)``,
  ``H(0.09,250)SGL(90)``, ``F(0.09,32,150)SGL(90)`` and ``LA(50)`` components

* ``A_GL lineshapes.csv`` - HOPG C 1s, three fits of ``A(0.4,0.38,20)GL(20)``
  plus a ``GL(30)`` component each
* ``GL SGL and T tail.csv`` - Au 4d doublet, eight fits: ``SGL(70)``,
  ``SGL(70)T(1|1.5|2)``, ``SGL(90)``, ``GL(70)T(1|1.5|2)``
* ``E with GL and SGL.csv`` - Au 4d, nine ``E(1,b,n)GL|SGL(m)`` fits
* ``more E and PHI.csv`` - three more ``E`` fits (a = 1, 0.5, 0.1) and six
  ``GL(m,n)PHI(j,w)`` fits (``SL`` makes the same call as ``PHI``)
* ``ls shapes.csv`` - Au 4d, ten ``LS(a,b,n)`` fits (a, b 0.5-1.5, n 10-500)
* ``qf.csv`` - Au 4d, 13 fits: nine ``QF(a,b,w,n[,k])`` and four ``QSF``
  (a 0.5-1, b 1-12, w 20-50, n 0-100, k 1-12)
* ``ds 3 param noe second cannot be 0.csv`` - C 1s, three ``DS(a,n,m)`` carbide
  components (a 0.01 and 0.1, n 2 and 3, m 400; the program refuses n of 0 and 1)
* ``TLA lineshape.csv`` - Si 2p, two ``TLA(1,9.46667,222)`` components

Only the shape is compared (the curve is scaled to the literal one); the
position is the table's centre, not its maximum. Skipped when the files are
not present (``XPS_LITERAL_DIR`` or ``D:\\Temp``).
"""
import csv
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import numpy as np
    HAVE_NP = True
except ImportError:                                  # pragma: no cover
    HAVE_NP = False

import lineshapes as ls  # noqa: E402

DIR = os.environ.get("XPS_LITERAL_DIR", r"D:\Temp")
HV = 1486.6


def _path(name):
    p = os.path.join(DIR, name)
    return p if os.path.isfile(p) else None


def _blocks(path):
    with open(path, encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    hdr = next(i for i, r in enumerate(rows) if r and r[0] == "B.E.")
    meta, head = rows[:hdr], rows[hdr]
    data = [r for r in rows[hdr + 1:] if r and r[0]]
    out, start = [], 0
    while start < len(head):
        if head[start] != "B.E.":
            start += 1
            continue
        end = start + 1
        while end < len(head) and head[end] != "B.E.":
            end += 1
        info = {m[0]: m[start + 2:end] for m in meta if m and m[0]}
        cols = np.array([[float(r[c]) if r[c] != "" else np.nan
                          for c in range(start, end)]
                         for r in data if len(r) >= end])
        out.append((info, head[start:end], cols))
        start = end
    return out


def _errors(path, wanted):
    """``{shape: max |mine - literal| / literal peak}`` of the components
    whose shape starts with one of ``wanted``."""
    res = {}
    for info, head, cols in _blocks(path):
        be = cols[:, 0]
        bg = cols[:, head.index("Background")]
        ke = HV - be
        order = np.argsort(ke)
        comp = [i for i, h in enumerate(head)
                if i > 1 and h not in ("Background", "Envelope")]
        for k, ci in enumerate(comp):
            shape = info["Lineshape"][k]
            if not shape.startswith(wanted):
                continue
            lit = cols[:, ci] - bg
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][k]),
                float(info["FWHM"][k]), 1.0)
            sc = float(np.dot(mine, lit) / np.dot(mine, mine))
            res[shape] = float(np.abs(sc * mine - lit).max() / lit.max())
    return res


@unittest.skipUnless(HAVE_NP and _path("LF shape.csv"), "LF shape.csv not found")
class TestLfLiteral(unittest.TestCase):
    def test_lf_matches_casaxps_curves(self):
        err = _errors(_path("LF shape.csv"), ("LF",))
        self.assertEqual(len(err), 2)
        for shape, e in err.items():
            self.assertLess(e, 0.02, shape)


@unittest.skipUnless(HAVE_NP and _path("DS F and H lineshapes.csv"),
                     "DS F and H lineshapes.csv not found")
class TestDsLaLiteral(unittest.TestCase):
    def test_ds_matches_casaxps_curve(self):
        err = _errors(_path("DS F and H lineshapes.csv"), ("DS",))
        self.assertEqual(list(err), ["DS(0.09,500)"])
        self.assertLess(err["DS(0.09,500)"], 0.01)

    def test_h_matches_casaxps_curve(self):
        err = _errors(_path("DS F and H lineshapes.csv"), ("H(",))
        self.assertEqual(list(err), ["H(0.09,250)SGL(90)"])
        self.assertLess(err["H(0.09,250)SGL(90)"], 0.01)

    def test_f_matches_casaxps_curve(self):
        err = _errors(_path("DS F and H lineshapes.csv"), ("F(",))
        self.assertEqual(list(err), ["F(0.09,32,150)SGL(90)"])
        self.assertLess(err["F(0.09,32,150)SGL(90)"], 0.015)

    def test_la_matches_casaxps_curves(self):
        err = _errors(_path("DS F and H lineshapes.csv"), ("LA",))
        self.assertEqual(list(err), ["LA(50)"])
        self.assertLess(err["LA(50)"], 0.005)


@unittest.skipUnless(HAVE_NP and _path("TLA lineshape.csv"),
                     "TLA lineshape.csv not found")
class TestTlaLiteral(unittest.TestCase):
    def test_tla_matches_casaxps_curves(self):
        # two components with the same shape string: the dict keeps the last,
        # so check them one by one
        for info, head, cols in _blocks(_path("TLA lineshape.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            for k, ci in enumerate((2, 3)):
                lit = cols[:, ci] - bg
                mine = np.empty_like(lit)
                mine[order] = ls.component_curve(
                    ke[order], info["Lineshape"][k],
                    HV - float(info["Position"][k]),
                    float(info["FWHM"][k]), float(info["Area"][k]))
                # the stated Area is the table's own: no free scale needed
                self.assertLess(float(np.abs(mine - lit).max() / lit.max()),
                                0.03, info["Position"][k])


@unittest.skipUnless(HAVE_NP and _path("A_GL lineshapes.csv"),
                     "A_GL lineshapes.csv not found")
class TestAGlLiteral(unittest.TestCase):
    """The stated Area is the table's own, so no scale is fitted: the A
    shape's Area is its symmetric half doubled (the whole peak is 1.35 times
    larger), and a GL(30)'s is its whole integral."""

    def _check(self, prefix, tol):
        seen = 0
        for info, head, cols in _blocks(_path("A_GL lineshapes.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            for k, ci in enumerate((2, 3)):
                shape = info["Lineshape"][k]
                if not shape.startswith(prefix):
                    continue
                lit = cols[:, ci] - bg
                mine = np.empty_like(lit)
                mine[order] = ls.component_curve(
                    ke[order], shape, HV - float(info["Position"][k]),
                    float(info["FWHM"][k]), float(info["Area"][k]))
                self.assertLess(float(np.abs(mine - lit).max() / lit.max()),
                                tol, f"{shape} at {info['Position'][k]}")
                seen += 1
        self.assertEqual(seen, 3)

    def test_a_matches_casaxps_curves(self):
        self._check("A(", 0.025)

    def test_gl_matches_casaxps_curves(self):
        self._check("GL(", 0.005)


@unittest.skipUnless(HAVE_NP and _path("GL SGL and T tail.csv"),
                     "GL SGL and T tail.csv not found")
class TestTailLiteral(unittest.TestCase):
    """``GL(m)T(k)`` / ``SGL(m)T(k)``, the stated Area used as it is."""

    def test_every_component_matches_casaxps(self):
        seen = {}
        for info, head, cols in _blocks(_path("GL SGL and T tail.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            for k, ci in enumerate((2, 3)):
                shape = info["Lineshape"][k]
                lit = cols[:, ci] - bg
                mine = np.empty_like(lit)
                mine[order] = ls.component_curve(
                    ke[order], shape, HV - float(info["Position"][k]),
                    float(info["FWHM"][k]), float(info["Area"][k]))
                err = float(np.abs(mine - lit).max() / lit.max())
                self.assertLess(err, 0.006, f"{shape} at {info['Position'][k]}")
                seen[shape] = seen.get(shape, 0) + 1
        self.assertEqual(sum(seen.values()), 16)
        self.assertEqual(sum(1 for s in seen if "T(" in s), 6)

    def test_a_larger_k_gives_a_shorter_tail(self):
        x = np.linspace(990.0, 1020.0, 30001)
        tails = []
        for k in (1.0, 1.5, 2.0):
            y = ls.component_curve(x, f"SGL(70)T({k})", 1005.0, 3.0, 1.0)
            tails.append(float(y[x < 1000.0].sum() / y.sum()))
        self.assertGreater(tails[0], tails[1])
        self.assertGreater(tails[1], tails[2])


@unittest.skipUnless(HAVE_NP and _path("E with GL and SGL.csv"),
                     "E with GL and SGL.csv not found")
class TestELiteral(unittest.TestCase):
    """``E(a,b,n)GL(m)`` / ``SGL(m)``, the stated Area used as it is. The first
    file has a = 1 (the ``(1 - a) / 2`` power is zero); the second adds
    a = 0.5 and 0.1, which test the power."""

    def test_every_fit_matches_casaxps(self):
        seen = 0
        for info, head, cols in _blocks(_path("E with GL and SGL.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            lit = cols[:, 2] - bg
            shape = info["Lineshape"][0]
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][0]),
                float(info["FWHM"][0]), float(info["Area"][0]))
            err = float(np.abs(mine - lit).max() / lit.max())
            self.assertLess(err, 0.004, shape)
            seen += 1
        self.assertEqual(seen, 9)

    @unittest.skipUnless(_path("more E and PHI.csv"), "more E and PHI.csv not found")
    def test_phi_matches_casaxps(self):
        seen = 0
        for info, head, cols in _blocks(_path("more E and PHI.csv")):
            shape = info["Lineshape"][0]
            if "PHI" not in shape:
                continue
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            lit = cols[:, 2] - bg
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][0]),
                float(info["FWHM"][0]), float(info["Area"][0]))
            err = float(np.abs(mine - lit).max() / lit.max())
            self.assertLess(err, 0.004, shape)
            self.assertTrue(ls.is_exact(shape), shape)
            seen += 1
        self.assertEqual(seen, 6)

    @unittest.skipUnless(_path("more E and PHI.csv"), "more E and PHI.csv not found")
    def test_a_below_one_matches_casaxps(self):
        seen = {}
        for info, head, cols in _blocks(_path("more E and PHI.csv")):
            shape = info["Lineshape"][0]
            if not shape.startswith("E("):
                continue
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            lit = cols[:, 2] - bg
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][0]),
                float(info["FWHM"][0]), float(info["Area"][0]))
            seen[shape] = float(np.abs(mine - lit).max() / lit.max())
        self.assertEqual(sorted(seen), ["E(0.1,5,30)SGL(50)",
                                        "E(0.5,2,50)SGL(50)",
                                        "E(1,1,50)SGL(90)"])
        for shape, err in seen.items():
            self.assertLess(err, 0.004, shape)


@unittest.skipUnless(HAVE_NP and _path("ls shapes.csv"), "ls shapes.csv not found")
class TestLsLiteral(unittest.TestCase):
    """``LS(a,b,n)``: the ``LA`` table with the "Moving Average" (boxcar)
    kernel, the stated Area used as it is. n above 80 takes the stretch /
    clamp path."""

    def test_every_fit_matches_casaxps(self):
        seen = 0
        for info, head, cols in _blocks(_path("ls shapes.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            lit = cols[:, 2] - bg
            shape = info["Lineshape"][0]
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][0]),
                float(info["FWHM"][0]), float(info["Area"][0]))
            err = float(np.abs(mine - lit).max() / lit.max())
            self.assertLess(err, 0.005, shape)
            self.assertTrue(ls.is_exact(shape), shape)
            seen += 1
        self.assertEqual(seen, 10)


@unittest.skipUnless(HAVE_NP and _path("qf.csv"), "qf.csv not found")
class TestQfLiteral(unittest.TestCase):
    """``QF`` (the ``LF`` table with a ``16 x^4`` base) and ``QSF`` (the same
    with the boxcar kernel), the stated Area used as it is."""

    def test_every_fit_matches_casaxps(self):
        seen = {"QF": 0, "QSF": 0}
        for info, head, cols in _blocks(_path("qf.csv")):
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            lit = cols[:, 2] - bg
            shape = info["Lineshape"][0]
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][0]),
                float(info["FWHM"][0]), float(info["Area"][0]))
            err = float(np.abs(mine - lit).max() / lit.max())
            self.assertLess(err, 0.006, shape)
            self.assertTrue(ls.is_exact(shape), shape)
            seen[shape.split("(")[0]] += 1
        self.assertEqual(seen, {"QF": 9, "QSF": 4})


@unittest.skipUnless(HAVE_NP and _path("ds 3 param noe second cannot be 0.csv"),
                     "ds 3 param noe second cannot be 0.csv not found")
class TestDsCutoffLiteral(unittest.TestCase):
    """``DS(a,n,m)``: the DS form with the program's finite-extent correction.
    The stated Area is used as it is. Agreement is 0.9-1.4 % of the peak (the
    older code, which read m as nothing, was 73-82 % out); the plain DS(a,n)
    under it is 0.6-1.3 % out against the same curves at these small n, so the
    remainder is not the cut-off."""

    def test_the_three_fits_are_close_to_casaxps(self):
        seen = 0
        for info, head, cols in _blocks(_path("ds 3 param noe second cannot be 0.csv")):
            k = 4                      # the carbide DS component of each block
            shape = info["Lineshape"][k]
            if shape.count(",") != 2:
                continue
            be = cols[:, 0]
            bg = cols[:, head.index("Background")]
            ke = HV - be
            order = np.argsort(ke)
            ci = [i for i, h in enumerate(head)
                  if i > 1 and h not in ("Background", "Envelope")][k]
            lit = cols[:, ci] - bg
            mine = np.empty_like(lit)
            mine[order] = ls.component_curve(
                ke[order], shape, HV - float(info["Position"][k]),
                float(info["FWHM"][k]), float(info["Area"][k]))
            err = float(np.abs(mine - lit).max() / lit.max())
            self.assertLess(err, 0.02, shape)
            # the far tails, where the cut-off correction acts, agree to a
            # few tenths of a per cent of the peak
            tail = float(np.abs(mine - lit)[(be > 290) | (be < 279)].max() / lit.max())
            self.assertLess(tail, 0.005, shape)
            seen += 1
        self.assertEqual(seen, 3)


if __name__ == "__main__":
    unittest.main()
