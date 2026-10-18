"""``LA`` / ``LF`` / ``DS`` against CasaXPS's own exported component curves.

Two CasaXPS "Export All to ASCII" files, one component per column (the
component column carries the background, so the peak is the column minus
``Background``), each with the Position / FWHM / Area / Lineshape of the
components in its header:

* ``LF shape.csv`` - Mo 3d, ``LF(0.7,1.5,75,100)`` and ``LF(0.5,1,45,200)``
* ``DS F and H lineshapes.csv`` - HOPG C 1s: ``DS(0.09,500)``,
  ``H(0.09,250)SGL(90)``, ``F(0.09,32,150)SGL(90)`` and ``LA(50)`` components

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


if __name__ == "__main__":
    unittest.main()
