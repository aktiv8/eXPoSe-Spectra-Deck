"""Optional: every PHI ``.spe`` in a folder loads, and agrees with PHI's own
VAMAS export of it where one sits beside it.

    set XPS_PHI_CORPUS=C:\\path\\to\\folder\\of\\spe_and_npl_files
    python -m unittest tests.test_phi_corpus

Checked on a PHI Model 5000 set (float64) and a Quantes set (float32, Al and
Cr sources, ``.spe`` + ``.npl`` pairs): the ``.npl`` is MultiPak's own VAMAS
export of the same scan, so energies, counts, photon energy, anode, pass
energy, dwell and the analyser mode must be the same. A ``.spe`` written by
another program (a big-endian CasaXPS output has no partner) is only
checked for being plausible.
"""

import glob
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from readers import load_file  # noqa: E402

CORPUS = os.environ.get("XPS_PHI_CORPUS", "")


def _partner(spe):
    base = os.path.splitext(spe)[0]
    for ext in (".npl", ".vms", ".NPL", ".VMS"):
        if os.path.exists(base + ext):
            return base + ext
    return None


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_PHI_CORPUS to a folder of PHI .spe files")
class TestPhiCorpus(unittest.TestCase):
    def test_every_file_loads_and_matches_its_export(self):
        files = sorted(p for p in glob.glob(os.path.join(CORPUS, "**", "*"),
                                            recursive=True)
                       if p.lower().endswith(".spe"))
        self.assertTrue(files, "no .spe files in the folder")
        paired = 0
        for spe in files:
            f = load_file(spe)
            self.assertTrue(f.regions, spe)
            for r in f.regions:
                self.assertEqual(len(r.counts), len(r.energy), spe)
                self.assertTrue(all(c == c and abs(c) < 1e12
                                    for c in r.counts), spe)
            other = _partner(spe)
            if not other:
                continue
            paired += 1
            g = load_file(other)
            self.assertEqual(len(f.regions), len(g.regions), spe)
            for a, b in zip(f.regions, g.regions):
                tag = f"{os.path.basename(spe)} {a.name}"
                if a.name != "Survey":      # an old export spells it "S UR"
                    self.assertEqual(a.name, b.name, tag)
                self.assertEqual(a.n_points, b.n_points, tag)
                for x, y in zip(a.energy, b.energy):
                    self.assertAlmostEqual(x, y, places=3, msg=tag)
                # MultiPak 2009's export is counts per step in "arb." units
                # (counts/s x time per step) rounded to one or two decimals;
                # the Quantes export is c/s to six figures
                per_step = b.count_units.lower() != "c/s"
                scale = a.dwell if per_step else 1.0
                tol = 6e-2 if per_step else 6e-4
                for x, y in zip(a.counts, b.counts):
                    self.assertAlmostEqual(x * scale, y,
                                           delta=tol + 1e-6 * abs(y), msg=tag)
                self.assertAlmostEqual(a.photon_energy, b.photon_energy,
                                       places=3, msg=tag)
                self.assertEqual(a.anode.split()[0], b.anode.split()[0], tag)
                self.assertAlmostEqual(a.pass_energy, b.pass_energy, msg=tag)
                self.assertAlmostEqual(a.dwell, b.dwell, places=5, msg=tag)
                self.assertEqual(a.extra.get("analyser_mode"),
                                 b.extra.get("analyser_mode"), tag)
                self.assertEqual(a.date, b.date[:10], tag)
        print(f"\nPHI corpus: {len(files)} .spe loaded, {paired} checked "
              "against their export")


if __name__ == "__main__":
    unittest.main()
