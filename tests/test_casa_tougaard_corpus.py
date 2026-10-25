"""CasaXPS's own Tougaard family (``lineshapes.casa_tougaard``) against real
fits: the residual of background + the file's own components against its raw
data. ``XPS_CASA_MODELS=<folder of .vms files>`` (the author's "CasaXPS
Models" folder) runs it; John Walton's al-cu-tougaard.vms and ZnSe.vms were
the first files with plain ``Tougaard`` (the literal-B/C/D reading was 15-27
times the data range out)."""
import os
import statistics
import unittest

import casafit
import readers

ROOT = os.environ.get("XPS_CASA_MODELS")


def _residuals(root, kinds):
    out = {}
    for dp, _dn, fn in os.walk(root):
        for f in sorted(fn):
            if not f.lower().endswith((".vms", ".vamas")):
                continue
            try:
                doc = readers.load_file(os.path.join(dp, f))
            except Exception:                              # noqa: BLE001
                continue
            for r in doc.regions:
                if not r.fit or not r.fit.regions:
                    continue
                dwell, scans = r.dwell_and_scans()
                for cv in casafit.curves(r.fit, r.energy, r.counts,
                                         r.photon_energy, dwell, scans):
                    kind = cv.fit_region.background
                    if kind in kinds and cv.components \
                            and cv.residual_rms is not None:
                        out.setdefault(kind, []).append(cv.residual_rms)
    return out


@unittest.skipUnless(ROOT and os.path.isdir(ROOT), "XPS_CASA_MODELS not set")
class TestCasaTougaardCorpus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = _residuals(ROOT, ("Tougaard", "U 3 Tougaard"))

    def test_plain_tougaard_is_within_a_few_percent_of_the_range(self):
        v = self.res.get("Tougaard")
        if not v:
            self.skipTest("no plain Tougaard fit in this folder")
        self.assertLess(statistics.median(v), 0.02)
        self.assertLess(max(v), 0.05)

    def test_u3_tougaard_is_drawn_and_close(self):
        v = self.res.get("U 3 Tougaard")
        if not v:
            self.skipTest("no U 3 Tougaard fit in this folder")
        self.assertLess(max(v), 0.05)


if __name__ == "__main__":
    unittest.main()
