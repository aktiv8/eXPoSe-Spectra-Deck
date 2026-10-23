"""Which transmission function a CasaXPS VAMAS ``Transmission`` column holds:
Avantage's polynomial, CasaXPS's recalculated VG function, or a measured file.

The three NaCl exports of one Nexsa G2 run (``XPS_TXFN_CORPUS=<folder>``) share
their counts and differ only in this column."""
import os
import unittest
from types import SimpleNamespace

import txfn
import quantview

AV = [4.70503, 0.606984, -0.600862, 0.071214]
C = [-0.116457, 0.683645, -0.625035, 0.0438625]
KE = [140.0 + 25.0 * i for i in range(55)]


def avantage(pe):
    return [txfn.transmission(e, pe, AV) for e in KE]


def casa(pe, coeffs=C):
    return [txfn.casa_formula(e, pe, coeffs) for e in KE]


LINES = ["PE = 50, Step = 0.1",
         "New VG Transmission Function: C0 = -0.116457, C1 = 0.683645, "
         "C2 = -0.625035, C3 = 0.0438625 ",
         "Transmission file: C:\\Users\\casa-\\Customers Data\\X-Ray001\\"
         "TransFunct-NEXSA_G2_Standard_RelativeToESCALAB250XiCorrectedPE40.vms"]


class TestParse(unittest.TestCase):
    def test_coefficients_and_file_name_only(self):
        got = txfn.parse_casa_lines(LINES)
        self.assertEqual(got["coeffs"], C)
        self.assertEqual(
            got["file"],
            "TransFunct-NEXSA_G2_Standard_RelativeToESCALAB250XiCorrectedPE40"
            ".vms")                       # no path leaves the app

    def test_slash_paths_and_missing_lines(self):
        got = txfn.parse_casa_lines(["Transmission file: /home/u/tf/a.vms"])
        self.assertEqual((got["coeffs"], got["file"]), (None, "a.vms"))
        self.assertEqual(txfn.parse_casa_lines(["nothing here"]),
                         {"coeffs": None, "file": ""})

    def test_gap_in_the_coefficients_is_not_used(self):
        got = txfn.parse_casa_lines(["New VG Transmission Function: C0 = 1, "
                                     "C2 = 3"])
        self.assertIsNone(got["coeffs"])


class TestClassify(unittest.TestCase):
    casa_lines = txfn.parse_casa_lines(LINES)

    def test_avantage_polynomial(self):
        got = txfn.classify(KE, 100.0, avantage(100.0), AV, {})
        self.assertEqual(got["kind"], "avantage")
        self.assertEqual(got["coeffs"], AV)

    def test_casa_recalculated_formula_has_no_hundredth(self):
        col = casa(100.0)
        got = txfn.classify(KE, 100.0, col, AV, txfn.parse_casa_lines(LINES[:2]))
        self.assertEqual(got["kind"], "casa_formula")
        self.assertEqual(got["coeffs"], C)
        # it is its own function, not Avantage's scaled
        self.assertGreater(max(abs(c / a - 1) for c, a in zip(
            col, avantage(100.0))), 0.05)

    def test_a_measured_file_follows_neither_formula(self):
        col = [t * (1 + 0.1 * (i % 7) / 7) for i, t in enumerate(avantage(150.0))]
        got = txfn.classify(KE, 150.0, col, AV, self.casa_lines)
        self.assertEqual(got["kind"], "casa_file")
        self.assertTrue(got["file"].startswith("TransFunct-"))
        lo, hi = got["vs_avantage"]
        self.assertLess(lo, 0)
        self.assertGreater(hi, 0)

    def test_a_stale_comment_cannot_label_the_column(self):
        # the comment names the formula, the numbers say it is something else
        col = [t * 1.01 for t in casa(50.0)]
        self.assertIsNone(txfn.classify(KE, 50.0, col, AV,
                                        txfn.parse_casa_lines(LINES[:2])))

    def test_nothing_is_claimed_without_a_statement(self):
        self.assertIsNone(txfn.classify(KE, 50.0, [2.0] * len(KE), None, {}))
        self.assertIsNone(txfn.classify(KE, None, casa(50.0), AV, {}))
        self.assertIsNone(txfn.classify(KE, 50.0, [], AV, {}))

    def test_the_avantage_dump_is_not_needed_for_the_casa_formula(self):
        got = txfn.classify(KE, 50.0, casa(50.0), None,
                            txfn.parse_casa_lines(LINES[:2]))
        self.assertEqual(got["kind"], "casa_formula")
        self.assertNotIn("vs_avantage", got)


class TestDescribe(unittest.TestCase):
    def test_wording(self):
        a = txfn.classify(KE, 50.0, avantage(50.0), AV, {})
        self.assertIn("Avantage polynomial", txfn.describe(a))
        f = txfn.classify(KE, 200.0, casa(200.0), AV,
                          txfn.parse_casa_lines(LINES[:2]))
        self.assertIn("recalculated", txfn.describe(f))
        self.assertIn("-0.116457", txfn.describe(f))
        col = [t * (1 + 0.3 * i / len(KE)) for i, t in enumerate(avantage(200.0))]
        m = txfn.classify(KE, 200.0, col, AV, txfn.parse_casa_lines(LINES))
        text = txfn.describe(m)
        self.assertIn("TransFunct-NEXSA_G2", text)
        self.assertIn("replacing Avantage", text)
        self.assertIn("departs from Avantage's", text)

    def test_a_tiny_spread_is_not_mentioned(self):
        col = [t * (1 + 1e-3 * (i % 2)) for i, t in enumerate(avantage(50.0))]
        m = txfn.classify(KE, 50.0, col, AV, txfn.parse_casa_lines(LINES))
        self.assertNotIn("departs", txfn.describe(m))

    def test_without_avantage_coefficients_nothing_is_said_to_be_replaced(self):
        m = txfn.classify(KE, 50.0, [2.0 + i for i in range(len(KE))], None,
                          txfn.parse_casa_lines(LINES[2:]))
        self.assertEqual(m["kind"], "casa_file")
        self.assertNotIn("replacing", txfn.describe(m))


def entry(kind=None, name="", area_t=1.0):
    row = {"area_t": area_t}
    if kind:
        row.update(tf_kind=kind, tf_file=name)
    return {"row": row}


def sample(*entries):
    return SimpleNamespace(levels=[SimpleNamespace(entries=list(entries))])


class TestQuantNotes(unittest.TestCase):
    def test_only_when_the_choice_is_on(self):
        s = sample(entry("casa_file", "tf.vms"))
        self.assertEqual(quantview.transmission_notes(
            s, SimpleNamespace(transmission=False)), [])

    def test_a_measured_file_is_named(self):
        s = sample(entry("casa_file", "tf.vms"), entry("casa_file", "tf.vms"))
        notes = quantview.transmission_notes(
            s, SimpleNamespace(transmission=True))
        self.assertEqual(len(notes), 1)
        self.assertIn("tf.vms", notes[0])

    def test_mixed_functions_are_flagged(self):
        s = sample(entry("avantage"), entry("casa_file", "tf.vms"),
                   entry(None), entry("casa_formula", area_t=None))
        notes = quantview.transmission_notes(
            s, SimpleNamespace(transmission=True))
        self.assertEqual(len(notes), 3)          # two functions + the warning
        self.assertIn("not on one scale", notes[-1])


@unittest.skipUnless(os.environ.get("XPS_TXFN_CORPUS"),
                     "set XPS_TXFN_CORPUS=<folder with the NaCl .vms files>")
class TestRealFiles(unittest.TestCase):
    def kinds(self, name):
        import readers
        f = readers.load_file(os.path.join(os.environ["XPS_TXFN_CORPUS"],
                                           name))
        return f, [(r.extra.get("tf_source") or {}).get("kind")
                   for r in f.regions]

    def find(self, name):
        root = os.environ["XPS_TXFN_CORPUS"]
        for dirpath, _d, files in os.walk(root):
            if name in files:
                return os.path.relpath(os.path.join(dirpath, name), root)
        self.skipTest(name + " not found")

    def test_the_three_exports_of_one_run(self):
        for name, kind in (("nacl-avg-transcurve.vms", "avantage"),
                           ("nacl-quant-recalculated-TF.vms", "casa_formula"),
                           ("nacl-quant-rel-TF-ESCALAB250Xi.vms",
                            "casa_file")):
            f, kinds = self.kinds(self.find(name))
            self.assertEqual(kinds, [kind] * 7, name)

    def test_the_metadata_says_so(self):
        f, _k = self.kinds(self.find("nacl-quant-rel-TF-ESCALAB250Xi.vms"))
        md = f.region_metadata(f.regions[-1])
        self.assertIn("TransFunct-NEXSA_G2", md["Transmission function"])
        self.assertNotIn("\\", md["Transmission function"])   # no path
        f, _k = self.kinds(self.find("nacl-avg-transcurve.vms"))
        md = f.region_metadata(f.regions[0])
        self.assertEqual(md["Transmission coefficients"],
                         "4.70503, 0.606984, -0.600862, 0.071214")

    def test_counts_are_untouched_and_not_flagged_for_rescaling(self):
        f, _k = self.kinds(self.find("nacl-avg-transcurve.vms"))
        for r in f.regions:
            self.assertNotIn("txfn", r.extra)      # casa_form must not apply


if __name__ == "__main__":
    unittest.main()
