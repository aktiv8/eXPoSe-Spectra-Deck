"""The repeat-group CSV in the hand-over ZIP (handover.repeat_group_parts and
Workspace.handover_parts).

Run:  python -m unittest discover tests
"""

import csv
import io
import os
import shutil
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import handover as ho  # noqa: E402
import quant  # noqa: E402
import repeats  # noqa: E402
import resultspages as rp  # noqa: E402
from test_repeats import sample  # noqa: E402


def groups():
    ss = [sample(f"S - {i}", [("C 1s", 1.0, 100.0 + 20 * i, ()),
                              ("O 1s", 1.0, 100.0, ())]) for i in range(3)]
    return repeats.compute(ss, {"S": {"kind": "locations",
                                      "samples": [s.key for s in ss]}})


class TestPart(unittest.TestCase):
    def test_the_csv(self):
        (part,) = ho.repeat_group_parts(groups())
        self.assertEqual(part.arc, ho.REPEAT_GROUPS_ARC)
        self.assertTrue(part.data.startswith(b"\xef\xbb\xbf"))
        rows = list(csv.reader(io.StringIO(part.data.decode("utf-8-sig"))))
        self.assertEqual(rows[0], list(repeats.HEADER))
        self.assertEqual(rows[1][:3], ["S", "different locations on the "
                                       "sample", "C 1s"])
        self.assertEqual(rows[1][6], "3")
        self.assertIn("S - 0=", rows[1][7])

    def test_nothing_when_there_are_no_groups(self):
        self.assertEqual(ho.repeat_group_parts([]), [])
        g = groups()[0]
        g.rows = []
        self.assertEqual(ho.repeat_group_parts([g]), [])

    def test_in_the_zip_with_a_checksum_and_a_readme_paragraph(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        out = os.path.join(d, "o.zip")
        names = ho.write_zip(out, "pkg", ho.repeat_group_parts(groups()),
                             {"title": "T"}, [], "")
        self.assertIn("pkg/quantification/quantification_repeat_groups.csv",
                      names)
        with zipfile.ZipFile(out) as zf:
            self.assertIn("quantification/quantification_repeat_groups.csv",
                          zf.read("pkg/SHA256SUMS.txt").decode())
            readme = zf.read("pkg/README.txt").decode()
        self.assertIn("About the repeat groups", readme)

    def test_no_paragraph_without_the_file(self):
        text = ho.readme_text({"title": "T"}, [ho.Part("a.txt", "a", data=b"x")])
        self.assertNotIn("repeat groups", text)


class TestQuantTable(unittest.TestCase):
    def results(self):
        ss = [sample(f"Film {i}", [("C 1s", 1.0, 100.0 + 20 * i, ()),
                                   ("Ti 2p", 2.0, 100.0,
                                    [("Ti-C", 60.0), ("Ti-O", 40.0)])])
              for i in range(2)]
        ss.append(rp.Sample("casaxps:X", "X", kind="casaxps"))   # no fit
        return rp.Results(samples=ss)

    def rows(self, part):
        return list(csv.reader(io.StringIO(part.data.decode("utf-8-sig"))))

    def test_the_table(self):
        (part,) = ho.quantification_parts(self.results())
        self.assertEqual(part.arc, ho.QUANT_ARC)
        self.assertEqual(part.data[:3], bytes([0xEF, 0xBB, 0xBF]))
        rows = self.rows(part)
        self.assertEqual(rows[0], list(quant.CSV_HEADER))
        regions = [r for r in rows[1:] if r[4]]
        self.assertEqual([(r[0], r[3]) for r in regions],
                         [("Film 0", "C 1s"), ("Film 0", "Ti 2p"),
                          ("Film 1", "C 1s"), ("Film 1", "Ti 2p")])
        # each sample is its own total: the at % of one sum to 100
        for name in ("Film 0", "Film 1"):
            self.assertAlmostEqual(
                sum(float(r[8]) for r in regions if r[0] == name), 100, 3)
        states = [r for r in rows[1:] if r[9]]
        self.assertEqual({r[9] for r in states}, {"Ti-C", "Ti-O"})

    def test_a_region_left_out_by_the_report_says_why(self):
        res = self.results()
        lv = res.samples[0].levels[0]
        lv.include = [True, False]
        rows = self.rows(ho.quantification_parts(res)[0])
        ti = next(r for r in rows if r[0] == "Film 0" and r[3] == "Ti 2p"
                  and r[4])
        self.assertEqual((ti[8], ti[11]), ("", "not included"))
        c = next(r for r in rows if r[0] == "Film 0" and r[3] == "C 1s")
        self.assertEqual(c[8], "100")

    def test_a_sample_left_out_of_the_report_is_left_out(self):
        res = self.results()
        rows = self.rows(ho.quantification_parts(
            res, skip=[res.samples[0].key])[0])
        self.assertEqual({r[0] for r in rows[1:]}, {"Film 1"})

    def test_nothing_without_a_fit(self):
        self.assertEqual(ho.quantification_parts(None), [])
        self.assertEqual(ho.quantification_parts(rp.Results(samples=[
            rp.Sample("casaxps:X", "X", kind="casaxps")])), [])

    def test_readme_paragraph_only_with_the_file(self):
        (part,) = ho.quantification_parts(self.results())
        text = ho.readme_text({"title": "T"}, [part])
        self.assertIn("About the quantification", text)
        self.assertIn("sensitivity factor", text)
        self.assertNotIn("About the quantification", ho.readme_text(
            {"title": "T"}, [ho.Part("a.txt", "a", data=b"x")]))


class TestInApp(unittest.TestCase):
    """Through the Workspace; the fixture is the one the tab tests use."""

    def test_handover_parts(self):
        import test_repeats_app as t
        if not (t.HAVE_MPL and t.HAVE_NP):
            self.skipTest("numpy / matplotlib / Tk not available")
        t.TestInApp.setUpClass()
        self.addCleanup(t.TestInApp.tearDownClass)
        case = t.TestInApp()
        ws, panel = case._ws()
        parts, notes = ws.handover_parts(("quantification",))
        # no groups made yet: the atomic % table only, and no complaint
        self.assertEqual([p.arc for p in parts], [ho.QUANT_ARC])
        self.assertEqual(notes, [])
        table = list(csv.reader(io.StringIO(
            parts[0].data.decode("utf-8-sig"))))
        self.assertEqual(table[0], list(quant.CSV_HEADER))
        self.assertEqual({r[0] for r in table[1:]},
                         {"S", "S - location 1", "S - location 2"})
        self.assertEqual(table[1][3], "Ti 2p")
        self.assertEqual(table[1][8], "100")
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.suggest()
        parts, notes = ws.handover_parts(("quantification",))
        self.assertEqual([p.arc for p in parts],
                         [ho.QUANT_ARC, ho.REPEAT_GROUPS_ARC])
        self.assertEqual(notes, [])
        rows = list(csv.reader(io.StringIO(
            parts[1].data.decode("utf-8-sig"))))
        self.assertEqual(rows[1][2], "Ti 2p")
        self.assertEqual(rows[1][6], "3")
        parts, _ = ws.handover_parts(("methods",))        # not asked for
        self.assertNotIn(ho.REPEAT_GROUPS_ARC, [p.arc for p in parts])
        self.assertNotIn(ho.QUANT_ARC, [p.arc for p in parts])


if __name__ == "__main__":
    unittest.main()
