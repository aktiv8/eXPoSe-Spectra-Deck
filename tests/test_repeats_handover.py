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
import repeats  # noqa: E402
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
        self.assertEqual(parts, [])
        self.assertTrue(any("no repeat groups" in n for n in notes))
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.suggest()
        parts, notes = ws.handover_parts(("quantification",))
        self.assertEqual([p.arc for p in parts], [ho.REPEAT_GROUPS_ARC])
        self.assertEqual(notes, [])
        rows = list(csv.reader(io.StringIO(parts[0].data.decode("utf-8-sig"))))
        self.assertEqual(rows[1][2], "Ti 2p")
        self.assertEqual(rows[1][6], "3")
        parts, _ = ws.handover_parts(("methods",))        # not asked for
        self.assertNotIn(ho.REPEAT_GROUPS_ARC, [p.arc for p in parts])


if __name__ == "__main__":
    unittest.main()
