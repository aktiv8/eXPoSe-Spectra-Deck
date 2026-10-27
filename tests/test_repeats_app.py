"""Repeat groups in the Quantification tab (quant_ui.QuantPanel,
RepeatGroupsDialog) and the pure helpers behind it (quantview.repeat_*).

Run:  python -m unittest discover tests
"""

import copy
import csv
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

import quantview
import repeats
from readers import SpectrumFile
from test_casafit import HAVE_NP, fitted_region
from test_repeats import sample


class TestHelpers(unittest.TestCase):
    def test_ticks_reach_the_group(self):
        a = sample("S", [("C 1s", 1.0, 100.0, ()), ("O 1s", 1.0, 100.0, ())])
        b = sample("S - two", [("C 1s", 1.0, 300.0, ()),
                               ("O 1s", 1.0, 100.0, ())])
        groups = {"g": {"kind": "locations", "samples": [a.key, b.key]}}
        view = quantview.ViewState()
        (g,) = quantview.repeat_groups([a, b], groups, view)
        c = next(r for r in g.rows if r.region == "C 1s")
        self.assertAlmostEqual(c.mean, (50 + 75) / 2)
        view.toggle(quantview.entry_key(b, 0, 1), True)      # drop O 1s of b
        (g,) = quantview.repeat_groups([a, b], groups, view)
        c = next(r for r in g.rows if r.region == "C 1s")
        self.assertAlmostEqual(c.mean, (50 + 100) / 2)

    def test_table_columns_and_partial_rows(self):
        a = sample("S", [("Al 2p", 1.0, 100.0,
                          [("Al2O3", 60.0), ("AlOH", 40.0)])])
        b = sample("S - two", [("Al 2p", 1.0, 100.0,
                                [("Al2O3", 70.0), ("Al2O3b", 30.0)])])
        (g,) = repeats.compute([a, b], {"g": {"kind": "locations",
                                              "samples": [a.key, b.key]}})
        header, rows = quantview.repeat_table(g)
        self.assertEqual(header[:5], ["Region", "State", "Mean at %",
                                      "SD (n-1)", "n"])
        self.assertEqual(len(header), 7)
        by = {(c[0].strip(), c[1]): (c, p) for _k, c, p in rows}
        self.assertFalse(by[("Al 2p", "")][1])
        cells, partial = by[("Al 2p", "AlOH")]
        self.assertTrue(partial)
        self.assertEqual(cells[3], "–")
        self.assertEqual(cells[5:], ["40.00", "–"])

    def test_many_members_share_one_column(self):
        ss = [sample(f"S{i}", [("C 1s", 1.0, 100.0, ())]) for i in range(8)]
        (g,) = repeats.compute(ss, {"g": {"kind": "locations",
                                          "samples": [s.key for s in ss]}})
        header, rows = quantview.repeat_table(g)
        self.assertEqual(header[-1], "Values")
        self.assertEqual(len(rows[0][1]), 6)


@unittest.skipUnless(HAVE_MPL and HAVE_NP,
                     "numpy / matplotlib / Tk not available")
class TestInApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.dir = tempfile.mkdtemp(prefix="repeats_app_")
        cls.data = os.path.join(cls.dir, "a.vms")
        with open(cls.data, "wb") as fh:
            fh.write(b"stand-in for an instrument file")
        cls._load = ee.load_file
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror)
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror) = cls._boxes
        ee.load_file = cls._load
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def _ws(self, names=("S", "S - location 1", "S - location 2"),
            profile=False):
        regions = []
        for n in names:
            r = fitted_region()
            r.sample = n
            regions.append(r)
        if profile:
            for lv in range(2):
                r = fitted_region()
                r.sample = "S - depth"
                r.etch_level = lv
                regions.append(r)
        f = SpectrumFile()
        f.path = self.data
        f.format_name = "Test"
        f.regions = regions
        f.instrument = {}
        f._finish()
        ee.load_file = lambda p: f
        ws = ee.Workspace(self.root)
        ws._add_file(self.data)
        ws.checked |= {id(r) for r in ws.docs[0].regions}
        ws._render()
        return ws, ws.quant_panel

    def test_suggest_then_show_the_group(self):
        ws, panel = self._ws()
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        serial = ws._ann_serial
        self.assertEqual(dlg.suggest(), ["S"])
        self.assertGreater(ws._ann_serial, serial)
        self.assertEqual(ws.ann.repeat_groups["S"]["kind"], "locations")
        label = "Group: S (3 samples)"
        self.assertIn(label, panel.sample_box["values"])
        panel.sample_var.set(label)
        panel._show_sample()
        self.assertEqual(panel.group, "S")
        header = [panel.tree.heading(c, "text") for c in panel.tree["columns"]]
        self.assertEqual(header[:5], ["Region", "State", "Mean at %",
                                      "SD (n-1)", "n"])
        self.assertEqual(len(header), 8)
        rows = [panel.tree.item(i, "values") for i in panel.tree.get_children()]
        self.assertEqual(rows[0][0], "Ti 2p")
        self.assertEqual(rows[0][2:5], ("100.00", "0.00", "3"))

    def test_the_group_stays_shown_after_a_refresh_and_follows_ticks(self):
        ws, panel = self._ws()
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.suggest()
        panel.sample_var.set("Group: S (3 samples)")
        panel._show_sample()
        panel.refresh()
        self.assertEqual(panel.group, "S")
        member = next(s for s in panel.samples if s.label == "S - location 1")
        quantview.toggle_entry(member, 0, 0, panel.view)    # untick its region
        panel._show_sample()
        rows = [panel.tree.item(i, "values") for i in panel.tree.get_children()]
        self.assertEqual(rows[0][4], "2")
        panel.reset_ticks()
        rows = [panel.tree.item(i, "values") for i in panel.tree.get_children()]
        self.assertEqual(rows[0][4], "3")

    def test_export_writes_the_group_table(self):
        ws, panel = self._ws()
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.suggest()
        panel.sample_var.set("Group: S (3 samples)")
        panel._show_sample()
        path = os.path.join(self.dir, "groups.csv")
        self.assertEqual(panel.export_csv(path), path)
        with open(path, "rb") as fh:
            self.assertTrue(fh.read(3) == b"\xef\xbb\xbf")
        with open(path, encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0][:7], list(repeats.HEADER[:7]))
        self.assertEqual(rows[1][2], "Ti 2p")
        self.assertEqual(rows[1][6], "3")
        panel.sample_var.set("S")
        panel._show_sample()
        path2 = os.path.join(self.dir, "samples.csv")
        panel.export_csv(path2)
        with open(path2, encoding="utf-8-sig", newline="") as fh:
            self.assertEqual(next(csv.reader(fh))[0], "Sample")

    def test_dialog_editing(self):
        ws, panel = self._ws(profile=True)
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        self.assertIs(panel.open_groups(), dlg)           # one window
        name = dlg.new_group()
        self.assertEqual(name, "Group")
        self.assertEqual(ws.ann.repeat_groups, {})        # empty: not stored
        keys = {s.label: s.key for s in panel.samples}
        self.assertTrue(dlg.toggle_sample(keys["S"]))
        self.assertTrue(dlg.toggle_sample(keys["S - location 1"]))
        self.assertFalse(dlg.toggle_sample(keys["S - depth"]))   # a profile
        self.assertEqual(ws.ann.repeat_groups["Group"]["samples"],
                         [keys["S"], keys["S - location 1"]])
        self.assertTrue(dlg.set_name("Spots"))
        self.assertEqual(list(ws.ann.repeat_groups), ["Spots"])
        self.assertEqual(dlg.new_group("Spots"), "Spots 2")
        self.assertFalse(dlg.set_name("Spots"))           # taken
        self.assertFalse(dlg.set_kind("nonsense"))
        dlg.delete_group("Spots 2")
        dlg.select("Spots")
        self.assertTrue(dlg.set_kind("repeats"))
        self.assertEqual(ws.ann.repeat_groups["Spots"]["kind"], "repeats")
        self.assertEqual(dlg.suggest(), ["S"])            # not "Spots"
        self.assertEqual(sorted(ws.ann.repeat_groups), ["S", "Spots"])
        dlg.delete_group("S")
        dlg.delete_group("Spots")
        self.assertEqual(ws.ann.repeat_groups, {})
        self.assertNotIn("Group: Spots", " ".join(panel.sample_box["values"]))

    def test_a_missing_sample_is_listed_and_not_counted(self):
        ws, panel = self._ws()
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.new_group("G")
        keys = {s.label: s.key for s in panel.samples}
        dlg.toggle_sample(keys["S"])
        dlg.toggle_sample("f/gone")                       # not a candidate
        self.assertEqual(ws.ann.repeat_groups["G"]["samples"], [keys["S"]])
        ws.ann.set_repeat_group("G", "locations", [keys["S"], "f/gone"])
        panel.refresh()
        dlg2 = copy.copy(dlg)
        dlg2.drafts = {"G": dict(ws.ann.repeat_groups["G"])}
        labels = [l for _k, l, _u in dlg2.candidates()]
        self.assertIn("f/gone (missing)", labels)
        panel.sample_var.set("Group: G (1 sample)")
        panel._show_sample()
        notes = panel.notes.get("1.0", "end")
        self.assertIn("f/gone", notes)

    def test_groups_are_saved_in_the_workbook_annotations(self):
        ws, panel = self._ws()
        dlg = panel.open_groups()
        self.addCleanup(dlg.destroy)
        dlg.suggest()
        import annotations
        back = annotations.Annotations.from_json(ws.ann.to_json())
        self.assertEqual(back.repeat_groups, ws.ann.repeat_groups)


if __name__ == "__main__":
    unittest.main()
