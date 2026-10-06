"""Importing plain column-text files in the real window: the mapping dialog,
what the Workspace does with its answer, and that the way a file was read
survives a workbook save and reopen (needs a display and matplotlib; skipped
otherwise).

Run:  python -m unittest discover tests
"""

import math
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

import columntext as ct  # noqa: E402
import workbook as wbk  # noqa: E402
from test_column_text import peak, text  # noqa: E402


def items_for(paths):
    out = []
    for p in paths:
        with open(p, "rb") as fh:
            t = ct.read_table(ct.decode(fh.read()))
        out.append((p, t, ct.guess_options(t, p)))
    return out


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestColumnImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror, ee.messagebox.askyesno)
        cls.shown = []
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: cls.shown.append(a)
        ee.messagebox.askyesno = lambda *a, **k: True
        cls.dir = tempfile.mkdtemp()
        be, y = peak(60, 280, 292, 285)
        o_be, o_y = peak(60, 521, 541, 531)
        cls.c1s = cls.write("c1s.asc", text(zip(be, y),
                                             ["BE_C1s", "CPS_C1s"]))
        cls.o1s = cls.write("o1s.txt", text(zip(o_y, o_be)))     # no header
        ke = [1486.6 - b for b in be][::-1]
        cls.ke = cls.write("ke.csv", text(zip(ke, y[::-1]),
                                          ["KE_C 1s", "CPS_C 1s"], sep=","))
        cls.multi = cls.write("multi.csv", text(
            [(a, b, b * 2) for a, b in zip(be, y)],
            ["BE", "Probe A", "Probe B"], sep=","))
        cls.be, cls.y = be, y

    @classmethod
    def write(cls, name, body):
        path = os.path.join(cls.dir, name)
        with open(path, "wb") as fh:
            fh.write(body.encode("utf-8"))
        return path

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror, ee.messagebox.askyesno) = cls._boxes
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def setUp(self):
        self.ws = ee.Workspace(self.root)
        self.shown.clear()

    def dialog(self, paths, defaults=None):
        from columnimport_ui import ColumnImportDialog
        dlg = ColumnImportDialog(self.root, self.ws, items_for(paths),
                                 defaults)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        return dlg

    # -- the dialog ----------------------------------------------------------------
    def test_it_lists_each_file_with_what_was_found(self):
        dlg = self.dialog([self.c1s, self.o1s, self.ke])
        rows = [dlg.tree.item(str(i)) for i in range(3)]
        self.assertEqual([r["text"] for r in rows],
                         ["c1s.asc", "o1s.txt", "ke.csv"])
        v = [r["values"] for r in rows]
        self.assertEqual([str(x) for x in v[0]],
                         ["60", "BE_C1s", "CPS_C1s", "Binding energy",
                          "counts/s", "C 1s"])             # from the header
        self.assertEqual(v[1][1], "Column 2")              # intensity first
        self.assertEqual(v[1][2], "Column 1")
        self.assertEqual(v[1][5], "O 1s ?")                # suggested, flagged
        self.assertEqual(v[2][3], "Kinetic energy")

    def test_a_suggested_name_stops_being_one_when_the_user_types(self):
        dlg = self.dialog([self.o1s])
        self.assertTrue(dlg.rows[0]["suggested"])
        dlg.set_row(0, name="O 1s")                        # same text: kept
        dlg.set_row(0, name="Oxygen")
        self.assertFalse(dlg.rows[0]["suggested"])
        self.assertEqual(dlg.tree.item("0")["values"][5], "Oxygen")

    def test_changing_the_columns_axis_and_unit(self):
        dlg = self.dialog([self.o1s])
        dlg.set_row(0, axis="KE", units="counts")
        o = dlg.rows[0]["opts"]
        self.assertEqual((o["axis"], o["units"]), ("KE", "counts"))
        dlg.set_row(0, energy_col=0, intensity=1)          # swap the columns
        o = dlg.rows[0]["opts"]
        self.assertEqual((o["energy_col"], o["intensity_cols"]), (0, [1]))
        dlg.set_row(0, energy_col=1)                       # the other way
        self.assertEqual(dlg.rows[0]["opts"]["intensity_cols"], [0])

    def test_every_other_column_can_be_its_own_spectrum(self):
        dlg = self.dialog([self.multi])
        o = dlg.rows[0]["opts"]
        self.assertEqual((o["intensity_cols"], o["names"]),
                         ([1, 2], ["Probe A", "Probe B"]))
        self.assertEqual(dlg.tree.item("0")["values"][5], "2 spectra")
        dlg.set_row(0, intensity=2)
        self.assertEqual(dlg.rows[0]["opts"]["intensity_cols"], [2])
        dlg.set_row(0, intensity="all")
        self.assertEqual(dlg.rows[0]["opts"]["intensity_cols"], [1, 2])

    def test_the_boxes_under_the_table_follow_the_selection(self):
        dlg = self.dialog([self.c1s, self.ke])
        dlg.tree.selection_set("1")
        dlg._on_select()
        self.assertEqual(dlg.v_axis.get(), "Kinetic energy")
        self.assertEqual(dlg.v_energy.get(), "1: KE_C 1s")
        dlg.v_units.set("counts")
        dlg._commit()
        self.assertEqual(dlg.rows[1]["opts"]["units"], "counts")
        self.assertEqual(dlg.rows[0]["opts"]["units"], "counts/s")

    def test_axis_and_unit_can_be_copied_to_every_file(self):
        dlg = self.dialog([self.c1s, self.o1s, self.ke])
        dlg.set_row(0, axis="KE", units="counts")
        dlg.tree.selection_set("0")
        dlg.axis_to_all()
        for row in dlg.rows:
            self.assertEqual((row["opts"]["axis"], row["opts"]["units"]),
                             ("KE", "counts"))

    def test_accept_returns_the_options_with_what_is_shared(self):
        dlg = self.dialog([self.c1s, self.o1s])
        dlg.set_shared(sample="Probe 7", photon="Al Kα (1486.6 eV)",
                       pass_energy=20)
        dlg.accept()
        res = dlg.result
        self.assertEqual(sorted(res), sorted([self.c1s, self.o1s]))
        for o in res.values():
            self.assertEqual((o["sample"], o["photon_energy"], o["anode"],
                              o["pass_energy"]),
                             ("Probe 7", 1486.6, "Al Kα", 20.0))
        self.assertEqual(res[self.o1s]["names"], ["O 1s"])   # the suggestion
        self.assertEqual(res[self.o1s]["energy_col"], 1)

    def test_unknown_photon_energy_stays_unset_and_cancel_returns_none(self):
        dlg = self.dialog([self.c1s])
        dlg.set_shared(photon="Unknown", pass_energy="")
        dlg.accept()
        o = dlg.result[self.c1s]
        self.assertIsNone(o["photon_energy"])
        self.assertIsNone(o["pass_energy"])
        self.assertEqual(o["anode"], "")
        dlg2 = self.dialog([self.c1s])
        dlg2.cancel()
        self.assertIsNone(dlg2.result)

    def test_the_last_answers_are_offered_next_time(self):
        dlg = self.dialog([self.c1s], {"sample": "S1", "photon_energy": 1253.6,
                                       "pass_energy": 40})
        self.assertEqual(dlg.v_sample.get(), "S1")
        self.assertEqual(dlg.v_photon.get(), "Mg Kα (1253.6 eV)")
        self.assertEqual(dlg.v_pass.get(), "40")

    def test_a_file_with_lines_under_it_says_so(self):
        be, y = peak(30)
        p = self.write("bad.asc", text(zip(be, y), ["BE", "Counts/s"],
                                       sep=",", post=["65406.7,", "65298.3,"]))
        dlg = self.dialog([p])
        self.assertIn("2 line(s) after the table", dlg.notes.cget("text"))

    # -- the Workspace ---------------------------------------------------------------
    def test_files_of_numbers_go_through_the_dialog_and_are_read_as_answered(self):
        seen = []

        def answer(items):
            seen.extend(p for p, _t, _o in items)
            return {p: dict(o, sample="Probe 1", photon_energy=1486.6,
                            anode="Al Kα") for p, _t, o in items}

        self.ws._ask_column_import = answer
        self.ws._add_files([self.c1s, self.ke])
        self.assertEqual(sorted(seen), sorted([self.c1s, self.ke]))
        by = {os.path.basename(d.path): d for d in self.ws.docs}
        r = by["c1s.asc"].regions[0]
        self.assertEqual((r.name, r.sample, r.photon_energy),
                         ("C 1s", "Probe 1", 1486.6))
        k = by["ke.csv"].regions[0]                       # converted with hν
        self.assertEqual(k.energy_label, "Binding Energy")
        for a, b in zip(k.energy, self.be):
            self.assertAlmostEqual(a, b, 6)
        self.assertEqual(len(self.ws.ann.imports), 2)

    def test_cancelling_leaves_those_files_out_and_the_others_in(self):
        self.ws._ask_column_import = lambda items: None
        self.ws._add_files([self.c1s, self.o1s])
        self.assertEqual(self.ws.docs, [])
        self.assertEqual(self.ws.ann.imports, {})

    def test_the_answers_are_remembered_for_the_next_dialog(self):
        self.ws._ask_column_import = lambda items: {
            p: dict(o, sample="Probe 9", photon_energy=1253.6,
                    pass_energy=40.0) for p, _t, o in items}
        self.ws._add_files([self.c1s])
        self.assertEqual(self.ws.cfg["column_import"],
                         {"sample": "Probe 9", "photon_energy": 1253.6,
                          "pass_energy": 40.0})

    def test_a_loaded_file_is_not_asked_about_again(self):
        asked = []
        self.ws._ask_column_import = lambda items: (
            asked.append(len(items)) or {p: o for p, _t, o in items})
        self.ws._add_files([self.c1s])
        self.ws._add_files([self.c1s])
        self.assertEqual(asked, [1])

    def test_the_way_a_file_was_read_survives_saving_and_reopening(self):
        self.ws._ask_column_import = lambda items: {
            p: dict(o, sample="Probe 1", photon_energy=1486.6,
                    anode="Al Kα", pass_energy=20.0, names=["C 1s"])
            for p, _t, o in items}
        self.ws._add_files([self.c1s])
        self.ws.checked = {id(r) for d in self.ws.docs for r in d.regions}
        path = os.path.join(self.dir, "book" + wbk.EXT)
        self.ws.wb_path = path
        self.assertTrue(self.ws.save_workbook())

        other = ee.Workspace(self.root)
        asked = []
        other._ask_column_import = lambda items: asked.append(items)
        other.open_workbook(path)
        self.assertEqual(asked, [])                       # not asked again
        (doc,) = other.docs
        r = doc.regions[0]
        self.assertEqual((r.name, r.sample, r.photon_energy, r.pass_energy,
                          r.anode), ("C 1s", "Probe 1", 1486.6, 20.0, "Al Kα"))
        self.assertEqual(r.energy, self.be)
        self.assertEqual(len(other.ann.imports), 1)

    def test_other_formats_never_reach_the_dialog(self):
        called = []
        self.ws._ask_column_import = lambda items: called.append(items)
        paths, imports = self.ws._column_import(
            [os.path.join(self.dir, "missing.vms"), self.dir])
        self.assertEqual(called, [])
        self.assertEqual(imports, {})
        self.assertEqual(len(paths), 2)

    def test_the_open_dialog_lists_the_format(self):
        names = [n for n, _p in ee.supported_patterns()]
        self.assertIn("Column text (CSV / ASC / TXT)", names)
        types = dict(self.ws._open_filetypes())
        self.assertIn("*.csv", types["Column text (CSV / ASC / TXT)"])
        self.assertIn("*.asc", types["All supported spectra files"])


if __name__ == "__main__":
    unittest.main()
