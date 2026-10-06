"""Opening a CasaXPS ASCII export (casacsv.py) the ordinary way: with *Open*, by
drag-and-drop or with a folder it is applied to the fits of the spectra loaded
with it instead of being refused as an "unrecognised file". Needs a display and
matplotlib; skipped otherwise.

Run:  python -m unittest discover tests
"""

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

from test_workbook_csvimport_app import (  # noqa: E402
    HAVE_NP, _be_axis, _make_doc, _write_csv)


def has_curves(ws):
    return [fr.csv_curves is not None for d in ws.docs for r in d.regions
            if r.fit for fr in r.fit.regions]


@unittest.skipUnless(HAVE_MPL and HAVE_NP, "numpy / matplotlib / Tk not available")
class TestOpeningAnExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.dir = tempfile.mkdtemp(prefix="casaopen_")
        cls.vms = os.path.join(cls.dir, "a.vms")
        with open(cls.vms, "wb") as fh:
            fh.write(b"stand-in for an instrument file")
        cls.csv = os.path.join(cls.dir, "fit.csv")
        _write_csv(cls.csv, _be_axis())
        cls.other = os.path.join(cls.dir, "notes.csv")        # not an export
        with open(cls.other, "w") as fh:
            fh.write("name,value\n" + "\n".join(f"k{i},v{i}" for i in range(20)))
        cls._load, cls._reader = ee.load_file, ee.reader_for
        ee.load_file = lambda p: _make_doc(p)
        real = ee.reader_for

        def reader_for(p):                   # a.vms stands in for a real file
            return object if p.endswith(".vms") else real(p)
        ee.reader_for = reader_for
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror)
        cls.shown = []
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: cls.shown.append(a)

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror) = cls._boxes
        ee.load_file, ee.reader_for = cls._load, cls._reader
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def setUp(self):
        self.ws = ee.Workspace(self.root)
        self.shown.clear()

    def summaries(self):
        return [a[1] for a in self.shown if "CSV block(s) matched" in a[1]]

    def test_the_export_is_recognised_as_one(self):
        self.assertTrue(self.ws._is_casa_export(self.csv))
        self.assertFalse(self.ws._is_casa_export(self.other))
        self.assertFalse(self.ws._is_casa_export(self.vms))
        self.assertFalse(self.ws._is_casa_export(self.dir))

    def test_with_its_spectra_selected_together_in_either_order(self):
        for order in ([self.vms, self.csv], [self.csv, self.vms]):
            ws = ee.Workspace(self.root)
            self.shown.clear()
            ws._add_files(order)
            self.assertEqual([os.path.basename(d.path) for d in ws.docs],
                             ["a.vms"])               # the export is no spectrum
            self.assertEqual(has_curves(ws), [True])
            self.assertEqual(len(self.summaries()), 1)
            self.assertIn("1 of 1 CSV block(s) matched", self.summaries()[0])

    def test_opened_alone_after_the_spectra_it_belongs_to(self):
        self.ws._add_files([self.vms])
        self.assertEqual(has_curves(self.ws), [False])
        self.ws._add_files([self.csv])
        self.assertEqual(has_curves(self.ws), [True])
        self.assertEqual(len(self.ws.docs), 1)

    def test_alone_with_nothing_open_says_what_to_do(self):
        self.ws._add_files([self.csv])
        self.assertEqual(self.ws.docs, [])
        self.assertTrue(any("Open the .vms" in a[1] for a in self.shown))

    def test_a_folder_applies_the_export_and_does_not_count_it_as_skipped(self):
        folder = tempfile.mkdtemp(dir=self.dir)
        for p in (self.vms, self.csv, self.other):
            shutil.copy(p, folder)
        self.ws._open_folder_path(folder)
        self.assertEqual(has_curves(self.ws), [True])
        self.assertEqual(len(self.ws.docs), 1)
        # only notes.csv was skipped, not the export
        self.assertIn("1 other file(s)", self.ws.status.cget("text"))

    def test_a_folder_with_only_the_export_uses_what_is_already_open(self):
        self.ws._add_files([self.vms])
        folder = tempfile.mkdtemp(dir=self.dir)
        shutil.copy(self.csv, folder)
        self.ws._open_folder_path(folder)
        self.assertEqual(has_curves(self.ws), [True])

    def test_it_is_remembered_so_a_workbook_reapplies_it(self):
        self.ws._add_files([self.vms, self.csv])
        self.assertEqual([f.name for f in self.ws.casa_csv_imports],
                         ["fit.csv"])


if __name__ == "__main__":
    unittest.main()
