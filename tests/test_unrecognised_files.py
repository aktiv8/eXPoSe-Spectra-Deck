"""A file no reader recognises: say what it looks like and what is supported;
opening a folder says how many files it skipped.

Run:  python -m unittest discover tests
"""

import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import readers  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_APP = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_APP = None, None, False


class TestMessage(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, True)

    def error(self, name, data):
        path = os.path.join(self.dir, name)
        with open(path, "wb") as fh:
            fh.write(data)
        with self.assertRaises(readers.UnsupportedFormat) as cm:
            readers.reader_for(path)
        return str(cm.exception)

    def test_it_names_the_file_what_it_looks_like_and_the_supported_formats(self):
        msg = self.error("a.csv", b"Energy,Counts\n1,2\n")
        self.assertTrue(msg.startswith("Unrecognised file format: a.csv"))
        self.assertIn("text file starting", msg)
        self.assertIn("Energy,Counts", msg)
        for name in readers.supported_names():
            self.assertIn(name, msg)

    def test_binary_and_empty_files_are_told_apart(self):
        self.assertIn("binary", self.error("b.bin", bytes(range(0, 200))))
        self.assertIn("empty", self.error("c.txt", b""))

    def test_a_long_first_line_is_trimmed(self):
        msg = self.error("d.dat", b"x" * 500 + b"\n")
        self.assertIn("…", msg)
        self.assertLess(len(msg), 800)

    def test_the_first_line_is_still_the_old_message(self):
        msg = self.error("e.csv", b"a,b\n")
        self.assertEqual(msg.splitlines()[0], "Unrecognised file format: e.csv")


@unittest.skipUnless(HAVE_APP, "matplotlib / Tk not available")
class TestFolder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, True)
        for name in ("a.vms", "b.vms", "notes.csv", "log.dat"):
            with open(os.path.join(self.dir, name), "wb") as fh:
                fh.write(b"x")
        self.added, self.shown = [], []
        self._patches = [
            (ee, "reader_for", ee.reader_for),
            (ee.messagebox, "showinfo", ee.messagebox.showinfo),
            (ee.Workspace, "_add_files", ee.Workspace._add_files),
            (ee.Workspace, "_scan_casa_quant", ee.Workspace._scan_casa_quant)]
        self.addCleanup(self._unpatch)

        def fake_reader(path):
            if not path.endswith(".vms"):
                raise ee.UnsupportedFormat("no")
            return object

        ee.reader_for = fake_reader
        ee.messagebox.showinfo = lambda *a, **k: self.shown.append(a)
        ee.Workspace._add_files = lambda ws, paths: self.added.extend(paths)
        ee.Workspace._scan_casa_quant = lambda ws, folder: None
        self.ws.status.config(text="2 files loaded")

    def _unpatch(self):
        for obj, name, value in self._patches:
            setattr(obj, name, value)

    def test_skipped_files_are_counted_in_the_status_bar(self):
        self.ws._open_folder_path(self.dir)
        self.assertEqual([os.path.basename(p) for p in self.added],
                         ["a.vms", "b.vms"])
        text = self.ws.status.cget("text")
        self.assertIn("2 other file(s)", text)
        self.assertTrue(text.startswith("2 files loaded"))

    def test_a_folder_with_nothing_readable_lists_what_was_skipped(self):
        for n in ("a.vms", "b.vms"):
            os.remove(os.path.join(self.dir, n))
        self.ws._open_folder_path(self.dir)
        self.assertEqual(self.added, [])
        self.assertEqual(len(self.shown), 1)
        text = self.shown[0][1]
        self.assertIn("No recognised spectra files", text)
        self.assertIn("2 file(s) were not recognised", text)
        self.assertIn("notes.csv", text)


if __name__ == "__main__":
    unittest.main()
