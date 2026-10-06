"""Menu entries that need spectra are dimmed while none are loaded, and the
file shortcuts exist and keep out of text boxes (real window; skipped without a
display).

Run:  python -m unittest discover tests
"""

import math
import os
import sys
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

from readers import Region, SpectrumFile  # noqa: E402


def doc(path):
    n = 41
    e = [292.0 - i * 12.0 / (n - 1) for i in range(n)]
    c = [100 + 900 * math.exp(-((x - 286) / 1.0) ** 2) for x in e]
    r = Region(name="C 1s", index=0, offset=0, energy=e, counts=c,
               decodable=True, sample="S", photon_energy=1486.6,
               pass_energy=20.0, dwell=0.1, step=0.2, source=path,
               count_units="counts/s")
    f = SpectrumFile()
    f.path = path
    f.format_name = "Test"
    f.regions = [r]
    f.instrument = {}
    f._finish()
    return f


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestMenuState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {"a.vms": doc("a.vms")}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def setUp(self):
        for d in list(self.ws.docs):
            self.ws.docs.remove(d)
        self.ws.file_ids.clear()
        self.ws.region_parser.clear()
        self.ws._populate_tree()

    def entries(self):
        """{label: state} of every command in the menu bar's menus."""
        out = {}
        for menu, i in self.ws._data_items:
            out[menu.entrycget(i, "label")] = menu.entrycget(i, "state")
        return out

    def menu_labels(self):
        bar = self.ws.menubar
        labels = {}
        for k in range(bar.index("end") + 1):
            if bar.type(k) != "cascade":
                continue
            sub = self.ws.root.nametowidget(bar.entrycget(k, "menu"))
            for i in range(sub.index("end") + 1):
                if sub.type(i) == "command":
                    labels[sub.entrycget(i, "label")] = (sub, i)
        return labels

    def test_every_listed_entry_exists_in_a_menu(self):
        labels = self.menu_labels()
        for group in (ee.FILE_NEEDS_DATA, ee.WORKBOOK_NEEDS_DATA,
                      ee.TOOLS_NEEDS_DATA, ee.VIEW_NEEDS_DATA):
            for name in group:
                self.assertIn(name, labels)
        self.assertEqual(len(self.ws._data_items),
                         len(ee.FILE_NEEDS_DATA | ee.WORKBOOK_NEEDS_DATA
                             | ee.TOOLS_NEEDS_DATA | ee.VIEW_NEEDS_DATA))

    def test_entries_that_need_data_are_dimmed_until_a_file_is_open(self):
        self.assertEqual(set(self.entries().values()), {"disabled"})
        self.ws._add_file("a.vms", file_id="a", refresh=False)
        self.ws._finish_adding()
        self.assertEqual(set(self.entries().values()), {"normal"})
        self.ws.close_all()
        self.assertEqual(set(self.entries().values()), {"disabled"})

    def test_opening_and_quitting_are_never_dimmed(self):
        labels = self.menu_labels()
        for name in ("Open spectra file(s)…", "Open folder…", "Quit",
                     "New workbook", "Open workbook…", ee.FORGET_DUP_LABEL,
                     "Plot style…", "About eXPoSe SpectraDeck…"):
            menu, i = labels[name]
            self.assertEqual(menu.entrycget(i, "state"), "normal", name)

    def test_the_file_shortcuts_are_bound_and_shown(self):
        bound = set(self.ws.root.bind())
        for seq in ("<Control-Key-o>", "<Control-Key-O>", "<Control-Key-e>",
                    "<Control-Key-q>", "<Control-Key-s>"):       # Tk's spelling
            self.assertIn(seq, bound)
        labels = self.menu_labels()
        for name, accel in (("Open spectra file(s)…", "Ctrl+O"),
                            ("Export spectra (choose regions/levels)…",
                             "Ctrl+E"), ("Quit", "Ctrl+Q"),
                            ("Save workbook", "Ctrl+S")):
            menu, i = labels[name]
            self.assertEqual(menu.entrycget(i, "accelerator"), accel, name)

    def test_a_shortcut_stands_aside_in_a_text_box_and_runs_elsewhere(self):
        calls = []
        run = self.ws._shortcut(lambda: calls.append(1))
        ws_root = self.ws.root
        entry = tk.Entry(ws_root)
        self.addCleanup(entry.destroy)
        old = ws_root.focus_get
        try:
            ws_root.focus_get = lambda: entry
            self.assertIsNone(run())
            self.assertEqual(calls, [])
            ws_root.focus_get = lambda: None
            self.assertEqual(run(), "break")
            self.assertEqual(calls, [1])
        finally:
            ws_root.focus_get = old


if __name__ == "__main__":
    unittest.main()
