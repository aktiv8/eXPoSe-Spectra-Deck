"""Notebook tabs must not be squeezed (tabfit.fit_tabs): the ISS / REELS and
Plot style dialogs, and a notebook in a pane that cannot grow.

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import tkinter as tk
    from tkinter import ttk
    import spectradeck as ee
    import tabfit
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

from test_iss_kratos import iss_object  # noqa: E402


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestTabFit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        import tempfile
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.dir = tempfile.mkdtemp()
        cls.kal = os.path.join(cls.dir, "iss.kal")
        with open(cls.kal, "w") as fh:
            fh.write("Dataset filename          = iss.dset\n"
                     + iss_object("ISS 1"))
        cls.ws = ee.Workspace(cls.root)
        cls.ws._add_files([cls.kal])
        cls.ws.checked = {id(cls.ws.docs[0].regions[0])}

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()
        import matplotlib
        import shutil
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def settle(self, win):
        for _ in range(3):
            win.update()

    def test_the_iss_reels_dialog_shows_every_tab_unsqueezed(self):
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(self.root, self.ws)
        self.addCleanup(dlg.destroy)
        self.settle(dlg)
        self.assertFalse(tabfit.is_clipped(dlg.nb))
        self.assertGreaterEqual(dlg.nb.winfo_width(), tabfit.needed_width(dlg.nb))

    def test_the_plot_style_dialog_shows_every_tab_unsqueezed(self):
        import plotstyle_ui
        dlg = plotstyle_ui.PlotStyleDialog(self.root, self.ws)
        self.addCleanup(dlg.destroy)
        self.settle(dlg)
        self.assertFalse(tabfit.is_clipped(dlg.nb))

    def test_a_dialog_grows_when_even_compact_tabs_do_not_fit(self):
        top = tk.Toplevel(self.root)
        self.addCleanup(top.destroy)
        nb = ttk.Notebook(top)
        nb.pack(fill="both", expand=True)
        for i in range(6):
            nb.add(ttk.Frame(nb), text="A fairly long tab name %d" % i)
        top.geometry("300x200")
        tabfit.fit_tabs(nb, top)
        self.settle(top)
        self.assertGreater(top.winfo_width(), 300)
        self.assertFalse(tabfit.is_clipped(nb))

    def test_hidden_tabs_take_no_room(self):
        top = tk.Toplevel(self.root)
        self.addCleanup(top.destroy)
        nb = ttk.Notebook(top)
        nb.pack(fill="both", expand=True)
        a, b = ttk.Frame(nb), ttk.Frame(nb)
        nb.add(a, text="One")
        nb.add(b, text="A very long hidden tab label indeed")
        shown = tabfit.needed_width(nb)
        nb.tab(b, state="hidden")
        self.assertLess(tabfit.needed_width(nb), shown)


if __name__ == "__main__":
    unittest.main()
