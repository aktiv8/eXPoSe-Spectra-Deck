"""The foldable view bar above the plot (needs a display; skipped otherwise).

Run:  python -m unittest tests.test_viewbar_app
"""

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestViewBar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls._load_config = ee.load_config
        ee.load_config = lambda: {}
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_config = cls._load_config
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def setUp(self):
        self.ws._set_viewbar_collapsed(False)

    def test_open_by_default(self):
        ws = self.ws
        self.assertFalse(ws.viewbar_collapsed)
        self.assertIn(ws.viewbody, ws.viewbar.pack_slaves())
        self.assertEqual(ws.viewbar_sum.cget("text"), "")

    def test_fold_hides_body_and_shows_summary(self):
        ws = self.ws
        ws.view_var.set("Heatmap")
        ws._set_viewbar_collapsed(True)
        self.assertNotIn(ws.viewbody, ws.viewbar.pack_slaves())
        self.assertIn("Heatmap", ws.viewbar_sum.cget("text"))
        self.assertFalse(ws.viewbar_var.get())
        ws.view_var.set("Stack")
        ws._on_view_changed()
        self.assertIn("Stack", ws.viewbar_sum.cget("text"))
        ws._toggle_viewbar()
        self.assertIn(ws.viewbody, ws.viewbar.pack_slaves())
        self.assertEqual(ws.viewbar_sum.cget("text"), "")
        self.assertTrue(ws.viewbar_var.get())

    def test_saved_in_config(self):
        ws = self.ws
        ws._set_viewbar_collapsed(True)
        self.assertTrue(ws.cfg["viewbar_collapsed"])
        ws._set_viewbar_collapsed(False)
        self.assertFalse(ws.cfg["viewbar_collapsed"])

    def test_not_part_of_saved_figure_state(self):
        self.assertNotIn("viewbar_collapsed", self.ws.capture_state())

    def test_repeated_toggles_do_not_hang(self):
        for _ in range(4):
            self.ws._toggle_viewbar()
            self.root.update_idletasks()

    def test_menu_entry(self):
        m = self.ws.view_menu
        labels = [m.entrycget(i, "label") for i in range(m.index("end") + 1)
                  if m.type(i) in ("command", "checkbutton")]
        self.assertIn("Show view options", labels)


if __name__ == "__main__":
    unittest.main()
