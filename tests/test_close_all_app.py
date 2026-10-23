"""Close all files leaves nothing of the old session behind: no stale display
copies (keyed by id(region), which Python reuses) and no old workbook name
(real window; skipped without a display).

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

from test_menu_state_app import doc  # noqa: E402


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestCloseAll(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {"a.vms": doc("a.vms"), "b.vms": doc("b.vms")}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def open(self, name, fid):
        self.ws._add_file(name, file_id=fid, refresh=False)
        self.ws._finish_adding()

    def test_a_stale_display_copy_is_not_handed_to_the_next_file(self):
        ws = self.ws
        self.open("a.vms", "a")
        ra = ws.docs[0].regions[0]
        ws._display(ra)
        self.assertTrue(ws._disp_cache)
        ws.close_all()
        self.assertEqual(ws._disp_cache, {})
        self.open("b.vms", "b")
        rb = ws.docs[0].regions[0]
        ws._disp_cache[id(ra)] = "stale"      # as if b reused a's id
        ws.close_all()
        self.assertNotIn(id(ra), ws._disp_cache)
        self.assertIsNotNone(rb)

    def test_the_old_workbook_name_is_gone_after_close_all(self):
        ws = self.ws
        self.open("a.vms", "a")
        ws.wb_path = os.path.join("x", "old.xpscontainer")
        ws._wb_sig = ws._signature()
        with mock.patch.object(ee.messagebox, "askyesnocancel",
                               return_value=False):
            ws.close_all()
        self.assertIsNone(ws.wb_path)
        self.assertFalse(ws._wb_active())
        self.open("b.vms", "b")
        asked = []
        with mock.patch.object(
                ee.messagebox, "askyesnocancel",
                side_effect=lambda *a, **k: asked.append(a) or False):
            ws._confirm_discard()
        for a in asked:
            self.assertNotIn("old.xpscontainer", " ".join(map(str, a)))
        ws.close_all()


if __name__ == "__main__":
    unittest.main()
