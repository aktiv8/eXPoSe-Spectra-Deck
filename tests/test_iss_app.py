"""The ISS / REELS dialog on Kratos ISS spectra in the real window: what it
offers, what "Use for this file" saves, and that the methods text states the
beam only after that (needs a display and matplotlib; skipped otherwise).

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

from test_iss_kratos import iss_object  # noqa: E402


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestIssDialog(unittest.TestCase):
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
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None
        ee.messagebox.askyesno = lambda *a, **k: True
        cls.dir = tempfile.mkdtemp()
        cls.kal = os.path.join(cls.dir, "iss.kal")
        with open(cls.kal, "w") as fh:
            fh.write("Dataset filename          = iss.dset\n"
                     + iss_object("ISS 1"))

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
        self.ws._add_files([self.kal])
        self.region = self.ws.docs[0].regions[0]
        self.ws.checked = {id(self.region)}

    def dialog(self):
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(self.root, self.ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        return dlg

    def test_it_offers_the_axis_ultra_geometry_not_the_123_degree_default(self):
        dlg = self.dialog()
        self.assertEqual((dlg.ion.get(), dlg.e0.get(), dlg.theta.get()),
                         ("He+", "1000", "135"))
        self.assertIn("Suggested", dlg.iss_note.cget("text"))
        self.assertIn("980 V", dlg.iss_note.cget("text"))   # as recorded
        self.assertEqual(self.ws.iss_saved(self.region), {})

    def test_nothing_is_stated_until_the_user_confirms(self):
        self.assertNotIn("scattering angle", self.ws.methods_generated())
        dlg = self.dialog()
        dlg._remember()
        self.assertEqual(self.ws.iss_saved(self.region),
                         {"ion": "He+", "e0": 1000.0, "theta": 135.0})
        self.assertIn("Saved with this file", dlg.iss_note.cget("text"))
        text = self.ws.methods_generated()
        self.assertIn("He+ ions, a beam energy of 1000 eV and a scattering "
                      "angle of 135°", text)
        self.assertNotIn("X-ray", text)

    def test_a_saved_value_comes_back_and_beats_the_offer(self):
        dlg = self.dialog()
        dlg.theta.set("136")
        dlg._remember()
        dlg.destroy()
        again = self.dialog()
        self.assertEqual(again.theta.get(), "136")

    def test_candidates_at_the_gold_peak(self):
        dlg = self.dialog()
        dlg._on_iss_click(933.0)
        syms = [c["symbol"] for c in dlg.cands]
        self.assertIn("Au", syms)
        self.assertEqual(syms[0], "Au")                    # nearest first


if __name__ == "__main__":
    unittest.main()
