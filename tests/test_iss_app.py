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
from test_iss_reels import plasmon_spectrum  # noqa: E402


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

    def test_measuring_a_plasmon_stores_it_and_clearing_removes_it(self):
        dlg = self.dialog()
        ke, y, _loss = plasmon_spectrum(level=100.0, ramp=5.0)
        dlg._ke_counts = lambda: (list(ke), list(y))
        dlg.elastic.set("1000")
        dlg.pl_feature.set("Custom")
        dlg.pl_lo.set("10")
        dlg.pl_hi.set("50")
        dlg.pl_bg.set("Shirley")
        dlg._measure_plasmon()
        got = self.ws.plasmon_get(self.region)
        self.assertAlmostEqual(got["ep"], 27.0, delta=0.5)
        self.assertEqual(got["background"], "Shirley")
        self.assertIn("Ep = 26.", dlg.pl_result.cget("text"))
        self.assertIsNone(self.ws.reels_get(self.region))   # band gap untouched
        # both methods: the Gaussian is stored, the second derivative shown
        dlg.pl_method.set(dlg.PL_METHODS[2])
        dlg._measure_plasmon()
        self.assertEqual(self.ws.plasmon_get(self.region)["method"], "gauss")
        self.assertIn("Second derivative", dlg.pl_detail.cget("text"))
        dlg._clear_plasmon()
        self.assertIsNone(self.ws.plasmon_get(self.region))

    def test_the_loss_fit_and_the_layer_estimate(self):
        from test_iss_reels import loss_spectrum, LAM_6_5
        dlg = self.dialog()
        ke, y, _loss = loss_spectrum()
        dlg._ke_counts = lambda: (list(ke), list(y))
        dlg.elastic.set("2000")
        dlg._fit_loss()
        got = self.ws.loss_fit_get(self.region)
        self.assertAlmostEqual(got["r"], 3.0, delta=0.5)
        self.assertIn("R = ", dlg.lf_result.cget("text"))
        self.assertIn("C sp2 bulk", dlg.lf_table.cget("text"))
        self.assertNotIn("layers", self.ws.methods_generated())
        # no geometry yet: nothing is assumed
        dlg._estimate_layers()
        self.assertIn("Enter", dlg.ly_result.cget("text"))
        self.assertNotIn("geometry", self.ws.loss_fit_get(self.region))
        for key, v in (("lam", f"{LAM_6_5:.4f}"), ("a_in", "0"),
                       ("a_out", "35")):
            dlg.ly[key].set(v)
        dlg._estimate_layers()
        got = self.ws.loss_fit_get(self.region)
        self.assertGreater(got["layers"]["n"], 1.0)
        self.assertIn("layers", dlg.ly_result.cget("text"))
        # a refit keeps the geometry and recomputes the estimate
        dlg._fit_loss()
        self.assertIn("geometry", self.ws.loss_fit_get(self.region))
        self.assertIn("layers", self.ws.loss_fit_get(self.region))
        dlg._clear_loss()
        self.assertIsNone(self.ws.loss_fit_get(self.region))

    def patch_spectrum(self, dlg, pi_area, sp_area=2000.0):
        from test_iss_reels import ratio_spectrum
        ke, y, _loss = ratio_spectrum(pi_area, sp_area)
        dlg._ke_counts = lambda: (list(ke), list(y))
        dlg._ke_counts_for = lambda r: (list(ke), list(y))
        dlg.elastic.set("")                   # use the spectrum's own peak

    def test_measuring_the_ratio_stores_it_with_the_windows_used(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)
        res = dlg._measure_ratio()
        got = self.ws.pi_ratio_get(self.region)
        self.assertAlmostEqual(got["ratio"], 0.16, delta=0.01)
        self.assertEqual((got["pi_window"], got["sp_window"], got["avg"]),
                         ([3.5, 8.0], [10.0, 35.0], 3))
        self.assertAlmostEqual(got["elastic"], 1000.0, delta=0.05)
        self.assertEqual(res["ratio"], got["ratio"])
        self.assertIn("ratio", dlg.rt_table.cget("text"))
        self.assertIn("0.1", dlg.rt_table.cget("text"))
        self.assertIn("No reference yet", dlg.rt_refs.cget("text"))
        md = self.ws.ann.apply_metadata(self.ws._marker_key(self.region)[0], 0,
                                       self.region, {"Sample": "S"})
        self.assertIn("REELS π / (σ+π) ratio", md)
        self.assertNotIn("REELS sp2 fraction", md)         # no reference yet
        dlg._clear_ratio()
        self.assertIsNone(self.ws.pi_ratio_get(self.region))

    def test_hopg_and_diamond_set_the_scale(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)                    # the HOPG
        dlg._set_reference("hopg")
        hopg = self.ws.ann.reels_refs["hopg"]
        self.assertEqual(hopg["kind"], "spectrum")
        self.assertEqual(hopg["pi_window"], [3.5, 8.0])
        self.assertAlmostEqual(hopg["e0"], 1000.0, delta=0.05)
        self.assertIn("1.00", dlg.rt_table.cget("text"))   # itself is 1
        self.patch_spectrum(dlg, 0.0)                      # the diamond
        dlg._set_reference("diamond")
        self.assertIn("diamond", dlg.rt_refs.cget("text"))
        self.patch_spectrum(dlg, 150.0)                    # a sample
        dlg._measure_ratio()
        sp = __import__("reels").sp2_for(self.ws.pi_ratio_get(self.region),
                                         self.ws.ann.reels_refs)
        self.assertAlmostEqual(sp["fraction"], 0.5, delta=0.05)
        self.assertIn("0.5", dlg.rt_table.cget("text"))
        md = self.ws.ann.apply_metadata(self.ws._marker_key(self.region)[0], 0,
                                       self.region, {"Sample": "S"})
        self.assertIn("HOPG", md["REELS sp2 fraction"])
        self.assertIn("diamond", md["REELS sp2 fraction"])
        dlg._clear_references()
        self.assertEqual(self.ws.ann.reels_refs, {})

    def test_a_typed_hopg_ratio_is_marked_as_typed(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)
        dlg._measure_ratio()
        dlg.rt_value.set("0.32")
        dlg.rt_source.set("a note from the lab")
        dlg._enter_reference()
        ref = self.ws.ann.reels_refs["hopg"]
        self.assertEqual((ref["kind"], ref["ratio"]), ("entered", 0.32))
        self.assertIn("typed", dlg.rt_refs.cget("text"))
        self.assertIn("a note from the lab", dlg.rt_refs.cget("text"))
        dlg.rt_value.set("x")
        dlg._enter_reference()                             # not a number
        self.assertEqual(self.ws.ann.reels_refs["hopg"]["ratio"], 0.32)

    def test_a_reference_measured_with_other_windows_is_not_used(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)
        dlg._set_reference("hopg")
        dlg.rt_pi_hi.set("9")                              # another window
        dlg._measure_ratio()
        self.assertIn("no sp² fraction", dlg.rt_detail.cget("text"))
        md = self.ws.ann.apply_metadata(self.ws._marker_key(self.region)[0], 0,
                                       self.region, {"Sample": "S"})
        self.assertIn("not calculated", md["REELS sp2 fraction"])

    def test_measuring_every_spectrum_and_the_bad_inputs(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)
        dlg._measure_ratio_all()
        self.assertIsNotNone(self.ws.pi_ratio_get(self.region))
        dlg._clear_ratio()
        dlg._ke_counts_for = lambda r: (None, None)        # no energies
        dlg._measure_ratio_all()
        self.assertIsNone(self.ws.pi_ratio_get(self.region))
        self.assertIn("left out", dlg.rt_detail.cget("text"))
        dlg.rt_sp_lo.set("40")                             # start after end
        dlg._measure_ratio_all()
        self.assertIn("two numbers", dlg.rt_detail.cget("text"))
        self.assertIsNone(self.ws.pi_ratio_get(self.region))

    def test_the_window_sensitivity_table(self):
        dlg = self.dialog()
        self.patch_spectrum(dlg, 300.0)
        dlg._ratio_sensitivity()
        text = dlg.rt_sens.cget("text")
        self.assertIn("σ+π start", text)
        self.assertEqual(len(text.splitlines()), 9)

    def test_tpp2m_fills_the_mean_free_path(self):
        dlg = self.dialog()
        dlg.elastic.set("2000")
        dlg._lam_from_tpp()
        self.assertAlmostEqual(float(dlg.ly["lam"].get()), 4.06, delta=0.01)

    def test_a_window_that_is_not_two_numbers_stores_nothing(self):
        dlg = self.dialog()
        dlg.pl_lo.set("x")
        dlg._measure_plasmon()
        self.assertIsNone(self.ws.plasmon_get(self.region))


if __name__ == "__main__":
    unittest.main()
