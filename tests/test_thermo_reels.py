"""Avantage REELS spectra: the reader (technique 22 / the REELS lens, a kinetic
axis, the electron beam energy kept as the primary energy and not read as a
photon energy), the methods text, and the calibration to the primary energy that
is applied when such a file is loaded (and can be removed, switched off, or
refused when the elastic peak is far from the primary energy).

The real files are `D:\\Temp\\hxps\\BM-20263874-25-55`: set ``XPS_REELS_VGD`` to
that project folder to compare the six `.VGD` spectra with the XPSView session
of the same data.

Run:  python -m unittest discover tests
"""

import gzip
import json
import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import annotations  # noqa: E402
import methods  # noqa: E402
import reels  # noqa: E402
from readers import SpectrumFile, load_file  # noqa: E402
from test_iss_reels import ratio_spectrum  # noqa: E402
from test_thermo_kinds import HEADER, Tmp, dump  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False


def reels_dump(technique=22, lens="REELS", ev_scale=0, energy=1000.0,
               title="REELS Survey", el_ke=1004.3, flood=True, n=601):
    """A REELS `.avg` dump: a kinetic axis from 950 eV in 0.1 eV steps with the
    elastic peak at ``el_ke`` (the spectrum of ``ratio_spectrum``)."""
    ke, y, _loss = ratio_spectrum(300.0, 2000.0, el_ke=el_ke)
    pairs = sorted(zip(ke, y))
    pairs = [(k, v) for k, v in pairs if 950.0 <= k <= 1010.0][:n]
    values = [v for _k, v in pairs]
    axes = [(pairs[0][0], 0.05, len(values), "ENERGY", "LINEAR", "E", "eV",
             "Energy")]
    text = dump(title, axes, [(0, len(values) - 1, 1)], [((), [], values)])
    text = text.replace("VT_R4   = 1486.680054", f"VT_R4   = {energy:f}")
    extra = []
    if technique is not None:
        extra.append(f"DS_GEPROPID_TECHNIQUE                       : VT_I4   "
                     f"= {technique}")
    if lens is not None:
        extra.append(f"DS_ANPROPID_LENS_MODE_NAME                  : VT_BSTR "
                     f"= '{lens}'")
    if ev_scale is not None:
        extra.append(f"DS_ACPROPID_EV_SCALE                        : VT_I4   "
                     f"= {ev_scale}")
    if flood:
        extra += ["DS_SOURCE_FLOODGUNPROPID_DESCRIPTION           : VT_BSTR "
                  "= 'FG07 in REELS mode'",
                  "DS_SOURCE_FLOODGUNPROPID_CURRENT               : VT_R4   "
                  "= 5.000000",
                  "DS_SOURCE_FLOODGUNPROPID_ENERGY                : VT_R4   "
                  "= 1000.000000"]
    return text.replace("$PROPERTIES=STD\n",
                        "$PROPERTIES=STD\n" + "\n".join(extra) + "\n", 1)


class TestReader(Tmp):
    def load(self, **kw):
        return load_file(self.write("REELS Survey.avg", reels_dump(**kw)))

    def test_the_technique_code_makes_it_reels_on_a_kinetic_axis(self):
        f = self.load()
        r = f.regions[0]
        self.assertEqual(r.technique, "REELS")
        self.assertTrue(r.is_reels)
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertIsNone(r.photon_energy)                # not a photon energy
        self.assertEqual(r.anode, "")
        self.assertEqual(r.lens_mode, "REELS")
        self.assertAlmostEqual(r.energy[0], 1004.3 - 50.0, delta=0.06)
        self.assertEqual(r.energy, sorted(r.energy))      # ascending, native

    def test_no_warning_for_the_flag_a_reels_file_has(self):
        self.assertEqual(self.load().warnings, [])

    def test_the_primary_energy_is_kept_not_read_as_a_photon_energy(self):
        f = self.load(energy=1500.0)
        r = f.regions[0]
        self.assertEqual(r.extra["primary_energy"], 1500.0)
        self.assertEqual(reels.primary_energy(r), 1500.0)
        md = f.region_metadata(r)
        self.assertEqual(md["Primary beam energy (eV)"], "1500")
        self.assertEqual(md["Technique"], "REELS")
        self.assertEqual(md.get("Photon energy (eV)", ""), "")
        self.assertEqual(md.get("X-ray source", ""), "")
        self.assertNotIn("X-ray source", f.instrument)

    def test_the_name_is_the_title_not_a_shortened_region_name(self):
        self.assertEqual(self.load().regions[0].name, "REELS Survey")
        self.assertEqual(self.load(title="Loss 2").regions[0].name, "Loss 2")

    def test_the_lens_name_alone_is_enough(self):
        r = self.load(technique=None).regions[0]
        self.assertTrue(r.is_reels)
        self.assertEqual(r.extra["primary_energy"], 1000.0)

    def test_the_technique_code_alone_is_enough(self):
        r = self.load(lens=None).regions[0]
        self.assertTrue(r.is_reels)

    def test_any_other_unknown_flag_still_warns(self):
        f = self.load(ev_scale=7)
        self.assertTrue(any("Unknown energy scale flag (7)" in w
                            for w in f.warnings))
        self.assertTrue(f.regions[0].is_reels)

    def test_the_flood_gun_energy_is_written_as_it_is(self):
        f = self.load()
        self.assertEqual(f.instrument["Charge neutraliser"],
                         "FG07 in REELS mode, 5 µA, 1000 eV")

    def test_a_reels_file_without_an_energy_has_no_primary_energy(self):
        f = self.load(energy=0.0)
        self.assertIsNone(reels.primary_energy(f.regions[0]))

    def test_an_xps_file_is_unchanged(self):
        text = dump("C1s Scan", [(1188.0, 0.5, 4, "ENERGY", "LINEAR", "E",
                                  "eV", "Energy")], [(0, 3, 1)],
                    [((), [], [1, 2, 3, 4])])
        f = load_file(self.write("C1s Scan.avg", text))
        r = f.regions[0]
        self.assertEqual((r.technique, r.energy_label), ("XPS",
                                                         "Binding Energy"))
        self.assertAlmostEqual(r.photon_energy, 1486.68, 2)
        self.assertFalse(r.is_reels)
        self.assertNotIn("primary_energy", r.extra)

    def test_an_iss_file_is_still_iss(self):
        f = self.load(technique=6, lens="ISS", energy=0.0)
        r = f.regions[0]
        self.assertEqual(r.technique, "ISS")
        self.assertFalse(r.is_reels)

    def test_the_lens_name_of_another_technique_is_not_reels(self):
        r = self.load(technique=1, lens="Standard", ev_scale=1).regions[0]
        self.assertEqual(r.technique, "XPS")


class TestMethods(unittest.TestCase):
    ROW = {"Technique": "REELS", "Region": "REELS Survey", "Sample": "S",
           "Instrument": "Nexsa", "Lens mode": "REELS",
           "Primary beam energy (eV)": "1000", "Pass energy (eV)": "10",
           "Analyser mode": "Constant analyser energy (CAE)",
           "Step (eV)": "0.1", "Dwell (s)": "0.05", "KE start (eV)": "950",
           "KE end (eV)": "1010"}

    def text(self, rows):
        return methods.generate(rows)

    def test_a_reels_measurement_is_described_as_one(self):
        t = self.text([dict(self.ROW), dict(self.ROW, Sample="T")])
        self.assertIn("Reflection electron energy loss spectroscopy (REELS)",
                      t)
        self.assertIn("primary electron beam of 1000 eV", t)
        self.assertIn("Energy losses are measured from the elastic peak", t)
        self.assertIn("pass energy of 10 eV", t)

    def test_no_xray_wording_and_no_charge_correction_sentence(self):
        t = self.text([dict(self.ROW)])
        self.assertNotIn("X-ray", t)
        self.assertNotIn("excited", t)
        self.assertNotIn("charge-corrected", t)
        self.assertNotIn("Survey spectra", t)      # a REELS "survey" is not one

    def test_the_calibration_is_stated_with_its_shifts(self):
        rows = [dict(self.ROW, **{"REELS energy calibration": "x",
                                  "REELS energy shift (eV)": s})
                for s in ("-4.42", "-4.30", "-4.24", "-4.10")]
        t = self.text(rows)
        self.assertIn("The energy scale of each spectrum was calibrated by "
                      "shifting its kinetic energies (-4.42 to -4.1 eV) so "
                      "that the elastic peak lay at the primary energy", t)
        few = self.text(rows[:2])                   # two values: both listed
        self.assertIn("(-4.42 and -4.3 eV)", few)

    def test_a_partly_calibrated_set_says_how_many(self):
        rows = [dict(self.ROW, **{"REELS energy calibration": "x",
                                  "REELS energy shift (eV)": "-4.3"}),
                dict(self.ROW), dict(self.ROW)]
        self.assertIn("The energy scale of 1 of 3 spectra was calibrated",
                      self.text(rows))

    def test_nothing_is_written_that_the_file_did_not_state(self):
        row = {k: v for k, v in self.ROW.items()
               if k != "Primary beam energy (eV)"}
        t = self.text([row])
        self.assertNotIn("primary electron beam", t)
        self.assertNotIn("calibrated", t)

    def test_xps_and_reels_together_keep_each_description(self):
        xps = {"Technique": "XPS", "Region": "C 1s", "Sample": "A",
               "Photon energy (eV)": "1486.6", "Anode": "Al Kα (mono)",
               "Pass energy (eV)": "20"}
        t = self.text([xps, dict(self.ROW)])
        self.assertIn("X-ray photoelectron spectroscopy (XPS) and reflection "
                      "electron energy loss spectroscopy (REELS)", t)
        self.assertIn("excited with Al Kα", t)           # the XPS source
        self.assertIn("primary electron beam", t)
        self.assertIn("Binding energies are not charge-corrected", t)

    def test_the_avantage_reels_lens_alone_identifies_it(self):
        row = {k: v for k, v in self.ROW.items() if k != "Technique"}
        self.assertTrue(methods.is_reels(row))
        self.assertFalse(methods.is_survey(row))
        self.assertFalse(methods.is_reels({"Technique": "XPS",
                                           "Lens mode": "Standard"}))


class TestAutoCalibrationFunction(unittest.TestCase):
    def spectrum(self, el_ke):
        ke, y, _l = ratio_spectrum(300.0, 2000.0, el_ke=el_ke)
        return ke[::-1], y[::-1]

    def test_it_is_the_calibration_marked_automatic(self):
        ke, y = self.spectrum(1004.3)
        cal, why = reels.auto_calibration(ke, y, 1000.0)
        self.assertEqual(why, "")
        self.assertTrue(cal["auto"])
        self.assertEqual(cal["source"], "file")
        self.assertAlmostEqual(cal["shift"], -4.3, delta=0.05)
        plain = reels.calibration(ke, y, 1000.0)
        self.assertEqual(plain["shift"], cal["shift"])
        self.assertNotIn("auto", plain)

    def test_a_peak_far_from_the_primary_energy_is_refused_with_a_reason(self):
        ke, y = self.spectrum(1080.0)
        cal, why = reels.auto_calibration(ke, y, 1000.0)
        self.assertIsNone(cal)
        self.assertIn("more than 30 eV", why)
        self.assertIn("1080.0", why)
        far, _ = reels.auto_calibration(ke, y, 1000.0, limit=100.0)
        self.assertIsNotNone(far)                      # the limit is a number

    def test_just_inside_the_limit(self):
        ke, y = self.spectrum(1029.0)
        self.assertIsNotNone(reels.auto_calibration(ke, y, 1000.0)[0])

    def test_no_data_and_a_bad_primary_energy(self):
        self.assertEqual(reels.auto_calibration([], [], 1000.0)[0], None)
        self.assertIn("no elastic peak",
                      reels.auto_calibration([], [], 1000.0)[1])
        ke, y = self.spectrum(1000.0)
        self.assertIsNone(reels.auto_calibration(ke, y, 0)[0])

    def test_the_text_says_it_was_automatic(self):
        cal = {"elastic": 1004.3, "nominal": 1000.0, "shift": -4.3,
               "source": "file", "auto": True}
        self.assertIn("applied automatically", reels.calibration_text(cal))
        self.assertNotIn("automatically",
                         reels.calibration_text({**cal, "auto": False}))


class TestAnnotationsForAuto(unittest.TestCase):
    CAL = {"elastic": 1004.3, "nominal": 1000.0, "shift": -4.3,
           "source": "file", "auto": True}

    def test_removing_can_be_remembered_and_a_new_calibration_forgets_it(self):
        a = annotations.Annotations()
        a.set_reels_cal("f", "S", "R", self.CAL)
        self.assertTrue(a.reels_cal_for("f", "S", "R")["auto"])
        a.set_reels_cal("f", "S", "R", None)                  # plain removal
        self.assertFalse(a.reels_cal_is_declined("f", "S", "R"))
        a.set_reels_cal("f", "S", "R", self.CAL)
        a.set_reels_cal("f", "S", "R", None, decline=True)    # the user's
        self.assertTrue(a.reels_cal_is_declined("f", "S", "R"))
        self.assertIsNone(a.reels_cal_for("f", "S", "R"))
        a.set_reels_cal("f", "S", "R", self.CAL)              # by hand again
        self.assertFalse(a.reels_cal_is_declined("f", "S", "R"))

    def test_declined_is_saved_and_counts_as_content(self):
        a = annotations.Annotations()
        self.assertTrue(a.is_empty())
        a.set_reels_cal("f", "S", "R", None, decline=True)
        self.assertFalse(a.is_empty())
        b = annotations.Annotations.from_json(a.to_json())
        self.assertTrue(b.reels_cal_is_declined("f", "S", "R"))
        a.set_reels_cal("f", "S", "R", None, decline=True)    # once only
        self.assertEqual(len(a.reels_cal_declined), 1)

    def test_the_switch_is_saved_only_when_it_is_off(self):
        a = annotations.Annotations()
        self.assertTrue(a.reels_auto)
        self.assertNotIn("reels_auto", a.to_json())
        a.reels_auto = False
        self.assertFalse(a.is_empty())
        self.assertIs(a.to_json()["reels_auto"], False)
        self.assertFalse(annotations.Annotations.from_json(
            a.to_json()).reels_auto)
        self.assertTrue(annotations.Annotations.from_json(
            {"reels_auto": "no"}).reels_auto)                 # tolerant

    def test_loading_is_tolerant(self):
        b = annotations.Annotations.from_json({"reels_cal_declined": [
            "a|S|R", "a|S|R", 5, None, "", "b|S|R"]})
        self.assertEqual(b.reels_cal_declined, ["a|S|R", "b|S|R"])
        self.assertEqual(annotations.Annotations.from_json(
            {"reels_cal_declined": "x"}).reels_cal_declined, [])

    def test_the_auto_flag_of_a_calibration_survives_storage(self):
        a = annotations.Annotations()
        a.set_reels_cal("f", "S", "R", self.CAL)
        b = annotations.Annotations.from_json(a.to_json())
        self.assertTrue(b.reels_cal_for("f", "S", "R")["auto"])
        a.set_reels_cal("f", "S", "R", {**self.CAL, "auto": False})
        self.assertFalse(a.reels_cal_for("f", "S", "R")["auto"])


def reels_file(el_ke=1004.3, primary="1000", name="REELS Survey", sample="S"):
    from test_reels_calibration import reels_region
    f = SpectrumFile()
    f.path = "t.vms"
    f.format_name = "Test"
    f.regions = [reels_region(el_ke=el_ke, primary=primary, name=name,
                              sample=sample)]
    f._finish()
    return f


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestOnLoad(unittest.TestCase):
    def workspace(self):
        import matplotlib
        saved = matplotlib.rcParams.copy()
        try:
            root = tk.Tk()
        except tk.TclError:
            self.skipTest("no display")
        root.withdraw()
        self.addCleanup(lambda: (root.destroy(),
                                 matplotlib.rcParams.update(saved)))
        boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                 ee.messagebox.showerror)
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None
        self.addCleanup(lambda: setattr(ee.messagebox, "showinfo", boxes[0])
                        or setattr(ee.messagebox, "showwarning", boxes[1])
                        or setattr(ee.messagebox, "showerror", boxes[2]))
        return ee.Workspace(root)

    def load(self, ws, f=None, fid="f1", **kw):
        f = f or reels_file(**kw)
        problems = ws._register_parser(f, "t.vms", file_id=fid)
        return f, problems

    def elastic_of_display(self, ws, r):
        d = ws._display(r)
        order = sorted(zip(d.energy, d.counts))
        return reels.elastic_peak_refined([k for k, _y in order],
                                          [y for _k, y in order])

    def test_a_reels_file_opens_at_the_primary_energy(self):
        ws = self.workspace()
        f, problems = self.load(ws)
        r = f.regions[0]
        cal = ws.reels_cal_get(r)
        self.assertEqual(problems, [])
        self.assertTrue(cal["auto"])
        self.assertEqual(cal["source"], "file")
        self.assertAlmostEqual(cal["shift"], -4.3, delta=0.05)
        self.assertAlmostEqual(self.elastic_of_display(ws, r), 1000.0, 9)
        self.assertAlmostEqual(r.energy[0], f.regions[0].energy[0], 12)  # file
        key = ws._marker_key(r)[0]
        md = ws.ann.apply_metadata(key, 0, r, {"Sample": "S"})
        self.assertIn("applied automatically",
                      md["REELS energy calibration"])

    def test_the_loss_axis_is_drawn_from_zero_on_load(self):
        ws = self.workspace()
        f, _p = self.load(ws)
        ws.checked = {id(f.regions[0])}
        ws._render()
        ax = next(iter(ws._axinfo))
        self.assertEqual(ax.get_xlabel(), "Energy loss (eV)")
        self.assertAlmostEqual(ws._axloss[ax], 1000.0, 9)

    def test_a_file_that_does_not_record_the_primary_energy_is_left_alone(self):
        ws = self.workspace()
        f, problems = self.load(ws, primary=None)
        self.assertIsNone(ws.reels_cal_get(f.regions[0]))
        self.assertEqual(problems, [])
        self.assertAlmostEqual(self.elastic_of_display(ws, f.regions[0]),
                               1004.3, delta=0.05)

    def test_a_spectrum_that_is_not_reels_is_left_alone(self):
        from test_reels_calibration import reels_region
        ws = self.workspace()
        f = SpectrumFile()
        f.path, f.format_name = "t.vms", "Test"
        f.regions = [reels_region(technique="XPS")]
        f._finish()
        self.load(ws, f)
        self.assertIsNone(ws.reels_cal_get(f.regions[0]))

    def test_the_switch_off_leaves_a_new_file_uncalibrated(self):
        ws = self.workspace()
        ws.reels_auto_set(False)
        f, _p = self.load(ws)
        self.assertIsNone(ws.reels_cal_get(f.regions[0]))
        ws.reels_auto_set(True)
        g, _p = self.load(ws, fid="f2")
        self.assertIsNotNone(ws.reels_cal_get(g.regions[0]))

    def test_a_peak_far_from_the_primary_energy_is_not_shifted_and_says_so(self):
        ws = self.workspace()
        f, problems = self.load(ws, el_ke=1080.0)
        self.assertIsNone(ws.reels_cal_get(f.regions[0]))
        self.assertEqual(len(problems), 1)
        self.assertIn("not calibrated automatically", problems[0])
        self.assertIn("more than 30 eV", problems[0])
        self.assertIn("t.vms", problems[0])                    # which file
        self.assertAlmostEqual(self.elastic_of_display(ws, f.regions[0]),
                               1080.0, delta=0.05)

    def test_stored_results_are_in_the_calibrated_frame(self):
        ws = self.workspace()
        f, _p = self.load(ws)
        r = f.regions[0]
        from iss_ui import IssReelsDialog
        ws.checked = {id(r)}
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        dlg.elastic.set("")
        res = dlg._measure_ratio()
        self.assertAlmostEqual(res["elastic"], 1000.0, delta=0.05)
        # the same ratio as on the uncalibrated spectrum
        ke, y, loss = ratio_spectrum(300.0, 2000.0, el_ke=1004.3)
        el = reels.elastic_peak_refined(sorted(ke), [
            v for _k, v in sorted(zip(ke, y))])
        raw = reels.pi_ratio(reels.loss_axis(sorted(ke), el),
                             [v for _k, v in sorted(zip(ke, y))])
        self.assertAlmostEqual(res["ratio"], raw["ratio"], 9)

    def test_remove_sticks_across_a_reload_and_a_saved_workbook(self):
        ws = self.workspace()
        f, _p = self.load(ws)
        r = f.regions[0]
        from iss_ui import IssReelsDialog
        ws.checked = {id(r)}
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        dlg._uncalibrate()
        self.assertIsNone(ws.reels_cal_get(r))
        self.assertAlmostEqual(self.elastic_of_display(ws, r), 1004.3,
                               delta=0.05)
        # the file is loaded again: not calibrated again
        ws._auto_calibrate_reels(f)
        self.assertIsNone(ws.reels_cal_get(r))
        # the workbook is saved and reopened: its annotations come back first
        ws.ann = annotations.Annotations.from_json(ws.ann.to_json())
        ws._auto_calibrate_reels(f)
        self.assertIsNone(ws.reels_cal_get(r))
        # the button still works, and forgets the removal
        dlg._calibrate()
        self.assertIsNotNone(ws.reels_cal_get(r))
        self.assertFalse(ws.reels_cal_get(r).get("auto"))      # made by hand
        self.assertFalse(ws.ann.reels_cal_is_declined(*ws._marker_key(r)))

    def test_a_stored_calibration_is_not_recomputed_on_reopen(self):
        ws = self.workspace()
        f, _p = self.load(ws)
        r = f.regions[0]
        ws.reels_cal_set(r, {"elastic": 1004.3, "nominal": 1000.0,
                             "shift": -4.0, "source": "entered"})
        ws._auto_calibrate_reels(f)
        self.assertEqual(ws.reels_cal_get(r)["shift"], -4.0)    # the user's

    def test_the_tab_and_the_switch(self):
        ws = self.workspace()
        f, _p = self.load(ws)
        r = f.regions[0]
        from iss_ui import IssReelsDialog
        ws.checked = {id(r)}
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        self.assertIn("Calibrated automatically",
                      dlg.cal_status.cget("text"))
        self.assertTrue(dlg.cal_auto.get())
        dlg.cal_auto.set(False)
        dlg._auto_toggled()
        self.assertFalse(ws.ann.reels_auto)
        dlg._uncalibrate()
        self.assertIn("calibration removed", dlg.cal_status.cget("text"))

    def test_several_files_in_one_go(self):
        ws = self.workspace()
        a, _p = self.load(ws, el_ke=1004.30, sample="A", fid="f1")
        b, _p = self.load(ws, el_ke=1004.42, sample="B", fid="f2")
        shifts = [ws.reels_cal_get(x.regions[0])["shift"] for x in (a, b)]
        self.assertAlmostEqual(shifts[0], -4.30, delta=0.05)
        self.assertAlmostEqual(shifts[1], -4.42, delta=0.05)

    def test_the_methods_text_of_the_workspace_describes_the_calibration(self):
        ws = self.workspace()
        self.load(ws)
        text = ws.methods_generated()
        self.assertIn("calibrated by shifting its kinetic energies", text)
        self.assertNotIn("X-ray", text)


CORPUS = os.environ.get("XPS_REELS_VGD", "")


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_REELS_VGD to the project folder holding "
                     "Experiment\\REELS Source and XPSView")
class TestAgainstTheRealFiles(unittest.TestCase):
    """The six Avantage REELS spectra against XPSView's session of the same
    data: acquisition, elastic peak, calibration and the pi / (sigma+pi)
    ratios it stored."""

    @classmethod
    def setUpClass(cls):
        src = os.path.join(CORPUS, "Experiment", "REELS Source")
        cls.files = {}
        for sample in sorted(os.listdir(src)):
            p = os.path.join(src, sample, "REELS Survey.VGD")
            if os.path.exists(p):
                cls.files[sample] = load_file(p)
        xpsv = None
        for dirpath, _dirs, names in os.walk(os.path.join(CORPUS, "XPSView")):
            xpsv = next((os.path.join(dirpath, n) for n in names
                         if n.endswith(".xpsv")), xpsv)
        if xpsv is None:
            raise unittest.SkipTest("no .xpsv beside the files")
        raw = json.loads(gzip.open(xpsv, "rb").read().decode("utf-8"))
        an = next(a for a in raw["content"]["analyses"]
                  if a["kind"] == "reels")
        cls.analysis = an
        cls.rows = {r["sample"]: r for r in an["rows"]}
        cls.session = load_file(xpsv)

    def test_all_six_are_there(self):
        self.assertEqual(len(self.files), 6)
        self.assertEqual(set(self.files), set(self.rows))

    def test_each_is_reels_on_a_kinetic_axis_with_no_warning(self):
        for sample, f in self.files.items():
            r = f.regions[0]
            self.assertEqual(f.warnings, [], sample)
            self.assertEqual((r.technique, r.energy_label),
                             ("REELS", "Kinetic Energy"), sample)
            self.assertIsNone(r.photon_energy, sample)
            self.assertEqual(reels.primary_energy(r), 1000.0, sample)
            self.assertEqual(len(r.energy), 601)
            self.assertAlmostEqual(r.energy[0], 950.0, 6)
            self.assertAlmostEqual(r.energy[-1], 1010.0, 6)

    def test_the_spectra_equal_the_xpsview_session(self):
        for sample, f in self.files.items():
            ours = f.regions[0]
            theirs = next(x for x in self.session.regions
                          if x.sample == sample and x.is_reels)
            self.assertEqual(len(ours.energy), len(theirs.energy))
            for a, b in zip(ours.energy, theirs.energy):
                self.assertAlmostEqual(a, b, 4)
            peak = max(theirs.counts)
            for a, b in zip(ours.counts, theirs.counts):
                self.assertAlmostEqual(a, b, delta=1e-3 * peak, msg=sample)

    def test_the_elastic_peaks_are_the_ones_xpsview_found(self):
        for sample, f in self.files.items():
            r = f.regions[0]
            peak = reels.elastic_peak_refined(list(r.energy), list(r.counts))
            self.assertAlmostEqual(peak, self.rows[sample]["elastic_ke"], 3,
                                   msg=sample)

    @unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
    def test_loaded_in_the_app_they_open_at_1000_ev_with_xpsviews_ratios(self):
        import matplotlib
        saved = matplotlib.rcParams.copy()
        try:
            root = tk.Tk()
        except tk.TclError:
            self.skipTest("no display")
        root.withdraw()
        self.addCleanup(lambda: (root.destroy(),
                                 matplotlib.rcParams.update(saved)))
        ws = ee.Workspace(root)
        exp = os.path.join(CORPUS, "Experiment")
        doc = load_file(exp)
        self.assertEqual(ws._register_parser(doc, exp, file_id="f1"), [])
        self.assertEqual(len(doc.regions), 6)
        from iss_ui import IssReelsDialog
        ws.checked = {id(r) for r in doc.regions}
        dlg = IssReelsDialog(root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        for r in doc.regions:
            cal = ws.reels_cal_get(r)
            self.assertTrue(cal["auto"], r.sample)
            d = ws._display(r)
            order = sorted(zip(d.energy, d.counts))
            now = reels.elastic_peak_refined([k for k, _y in order],
                                             [y for _k, y in order])
            self.assertAlmostEqual(now, 1000.0, 9, r.sample)
            dlg.reg_cb.current(dlg.regions.index(r))
            dlg._region_changed()
            dlg.elastic.set("")
            res = dlg._measure_ratio()
            theirs = self.rows[r.sample]
            self.assertAlmostEqual(res["ratio"], theirs["ratio"], 6, r.sample)
            self.assertAlmostEqual(res["i_pi"], theirs["pi"]["integral"],
                                   delta=2e-5 * max(1, abs(
                                       theirs["pi"]["integral"])))
            self.assertAlmostEqual(res["elastic"], 1000.0, 9)
        ws.checked = {id(doc.regions[0])}
        ws._render()
        ax = next(iter(ws._axinfo))
        self.assertEqual(ax.get_xlabel(), "Energy loss (eV)")
        self.assertAlmostEqual(ws._axloss[ax], 1000.0, 9)
        text = ws.methods_generated()
        self.assertIn("primary electron beam of 1000 eV", text)
        self.assertNotIn("X-ray", text)


if __name__ == "__main__":
    unittest.main()
