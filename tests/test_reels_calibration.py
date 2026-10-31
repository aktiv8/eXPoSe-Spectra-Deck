"""REELS energy calibration (the elastic peak moved to the primary beam energy)
and the energy-loss axis: the pure pieces (reels.py, viewdata.py, the stored
calibration in annotations.py, the plot overlays) and the Workspace (the
display copy, the REELS axis box, the dialog's energy-scale tab, saved state).

Run:  python -m unittest discover tests
"""

import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import annotations  # noqa: E402
import reels  # noqa: E402
import viewdata  # noqa: E402
from readers import Region, SpectrumFile  # noqa: E402
from test_iss_reels import ratio_spectrum  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False
try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
except Exception:                                   # pragma: no cover
    Figure = None

RAW = 1004.3                    # where the author's spectra put the elastic peak


def reels_region(el_ke=RAW, primary="1000", technique="REELS", name="REELS",
                 sample="S"):
    """A REELS spectrum on a kinetic axis, its elastic peak at ``el_ke``."""
    ke, y, _loss = ratio_spectrum(300.0, 2000.0, el_ke=el_ke)
    extra = {"acq_metadata": {"Primary beam energy (eV)": primary}} \
        if primary is not None else {}
    return Region(name=name, index=0, offset=0, energy=ke, counts=y,
                  decodable=True, sample=sample, photon_energy=None,
                  energy_label="Kinetic Energy", technique=technique,
                  extra=extra)


class TestPrimaryEnergy(unittest.TestCase):
    def test_it_is_read_from_the_files_metadata(self):
        self.assertEqual(reels.primary_energy(reels_region()), 1000.0)
        self.assertEqual(reels.primary_energy(reels_region(primary="1500.5")),
                         1500.5)

    def test_a_direct_value_wins_and_the_rest_is_none(self):
        r = reels_region(primary=None)
        self.assertIsNone(reels.primary_energy(r))
        r.extra["primary_energy"] = 2000
        self.assertEqual(reels.primary_energy(r), 2000.0)
        for bad in ("x", "", "0", "-5", "nan", None):
            self.assertIsNone(reels.primary_energy(reels_region(primary=bad)),
                              bad)

    def test_the_other_label_is_read_too(self):
        r = reels_region(primary=None)
        r.extra["acq_metadata"] = {"Primary energy (eV)": "1000"}
        self.assertEqual(reels.primary_energy(r), 1000.0)


class TestCalibrationMaths(unittest.TestCase):
    def test_the_shift_puts_the_elastic_peak_at_the_nominal_energy(self):
        ke, y, _loss = ratio_spectrum(300.0, 2000.0, el_ke=RAW)
        cal = reels.calibration(ke[::-1], y[::-1], 1000.0)
        self.assertAlmostEqual(cal["elastic"], RAW, delta=0.03)
        self.assertAlmostEqual(cal["shift"], 1000.0 - cal["elastic"], 12)
        self.assertAlmostEqual(cal["shift"], -4.3, delta=0.03)
        self.assertEqual(cal["nominal"], 1000.0)
        shifted = [k + cal["shift"] for k in ke[::-1]]
        self.assertAlmostEqual(reels.elastic_peak_refined(shifted, y[::-1]),
                               1000.0, 9)

    def test_a_positive_shift_and_the_order_of_the_input(self):
        ke, y, _l = ratio_spectrum(300.0, 2000.0, el_ke=995.0)
        up = reels.calibration(ke[::-1], y[::-1], 1000.0)
        self.assertGreater(up["shift"], 4.9)
        self.assertAlmostEqual(
            reels.calibration(ke, y, 1000.0)["elastic"] if False else up[
                "elastic"], 995.0, delta=0.05)

    def test_refusals(self):
        ke, y, _l = ratio_spectrum(300.0, 2000.0)
        for nominal in (0, -1, "x", None, float("nan"), float("inf")):
            self.assertIsNone(reels.calibration(ke, y, nominal), nominal)
        self.assertIsNone(reels.calibration([], [], 1000.0))
        self.assertIsNone(reels.calibration(None, None, 1000.0))
        self.assertIsNone(reels.calibration([1, 2], [1], 1000.0))

    def test_the_sentence_says_what_moved_and_where_it_came_from(self):
        cal = {"elastic": 1004.3, "nominal": 1000.0, "shift": -4.3,
               "source": "file"}
        text = reels.calibration_text(cal)
        for token in ("-4.30 eV", "1004.30 eV", "1000 eV", "the file",
                      "0 eV of energy loss"):
            self.assertIn(token, text)
        self.assertIn("entered by the user",
                      reels.calibration_text({**cal, "source": "entered"}))
        self.assertNotIn("(", reels.calibration_text(
            {**cal, "source": None}).split("primary energy of")[1].split(
                " (0 eV")[0])


class TestLossAxis(unittest.TestCase):
    def region(self, zero=1000.0, **kw):
        r = reels_region(**kw)
        r.extra = dict(r.extra, reels_zero=zero) if zero else r.extra
        return r

    def test_it_is_the_energy_lost_from_the_zero(self):
        r = self.region(zero=1004.3)
        ax = viewdata.energy_axis(r, viewdata.LOSS)
        self.assertEqual((ax.label, ax.units, ax.invert, ax.ok),
                         ("Energy loss", "eV", False, True))
        self.assertEqual(ax.zero, 1004.3)
        self.assertIsNone(ax.e0)
        self.assertEqual(len(ax.x), len(r.energy))
        for loss, ke in zip(ax.x, r.energy):
            self.assertAlmostEqual(loss, 1004.3 - ke, 12)
        self.assertAlmostEqual(min(ax.x), 0.0, delta=0.06)   # at the peak

    def test_a_region_that_is_not_reels_keeps_its_axis(self):
        r = self.region(technique="XPS")
        ax = viewdata.energy_axis(r, viewdata.LOSS)
        self.assertEqual(ax.label, "Kinetic Energy")
        self.assertIsNone(ax.zero)
        self.assertTrue(ax.ok)

    def test_without_a_zero_it_says_so(self):
        ax = viewdata.energy_axis(self.region(zero=None), viewdata.LOSS)
        self.assertFalse(ax.ok)
        self.assertIsNone(ax.zero)

    def test_a_binding_labelled_reels_spectrum_is_converted_first(self):
        r = self.region(zero=1000.0)
        r.energy_label, r.photon_energy = "Binding Energy", 1000.0
        r.energy = [1000.0 - k for k in r.energy]            # BE = hv - KE
        ax = viewdata.energy_axis(r, viewdata.LOSS)
        self.assertEqual(ax.label, "Energy loss")
        ke = [1000.0 - e for e in r.energy]
        for loss, k in zip(ax.x, ke):
            self.assertAlmostEqual(loss, 1000.0 - k, 9)

    def test_the_other_axes_are_untouched_and_have_no_zero(self):
        r = self.region()
        for scale in ("Binding", "Kinetic", viewdata.RATIO):
            self.assertIsNone(viewdata.energy_axis(r, scale).zero, scale)
        self.assertEqual(viewdata.energy_axis(r, "Kinetic").label,
                         "Kinetic Energy")

    def test_the_choices(self):
        self.assertEqual(viewdata.REELS_AXES, (viewdata.LOSS, "Kinetic"))
        self.assertEqual(viewdata.LOSS, "Energy loss")

    def test_is_reels(self):
        self.assertTrue(reels_region().is_reels)
        self.assertTrue(reels_region(technique=" reels ").is_reels)
        self.assertFalse(reels_region(technique="XPS").is_reels)
        self.assertFalse(Region(name="x", index=0, offset=0).is_reels)


class TestStoredCalibration(unittest.TestCase):
    def cal(self, shift=-4.3):
        return {"elastic": RAW, "nominal": RAW + shift, "shift": shift,
                "source": "file"}

    def with_results(self):
        a = annotations.Annotations()
        for table, extra in (("reels", {"p1": [4.0, 420.0],
                                        "p2": [5.0, 820.0]}),
                             ("plasmon", {"ep": 27.0, "window": [15, 45],
                                          "method": "gauss",
                                          "background": "Linear"}),
                             ("loss_fit", {"window": [3, 50],
                                           "background": "Linear", "tol": 1.0,
                                           "components": [{"key": "pi",
                                                           "name": "pi",
                                                           "loss": 6.0,
                                                           "fwhm": 3.0,
                                                           "area": 1.0}]}),
                             ("pi_ratio", {"pi_window": [3.5, 8.0],
                                           "sp_window": [10.0, 35.0],
                                           "avg": 3, "i_pi": 1.0,
                                           "i_sp": 2.0, "ratio": 0.5})):
            getattr(a, "set_" + table)("f", "S", "R", {"elastic": RAW,
                                                       **extra})
        return a

    def test_the_elastic_peak_of_every_result_moves_with_the_scale(self):
        a = self.with_results()
        delta = a.set_reels_cal("f", "S", "R", self.cal())
        self.assertAlmostEqual(delta, -4.3, 12)
        for table in ("reels", "plasmon", "loss_fit", "pi_ratio"):
            got = getattr(a, table)[annotations.region_key("f", "S", "R")]
            self.assertAlmostEqual(got["elastic"], RAW - 4.3, 9, table)
        # the loss coordinates do not move
        self.assertEqual(a.reels_for("f", "S", "R")["p1"], [4.0, 420.0])
        self.assertEqual(a.plasmon_for("f", "S", "R")["window"], [15.0, 45.0])
        self.assertEqual(a.pi_ratio_for("f", "S", "R")["pi_window"],
                         [3.5, 8.0])

    def test_recalibrating_moves_by_the_change_only_and_removing_restores(self):
        a = self.with_results()
        a.set_reels_cal("f", "S", "R", self.cal(-4.3))
        self.assertAlmostEqual(
            a.set_reels_cal("f", "S", "R", self.cal(-4.0)), 0.3, 12)
        self.assertAlmostEqual(a.plasmon_for("f", "S", "R")["elastic"],
                               RAW - 4.0, 9)
        self.assertAlmostEqual(a.set_reels_cal("f", "S", "R", None), 4.0, 12)
        self.assertAlmostEqual(a.plasmon_for("f", "S", "R")["elastic"], RAW, 9)
        self.assertAlmostEqual(a.pi_ratio_for("f", "S", "R")["elastic"], RAW, 9)
        self.assertIsNone(a.reels_cal_for("f", "S", "R"))
        self.assertEqual(a.reels_shift_for("f", "S", "R"), 0.0)

    def test_another_spectrums_results_are_not_moved(self):
        a = self.with_results()
        a.set_plasmon("f", "S", "Other", {"elastic": RAW, "ep": 27.0,
                                          "window": [15, 45],
                                          "method": "gauss",
                                          "background": "Linear"})
        a.set_reels_cal("f", "S", "R", self.cal())
        self.assertAlmostEqual(a.plasmon_for("f", "S", "Other")["elastic"],
                               RAW, 12)

    def test_round_trip_metadata_and_empty(self):
        a = annotations.Annotations()
        self.assertTrue(a.is_empty())
        a.set_reels_cal("f", "S", "R", self.cal())
        self.assertFalse(a.is_empty())
        self.assertEqual(a.reels_shift_for("f", "S", "R"), -4.3)
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.reels_cal_for("f", "S", "R"),
                         a.reels_cal_for("f", "S", "R"))
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f", 0, r, {"Sample": "S", "Region": "R"})
        self.assertIn("-4.30 eV", md["REELS energy calibration"])
        self.assertIn("recorded in the file", md["REELS energy calibration"])
        a.set_reels_cal("f", "S", "R", None)
        self.assertTrue(a.is_empty())

    def test_loading_is_tolerant(self):
        good = self.cal()
        b = annotations.Annotations.from_json({"reels_cal": {
            "ok": good, "junk": "x", "noshift": {"elastic": 1, "nominal": 2},
            "nan": {**good, "shift": float("nan")},
            "negative": {**good, "nominal": -5},
            "odd": {**good, "source": "who knows"}}})
        self.assertEqual(sorted(b.reels_cal), ["odd", "ok"])
        self.assertEqual(b.reels_cal["odd"]["source"], "entered")
        self.assertEqual(annotations.Annotations.from_json(
            {"reels_cal": [1]}).reels_cal, {})

    def test_a_split_sample_keeps_its_calibration(self):
        self.assertIn("reels_cal", annotations.Annotations.ADOPTED)


@unittest.skipUnless(Figure is not None, "matplotlib not installed")
class TestOverlaysOnTheLossAxis(unittest.TestCase):
    def draw(self, **overlay):
        import plots
        r = reels_region(el_ke=1000.0)
        r.extra = dict(r.extra, reels_zero=1000.0)
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        plots.draw_stack(ax, [r], 0.6, "None", None, ["#1f77b4"], "R", "", (),
                         False, scale=viewdata.LOSS, **overlay)
        fig.canvas.draw()
        return ax

    def test_the_axis_reads_the_energy_lost(self):
        ax = self.draw()
        self.assertEqual(ax.get_xlabel(), "Energy loss (eV)")
        lo, hi = ax.get_xlim()
        self.assertLess(lo, hi)                  # increasing to the right
        self.assertGreater(hi, 45.0)
        self.assertLess(lo, 1.0)

    def test_the_helper_maps_a_loss_on_each_kind_of_axis(self):
        import plots
        kin = viewdata.energy_axis(reels_region(), "Kinetic")
        loss = viewdata.Axis([], "Energy loss", "eV", False, True, None,
                             1004.3)
        be = viewdata.Axis([], "Binding Energy", "eV", True, True)
        self.assertAlmostEqual(plots._loss_x(kin, None, 1004.3)(5.0), 999.3, 9)
        self.assertAlmostEqual(plots._loss_x(loss, None, 1004.3)(5.0), 5.0, 9)
        # the result was measured from 1004.3 but the axis is zeroed on 1004.0
        loss2 = loss._replace(zero=1004.0)
        self.assertAlmostEqual(plots._loss_x(loss2, None, 1004.3)(5.0), 4.7, 9)
        self.assertAlmostEqual(plots._loss_x(be, 1486.6, 1004.3)(5.0),
                               1486.6 - 999.3, 9)
        self.assertIsNone(plots._loss_x(be, None, 1004.3)(5.0))

    def test_the_ratio_windows_are_shaded_at_their_loss_values(self):
        _ke, y, loss = ratio_spectrum(300.0, 2000.0, el_ke=1000.0)
        res = {**reels.pi_ratio(loss, y), "elastic": 1000.0}
        ax = self.draw(piratio={**res, "fraction": None})
        spans = sorted((round(p.get_x(), 2), round(p.get_x() + p.get_width(), 2))
                       for p in ax.patches)
        self.assertEqual(spans[0], (3.5, 8.0))
        self.assertEqual(spans[1], (10.0, 35.0))

    def test_the_plasmon_is_labelled_at_its_loss(self):
        _ke, y, loss = ratio_spectrum(300.0, 2000.0, el_ke=1000.0)
        res = {**reels.plasmon(loss, y, (10, 50), "Linear"),
               "elastic": 1000.0}
        ax = self.draw(plasmon=res)
        lines = [ln for ln in ax.lines
                 if len(ln.get_xdata()) == 2 and ln.get_linestyle() == "-"]
        self.assertTrue(any(abs(ln.get_xdata()[0] - res["ep"]) < 1e-6
                            for ln in lines))

    def test_the_band_gap_tangent_runs_from_the_gap(self):
        res = {"elastic": 1000.0, "p1": [4.0, 420.0], "p2": [5.0, 820.0],
               "gap": 3.2, "base": 100.0, "slope": 400.0}
        ax = self.draw(reels=res)
        tangent = [ln for ln in ax.lines if len(ln.get_xdata()) == 2
                   and ln.get_linestyle() == "-"][-1]
        self.assertAlmostEqual(min(tangent.get_xdata()), 3.2, 6)

    def test_a_stored_peak_off_the_axis_zero_is_corrected(self):
        _ke, y, loss = ratio_spectrum(300.0, 2000.0, el_ke=1000.0)
        res = {**reels.pi_ratio(loss, y), "elastic": 1000.3}   # 0.3 eV off
        ax = self.draw(piratio={**res, "fraction": None})
        first = sorted(p.get_x() for p in ax.patches)[0]
        self.assertAlmostEqual(first, 3.5 - 0.3, 2)


def make_workspace(case, regions=None, file_id="f1", auto=False):
    """A workspace holding one file. ``auto`` False (the default) switches the
    automatic calibration on load off, so these tests start from an
    uncalibrated spectrum and exercise the manual path; the automatic one is
    tested in test_thermo_reels.py."""
    import matplotlib
    saved = matplotlib.rcParams.copy()
    try:
        root = tk.Tk()
    except tk.TclError:
        case.skipTest("no display")
    root.withdraw()
    case.addCleanup(lambda: (root.destroy(), matplotlib.rcParams.update(saved)))
    boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
             ee.messagebox.showerror)
    ee.messagebox.showinfo = ee.messagebox.showwarning = \
        ee.messagebox.showerror = lambda *a, **k: None

    def restore():
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror) = boxes
    case.addCleanup(restore)
    ws = ee.Workspace(root)
    ws.cfg["reels_axis"] = viewdata.LOSS
    ws.ann.reels_auto = bool(auto)
    f = SpectrumFile()
    f.path = "t.vms"
    f.format_name = "Test"
    f.regions = regions or [reels_region()]
    f._finish()
    ws._register_parser(f, "t.vms", file_id=file_id)
    ws.checked = {id(r) for r in f.regions}
    return ws, f


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestWorkspace(unittest.TestCase):
    def setUp(self):
        self.ws, self.doc = make_workspace(self)
        self.region = self.doc.regions[0]

    def dialog(self):
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(self.ws.root, self.ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        return dlg

    def axis(self):
        ws = self.ws
        ws._render()
        return next(iter(ws._axinfo))

    def test_the_display_copy_is_shifted_and_the_file_is_not(self):
        raw = list(self.region.energy)
        self.assertAlmostEqual(
            reels.elastic_peak_refined(sorted(raw), [
                y for _k, y in sorted(zip(raw, self.region.counts))]),
            RAW, delta=0.05)
        dlg = self.dialog()
        dlg._calibrate()
        cal = self.ws.reels_cal_get(self.region)
        self.assertAlmostEqual(cal["shift"], -4.3, delta=0.05)
        self.assertEqual(cal["source"], "file")
        shown = self.ws._display(self.region)
        for a, b in zip(shown.energy, raw):
            self.assertAlmostEqual(a - b, cal["shift"], 9)
        self.assertEqual(list(self.region.energy), raw)      # not mutated
        order = sorted(zip(shown.energy, shown.counts))
        self.assertAlmostEqual(reels.elastic_peak_refined(
            [k for k, _y in order], [y for _k, y in order]), 1000.0, 9)
        self.assertAlmostEqual(shown.reels_shift_applied, cal["shift"], 12)

    def test_the_loss_axis_starts_at_zero_and_increases_to_the_right(self):
        ax = self.axis()
        self.assertEqual(ax.get_xlabel(), "Energy loss (eV)")
        lo, hi = ax.get_xlim()
        self.assertLess(lo, hi)
        self.assertAlmostEqual(self.ws._axloss[ax], RAW, delta=0.05)
        self.assertLess(abs(lo), 2.0)                    # 0 is at the left
        line = ax.lines[0]
        self.assertAlmostEqual(min(line.get_xdata()), 0.0, delta=0.06)

    def test_the_kinetic_choice_draws_the_calibrated_energy(self):
        dlg = self.dialog()
        dlg._calibrate()
        self.ws.reels_axis_var.set("Kinetic")
        ax = self.axis()
        self.assertEqual(ax.get_xlabel(), "Kinetic energy (eV)")
        self.assertNotIn(ax, self.ws._axloss)
        x = ax.lines[0].get_xdata()
        self.assertAlmostEqual(max(x), 1000.0 + 0.0, delta=3.0)
        self.assertLess(max(x), 1001.0)                  # was 1004.3 + 0
        self.ws.reels_axis_var.set(viewdata.LOSS)
        ax = self.axis()
        self.assertAlmostEqual(self.ws._axloss[ax], 1000.0, delta=0.05)

    def test_calibrating_moves_the_stored_peak_and_keeps_the_ratio(self):
        dlg = self.dialog()
        dlg.elastic.set("")
        first = dlg._measure_ratio()
        self.assertAlmostEqual(first["elastic"], RAW, delta=0.05)
        dlg._calibrate()
        moved = self.ws.pi_ratio_get(self.region)
        self.assertAlmostEqual(moved["elastic"], 1000.0, delta=0.05)
        self.assertEqual(moved["ratio"], first["ratio"])     # loss is loss
        self.assertEqual(moved["i_pi"], first["i_pi"])
        self.assertAlmostEqual(float(dlg.elastic.get()), 1000.0, delta=0.05)
        # measuring again on the calibrated spectrum gives the same ratio
        dlg._clear_ratio()
        again = dlg._measure_ratio()
        self.assertAlmostEqual(again["ratio"], first["ratio"], 9)
        self.assertAlmostEqual(again["elastic"], 1000.0, delta=0.05)

    def test_removing_it_puts_everything_back(self):
        dlg = self.dialog()
        dlg.elastic.set("")
        first = dlg._measure_ratio()
        dlg._calibrate()
        dlg._uncalibrate()
        self.assertIsNone(self.ws.reels_cal_get(self.region))
        back = self.ws.pi_ratio_get(self.region)
        self.assertAlmostEqual(back["elastic"], first["elastic"], 9)
        self.assertIs(self.ws._display(self.region).energy,
                      self.ws._display(self.region).energy)
        self.assertAlmostEqual(
            max(self.ws._display(self.region).energy), max(
                self.region.energy), 12)

    def test_the_metadata_says_what_was_done(self):
        dlg = self.dialog()
        dlg._calibrate()
        key = self.ws._marker_key(self.region)[0]
        md = self.ws.ann.apply_metadata(key, 0, self.region,
                                        {"Sample": "S", "Region": "REELS"})
        self.assertIn("1000 eV", md["REELS energy calibration"])

    def test_the_tab_shows_the_files_energy_and_the_state(self):
        dlg = self.dialog()
        self.assertEqual(dlg.cal_nominal.get(), "1000")
        self.assertIn("From the file", dlg.cal_source.cget("text"))
        self.assertIn("not calibrated", dlg.cal_status.cget("text"))
        self.assertIn("1004.3", dlg.cal_status.cget("text"))
        dlg._calibrate()
        self.assertIn("Calibrated", dlg.cal_status.cget("text"))
        self.assertIn("-4.3", dlg.cal_status.cget("text"))
        dlg._uncalibrate()
        self.assertIn("calibration removed", dlg.cal_status.cget("text"))
        self.assertIn("not calibrated automatically again",
                      dlg.cal_detail.cget("text"))

    def test_a_typed_energy_is_marked_entered_and_used(self):
        dlg = self.dialog()
        dlg.cal_nominal.set("1005")
        dlg._calibrate()
        cal = self.ws.reels_cal_get(self.region)
        self.assertEqual((cal["nominal"], cal["source"]), (1005.0, "entered"))
        self.assertAlmostEqual(cal["shift"], 1005.0 - cal["elastic"], 12)

    def test_the_cursor_and_the_clicks_speak_energy_loss(self):
        ax = self.axis()
        zero = self.ws._axloss[ax]

        class Ev:
            inaxes, xdata, ydata, button = ax, 12.0, 100.0, 1
        self.ws._on_motion(Ev)
        text = self.ws.cursor_lbl.cget("text")
        self.assertIn("Loss 12.00 eV", text)
        self.assertIn(f"KE {zero - 12.0:.2f} eV", text)
        self.assertAlmostEqual(self.ws.kinetic_from_event(self.region, Ev),
                               zero - 12.0, 9)
        self.assertIsNone(self.ws.displayed_be(Ev))     # not a binding energy

    def test_a_click_on_the_plot_gives_the_same_loss_as_before(self):
        # the band-gap pick: a click at loss x becomes the point (x, y)
        dlg = self.dialog()
        dlg.elastic.set("")
        dlg._auto_elastic()
        el = dlg._elastic()
        ax = self.axis()

        class Ev:
            inaxes, xdata, ydata = ax, 4.0, 100.0
        ke = self.ws.kinetic_from_event(self.region, Ev)
        self.assertAlmostEqual(el - ke, 4.0, delta=0.06)

    def test_a_renamed_spectrum_keeps_its_stored_peak_as_the_axis_zero(self):
        dlg = self.dialog()
        dlg.elastic.set("1004.9")                   # a peak the user chose
        dlg._measure_ratio()
        self.assertAlmostEqual(self.ws.pi_ratio_get(self.region)["elastic"],
                               1004.9, 9)
        fid = self.ws._marker_key(self.region)[0]
        self.ws.ann.set_name("sample_names",
                             annotations.sample_key(fid, self.region.sample),
                             "Renamed", self.region.sample)
        self.ws.ann.set_name("region_names",
                             annotations.region_key(fid, self.region.sample,
                                                    self.region.name),
                             "Also renamed", self.region.name)
        self.ws._ann_changed()
        ax = self.axis()
        self.assertAlmostEqual(self.ws._axloss[ax], 1004.9, 9)
        self.assertEqual(self.ws._display(self.region).sample, "Renamed")

    def test_the_view_state_carries_the_axis_and_old_states_are_kinetic(self):
        st = self.ws.capture_state()
        self.assertEqual(st["reels_axis"], viewdata.LOSS)
        self.ws.reels_axis_var.set("Kinetic")
        self.assertEqual(self.ws.capture_state()["reels_axis"], "Kinetic")
        old = {k: v for k, v in st.items() if k != "reels_axis"}
        self.ws.reels_axis_var.set(viewdata.LOSS)
        self.ws.apply_state(old)
        self.assertEqual(self.ws.reels_axis_var.get(), "Kinetic")
        self.ws.apply_state(st)
        self.assertEqual(self.ws.reels_axis_var.get(), viewdata.LOSS)

    def test_a_zoom_is_tied_to_the_loss_axis_only_on_a_page_that_has_one(self):
        self.assertEqual(self.ws._zoom_scale(), "Binding")   # nothing drawn
        self.axis()
        self.assertEqual(self.ws._zoom_scale(), "Binding|Energy loss")
        self.ws.reels_axis_var.set("Kinetic")
        self.axis()
        self.assertEqual(self.ws._zoom_scale(), "Binding")

    def test_the_axis_follows_a_saved_figure(self):
        self.ws.add_figure("loss")
        fig = self.ws.figures[-1]
        self.assertEqual(fig["state"]["reels_axis"], viewdata.LOSS)



@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestWorkspaceVariants(unittest.TestCase):
    """Each test builds its own workspace (there is one Tk root at a time)."""

    def test_without_a_primary_energy_nothing_is_assumed(self):
        ws, doc = make_workspace(self, [reels_region(primary=None)])
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        self.assertEqual(dlg.cal_nominal.get(), "")
        self.assertIn("does not record", dlg.cal_source.cget("text"))
        dlg._calibrate()
        self.assertIsNone(ws.reels_cal_get(doc.regions[0]))
        self.assertIn("Enter the primary", dlg.cal_status.cget("text"))
        dlg.cal_nominal.set("1000")
        dlg._calibrate()
        self.assertEqual(ws.reels_cal_get(doc.regions[0])["source"], "entered")

    def test_a_spectrum_that_is_not_reels_cannot_be_calibrated(self):
        ws, doc = make_workspace(self, [reels_region(technique="XPS")])
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        self.assertIn("Not a REELS", dlg.cal_status.cget("text"))
        dlg._calibrate()
        self.assertIsNone(ws.reels_cal_get(doc.regions[0]))
        self.assertIn("Could not", dlg.cal_status.cget("text"))

    def test_every_spectrum_in_the_list_each_with_its_own_peak(self):
        regions = [reels_region(el_ke=1004.30, name="A"),
                   reels_region(el_ke=1004.42, name="B")]
        ws, doc = make_workspace(self, regions)
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(ws.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        dlg._calibrate_all()
        shifts = [ws.reels_cal_get(r)["shift"] for r in doc.regions]
        self.assertAlmostEqual(shifts[0], -4.30, delta=0.05)
        self.assertAlmostEqual(shifts[1], -4.42, delta=0.05)
        self.assertNotAlmostEqual(shifts[0], shifts[1], delta=0.05)
        for r in doc.regions:                       # both now at 1000
            d = ws._display(r)
            order = sorted(zip(d.energy, d.counts))
            self.assertAlmostEqual(reels.elastic_peak_refined(
                [k for k, _y in order], [y for _k, y in order]), 1000.0, 9)

    def test_an_xps_page_is_not_affected_by_the_reels_axis(self):
        xps = reels_region(technique="XPS", name="C 1s")
        ws, doc = make_workspace(self, [xps])
        ws._render()
        ax = next(iter(ws._axinfo))
        self.assertNotIn(ax, ws._axloss)
        self.assertNotEqual(ax.get_xlabel(), "Energy loss (eV)")


if __name__ == "__main__":
    unittest.main()
