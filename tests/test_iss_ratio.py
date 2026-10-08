"""The energy-ratio axis (E/E0) for ion scattering spectra: the pure axis, the
drawing of markers and cursors on it, and the real window (needs a display
and matplotlib for the last part; skipped otherwise).

The beam energy comes from what the user saved with the file or what the file
records. The Kratos suggestion (1000 eV) is never used for the axis, so a
Kratos spectrum stays on kinetic energy, with a note, until it is confirmed.

Run:  python -m unittest discover tests
"""

import os
import sys
import tempfile
import types
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import viewdata  # noqa: E402
from readers import Region  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

from test_iss_kratos import iss_object  # noqa: E402
from test_iss_thermo import iss_avg  # noqa: E402
from test_readers import AVG_TEMPLATE  # noqa: E402


def iss(e0=None, **kw):
    r = Region("ISS", 0, 0, technique="ISS", energy=[100.0, 500.0, 933.0],
               counts=[1, 2, 3], energy_label="Kinetic Energy",
               decodable=True, **kw)
    if e0:
        r.extra["iss_e0"] = e0
    return r


class TestRatioAxis(unittest.TestCase):
    def test_it_is_an_iss_axis_not_an_xps_scale(self):
        self.assertEqual(viewdata.ENERGY_SCALES, ("Binding", "Kinetic"))
        self.assertEqual(viewdata.ISS_AXES, ("Kinetic", "Energy ratio"))

    def test_it_divides_by_the_beam_energy(self):
        a = viewdata.energy_axis(iss(1000.0), viewdata.RATIO)
        self.assertEqual(a.x, [0.1, 0.5, 0.933])
        self.assertEqual((a.label, a.units, a.invert, a.ok, a.e0),
                         ("Energy ratio", "E/E$_0$", False, True, 1000.0))

    def test_the_calibrated_energy_is_what_is_divided_by(self):
        a = viewdata.energy_axis(iss(966.313), viewdata.RATIO)
        self.assertAlmostEqual(a.x[-1], 933.0 / 966.313)

    def test_without_a_beam_energy_the_kinetic_axis_is_kept(self):
        a = viewdata.energy_axis(iss(), viewdata.RATIO)
        self.assertEqual((a.label, a.ok, a.e0), ("Kinetic Energy", False, None))
        self.assertEqual(a.x, [100.0, 500.0, 933.0])

    def test_other_spectra_keep_their_own_axis(self):
        xps = Region("C 1s", 0, 0, energy=[290.0, 285.0, 280.0],
                     counts=[1, 2, 3], photon_energy=1486.6, decodable=True)
        a = viewdata.energy_axis(xps, viewdata.RATIO)
        self.assertEqual((a.label, a.ok, a.invert),
                         ("Binding Energy", True, True))
        ke = Region("T", 0, 0, energy=[100.0, 200.0], counts=[1, 2],
                    energy_label="Kinetic Energy", decodable=True)
        self.assertEqual(viewdata.energy_axis(ke, viewdata.RATIO).label,
                         "Kinetic Energy")
        self.assertTrue(viewdata.energy_axis(ke, viewdata.RATIO).ok)

    def test_the_other_scales_are_unchanged(self):
        r = iss(1000.0)
        for scale in ("Binding", "Kinetic"):
            a = viewdata.energy_axis(r, scale)
            self.assertEqual((a.label, a.e0), ("Kinetic Energy", None))


class TestDrawing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib.figure import Figure
        except Exception:                           # pragma: no cover
            raise unittest.SkipTest("no matplotlib")
        cls._rc = matplotlib.rcParams.copy()
        cls.Figure = Figure

    @classmethod
    def tearDownClass(cls):
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def panel(self, **kw):
        import plots
        fig = self.Figure()
        ax = fig.add_subplot(111)
        plots.draw_stack(ax, [iss(1000.0)], scale=viewdata.RATIO, **kw)
        fig.canvas.draw()
        return ax

    def test_label_and_direction(self):
        ax = self.panel()
        self.assertEqual(ax.get_xlabel(), "Energy ratio (E/E$_0$)")
        lo, hi = ax.get_xlim()
        self.assertLess(lo, hi)                    # ISS reads left to right
        self.assertAlmostEqual(ax.lines[0].get_xdata()[-1], 0.933)

    def test_a_kinetic_marker_sits_at_its_ratio(self):
        ax = self.panel(markers=[(933.0, "Au", True)])
        xs = [l.get_xdata()[0] for l in ax.lines if len(l.get_xdata()) == 2
              and l.get_xdata()[0] == l.get_xdata()[1]]
        self.assertTrue(any(abs(x - 0.933) < 1e-9 for x in xs), xs)

    def test_the_cursor_line_sits_at_its_ratio(self):
        ax = self.panel(norm="At cursor", cursor=500.0)
        vlines = [l for l in ax.lines
                  if len(set(l.get_xdata())) == 1 and len(l.get_xdata()) == 2]
        self.assertTrue(any(abs(l.get_xdata()[0] - 0.5) < 1e-9
                            for l in vlines))


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.dir = tempfile.mkdtemp()
        cls.avg = os.path.join(cls.dir, "ISS Survey.avg")
        with open(cls.avg, "w") as fh:
            fh.write(iss_avg())
        cls.xps = os.path.join(cls.dir, "C1s Scan.avg")
        with open(cls.xps, "w") as fh:
            fh.write(AVG_TEMPLATE.format(
                title="C1s Scan",
                rows="LIST@   0=  1.0,  2.0,  3.0,  4.0\n"
                     "LIST@   4=  5.0,  6.0,  7.0,  8.0"))
        cls.kal = os.path.join(cls.dir, "iss.kal")
        with open(cls.kal, "w") as fh:
            fh.write("Dataset filename          = iss.dset\n"
                     + iss_object("ISS 1"))

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        import shutil
        shutil.rmtree(cls.dir, ignore_errors=True)

    def open(self, path):
        ws = ee.Workspace(self.root)
        ws._add_files([path])
        region = ws.docs[0].regions[0]
        ws.checked = {id(region)}
        ws.iss_axis_var.set("Energy ratio")
        ws._render()
        return ws, region

    def axis(self, ws):
        return next(a for a in ws.fig.axes if a.get_xlabel())

    def event(self, ws, x):
        ax = self.axis(ws)
        return types.SimpleNamespace(inaxes=ax, xdata=x, ydata=1.0)

    def test_a_thermo_file_plots_on_its_calibrated_ratio(self):
        ws, r = self.open(self.avg)
        ax = self.axis(ws)
        self.assertEqual(ax.get_xlabel(), "Energy ratio (E/E$_0$)")
        line = next(l for l in ax.lines if len(l.get_xdata()) == len(r.energy))
        self.assertAlmostEqual(line.get_xdata()[-1], r.energy[-1] / 966.313)
        self.assertEqual(ws._axe0[ax], 966.313)

    def test_clicks_and_the_read_out_use_the_beam_energy(self):
        ws, r = self.open(self.avg)
        ev = self.event(ws, 0.9)
        self.assertAlmostEqual(ws.kinetic_from_event(r, ev), 0.9 * 966.313)
        self.assertIsNone(ws.displayed_be(ev))      # not a binding energy
        ws._on_motion(ev)
        self.assertIn("E/E₀ 0.9000", ws.cursor_lbl.cget("text"))
        self.assertIn("KE 869.7 eV", ws.cursor_lbl.cget("text"))

    def test_a_kratos_file_waits_for_a_confirmed_energy(self):
        ws, r = self.open(self.kal)
        self.assertEqual(self.axis(ws).get_xlabel(),
                         "Kinetic energy (eV)")
        self.assertTrue(any("beam energy unknown" in n
                            for n in ws._view_notes))
        ws.iss_save(r, {"ion": "He+", "e0": 1000.0, "theta": 135.0})
        ws._render()
        self.assertEqual(self.axis(ws).get_xlabel(), "Energy ratio (E/E$_0$)")
        self.assertEqual(ws._axe0[self.axis(ws)], 1000.0)

    def test_kinetic_is_the_default_and_the_switch_is_one_box(self):
        ws = ee.Workspace(self.root)
        self.assertEqual(ws.iss_axis_var.get(), "Kinetic")
        ws._add_files([self.avg])
        r = ws.docs[0].regions[0]
        ws.checked = {id(r)}
        ws._render()
        self.assertEqual(self.axis(ws).get_xlabel(), "Kinetic energy (eV)")
        ws.iss_axis_var.set("Energy ratio")
        ws._on_view_changed()
        ws._render()
        self.assertEqual(self.axis(ws).get_xlabel(), "Energy ratio (E/E$_0$)")
        self.assertEqual(ws.scale_var.get(), "Binding")     # XPS box untouched

    def test_xps_keeps_its_own_scale_beside_ratio_iss(self):
        ws = ee.Workspace(self.root)
        ws._add_files([self.avg, self.xps])
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        ws.scale_var.set("Kinetic")
        ws.iss_axis_var.set("Energy ratio")
        ws._render()
        labels = sorted(a.get_xlabel() for a in ws.fig.axes if a.get_xlabel())
        self.assertEqual(labels, ["Energy ratio (E/E$_0$)",
                                  "Kinetic energy (eV)"])

    def test_the_axis_is_saved_with_the_view_and_old_states_stay_kinetic(self):
        ws, _r = self.open(self.avg)
        st = ws.capture_state()
        self.assertEqual((st["energy_scale"], st["iss_axis"]),
                         ("Binding", "Energy ratio"))
        self.assertIn("energy-ratio axis",
                      __import__("pptx_export").look_notes(st))
        old = dict(st)
        del old["iss_axis"]                  # a state from before the choice
        ws.apply_state(old, render=False)
        self.assertEqual(ws.iss_axis_var.get(), "Kinetic")
        first = dict(st, energy_scale="Energy ratio")  # the first version
        del first["iss_axis"]
        ws.apply_state(first, render=False)
        self.assertEqual((ws.scale_var.get(), ws.iss_axis_var.get()),
                         ("Binding", "Energy ratio"))

    def test_the_iss_dialog_has_the_same_switch(self):
        from iss_ui import IssReelsDialog
        ws, _r = self.open(self.avg)
        ws.iss_axis_var.set("Kinetic")
        dlg = IssReelsDialog(self.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        self.assertFalse(dlg.ratio_var.get())
        dlg.ratio_var.set(True)
        dlg._ratio_toggled()
        self.assertEqual(ws.iss_axis_var.get(), "Energy ratio")
        ws.iss_axis_var.set("Kinetic")             # changed in the main window
        self.assertFalse(dlg.ratio_var.get())


if __name__ == "__main__":
    unittest.main()
