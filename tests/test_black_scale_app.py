"""The single-colour ("Black") trace scale and the bundled plot fonts in the
real window (needs a display and matplotlib; skipped otherwise).

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
import themes  # noqa: E402


def reg(name, sample, peak, level=None, n=61):
    e = [292.0 - i * 12.0 / (n - 1) for i in range(n)]
    c = [100 + 900 * math.exp(-((x - peak) / 1.0) ** 2) for x in e]
    return Region(name=name, index=0, offset=0, energy=e, counts=c,
                  decodable=True, sample=sample, photon_energy=1486.6,
                  pass_energy=20.0, dwell=0.1, step=0.2, etch_level=level,
                  etch_time=None if level is None else 30.0 * level,
                  source="a.vms", count_units="counts/s")


def doc(path, regions):
    f = SpectrumFile()
    f.path = path
    f.format_name = "Test"
    f.regions = regions
    f.instrument = {"Instrument": "Test Spec"}
    f._finish()
    return f


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestBlackScaleInApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {"depth.vms": doc("depth.vms", [
            reg("C 1s", "Depth", 286 - i / 5, level=i + 1) for i in range(4)])}
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
        ws = self.ws
        for d in list(ws.docs):
            ws.docs.remove(d)
        ws.file_ids.clear()
        ws.region_parser.clear()
        ws.group_var.set("Element, per sample")
        ws.view_var.set("Stack")
        ws.colscale_var.set("Theme default")
        ws.colrev_var.set(False)
        ws.panel_views = {}
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}

    def line_colours(self):
        from matplotlib.colors import to_hex
        self.ws._render()
        cols = []
        for ax in self.ws.fig.axes:
            for ln in ax.get_lines():
                if len(ln.get_xdata()) > 10:
                    cols.append(to_hex(ln.get_color()))
        return cols

    def test_the_choice_is_listed_and_offers_no_reverse(self):
        self.assertIn(themes.BLACK_SCALE, themes.TRACE_SCALE_NAMES)
        self.ws.colscale_var.set(themes.BLACK_SCALE)
        self.ws._colour_scale_changed()
        self.assertIn("disabled", self.ws.colrev_btn.state())
        self.ws.colscale_var.set("Viridis")
        self.ws._colour_scale_changed()
        self.assertNotIn("disabled", self.ws.colrev_btn.state())

    def test_every_trace_is_one_colour_and_the_tree_agrees(self):
        ws = self.ws
        ws.colscale_var.set(themes.BLACK_SCALE)
        cols = self.line_colours()
        self.assertEqual(len(cols), 4)
        self.assertEqual(len(set(cols)), 1)
        ink = themes.scale_colours(themes.BLACK_SCALE, False, 1,
                                   ws.palette)[0]
        self.assertEqual(cols[0].upper(), ink.upper())
        self.assertEqual({c.upper() for c in ws.trace_color.values()},
                         {ink.upper()})

    def test_it_is_part_of_the_saved_look_and_a_bad_name_is_ignored(self):
        ws = self.ws
        ws.colscale_var.set(themes.BLACK_SCALE)
        st = ws.capture_state()
        self.assertEqual(st["colour_scale"], themes.BLACK_SCALE)
        ws.colscale_var.set("Viridis")
        ws.apply_state(st)
        self.assertEqual(ws.colscale_var.get(), themes.BLACK_SCALE)
        self.assertIn("disabled", ws.colrev_btn.state())
        ws.apply_state({"colour_scale": "Not a scale"})
        self.assertEqual(ws.colscale_var.get(), themes.BLACK_SCALE)

    def test_a_heatmap_still_draws_with_it(self):
        ws = self.ws
        ws.colscale_var.set(themes.BLACK_SCALE)
        ws.view_var.set("Heatmap")
        ws._render()
        self.assertTrue(ws.fig.axes)

    def test_a_bundled_plot_font_reaches_the_drawn_text(self):
        import plotstyle
        import fonts
        ws = self.ws
        old = ws.plot_style
        self.addCleanup(lambda: ws.set_plot_style(old, save=False))
        for fam in ("Inter", "STIX Two Text"):
            self.assertIn(fam, fonts.FAMILIES)
            ws.set_plot_style(plotstyle.sanitise({"font": fam}), save=False)
            ws._render()
            ax = next(a for a in ws.fig.axes if a.get_xlabel())
            self.assertEqual(ax.xaxis.label.get_fontname(), fam)


if __name__ == "__main__":
    unittest.main()
