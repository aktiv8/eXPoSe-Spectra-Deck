"""Columns and panel shape of the plot page (the grid helper, the real window,
saved looks). The window tests need a display and matplotlib; skipped
otherwise.

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


def reg(name, sample, peak, n=61):
    e = [292.0 - i * 12.0 / (n - 1) for i in range(n)]
    c = [100 + 900 * math.exp(-((x - peak) / 1.0) ** 2) for x in e]
    return Region(name=name, index=0, offset=0, energy=e, counts=c,
                  decodable=True, sample=sample, photon_energy=1486.6,
                  pass_energy=20.0, dwell=0.1, step=0.2, source="a.vms",
                  count_units="counts/s")


def doc(path, regions):
    f = SpectrumFile()
    f.path = path
    f.format_name = "Test"
    f.regions = regions
    f.instrument = {"Instrument": "Test Spec"}
    f._finish()
    return f


@unittest.skipUnless(ee is not None, "Tk not available")
class TestGridDims(unittest.TestCase):
    def test_auto_unchanged(self):
        self.assertEqual(ee._grid_dims(4), (2, 2))
        self.assertEqual(ee._grid_dims(1), (1, 1))

    def test_forced_columns(self):
        self.assertEqual(ee._grid_dims(4, 2), (2, 2))
        self.assertEqual(ee._grid_dims(6, 2), (3, 2))
        self.assertEqual(ee._grid_dims(5, 2), (3, 2))
        self.assertEqual(ee._grid_dims(2, 1), (2, 1))

    def test_columns_never_exceed_panels(self):
        self.assertEqual(ee._grid_dims(1, 4), (1, 1))
        self.assertEqual(ee._grid_dims(3, 4), (1, 3))


class TestSizePresets(unittest.TestCase):
    def test_presets_inside_dialog_limits(self):
        import plotstyle
        for name, (w, h) in plotstyle.SIZE_PRESETS.items():
            self.assertTrue(1 <= w <= 40 and 1 <= h <= 40, name)

    def test_fit_height(self):
        import plotstyle
        # two columns of 3:4 portrait panels, 7 in wide: cells 3.5 x 4.67 in
        self.assertAlmostEqual(plotstyle.fit_height(7.0, 1, 2, 4 / 3), 4.667,
                               places=2)
        self.assertAlmostEqual(plotstyle.fit_height(7.0, 2, 2, 4 / 3), 9.333,
                               places=2)
        self.assertIsNone(plotstyle.fit_height(7.0, 1, 2, None))
        self.assertEqual(plotstyle.fit_height(40.0, 8, 1, 1.5), 40.0)


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestPanelGridInApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {
            "a.vms": doc("a.vms", [reg("C 1s", "A", 286.0),
                                   reg("O 1s", "A", 288.0)]),
            "b.vms": doc("b.vms", [reg("C 1s", "B", 287.0),
                                   reg("O 1s", "B", 289.0)])}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        cls._load_config = ee.load_config
        ee.load_config = lambda: {}
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.load_config = cls._load_config
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
        ws.norm_var.set("None")
        ws.scale_var.set("Binding")
        ws.panels_var.set("Auto")
        ws.cols_var.set("Auto")
        ws.shape_var.set("Auto")
        ws.traces_var.set("All")
        ws.panel_views = {}
        ws.panel_order = []
        ws.trace_start = 0
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}

    def axes(self):
        self.ws._render()
        self.ws.fig.canvas.draw()
        return [a for a in self.ws.fig.axes if a in self.ws._axmap]

    def test_four_panels_one_column(self):
        self.ws.cols_var.set("1")
        axs = self.axes()
        self.assertEqual(len(axs), 4)
        x0s = {round(a.get_position().x0, 3) for a in axs}
        self.assertEqual(len(x0s), 1)

    def test_four_panels_two_columns(self):
        self.ws.cols_var.set("2")
        axs = self.axes()
        x0s = {round(a.get_position().x0, 3) for a in axs}
        y0s = {round(a.get_position().y0, 3) for a in axs}
        self.assertEqual((len(x0s), len(y0s)), (2, 2))

    def test_portrait_shape(self):
        self.ws.cols_var.set("2")
        self.ws.shape_var.set("Portrait 3:4")
        axs = self.axes()
        fw, fh = self.ws.fig.get_size_inches()
        for a in axs:
            p = a.get_position()
            self.assertAlmostEqual((p.height * fh) / (p.width * fw),
                                   4 / 3, places=2)

    def test_report_per_page_follows_grid(self):
        ws = self.ws
        self.assertEqual(ws._pdf_per_page((11.7, 8.3), (0, .17, 1, .94)), 6)
        ws.cols_var.set("4")
        ws.shape_var.set("Portrait 3:4")
        # slide strip 12.1 x 4.95: one row of four 3.0 x 4.0 in panels
        self.assertEqual(ws._pdf_per_page((12.1, 4.95), None), 4)
        # A4 portrait page: 4 columns of 2.08 x 2.77 in cells, 9.0 in usable
        self.assertEqual(ws._pdf_per_page((8.3, 11.7), (0, .17, 1, .94)), 12)
        ws.panels_var.set("2")                    # a number always wins
        self.assertEqual(ws._pdf_per_page((12.1, 4.95), None), 2)

    def test_figure_pages_use_columns(self):
        ws = self.ws
        ws.cols_var.set("4")
        ws.shape_var.set("Portrait 3:4")
        st = ws.capture_state()
        pages = []
        ws._render_figure_pages({"name": "F", "caption": "", "state": st},
                                pages.append, size=(12.1, 4.95), rect=None,
                                decorate=False)
        self.assertEqual(len(pages), 1)          # all four on one slide
        axs = [a for a in pages[0].axes][:4]
        self.assertEqual(len({round(a.get_position().y0, 3)
                              for a in axs}), 1)

    def test_save_dialog_preset_and_fit(self):
        self.ws.cols_var.set("2")
        self.ws.shape_var.set("Portrait 3:4")
        self.ws.save_grid = (2, 2, 4 / 3)
        dlg = ee.workbook_ui.ImageOptionsDialog(self.root, self.ws,
                                                lambda *a: None)
        try:
            dlg.apply_preset("Double column (180 mm)")
            self.assertEqual(dlg.vars["w"].get(), "7.09")
            dlg.fit_to_grid()
            self.assertAlmostEqual(float(dlg.vars["h"].get()),
                                   2 * 7.09 / 2 * 4 / 3, places=1)
            self.assertEqual(dlg.preset.get(), "Custom")
        finally:
            dlg.destroy()

    def test_state_round_trip_and_old_state(self):
        self.ws.cols_var.set("2")
        self.ws.shape_var.set("Square")
        st = self.ws.capture_state()
        self.assertEqual((st["panel_columns"], st["panel_shape"]),
                         ("2", "Square"))
        self.ws.cols_var.set("Auto")
        self.ws.shape_var.set("Auto")
        self.ws.apply_state(st, render=False)
        self.assertEqual((self.ws.cols_var.get(), self.ws.shape_var.get()),
                         ("2", "Square"))
        st.pop("panel_columns")
        st.pop("panel_shape")
        self.ws.apply_state(st, render=False)         # not the live value
        self.assertEqual((self.ws.cols_var.get(), self.ws.shape_var.get()),
                         ("Auto", "Auto"))


if __name__ == "__main__":
    unittest.main()
