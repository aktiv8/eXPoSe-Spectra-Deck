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


class TestPanelShapeOverride(unittest.TestCase):
    def test_sanitise_resolve_describe(self):
        import panelview
        self.assertEqual(panelview.sanitise({"shape": "Square"}),
                         {"shape": "Square"})
        self.assertEqual(panelview.sanitise({"shape": "Round"}), {})
        defaults = {"view": "Stack", "norm": "None", "offset": 0.6,
                    "z_axis": "Auto", "reverse": False, "fit_show": {},
                    "ident_show": {}, "shape": "Auto"}
        self.assertEqual(panelview.resolve(defaults, {"shape": "Portrait 2:3"}
                                           )["shape"], "Portrait 2:3")
        self.assertEqual(panelview.resolve(defaults)["shape"], "Auto")
        self.assertIn("panel shape Portrait 2:3",
                      panelview.describe({"shape": "Portrait 2:3"}, defaults))
        views = panelview.with_value({}, "C 1s", "shape", "Square")
        self.assertEqual(views, {"C 1s": {"shape": "Square"}})
        self.assertEqual(panelview.reset(views, "C 1s"), {})


class TestPlace(unittest.TestCase):
    def test_plain_grid_is_row_major(self):
        import panelview
        rows, cells = panelview.place([(1, 1)] * 5, 2)
        self.assertEqual(rows, 3)
        self.assertEqual([(c.r0, c.c0) for c in cells],
                         [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)])
        self.assertEqual([c.bottom_row for c in cells],
                         [False, False, False, True, True])
        self.assertEqual([c.first_col for c in cells],
                         [True, False, True, False, True])

    def test_tall_panel_beside_two_stacked(self):
        import panelview
        rows, cells = panelview.place([(2, 1), (1, 1), (1, 1)], 2)
        self.assertEqual(rows, 2)
        tall, top, bottom = cells
        self.assertEqual((tall.r0, tall.r1, tall.c0, tall.c1), (0, 2, 0, 1))
        self.assertEqual((top.r0, top.c0), (0, 1))
        self.assertEqual((bottom.r0, bottom.c0), (1, 1))
        self.assertTrue(tall.bottom_row and bottom.bottom_row)
        self.assertFalse(top.bottom_row)

    def test_wide_panel_and_clamping(self):
        import panelview
        rows, cells = panelview.place([(1, 1), (1, 2), (1, 1)], 2)
        # the wide one cannot start in column 1, so it goes to the next row
        self.assertEqual((cells[1].r0, cells[1].c0, cells[1].c1), (1, 0, 2))
        self.assertEqual(cells[2].r0, 2)
        rows, cells = panelview.place([(1, 2)], 1)      # wider than the grid
        self.assertEqual((cells[0].c0, cells[0].c1), (0, 1))
        self.assertEqual(panelview.place([], 3)[0], 1)

    def test_span_key(self):
        import panelview
        self.assertEqual(panelview.sanitise({"span": "Tall (2 rows)"}),
                         {"span": "Tall (2 rows)"})
        self.assertEqual(panelview.sanitise({"span": "1 × 1"}), {})
        self.assertEqual(panelview.sanitise({"span": "huge"}), {})


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

    def test_every_panel_can_be_reached(self):
        """3 panels, 2 a page in 2 columns: the third was never shown."""
        ws = self.ws
        regs = [r for d in ws.docs for r in d.regions]
        ws.checked = {id(r) for r in regs[:3]}
        ws.panels_var.set("2")
        ws.cols_var.set("2")
        keys = [k for k, _ in ws._groups()]
        self.assertEqual(len(keys), 3)
        ws.panel_start = 99                     # as far as it will go
        ws._render()
        self.assertEqual(ws.panel_start, 2)
        self.assertEqual(list(ws._axmap.values()), keys[2:])
        ws.prev_page()
        ws._render()
        self.assertEqual(ws.panel_start, 0)
        ws.next_page()
        ws._render()
        self.assertEqual(ws.panel_start, 2)
        self.assertEqual(list(ws._axmap.values()), keys[2:])

    def test_page_turns_never_skip_a_panel(self):
        ws = self.ws
        ws.panels_var.set("6")
        ws.cols_var.set("4")                    # two ragged rows: 4 + 2
        self.assertEqual(ws._page_step(10), 4)  # not 6: panels 5-6 would skip
        ws.cols_var.set("2")
        self.assertEqual(ws._page_step(10), 6)

    def test_tall_panel_layout_on_the_page(self):
        ws = self.ws
        regs = [r for d in ws.docs for r in d.regions]
        ws.checked = {id(r) for r in regs[:3]}
        ws.cols_var.set("2")
        keys = [k for k, _ in ws._groups()]
        ws._set_panel_view(keys[0], "span", "Tall (2 rows)")
        self.axes()
        by_key = {ws._axmap[a]: a.get_position() for a in ws.fig.axes
                  if a in ws._axmap}
        tall, top, bottom = (by_key[k] for k in keys)
        self.assertGreater(tall.height, 1.6 * top.height)   # spans both rows
        self.assertAlmostEqual(top.x0, bottom.x0, places=3)  # stacked
        self.assertGreater(top.x0, tall.x1)                  # to its right
        self.assertGreater(top.y0, bottom.y0)
        st = ws.capture_state()
        self.assertEqual(st["panel_views"][keys[0]]["span"], "Tall (2 rows)")

    def test_one_panel_with_its_own_shape(self):
        ws = self.ws
        ws.cols_var.set("2")                      # page shape stays Auto
        keys = [k for k, _ in ws._groups()]
        ws._set_panel_view(keys[0], "shape", "Portrait 2:3")
        self.axes()
        by_key = {ws._axmap[a]: a for a in ws.fig.axes if a in ws._axmap}
        self.assertAlmostEqual(by_key[keys[0]].get_box_aspect(), 1.5)
        for k in keys[1:]:
            self.assertIsNone(by_key[k].get_box_aspect())
        # saved with the look, and gone with the reset
        st = ws.capture_state()
        self.assertEqual(st["panel_views"][keys[0]]["shape"], "Portrait 2:3")
        ws.reset_panel_views()
        self.axes()
        self.assertTrue(all(a.get_box_aspect() is None for a in ws.fig.axes
                            if a in ws._axmap))

    def test_panel_override_beats_the_page_shape(self):
        ws = self.ws
        ws.shape_var.set("Square")
        keys = [k for k, _ in ws._groups()]
        ws._set_panel_view(keys[0], "shape", "Auto")     # fill the cell
        self.axes()
        by_key = {ws._axmap[a]: a for a in ws.fig.axes if a in ws._axmap}
        self.assertIsNone(by_key[keys[0]].get_box_aspect())
        self.assertAlmostEqual(by_key[keys[1]].get_box_aspect(), 1.0)

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
