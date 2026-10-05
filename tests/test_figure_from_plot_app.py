"""Right-click and Identify while the zoom tool is on, and saving the current
view (with its zoom) as a figure from the plot (needs a display and
matplotlib / Tk; skipped otherwise).

Run:  python -m unittest discover tests
"""

import math
import os
import sys
import unittest
from types import SimpleNamespace
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

from readers import Region, SpectrumFile  # noqa: E402


def survey(name="C 1s", lo=280.0, hi=292.0, peak=285.0, n=121):
    e = [hi - i * (hi - lo) / (n - 1) for i in range(n)]
    c = [100.0 + 900.0 * math.exp(-0.5 * ((x - peak) / 0.8) ** 2) for x in e]
    return Region(name=name, index=0, offset=0, energy=e, counts=c,
                  decodable=True, sample="S", photon_energy=1486.6,
                  source="s.vms", count_units="counts/s")


def doc(path, regions):
    f = SpectrumFile()
    f.path = path
    f.format_name = "Test"
    f.regions = regions
    f.instrument = {"Instrument": "Test Spec"}
    f._finish()
    return f


def mouse(button, x, y, ax=None, xdata=None):
    """A matplotlib-like mouse event."""
    return SimpleNamespace(button=button, x=x, y=y, inaxes=ax, xdata=xdata,
                           ydata=None,
                           guiEvent=SimpleNamespace(x_root=5, y_root=5))


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestFigureFromThePlot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.docs = {"s.vms": doc("s.vms", [survey("C 1s"),
                                           survey("O 1s", 525.0, 540.0, 531.0)])}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()

    def setUp(self):
        ws = self.ws
        for d in list(ws.docs):
            ws.docs.remove(d)
        ws.file_ids.clear()
        ws.region_parser.clear()
        ws.figures.clear()
        ws.group_var.set("Element, per sample")
        ws.view_var.set("Stack")
        ws.norm_var.set("None")
        ws.scale_var.set("Binding")
        ws.panels_var.set("Auto")
        ws.traces_var.set("All")
        ws.panel_views = {}
        ws._click_cb = ws._pick_cb = None
        ws._press = None
        ws._zoom_request = None
        ws._disp_cache.clear()
        ws._axmap, ws._zoom_sig = {}, {}         # nothing zoomed carries over
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        ws._render()
        self.addCleanup(self._tool, "")

    def _tool(self, mode):
        """Pretend the toolbar's zoom / pan tool is (not) switched on."""
        self.ws.toolbar.mode = mode

    def _panel(self, key="C 1s"):
        ws = self.ws
        for ax, k in ws._axmap.items():
            if k.startswith(key):
                return ax, k
        self.fail(f"no panel for {key}: {list(ws._axmap.values())}")

    # -- clicks with the zoom tool on -----------------------------------------
    def test_right_click_opens_the_menu_with_the_zoom_tool_on(self):
        ax, _k = self._panel()
        self._tool("zoom rect")
        with mock.patch.object(self.ws, "_panel_menu") as menu:
            self.ws._on_plot_click(mouse(3, 100, 100, ax, 285.0))
            menu.assert_not_called()                  # not yet: it may be a drag
            self.ws._on_plot_release(mouse(3, 102, 99, ax, 285.0))
            menu.assert_called_once()

    def test_a_drag_with_the_tool_on_is_left_to_the_tool(self):
        ax, _k = self._panel()
        self._tool("zoom rect")
        with mock.patch.object(self.ws, "_panel_menu") as menu:
            self.ws._on_plot_click(mouse(3, 100, 100, ax, 285.0))
            self.ws._on_plot_release(mouse(3, 160, 140, ax, 288.0))
            menu.assert_not_called()

    def test_right_click_still_opens_the_menu_at_once_without_a_tool(self):
        ax, _k = self._panel()
        with mock.patch.object(self.ws, "_panel_menu") as menu:
            self.ws._on_plot_click(mouse(3, 100, 100, ax, 285.0))
            menu.assert_called_once()

    def test_right_click_outside_every_panel_offers_the_figure_command(self):
        with mock.patch.object(self.ws, "_figure_menu") as menu:
            self.ws._on_plot_click(mouse(3, 5, 5, None, None))
            menu.assert_called_once()

    def test_identify_gets_a_left_click_with_the_zoom_tool_on(self):
        ax, _k = self._panel()
        got = []
        self.ws._click_cb = lambda be, ev=None: got.append(be)
        self._tool("zoom rect")
        self.ws._on_plot_click(mouse(1, 100, 100, ax, 285.0))
        self.assertEqual(got, [])                    # waits for the release
        self.ws._on_plot_release(mouse(1, 101, 100, ax, 285.0))
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got[0], 285.0, delta=0.5)

    def test_a_right_click_is_not_an_identify_pick(self):
        ax, _k = self._panel()
        got = []
        self.ws._click_cb = lambda be, ev=None: got.append(be)
        with mock.patch.object(self.ws, "_panel_menu"):
            self.ws._on_plot_click(mouse(3, 100, 100, ax, 285.0))
        self.assertEqual(got, [])

    # -- the figure -------------------------------------------------------------
    def _zoom_c1s(self):
        ax, key = self._panel("C 1s")
        ax.set_xlim(288.0, 283.0)
        ax.set_ylim(0.0, 500.0)
        return ax, key

    def test_a_figure_keeps_only_the_panels_that_are_zoomed(self):
        _ax, key = self._zoom_c1s()
        fig = self.ws.add_figure("Zoomed")
        z = fig["state"]["zoom"]
        self.assertEqual(z["scale"], "Binding")
        self.assertEqual(list(z["panels"]), [key])
        self.assertEqual(z["panels"][key]["x"], [288.0, 283.0])

    def test_an_unzoomed_view_stores_an_empty_zoom(self):
        fig = self.ws.add_figure("Plain")
        self.assertEqual(fig["state"]["zoom"]["panels"], {})

    def test_a_zoom_does_not_change_the_workbook_signature(self):
        before = self.ws._signature()
        self._zoom_c1s()
        self.assertEqual(self.ws._signature(), before)
        self.assertNotIn("zoom", self.ws.capture_state())

    def test_save_view_asks_for_a_name_and_says_where_to_caption(self):
        self._zoom_c1s()
        with mock.patch.object(ee.simpledialog, "askstring",
                               return_value="  Carbon  "):
            self.ws.save_view_as_figure()
        self.assertEqual([f["name"] for f in self.ws.figures], ["Carbon"])
        self.assertIn("Carbon", self.ws.status.cget("text"))

    def test_cancelling_the_name_saves_nothing(self):
        with mock.patch.object(ee.simpledialog, "askstring", return_value=None):
            self.ws.save_view_as_figure()
        self.assertEqual(self.ws.figures, [])

    def test_nothing_ticked_is_said_not_saved(self):
        self.ws.checked = set()
        self.ws._render()
        with mock.patch.object(ee.messagebox, "showinfo") as info:
            self.ws.save_view_as_figure()
        info.assert_called_once()
        self.assertEqual(self.ws.figures, [])

    def test_recall_puts_the_zoom_back_and_clears_a_newer_one(self):
        _ax, key = self._zoom_c1s()
        fig = self.ws.add_figure("Zoomed")
        # change the live zoom, then recall
        ax, _k = self._panel("C 1s")
        ax.set_xlim(292.0, 280.0)
        self.ws.apply_state(fig["state"])
        self.ws._render()
        ax, _k = self._panel("C 1s")
        self.assertEqual(tuple(ax.get_xlim()), (288.0, 283.0))
        self.assertEqual(tuple(ax.get_ylim()), (0.0, 500.0))
        o, _k = self._panel("O 1s")                       # never zoomed
        self.assertAlmostEqual(o.get_xlim()[0], 540.0, delta=1.0)

    def test_an_older_figure_without_a_zoom_leaves_the_live_zoom_alone(self):
        self._zoom_c1s()
        old = self.ws.capture_state()                       # no "zoom" key
        self.ws.apply_state(old)
        self.assertIsNone(self.ws._zoom_request)
        self.ws._render()
        ax, _k = self._panel("C 1s")
        self.assertEqual(tuple(ax.get_xlim()), (288.0, 283.0))

    def test_report_pages_show_the_zoom_and_the_live_view_is_untouched(self):
        _ax, _key = self._zoom_c1s()
        fig = self.ws.add_figure("Zoomed")
        before = (self.ws.capture_state(), self.ws.current_zoom())
        pages = []
        self.ws._render_figure_pages(fig, pages.append, decorate=False)
        self.assertTrue(pages)
        limits = {}
        for page in pages:
            for ax in page.axes:
                limits[ax.get_title() or ax.get_ylabel() or id(ax)] = (
                    tuple(ax.get_xlim()))
        # one of the page's panels is the zoomed C 1s
        self.assertIn((288.0, 283.0), limits.values())
        # and the live view is as it was
        self.assertEqual(self.ws.capture_state(), before[0])
        self.assertEqual(self.ws.current_zoom(), before[1])
        self.assertIsNone(self.ws._zoom_request)

    def test_a_saved_range_outside_the_data_is_ignored(self):
        _ax, key = self._zoom_c1s()
        fig = self.ws.add_figure("Zoomed")
        fig["state"]["zoom"]["panels"][key]["x"] = [900.0, 880.0]
        pages = []
        self.ws._render_figure_pages(fig, pages.append, decorate=False)
        xs = [tuple(ax.get_xlim()) for page in pages for ax in page.axes]
        self.assertNotIn((900.0, 880.0), xs)

    def test_the_figures_dialog_routes_use_the_same_state(self):
        _ax, key = self._zoom_c1s()
        fig = self.ws.add_figure("One")
        ax, _k = self._panel("C 1s")
        ax.set_xlim(290.0, 282.0)
        with mock.patch.object(ee.workbook_ui.messagebox, "askyesno",
                               return_value=True):
            dlg = ee.workbook_ui.FiguresDialog(self.root, self.ws)
            self.addCleanup(dlg.destroy)
            dlg.list.selection_set(0)
            dlg._update()
        self.assertEqual(fig["state"]["zoom"]["panels"][key]["x"],
                         [290.0, 282.0])


if __name__ == "__main__":
    unittest.main()
