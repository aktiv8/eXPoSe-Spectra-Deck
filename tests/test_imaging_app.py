"""Kratos imaging maps in the app: the tree row, opening the series viewer,
stepping through frames, the filter, an area, saving, and that the report
builders leave these maps alone (needs a display and matplotlib / Tk; skipped
otherwise).

Run:  python -m unittest discover tests
"""

import os
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

import imagepages  # noqa: E402
from readers import load_file  # noqa: E402
from test_kratosmap import NX, NY, kal_map, kal_position, kal_text  # noqa: E402


def imaging_first(node):
    from spectradeck import regions_under
    return regions_under(node)[0]


def focus_series_file(directory):
    """Three maps at one position with the stage height stepping 30 um, then
    one at another position, as a Kratos focus run would write them."""
    objs = [kal_position("Grid on Tape", 1)]
    for k, z in ((2, 900), (3, 930), (4, 960)):
        objs.append(kal_map("Au 4f", k, z_m=z / 1e6,
                            when=f"17/09/08 09:0{k}:00"))
    objs.append(kal_position("Elsewhere", 5))
    objs.append(kal_map("Au 4f", 6, x_m=0.05, z_m=960 / 1e6,
                        when="17/09/08 09:30:00"))
    p = os.path.join(directory, "focus.kal")
    with open(p, "w") as fh:
        fh.write(kal_text(*objs))
    return p


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestImagingApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = focus_series_file(cls.tmp.name)
        cls.ws = ee.Workspace(cls.root)
        cls.parser = load_file(cls.path)
        cls.problems = cls.ws._register_parser(cls.parser, cls.path)

    @classmethod
    def tearDownClass(cls):
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()
        cls.tmp.cleanup()

    def open_dialog(self, region):
        made = []
        real = ee.imaging_ui.MapSeriesDialog

        class Recorded(real):
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                made.append(self)

        old = ee.imaging_ui.MapSeriesDialog
        ee.imaging_ui.MapSeriesDialog = Recorded
        try:
            self.ws.open_snapmap(region)
        finally:
            ee.imaging_ui.MapSeriesDialog = old
        self.assertEqual(len(made), 1)
        self.addCleanup(made[0].destroy)
        return made[0]

    @property
    def maps(self):
        return [r for r in self.parser.regions if r.extra.get("cube") is not None]

    def test_loads_without_a_problem_and_rows_are_not_tickable(self):
        self.assertEqual(self.problems, [])
        self.assertEqual(len(self.maps), 4)
        labels = [self.ws.tree.item(i, "text") for i in self.ws.node_map]
        self.assertTrue(any("(Map)" in t for t in labels))
        # nothing to plot, so no row of this file carries a tick box
        self.assertFalse(self.ws.leaf_ids)

    def test_a_focus_series_opens_on_its_position_and_plots_against_z(self):
        dlg = self.open_dialog(self.maps[1])
        self.assertEqual(dlg.how.get(), "Same stage position")
        self.assertEqual(len(dlg.shown), 3)
        self.assertEqual(dlg.frame.region, self.maps[1])
        xl, xs, yl, ys = dlg.series()
        self.assertEqual(xl, "Stage Z (um)")
        self.assertEqual(xs, [900.0, 930.0, 960.0])
        self.assertEqual(len(ys), 3)
        self.assertIn("whole image", yl)

    def test_stepping_filtering_and_the_area(self):
        import numpy as np
        dlg = self.open_dialog(self.maps[0])
        dlg.step(1)
        self.assertEqual(dlg.i, 1)
        dlg.step(10)                                    # stops at the last frame
        self.assertEqual(dlg.i, 2)
        dlg.set_frame(0)
        dlg.set_filter("All maps")
        self.assertEqual(len(dlg.shown), 4)
        self.assertEqual(dlg.frame.region, self.maps[0])   # stays on the frame
        mask = np.zeros((NY, NX), bool)
        mask[:2, :2] = True
        dlg.set_roi(mask)
        self.assertIn("area", dlg.series()[2])
        means = dlg.series()[3]
        want = [float(np.mean(
            np.asarray(r.extra["cube"].array3d()[:, :, 0])[:2, :2]))
            for r in self.maps]
        for a, b in zip(means, want):
            self.assertAlmostEqual(a, b, 4)
        dlg.clear_roi()
        self.assertIn("whole image", dlg.series()[2])

    def test_saved_values_are_the_recorded_counts_even_when_smoothed(self):
        dlg = self.open_dialog(self.maps[0])
        dlg.sigma.set(2.0)
        dlg._rebuild()
        text = dlg.map_csv()
        first = text.splitlines()[2].split(",")           # note, X row, Y row 0
        raw =self.maps[0].extra["cube"].array3d()[0, :, 0].tolist()
        self.assertEqual([float(v) for v in first[1:]], raw)
        rows = dlg.series_table()
        self.assertEqual(rows[0]["position"], "Grid on Tape")
        self.assertEqual(rows[0]["stage_z_um"], 900.0)
        self.assertIn("sharpness", rows[0])

    def test_the_picture_draws_with_a_scale_bar_and_a_colour_bar(self):
        dlg = self.open_dialog(self.maps[0])
        dlg.update_idletasks()
        texts = [t.get_text() for t in dlg.ax_map.texts]
        self.assertTrue(any("(approx.)" in t for t in texts))
        self.assertEqual(len(dlg.fig.axes), 3)            # image, series, colours

    def test_the_viewers_are_owned_by_the_main_window(self):
        # an owned window cannot sink behind its owner; stacking order itself
        # cannot be seen headlessly
        dlg = self.open_dialog(self.maps[0])
        self.assertEqual(str(dlg.wm_transient()), str(self.root))
        import snapmap
        import snapmap_ui
        cube = snapmap.build([284.0, 285.0, 286.0], 2, 2, 0.0, 1.0, 0.0, 1.0,
                             [((x, y), [1, 2, 3]) for x in range(2)
                              for y in range(2)])
        r = self.maps[0]
        sm = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                      type(r)(name="C 1s", index=9, offset=9,
                                              energy=[286.0, 285.0, 284.0],
                                              counts=[1.0, 2.0, 3.0],
                                              decodable=True, sample=r.sample,
                                              extra={"cube": cube}))
        self.addCleanup(sm.destroy)
        self.assertEqual(str(sm.wm_transient()), str(self.root))

    def _row_of(self, region):
        for iid, (_p, node) in self.ws.node_map.items():
            if node.region is region:
                return iid
        raise AssertionError("row not found")

    def test_enter_opens_a_selected_map_and_is_left_alone_otherwise(self):
        made = []
        orig = self.ws.open_snapmap
        self.ws.open_snapmap = lambda region=None: made.append(region)
        try:
            self.ws.tree.selection_set(self._row_of(self.maps[1]))
            self.ws._on_select()
            self.assertEqual(self.ws._open_selected_map(), "break")
            self.assertEqual(len(made), 1)
            self.ws.tree.selection_set(())
            self.ws._on_select()
            self.assertIsNone(self.ws._open_selected_map())
            self.assertEqual(len(made), 1)
        finally:
            self.ws.open_snapmap = orig

    def test_double_click_on_a_sample_row_of_maps_opens_the_first(self):
        from types import SimpleNamespace
        made = []
        orig = self.ws.open_snapmap
        self.ws.open_snapmap = lambda region=None: made.append(region)
        try:
            sample_rows = [iid for iid, (_p, n) in self.ws.node_map.items()
                           if n.region is None and n.type_name == "sample"]
            self.assertTrue(sample_rows)
            tree = self.ws.tree
            real = tree.identify_row
            tree.identify_row = lambda y: sample_rows[0]
            try:
                out = self.ws._on_tree_double(SimpleNamespace(x=5, y=5))
            finally:
                tree.identify_row = real
            self.assertEqual(out, "break")
            self.assertEqual(len(made), 1)
            under = self.ws.node_map[sample_rows[0]][1]
            self.assertIs(made[0], imaging_first(under))
            # a row that is not a map keeps its normal double-click behaviour
            tree.identify_row = lambda y: ""
            try:
                self.assertIsNone(self.ws._on_tree_double(
                    SimpleNamespace(x=5, y=5)))
            finally:
                tree.identify_row = real
        finally:
            self.ws.open_snapmap = orig

    def test_status_bar_tells_how_to_open_a_selected_map(self):
        self.ws.tree.selection_set(self._row_of(self.maps[0]))
        self.ws._on_select()
        self.assertIn("press Enter to open", self.ws.status.cget("text"))
        self.ws.tree.selection_set(())
        self.ws._on_select()
        self.assertNotIn("press Enter", self.ws.status.cget("text"))

    def test_report_pages_do_not_pick_up_imaging_maps(self):
        self.assertFalse(imagepages.available([self.parser]))
        self.assertEqual(imagepages.items([self.parser]), [])
        self.assertEqual(imagepages.plan([self.parser]), [])


if __name__ == "__main__":
    unittest.main()
