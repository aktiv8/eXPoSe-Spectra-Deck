"""A SnapMap view as data and as a figure: the spec, the drawing, a picture
file, the saved figure in the report paths, and the viewer reopening it (the
app part needs a display and matplotlib / Tk; skipped otherwise).

Run:  python -m unittest discover tests
"""

import json
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import numpy  # noqa: F401
    import matplotlib
    HAVE_PLOT = True
except Exception:                                   # pragma: no cover
    HAVE_PLOT = False

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

import snapmapfig  # noqa: E402
import themes  # noqa: E402

if HAVE_PLOT:
    from test_htmlbrowser import cube_region, doc, map_cube  # noqa: E402


class TestSpec(unittest.TestCase):
    def test_a_bad_view_becomes_a_valid_one(self):
        s = snapmapfig.sanitise({"window": [1, "x"], "overlay": {"alpha": 9},
                                 "roi": [1, 2, 3], "smooth": 5, "junk": 1})
        self.assertEqual(s["kind"], "snapmap")
        self.assertIsNone(s["window"])
        self.assertEqual(s["rois"], [])
        self.assertEqual(s["overlay"]["alpha"], 1.0)
        self.assertEqual(s["smooth"]["method"], "None")
        self.assertNotIn("junk", s)
        self.assertEqual(snapmapfig.sanitise(None)["scale"], "Viridis")

    def test_it_survives_json_unchanged(self):
        s = snapmapfig.sanitise({"ref": {"file": "f1", "pos": 0},
                                 "window": [284.0, 288.5], "background": True,
                                 "roi": [-10, -10, 5, 5], "scale": "Magma"})
        self.assertEqual(snapmapfig.sanitise(json.loads(json.dumps(s))), s)
        self.assertTrue(snapmapfig.is_view(s))
        self.assertFalse(snapmapfig.is_view({"group_by": "x"}))

    def test_scale_bar_lengths_are_round(self):
        for width, want in ((100.0, 20.0), (250.0, 50.0), (60.0, 10.0),
                            (1000.0, 200.0), (7.0, 1.0)):
            self.assertEqual(snapmapfig.nice_length(width), want)


@unittest.skipUnless(HAVE_PLOT, "numpy / matplotlib not available")
class TestDrawing(unittest.TestCase):
    def draw(self, spec, cube=None):
        from matplotlib.figure import Figure
        cube = cube or map_cube()
        r = cube_region(cube)
        fig = Figure(figsize=(10, 5))
        axes = snapmapfig.draw_view(
            fig, None, spec, r, themes.PRINT,
            themes.scale_colourmap("Viridis", False, themes.PRINT))
        return fig, axes, cube

    def test_map_and_spectrum_are_drawn_with_the_chosen_window(self):
        fig, (ax_map, ax_spec), cube = self.draw({"window": [286.0, 290.0]})
        self.assertIsNotNone(ax_spec)
        want = cube.image(286.0, 290.0)
        got = ax_map.images[0].get_array()
        self.assertTrue((got == want).all())
        self.assertIn("286.0–290.0", ax_map.get_title())
        self.assertEqual(len(fig.axes), 3)               # map, spectrum, colour bar

    def test_map_only_has_no_spectrum_panel(self):
        fig, (ax_map, ax_spec), _ = self.draw({"show_spectrum": False})
        self.assertIsNone(ax_spec)
        self.assertEqual(len(fig.axes), 2)

    def test_the_area_is_outlined_and_gives_its_own_spectrum(self):
        cube = map_cube()
        roi = [cube.x_of(0) - 1, cube.y_of(0) - 1, cube.x_of(1) + 1,
               cube.y_of(1) + 1]
        _fig, (ax_map, ax_spec), _ = self.draw({"roi": roi}, cube)   # legacy key
        self.assertEqual(len(ax_map.patches), 1)
        labels = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        self.assertIn("Area 1, 4 px", labels)
        mask = snapmapfig.roi_mask(cube, roi)
        self.assertEqual(int(mask.sum()), 4)

    def test_several_areas_each_get_a_colour_outline_and_spectrum(self):
        cube = map_cube()
        rois = [snapmapfig.new_roi(snapmapfig.roi_rect_of_pixel(cube, 0, 0)),
                snapmapfig.new_roi([cube.x_of(2) - 1, cube.y_of(1) - 1,
                                    cube.x_of(4) + 1, cube.y_of(3) + 1],
                                   "ellipse", "Bright side")]
        _fig, (ax_map, ax_spec), _ = self.draw({"rois": rois}, cube)
        self.assertEqual(len(ax_map.patches), 2)
        names = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        self.assertEqual(names[0], "whole map")
        self.assertTrue(names[1].startswith("Area 1"))
        self.assertTrue(names[2].startswith("Bright side"))
        from matplotlib.patches import Ellipse
        self.assertTrue(any(isinstance(p, Ellipse) for p in ax_map.patches))
        colours = {l.get_color() for l in ax_spec.lines[1:3]}
        self.assertEqual(len(colours), 2)

    def test_an_ellipse_holds_fewer_pixels_than_its_rectangle(self):
        cube = map_cube(nx=12, ny=12)
        rect = [cube.x_of(0), cube.y_of(0), cube.x_of(11), cube.y_of(11)]
        r = int(snapmapfig.roi_mask(cube, snapmapfig.new_roi(rect)).sum())
        e = int(snapmapfig.roi_mask(cube,
                                    snapmapfig.new_roi(rect, "ellipse")).sum())
        self.assertEqual(r, 144)
        self.assertLess(e, r)
        self.assertGreater(e, 0.6 * r)                 # about pi/4 of it
        tiny = snapmapfig.new_roi(snapmapfig.roi_rect_of_pixel(cube, 3, 3),
                                  "ellipse")
        self.assertEqual(int(snapmapfig.roi_mask(cube, tiny).sum()), 1)

    def test_area_numbers_are_those_of_the_map_on_show(self):
        cube = map_cube()
        roi = snapmapfig.new_roi([cube.x_of(0) - 1, cube.y_of(0) - 1,
                                  cube.x_of(1) + 1, cube.y_of(1) + 1],
                                 label="corner")
        row, = snapmapfig.roi_stats(cube, (284.0, 296.0), False, [roi])
        img = cube.image(284.0, 296.0)
        want = img[:2, :2]
        self.assertEqual(row["px"], 4)
        self.assertAlmostEqual(row["area_um2"], 4 * 10.0 * 10.0)
        self.assertAlmostEqual(row["mean"], float(want.mean()))
        self.assertAlmostEqual(row["sum"], float(want.sum()))
        self.assertAlmostEqual(row["sd"], float(want.std(ddof=1)))
        self.assertAlmostEqual(row["ratio"], float(want.mean() / img.mean()))
        self.assertEqual(row["label"], "corner")

    def test_areas_are_capped_labelled_and_survive_json(self):
        many = [{"rect": [0, 0, i + 1, 1]} for i in range(10)]
        s = snapmapfig.sanitise({"rois": many + [{"rect": "bad"}]})
        self.assertEqual(len(s["rois"]), snapmapfig.MAX_ROIS)
        self.assertEqual(s["rois"][2]["label"], "Area 3")
        self.assertEqual(snapmapfig.sanitise(json.loads(json.dumps(s))), s)
        legacy = snapmapfig.sanitise({"roi": [1, 2, 3, 4]})
        self.assertEqual(legacy["rois"][0]["rect"], [1.0, 2.0, 3.0, 4.0])
        self.assertNotIn("roi", legacy)

    def test_one_clicked_pixel_makes_a_one_pixel_area(self):
        cube = map_cube()
        roi = snapmapfig.roi_rect_of_pixel(cube, 3, 2)
        mask = snapmapfig.roi_mask(cube, roi)
        self.assertEqual(int(mask.sum()), 1)
        self.assertTrue(mask[2, 3])

    def test_the_scale_bar_is_labelled_in_micrometres(self):
        _fig, (ax_map, _a), cube = self.draw({})
        self.assertTrue(any("µm" in t.get_text() for t in ax_map.texts))
        _fig, (ax_map, _a), cube = self.draw({"scalebar": False})
        self.assertFalse(ax_map.texts)

    def test_drawing_never_touches_the_cube(self):
        cube = map_cube()
        before = bytes(cube.data)
        self.draw({"background": True, "smooth": {"method": "None"}}, cube)
        self.assertEqual(bytes(cube.data), before)

    def test_describe_says_what_the_view_shows(self):
        text = snapmapfig.describe({"window": [284, 290], "background": True,
                                    "rois": [{"rect": [0, 0, 1, 1]}, {"rect": [1, 1, 2, 2]}]})
        self.assertIn("284.0-290.0", text)
        self.assertIn("background removed", text)
        self.assertIn("2 areas of interest", text)
        import pptx_export
        self.assertEqual(
            pptx_export.look_notes({"kind": "snapmap", "window": [1, 2]}),
            snapmapfig.describe({"kind": "snapmap", "window": [1, 2]}))


def flat_cube(value, stage, nx=4, ny=3, dx=10.0, ne=12):
    """A map whose every pixel counts ``value`` in every channel."""
    import snapmap
    blocks = [((ix, iy), [value] * ne) for iy in range(ny) for ix in range(nx)]
    c = snapmap.build([float(300 - i) for i in range(ne)], nx, ny,
                      -(nx - 1) * dx / 2, dx, -(ny - 1) * dx / 2, dx, blocks)
    c.stage_x_mm, c.stage_y_mm = stage
    return c


@unittest.skipUnless(HAVE_PLOT, "numpy / matplotlib not available")
class TestMosaic(unittest.TestCase):
    W = (290.0, 300.0)

    def test_a_map_to_the_right_on_the_stage_lies_to_the_left_in_the_image(self):
        # image x runs against stage X (snapshot.X_DIR), as for camera pictures
        a = flat_cube(1.0, (0.0, 0.0))
        b = flat_cube(5.0, (0.1, 0.0))              # 100 um further along stage X
        m = snapmapfig.mosaic_image([a, b], self.W)
        x0, x1, y1, y0 = m.extent
        self.assertAlmostEqual(x0, -100.0 - 20.0)   # b occupies -120..-80
        self.assertAlmostEqual(x1, 20.0)
        left_half = m.image[:, : m.image.shape[1] // 2]
        self.assertTrue(np_all(left_half[~np_isnan(left_half)] > 0))
        self.assertGreater(np_nanmax(m.image[:, 0]), np_nanmax(m.image[:, -1]))

    def test_stage_y_runs_down_the_image(self):
        a = flat_cube(1.0, (0.0, 0.0))
        b = flat_cube(5.0, (0.0, 0.1))              # +Y: further down the image
        m = snapmapfig.mosaic_image([a, b], self.W)
        self.assertGreater(np_nanmax(m.image[-1]), np_nanmax(m.image[0]))

    def test_overlap_is_averaged_and_gaps_are_blank(self):
        a = flat_cube(2.0, (0.0, 0.0))
        b = flat_cube(4.0, (0.02, 0.0))             # 20 um over: 2 px overlap
        m = snapmapfig.mosaic_image([a, b], self.W)
        vals = sorted({round(v, 6) for v in m.image[~np_isnan(m.image)]})
        self.assertEqual(vals, sorted({2.0 * 11, 3.0 * 11, 4.0 * 11}))
        self.assertTrue(any("overlap" in n for n in m.notes))
        c = flat_cube(1.0, (0.0, 0.0))
        d = flat_cube(1.0, (0.0, 0.2))              # well apart
        gap = snapmapfig.mosaic_image([c, d], self.W)
        self.assertTrue(np_isnan(gap.image).any())
        self.assertTrue(any("blank" in n for n in gap.notes))

    def test_maps_of_different_pixel_size_share_the_finest_grid(self):
        a = flat_cube(1.0, (0.0, 0.0), dx=10.0)
        b = flat_cube(1.0, (0.0, 0.0), nx=2, ny=2, dx=20.0)
        m = snapmapfig.mosaic_image([a, b], self.W)
        self.assertAlmostEqual((m.extent[1] - m.extent[0]) / m.image.shape[1],
                               10.0)

    def test_a_map_without_a_stage_position_is_left_out_with_a_note(self):
        a = flat_cube(1.0, (0.0, 0.0))
        b = flat_cube(1.0, (0.0, 0.0))
        b.stage_x_mm = None
        m = snapmapfig.mosaic_image([a, b], self.W)
        self.assertEqual(len(m.used), 1)
        self.assertTrue(any("stage position" in n for n in m.notes))

    def test_mosaic_drawing_uses_one_window_and_one_range_for_all(self):
        from matplotlib.figure import Figure
        a, b = flat_cube(1.0, (0.0, 0.0)), flat_cube(9.0, (0.1, 0.0))
        ra, rb = cube_region(a), cube_region(b)
        fig = Figure(figsize=(10, 5))
        ax_map, ax_spec = snapmapfig.draw_view(
            fig, None, {"window": list(self.W), "scalebar": False}, ra,
            themes.PRINT, themes.scale_colourmap("Viridis", False, themes.PRINT),
            tiles=[rb])
        im = ax_map.images[0]
        lo, hi = im.get_clim()
        self.assertLess(lo, 12.0)
        self.assertGreater(hi, 12.0)
        self.assertEqual(len(ax_map.patches), 2)       # a footprint per map
        self.assertTrue(snapmapfig.sanitise({"tiles": [{"file": "f"}]})["tiles"])


def np_all(a):
    import numpy
    return bool(numpy.all(a))


def np_isnan(a):
    import numpy
    return numpy.isnan(a)


def np_nanmax(a):
    import numpy
    return float(numpy.nanmax(a))


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestFigureInTheApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.tmp = tempfile.TemporaryDirectory()
        cls.region = cube_region(map_cube(), "C 1s", "S1")
        cls.parser = doc("m.vms", [cls.region])
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.parser
        cls.ws = ee.Workspace(cls.root)
        cls.ws._add_file("m.vms", file_id="m", refresh=False)
        cls.ws._finish_adding()

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()
        cls.tmp.cleanup()

    def setUp(self):
        self.ws.figures.clear()

    def dialog(self, **kw):
        import snapmap_ui
        dlg = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                       self.region, **kw)
        self.addCleanup(dlg.destroy)
        return dlg

    def test_the_viewer_saves_a_picture_file(self):
        dlg = self.dialog()
        for ext in (".png", ".pdf", ".svg"):
            path = os.path.join(self.tmp.name, "map" + ext)
            self.assertTrue(dlg.save_picture(path))
            self.assertGreater(os.path.getsize(path), 1000)
        with open(os.path.join(self.tmp.name, "map.png"), "rb") as fh:
            self.assertEqual(fh.read(8), b"\x89PNG\r\n\x1a\n")

    def test_saved_as_a_figure_it_holds_the_view_and_the_map_reference(self):
        dlg = self.dialog()
        dlg.window = (285.0, 289.0)
        dlg.background.set(True)
        fig = dlg.save_as_figure("C 1s map")
        st = fig["state"]
        self.assertTrue(snapmapfig.is_view(st))
        self.assertEqual(st["window"], [285.0, 289.0])
        self.assertTrue(st["background"])
        self.assertEqual(st["ref"]["name"], "C 1s")
        self.assertEqual(st["ref"]["file"], "m")
        self.assertEqual(self.ws.figures, [fig])
        json.dumps(self.ws.figures)                      # storable in the workbook

    def test_every_report_path_renders_the_figure(self):
        dlg = self.dialog()
        dlg.save_as_figure("C 1s map")
        fig = self.ws.figures[0]
        fig["caption"] = "Carbon distribution."
        pngs = self.ws._figure_pngs(fig, dpi=60)         # HTML browser, hand-over
        self.assertEqual(len(pngs), 1)
        self.assertEqual(pngs[0][:4], b"\x89PNG")
        slide = self.ws._deck_images(1, fig)             # slides and Word
        self.assertEqual(len(slide), 1)
        pages = []
        self.ws._render_figure_pages(fig, pages.append, number=1, dpi=60)
        self.assertEqual(len(pages), 1)
        # the report's own figure list includes it
        self.assertIn(fig, self.ws._report_figures())

    def test_a_figure_whose_map_is_gone_gives_no_page(self):
        dlg = self.dialog()
        fig = dlg.save_as_figure("gone")
        fig["state"]["ref"]["file"] = "other"
        self.assertEqual(self.ws._figure_pngs(fig, dpi=50), [])

    def test_reopening_restores_the_view_and_updating_changes_the_figure(self):
        dlg = self.dialog()
        dlg.window = (285.0, 289.0)
        cube = self.region.extra["cube"]
        dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 1, 1))
        dlg.scale.set("Magma")
        dlg.save_as_figure("again")
        saved = self.ws.figures[0]["state"]
        again = self.dialog(view=saved, figure=0)
        self.assertEqual(again.window, (285.0, 289.0))
        self.assertEqual(again.scale.get(), "Magma")
        self.assertEqual(len(again.rois), 1)
        self.assertEqual(int(again.masks()[0].sum()), 1)
        again.window = (286.0, 288.0)
        self.assertTrue(again.update_figure())
        self.assertEqual(self.ws.figures[0]["state"]["window"], [286.0, 288.0])

    def test_recall_opens_the_viewer_on_the_saved_view(self):
        dlg = self.dialog()
        dlg.window = (285.0, 289.0)
        dlg.save_as_figure("recall me")
        made = []
        real = ee.snapmap_ui.SnapMapDialog

        class Recorded(real):
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                made.append(self)

        ee.snapmap_ui.SnapMapDialog = Recorded
        try:
            self.assertTrue(self.ws.open_snapmap_view(
                self.ws.figures[0]["state"], figure=0))
        finally:
            ee.snapmap_ui.SnapMapDialog = real
        self.addCleanup(made[0].destroy)
        self.assertEqual(made[0].window, (285.0, 289.0))
        self.assertEqual(made[0].figure, 0)
        self.assertFalse(self.ws.open_snapmap_view({"kind": "snapmap",
                                                    "ref": {"file": "x"}}))

    def test_a_mosaic_is_offered_saved_and_reopened(self):
        other = cube_region(map_cube(), "C 1s", "S2")
        other.extra["cube"].stage_x_mm = 10.05       # 50 um along from S1
        other.extra["cube"].stage_y_mm = 20.0
        extra = doc("n.vms", [other])
        ee.load_file = lambda p: {"m.vms": self.parser, "n.vms": extra}[
            os.path.basename(p)]
        self.ws._add_file("n.vms", file_id="n", refresh=False)
        self.ws._finish_adding()
        self.addCleanup(lambda: (self.ws.docs.remove(extra),
                                 self.ws.file_ids.pop(id(extra), None)))
        cands = self.ws.snapmap_candidates(self.parser, self.region)
        self.assertEqual([r for _p, r, _d in cands], [other])
        self.assertAlmostEqual(cands[0][2], 0.05)
        dlg = self.dialog()
        dlg.choose_tiles([(extra, other)])
        wide = dlg._extent()
        self.assertGreater(wide[1] - wide[0],
                           self.region.extra["cube"].extent()[1]
                           - self.region.extra["cube"].extent()[0])
        fig = dlg.save_as_figure("two maps")
        self.assertEqual(len(fig["state"]["tiles"]), 1)
        self.assertEqual(len(self.ws._figure_pngs(fig, dpi=50)), 1)
        again = self.dialog(view=fig["state"], figure=0)
        self.assertEqual(len(again.tiles), 1)
        dlg.choose_tiles([])
        self.assertEqual(dlg.view_spec()["tiles"], [])

    def test_the_window_is_resizable_with_a_toolbar_and_an_area_panel(self):
        dlg = self.dialog()
        dlg.update_idletasks()
        self.assertEqual(dlg.resizable(), (True, True))
        self.assertTrue(hasattr(dlg, "toolbar"))
        # a zoom tool on: the selectors stand aside; off: they act again
        dlg.toolbar.mode = "zoom rect"
        dlg._sync_tools()
        self.assertFalse(dlg.rect.active)
        self.assertFalse(dlg.span.active)
        dlg.toolbar.mode = ""
        dlg._sync_tools()
        self.assertTrue(dlg.rect.active)

    def test_areas_are_added_listed_kept_through_a_zoom_and_removed(self):
        dlg = self.dialog()
        cube = self.region.extra["cube"]
        dlg.ax_map.set_xlim(-10, 10)                   # zoomed in
        a = dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 0, 0))
        b = dlg.add_roi([cube.x_of(2), cube.y_of(1), cube.x_of(4), cube.y_of(3)],
                        "ellipse", "Edge")
        self.assertEqual((a, b), (0, 1))
        self.assertEqual(dlg.ax_map.get_xlim(), (-10, 10))   # zoom kept
        self.assertEqual(len(dlg.area_tree.get_children()), 2)
        self.assertEqual(dlg.area_tree.item("1", "values")[0], "Edge")
        rows = dlg.stats()
        self.assertEqual([r["label"] for r in rows], ["Area 1", "Edge"])
        csv_text = dlg.stats_csv()
        self.assertIn('"Edge",ellipse', csv_text)
        self.assertEqual(len(csv_text.strip().splitlines()), 4)
        head, srows = dlg.spectra_csv_rows()
        self.assertEqual(len(head), 4)                 # energy, map, two areas
        self.assertEqual(len(srows), len(cube.energy))
        dlg.rename_roi(0, "Spot")
        self.assertEqual(dlg.view_spec()["rois"][0]["label"], "Spot")
        dlg.remove_roi(0)
        self.assertEqual([r["label"] for r in dlg.rois], ["Edge"])
        self.assertEqual(len(dlg.area_tree.get_children()), 1)
        self.assertIsNone(dlg.add_roi([1e6, 1e6, 1e6 + 1, 1e6 + 1]))   # off the map

    def test_changing_the_window_updates_the_numbers(self):
        dlg = self.dialog()
        cube = self.region.extra["cube"]
        dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 1, 1))
        before = dlg.stats()[0]["mean"]
        dlg.window = (300.0, 301.0)
        dlg._update_map()
        shown = float(dlg.area_tree.item("0", "values")[3])
        self.assertNotAlmostEqual(before, dlg.stats()[0]["mean"])
        self.assertAlmostEqual(shown, dlg.stats()[0]["mean"], delta=abs(shown) * 1e-3)

    def test_a_saved_figure_keeps_all_its_areas(self):
        dlg = self.dialog()
        cube = self.region.extra["cube"]
        dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 0, 0))
        dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 2, 2), "ellipse", "B")
        fig = dlg.save_as_figure("areas")
        self.assertEqual(len(fig["state"]["rois"]), 2)
        self.assertEqual(len(self.ws._figure_pngs(fig, dpi=50)), 1)
        again = self.dialog(view=fig["state"])
        self.assertEqual([r["label"] for r in again.rois], ["Area 1", "B"])

    def test_spectrum_figures_are_untouched(self):
        fig = self.ws.add_figure("plain")
        self.assertFalse(snapmapfig.is_view(fig["state"]))
        self.assertIn("group_by", fig["state"])


if __name__ == "__main__":
    unittest.main()
