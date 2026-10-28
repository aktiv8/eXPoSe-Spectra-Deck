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
        self.assertIsNone(s["roi"])
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
        _fig, (ax_map, ax_spec), _ = self.draw({"roi": roi}, cube)
        self.assertEqual(len(ax_map.patches), 1)
        labels = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        self.assertIn("area, 4 px", labels)
        mask = snapmapfig.roi_mask(cube, roi)
        self.assertEqual(int(mask.sum()), 4)

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
                                    "roi": [0, 0, 1, 1]})
        self.assertIn("284.0-290.0", text)
        self.assertIn("background removed", text)
        self.assertIn("area of interest", text)
        import pptx_export
        self.assertEqual(
            pptx_export.look_notes({"kind": "snapmap", "window": [1, 2]}),
            snapmapfig.describe({"kind": "snapmap", "window": [1, 2]}))


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
        dlg.roi = snapmapfig.roi_rect_of_pixel(cube, 1, 1)
        dlg.mask = snapmapfig.roi_mask(cube, dlg.roi)
        dlg.scale.set("Magma")
        dlg.save_as_figure("again")
        saved = self.ws.figures[0]["state"]
        again = self.dialog(view=saved, figure=0)
        self.assertEqual(again.window, (285.0, 289.0))
        self.assertEqual(again.scale.get(), "Magma")
        self.assertEqual(int(again.mask.sum()), 1)
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

    def test_spectrum_figures_are_untouched(self):
        fig = self.ws.add_figure("plain")
        self.assertFalse(snapmapfig.is_view(fig["state"]))
        self.assertIn("group_by", fig["state"])


if __name__ == "__main__":
    unittest.main()
