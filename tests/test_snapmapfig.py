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

    def test_a_smoothed_spectrum_is_never_drawn_unlabelled(self):
        import smoothing
        method = next(m for m in smoothing.METHODS if m != "None")
        cube = map_cube()
        roi = snapmapfig.new_roi(snapmapfig.roi_rect_of_pixel(cube, 1, 1))
        _f, (_m, ax_spec), _ = self.draw(
            {"rois": [roi], "smooth": {"method": method, "strength": 0.5}}, cube)
        names = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        self.assertEqual(names[0], f"whole map ({method})")
        self.assertTrue(names[1].endswith(f"px, {method}"), names)
        _f, (_m, ax_spec), _ = self.draw({"rois": [roi]}, cube)
        self.assertEqual(ax_spec.get_legend().get_texts()[0].get_text(),
                         "whole map")

    def test_area_outlines_have_a_dark_halo_for_contrast(self):
        cube = map_cube()
        roi = snapmapfig.new_roi(snapmapfig.roi_rect_of_pixel(cube, 1, 1))
        _f, (ax_map, _s), _ = self.draw({"rois": [roi]}, cube)
        self.assertTrue(ax_map.patches[0].get_path_effects())

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

    def test_a_reference_without_a_stage_position_gives_no_mosaic(self):
        a, b = flat_cube(1.0, (0.0, 0.0)), flat_cube(2.0, (0.1, 0.0))
        a.stage_x_mm = None
        self.assertIsNone(snapmapfig.mosaic_image([a, b], self.W))

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


def title_of(ax):
    """An axes' title whichever side the active style puts it."""
    return ax.get_title(loc="left") or ax.get_title()


def element_cubes():
    """Three elements of one map file on one pixel grid: C 1s bright on the
    left, O 1s bright on the right, N 1s weak and flat."""
    import snapmap
    nx, ny, ne = 6, 4, 30
    out = []
    for name, centre, left in (("C 1s", 12, 1.0), ("O 1s", 18, -1.0),
                               ("N 1s", 15, 0.0)):
        blocks = []
        for iy in range(ny):
            for ix in range(nx):
                f = (1.0 - 0.25 * left * (ix - 2.5)) if left else 0.4
                y = [20 + 80 * f * 2.718 ** (-((c - centre) / 2.0) ** 2)
                     for c in range(ne)]
                blocks.append(((ix, iy), y))
        c = snapmap.build([300.0 - 0.5 * i for i in range(ne)], nx, ny,
                          -25.0, 10.0, -15.0, 10.0, blocks)
        c.stage_x_mm, c.stage_y_mm = 10.0, 20.0
        out.append((name, c))
    return out


@unittest.skipUnless(HAVE_PLOT, "numpy / matplotlib not available")
class TestLayerHelpers(unittest.TestCase):
    def test_a_ratio_leaves_blank_where_the_denominator_is_noise(self):
        import numpy as np
        num = np.array([[100.0, 400.0], [6.0, 8.0]])
        den = np.array([[50.0, 100.0], [0.0, 100.0]])
        out = snapmapfig.ratio_image(num, den)
        self.assertEqual(out[0, 0], 2.0)
        self.assertEqual(out[0, 1], 4.0)
        self.assertTrue(np.isnan(out[1, 0]))                 # zero
        self.assertAlmostEqual(out[1, 1], 0.08)
        tiny = snapmapfig.ratio_image(num, np.array([[100.0, 100.0],
                                                     [100.0, 1.0]]))
        self.assertTrue(np.isnan(tiny[1, 1]))                # under 5 % of p99

    def test_rgb_channels_are_scaled_on_their_own_range(self):
        import numpy as np
        a = np.linspace(0, 1000, 12).reshape(3, 4)
        b = np.linspace(0, 1, 12).reshape(3, 4)
        rgb = snapmapfig.rgb_image([a, b])
        self.assertEqual(rgb.shape, (3, 4, 3))
        self.assertAlmostEqual(float(rgb[..., 0].max()), 1.0, places=6)
        self.assertAlmostEqual(float(rgb[..., 1].max()), 1.0, places=6)
        self.assertEqual(float(rgb[..., 2].max()), 0.0)      # no blue given

    def test_only_maps_on_the_same_pixels_combine(self):
        (_, a), (_, b), _ = element_cubes()
        self.assertTrue(snapmapfig.same_grid(a, b))
        b.dx = 20.0
        self.assertFalse(snapmapfig.same_grid(a, b))

    def test_modes_and_layers_are_sanitised(self):
        s = snapmapfig.sanitise({"mode": "bogus"})
        self.assertEqual((s["mode"], s["layers"]), ("single", []))
        s = snapmapfig.sanitise({"mode": "ratio", "layers": [
            {"ref": {"file": "f"}, "window": [1, 2]}, {"ref": 3}, {"x": 1},
            {"ref": {"file": "g"}, "window": "no"}]})
        self.assertEqual(s["mode"], "ratio")
        self.assertEqual([l["window"] for l in s["layers"]], [[1.0, 2.0], None])
        self.assertEqual(snapmapfig.sanitise(json.loads(json.dumps(s))), s)
        self.assertIn("ratio map of 2 elements", snapmapfig.describe(s))


@unittest.skipUnless(HAVE_PLOT, "numpy / matplotlib not available")
class TestLayerDrawing(unittest.TestCase):
    def setUp(self):
        self.pairs = element_cubes()
        self.regions = [cube_region(c, n, "S1") for n, c in self.pairs]
        self.layers = [(r, None) for r in self.regions]

    def draw(self, mode, layers=None, **spec):
        from matplotlib.figure import Figure
        fig = Figure(figsize=(10, 5))
        axes = snapmapfig.draw_view(
            fig, None, dict(spec, mode=mode), self.regions[0], themes.PRINT,
            themes.scale_colourmap("Viridis", False, themes.PRINT),
            layers=self.layers if layers is None else layers)
        return fig, axes

    def test_a_grid_has_one_map_and_one_colour_bar_per_element(self):
        fig, (ax, spec_ax) = self.draw("grid")
        self.assertIsNone(spec_ax)
        titles = [title_of(a) for a in fig.axes if title_of(a)]
        self.assertEqual(len(titles), 3)
        self.assertTrue(titles[0].startswith("C 1s"))
        self.assertEqual(len(fig.axes), 6)           # three maps, three colour bars

    def test_a_composite_says_which_element_is_which_colour(self):
        fig, (ax, _s) = self.draw("composite")
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        self.assertEqual(len(labels), 3)
        self.assertTrue(labels[0].startswith("C 1s"))
        img = ax.images[0].get_array()
        self.assertEqual(img.shape[-1], 3)
        # C is brighter on the left, O on the right: red wins left, green right
        self.assertGreater(img[0, 0, 0], img[0, 0, 1])
        self.assertGreater(img[0, -1, 1], img[0, -1, 0])

    def test_a_ratio_map_is_the_first_over_the_second(self):
        import numpy as np
        fig, (ax, _s) = self.draw("ratio")
        want = snapmapfig.ratio_image(
            self.pairs[0][1].image(*snapmapfig.snapmap.default_window(
                self.pairs[0][1])),
            self.pairs[1][1].image(*snapmapfig.snapmap.default_window(
                self.pairs[1][1])))
        got = np.asarray(ax.images[0].get_array(), dtype=float)
        self.assertTrue(np.allclose(got, want, equal_nan=True))
        self.assertEqual(title_of(ax), "C 1s / O 1s")
        self.assertGreater(got[0, 0], got[0, -1])    # C/O falls to the right

    def test_areas_are_outlined_on_the_new_views_too(self):
        cube = self.pairs[0][1]
        rois = [snapmapfig.new_roi(snapmapfig.roi_rect_of_pixel(cube, 1, 1))]
        for mode in ("grid", "composite", "ratio"):
            fig, (ax, _s) = self.draw(mode, rois=rois)
            self.assertEqual(len(ax.patches), 1, mode)

    def test_too_few_elements_fall_back_to_the_single_map(self):
        fig, (ax, spec_ax) = self.draw("composite", layers=self.layers[:1])
        self.assertIsNotNone(spec_ax)                # the map and spectrum view

    def test_elements_on_other_pixels_are_left_out(self):
        self.pairs[1][1].dx = 25.0
        fig, (ax, spec_ax) = self.draw("ratio", layers=self.layers[:2])
        self.assertIsNotNone(spec_ax)                # only one usable: single map


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


def si_view(**over):
    """A fit-components view of the two-state Si 2p test map."""
    from test_snapfit import two_state_map
    cube, truth = two_state_map()
    region = cube_region(cube, "Si 2p", "S1")
    region.energy = list(truth["be"])
    fit = {"spec": truth["spec"]}
    fit.update(over.pop("fit", {}))
    return region, dict({"mode": "fit", "fit": fit}, **over), truth


@unittest.skipUnless(HAVE_PLOT, "numpy / matplotlib not available")
class TestFitView(unittest.TestCase):
    def draw(self, spec, region):
        from matplotlib.figure import Figure
        fig = Figure(figsize=(12, 5))
        axes = snapmapfig.draw_view(
            fig, None, spec, region, themes.PRINT,
            themes.scale_colourmap("Viridis", False, themes.PRINT))
        return fig, axes

    def test_the_fit_part_is_sanitised_and_survives_json(self):
        s = snapmapfig.sanitise({"mode": "fit", "fit": {
            "spec": {"components": [{"name": "Si", "be": 99.4},
                                    {"name": "SiO2", "be": 103.4}]},
            "show": ["fraction", "junk", "spectrum"], "on": -2}})
        f = s["fit"]
        self.assertEqual(s["mode"], "fit")
        self.assertEqual(f["show"], ["fraction", "spectrum"])
        self.assertEqual(f["numerator"], ["SiO2"])      # the last of two or more
        self.assertIsNone(f["on"])
        self.assertEqual(snapmapfig.sanitise(json.loads(json.dumps(s))), s)
        one = snapmapfig.sanitise({"fit": {"spec": {"components": [
            {"be": 100.0}]}}})["fit"]
        self.assertEqual(one["numerator"], [])
        self.assertEqual(snapmapfig.sanitise({"fit": {"spec": {"components": [
            {"name": "A", "be": 1.0}, {"name": "B", "be": 2.0}]},
            "numerator": ["B", "zzz"]}})["fit"]["numerator"], ["B"])
        self.assertIn("fit of Si, SiO2", snapmapfig.describe(
            {"mode": "fit", "fit": {"spec": {"components": [
                {"name": "Si", "be": 99.4}, {"name": "SiO2", "be": 103.4}]}}}))

    def test_component_fraction_residual_maps_and_the_spectrum(self):
        region, spec, truth = si_view()
        fig, (ax, ax_spec) = self.draw(spec, region)
        titles = [title_of(a) for a in fig.axes if title_of(a)]
        self.assertIn("Si", titles)
        self.assertIn("SiO2", titles)
        self.assertIn("SiO2 fraction", titles)
        self.assertTrue(any("\u03c7" in t for t in titles))   # residual map
        legend = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        self.assertEqual(legend[0], "data")
        self.assertTrue(legend[1].startswith("Si  99.40 eV"))
        self.assertTrue(legend[2].startswith("SiO2  103.40 eV"))
        self.assertEqual(legend[-1], "fit")
        # the fraction map is the oxide share, drawn from 0 to 1
        import numpy as np
        frac_ax = next(a for a in fig.axes if title_of(a) == "SiO2 fraction")
        img = frac_ax.images[0]
        self.assertEqual(img.get_clim(), (0.0, 1.0))
        got = np.asarray(img.get_array(), dtype=float)
        self.assertGreater(np.nanmean(got[:, -5:]), np.nanmean(got[:, :5]) + 0.4)

    def test_the_percentages_in_the_legend_add_up_to_100(self):
        region, spec, _t = si_view()
        _fig, (_ax, ax_spec) = self.draw(spec, region)
        legend = [t.get_text() for t in ax_spec.get_legend().get_texts()]
        pcts = [int(x.split()[-2]) for x in legend[1:3]]
        self.assertEqual(sum(pcts), 100)

    def test_panels_can_be_chosen(self):
        region, spec, _t = si_view(fit={"show": ["fraction"]})
        fig, (ax, ax_spec) = self.draw(spec, region)
        self.assertIsNone(ax_spec)
        self.assertEqual([title_of(a) for a in fig.axes if title_of(a)],
                         ["SiO2 fraction"])
        region, spec, _t = si_view(fit={"show": ["spectrum"]})
        fig, (ax, ax_spec) = self.draw(spec, region)
        self.assertIsNotNone(ax_spec)
        self.assertEqual(len(fig.axes), 1)

    def test_the_spectrum_can_be_that_of_an_area(self):
        region, spec, _t = si_view()
        cube = region.extra["cube"]
        spec["rois"] = [snapmapfig.new_roi(
            [cube.x_of(0), cube.y_of(0), cube.x_of(5), cube.y_of(29)],
            label="left edge")]
        spec["fit"]["on"] = 0
        fig, (ax, ax_spec) = self.draw(spec, region)
        self.assertIn("left edge", title_of(ax_spec))
        # the area is outlined on the maps too
        self.assertTrue(any(a.patches for a in fig.axes if title_of(a) == "Si"))

    def test_without_a_binding_axis_or_components_it_says_so(self):
        region, spec, _t = si_view()
        region.energy_label = "Kinetic Energy"
        fig, _ = self.draw(spec, region)
        self.assertTrue(any("binding-energy" in t.get_text()
                            for a in fig.axes for t in a.texts))
        region, spec, _t = si_view(fit={"spec": {"components": []}})
        fig, _ = self.draw(spec, region)
        self.assertTrue(any("No components" in t.get_text()
                            for a in fig.axes for t in a.texts))

    def test_the_energy_shift_of_the_shown_region_is_followed(self):
        # the user shifted the spectrum by +0.7 eV: the fit looks there
        import numpy as np
        region, spec, truth = si_view()
        region.energy = [e + 0.7 for e in region.energy]
        fig, (ax, ax_spec) = self.draw(spec, region)
        x = np.asarray(ax_spec.lines[0].get_xdata())    # the data, on the shown axis
        self.assertAlmostEqual(float(x.min()), 95.7, delta=0.05)
        self.assertLess(float(x.max()), 110.0)          # still cut to the window

    def test_drawing_does_not_touch_the_cube(self):
        region, spec, _t = si_view()
        before = bytes(region.extra["cube"].data)
        self.draw(spec, region)
        self.assertEqual(bytes(region.extra["cube"].data), before)


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestFitFigureInTheApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.region, cls.spec, cls.truth = si_view()
        cls.parser = doc("si.vms", [cls.region])
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.parser
        cls.ws = ee.Workspace(cls.root)
        cls.ws._add_file("si.vms", file_id="si", refresh=False)
        cls.ws._finish_adding()

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()

    def test_a_fit_view_is_a_picture_and_a_report_figure(self):
        import snapmap_ui
        dlg = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                       self.region)
        self.addCleanup(dlg.destroy)
        view = dict(dlg.view_spec(), mode="fit",
                    fit=snapmapfig.sanitise_fit(self.spec["fit"]))
        path = os.path.join(tempfile.gettempdir(), "snapmap_fit.png")
        self.assertTrue(self.ws.snapmap_picture(view, path, size=(12, 5)))
        self.addCleanup(os.remove, path)
        self.assertGreater(os.path.getsize(path), 5000)
        fig = self.ws.add_figure("Si fit", state=view)
        self.ws.figures.append(fig) if fig not in self.ws.figures else None
        self.assertEqual(len(self.ws._figure_pngs(fig, dpi=50)), 1)
        self.assertEqual(len(self.ws._deck_images(1, fig)), 1)
        import pptx_export
        self.assertIn("fit of Si, SiO2", pptx_export.look_notes(fig["state"]))
        self.ws.figures.remove(fig)


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestFitInTheDialog(unittest.TestCase):
    """The Fit tab: seeding, refinement within +/-0.5 eV, components, the
    fraction, areas, saving and reopening (the Si 2p test map, shown 0.3 eV
    below where the library puts its states)."""

    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls._load = ee.load_file
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()

    def setUp(self):
        from test_snapfit import two_state_map
        cube, truth = two_state_map()
        self.region = cube_region(cube, "Si 2p", "S1")
        self.region.energy = [e - 0.3 for e in truth["be"]]     # charging / calibration
        self.truth = truth
        self.parser = doc("si.vms", [self.region])
        ee.load_file = lambda p: self.parser
        for d in list(self.ws.docs):
            self.ws.docs.remove(d)
        self.ws.file_ids.clear()
        self.ws.region_parser.clear()
        self.ws.figures.clear()
        self.ws._add_file("si.vms", file_id="si", refresh=False)
        self.ws._finish_adding()
        import snapmap_ui
        self.dlg = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                            self.region)
        self.addCleanup(self.dlg.destroy)

    def fit_mode(self):
        self.dlg.set_mode("fit")
        return self.dlg

    def test_the_fit_view_seeds_si_and_sio2_and_shows_its_tab(self):
        d = self.fit_mode()
        self.assertEqual(d.effective_mode(), "fit")
        comps = d.fit_spec["components"]
        self.assertEqual([(c["name"], c["be"]) for c in comps],
                         [("Si", 99.4), ("SiO2", 103.4)])
        self.assertTrue(all(c["doublet"] for c in comps))
        self.assertEqual(str(d.side_nb.tab(d.fit_tab, "state")), "normal")
        self.assertEqual(d.side_nb.select(), str(d.fit_tab))
        lo, hi = float(d.v_fwin_lo.get()), float(d.v_fwin_hi.get())
        self.assertLess(lo, 98.9)
        self.assertGreater(hi, 104.5)
        titles = {title_of(a) for a in d.fig.axes}
        self.assertTrue({"Si", "SiO2", "SiO2 fraction"} <= titles)
        d.set_mode("single")
        self.assertEqual(str(d.side_nb.tab(d.fit_tab, "state")), "hidden")

    def test_the_refinement_follows_the_shift_and_reports_what_it_did(self):
        d = self.fit_mode()
        res = d.fit_spectrum()
        self.assertIsNotNone(res)
        got = {c["name"]: c["be"] for c in d.fit_spec["components"]}
        self.assertAlmostEqual(got["Si"], 99.1, delta=0.08)       # the seed was 0.3 off
        self.assertAlmostEqual(got["SiO2"], 103.1, delta=0.08)
        rows = d.fit_results_rows()
        self.assertEqual([r[0] for r in rows], ["Si", "SiO2"])
        self.assertAlmostEqual(sum(r[-1] for r in rows), 100.0, places=6)
        self.assertIsNotNone(rows[0][2])                          # a BE uncertainty
        self.assertIn("χ²", d.fit_note.cget("text"))
        self.assertEqual(len(d.fit_res.get_children()), 2)
        self.assertEqual(d.fit_tree.item("0", "values")[0], "Si")

    def test_a_shift_beyond_the_allowed_room_is_said(self):
        d = self.fit_mode()
        d.fit_edit(0, be_range=0.1)
        d.fit_edit(1, be_range=0.1)
        d.fit_spectrum()
        self.assertIn("limit", d.fit_note.cget("text"))

    def test_the_fraction_and_the_areas_follow_the_fit(self):
        d = self.fit_mode()
        d.fit_spectrum()
        cube = self.region.extra["cube"]
        self.assertEqual(d.fit_view()["numerator"], None)
        d.add_roi([cube.x_of(0), cube.y_of(2), cube.x_of(5), cube.y_of(27)],
                  label="Si side")
        d.add_roi([cube.x_of(34), cube.y_of(2), cube.x_of(39), cube.y_of(27)],
                  label="oxide side")
        rows = {r["label"]: r for r in d.stats()}
        self.assertLess(rows["Si side"]["mean"], 0.3)
        self.assertGreater(rows["oxide side"]["mean"], 0.7)
        self.assertEqual(len(d.area_tree.get_children()), 2)
        self.assertIn("fraction", d.area_note.cget("text"))
        # take the fraction of Si instead: it flips
        d._numer_vars["SiO2"].set(False)
        d._numer_vars["Si"].set(True)
        d._numer_changed()
        self.assertEqual(d.fit_numer, ["Si"])
        rows = {r["label"]: r for r in d.stats()}
        self.assertGreater(rows["Si side"]["mean"], 0.7)

    def test_the_spectrum_can_be_fitted_on_an_area(self):
        d = self.fit_mode()
        cube = self.region.extra["cube"]
        d.add_roi([cube.x_of(0), cube.y_of(2), cube.x_of(5), cube.y_of(27)],
                  label="Si side")
        whole = [r[-1] for r in d.fit_results_rows()]
        self.assertIn("Si side", d.fon_cb.cget("values"))
        d.v_fon.set("Si side")
        d._fit_on_changed()
        self.assertEqual(d.fit_on, 0)
        area = [r[-1] for r in d.fit_results_rows()]
        self.assertGreater(area[0], whole[0] + 15)       # more elemental Si there
        self.assertIsNotNone(d.fit_spectrum())            # refines on that area
        d.remove_roi(0)
        self.assertIsNone(d.fit_on)

    def test_components_are_added_edited_and_removed(self):
        d = self.fit_mode()
        i = d.fit_add_custom("Sub", 101.4, 1.2)
        self.assertEqual(i, 2)
        d.fit_edit(2, be_range=0.3, free_fwhm=False, doublet=None,
                   name="Suboxide")
        c = d.fit_spec["components"][2]
        self.assertEqual((c["name"], c["be_range"], c["free_fwhm"], c["doublet"]),
                         ("Suboxide", 0.3, False, None))
        self.assertEqual(len(d.fit_tree.get_children()), 3)
        d.fit_remove(2)
        self.assertEqual([c["name"] for c in d.fit_spec["components"]],
                         ["Si", "SiO2"])
        d.fit_remove(1)
        self.assertEqual(d.fit_view()["numerator"], None)     # back to default
        j = d.fit_add_library("SiO2")
        self.assertEqual(d.fit_spec["components"][j]["be"], 103.4)
        for k in range(10):
            d.fit_add_custom(f"X{k}")
        self.assertEqual(len(d.fit_spec["components"]), snapfit_max())

    def test_options_and_panels_reach_the_view_and_the_picture(self):
        d = self.fit_mode()
        d.v_fbg.set("linear")
        d.v_fbin.set("2")
        d.v_fshared.set(True)
        d.v_fshow["residual"].set(False)
        d._fit_changed()
        v = d.view_spec()["fit"]
        self.assertEqual((v["spec"]["background"], v["spec"]["bin"],
                          v["spec"]["shared_fwhm"]), ("linear", 2, True))
        self.assertNotIn("residual", v["show"])
        titles = {title_of(a) for a in d.fig.axes}
        self.assertFalse(any("\u03c7" in t and "whole" not in t for t in titles))
        path = os.path.join(tempfile.gettempdir(), "snapmap_dlg_fit.png")
        self.assertTrue(d.save_picture(path))
        self.addCleanup(os.remove, path)

    def test_a_saved_fit_figure_comes_back_with_its_fitted_components(self):
        d = self.fit_mode()
        d.fit_spectrum()
        d.v_fshow["components"].set(False)
        d._fit_changed()
        fig = d.save_as_figure("Si partition")
        st = fig["state"]
        self.assertEqual(st["mode"], "fit")
        self.assertAlmostEqual(st["fit"]["spec"]["components"][0]["be"], 99.1,
                               delta=0.08)
        self.assertEqual(len(self.ws._figure_pngs(fig, dpi=50)), 1)
        import snapmap_ui
        again = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                         self.region, view=st, figure=0)
        self.addCleanup(again.destroy)
        self.assertEqual(again.effective_mode(), "fit")
        self.assertAlmostEqual(again.fit_spec["components"][0]["be"],
                               d.fit_spec["components"][0]["be"], places=6)
        self.assertFalse(again.v_fshow["components"].get())
        self.assertEqual(again.view_spec()["fit"]["show"],
                         ["fraction", "residual", "spectrum"])

    def test_a_mosaic_ends_the_fit_view_and_the_level_name_is_found(self):
        d = self.fit_mode()
        other = cube_region(map_cube(), "Si 2p", "S2")
        d.choose_tiles([(self.parser, other)])         # this map has no position
        self.assertEqual(d.mode.get(), "fit")
        self.assertIn("stage position", d.hint.cget("text"))
        cube = self.region.extra["cube"]
        cube.stage_x_mm, cube.stage_y_mm = 10.0, 20.0
        other.extra["cube"].stage_x_mm, other.extra["cube"].stage_y_mm = 10.05, 20.0
        d.choose_tiles([(self.parser, other)])
        self.assertEqual(d.mode.get(), "single")
        self.assertEqual(d.effective_mode(), "single")
        d.choose_tiles([])
        for name, want in (("Si 2p", "Si 2p"), ("Cu 2p3/2", "Cu 2p"),
                           ("C 1s", "C 1s"), ("Survey", "Survey")):
            self.region.name = name
            self.assertEqual(d._level_name(), want)
        self.region.name = "Si 2p"

    def test_without_a_binding_axis_the_tab_says_so(self):
        self.region.energy_label = "Kinetic Energy"
        d = self.fit_mode()
        self.assertEqual(d.fit_spec["components"], [])
        self.assertIn("binding-energy", d.hint.cget("text"))
        self.assertIsNone(d.fit_spectrum())
        self.region.energy_label = "Binding Energy"


def snapfit_max():
    import snapfit
    return snapfit.MAX_COMPONENTS


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestModesInTheApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.regions = [cube_region(c, n, "S1") for n, c in element_cubes()]
        cls.parser = doc("e.vms", cls.regions)
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.parser
        cls.ws = ee.Workspace(cls.root)
        cls.ws._add_file("e.vms", file_id="e", refresh=False)
        cls.ws._finish_adding()

    @classmethod
    def tearDownClass(cls):
        ee.load_file = cls._load
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()

    def setUp(self):
        self.ws.figures.clear()

    def dialog(self, **kw):
        import snapmap_ui
        dlg = snapmap_ui.SnapMapDialog(self.root, self.ws, self.parser,
                                       self.regions[0], **kw)
        self.addCleanup(dlg.destroy)
        return dlg

    def test_each_mode_draws_and_the_rows_follow_the_mode(self):
        dlg = self.dialog()
        self.assertEqual(dlg.effective_mode(), "single")
        dlg.set_mode("composite")
        self.assertEqual(dlg.effective_mode(), "composite")
        self.assertEqual([r[0].get() for r in dlg.layer_rows],
                         ["C 1s", "O 1s", "N 1s"])
        self.assertIsNone(dlg.ax_spec)               # no spectrum panel
        dlg.set_mode("ratio")
        self.assertEqual(dlg._row_count(), 2)
        self.assertEqual(title_of(dlg.ax_map), "C 1s / O 1s")
        dlg.set_mode("grid")
        self.assertEqual(len([a for a in dlg.fig.axes if title_of(a)]), 3)
        dlg.set_mode("single")
        self.assertIsNotNone(dlg.ax_spec)

    def test_a_ratio_gets_its_own_area_numbers_and_a_composite_none(self):
        import numpy as np
        dlg = self.dialog()
        dlg.set_mode("ratio")
        cube = self.regions[0].extra["cube"]
        dlg.add_roi(snapmapfig.roi_rect_of_pixel(cube, 0, 0))
        row, = dlg.stats()
        img = dlg._image()
        self.assertAlmostEqual(row["mean"], float(img[0, 0]))
        self.assertEqual(len(dlg.area_tree.get_children()), 1)
        dlg.set_mode("composite")
        self.assertEqual(dlg.stats(), [])
        self.assertEqual(len(dlg.area_tree.get_children()), 0)
        self.assertIn("single map and the ratio map", dlg.area_note.cget("text"))

    def test_changing_an_element_and_its_window_redraws(self):
        dlg = self.dialog()
        dlg.set_mode("ratio")
        dlg.layer_rows[1][0].set("N 1s")
        dlg._layer_changed(1)
        self.assertEqual(title_of(dlg.ax_map), "C 1s / N 1s")
        before = dlg._image().copy()
        dlg.layer_rows[0][1].set("296.0")
        dlg.layer_rows[0][2].set("298.0")
        dlg._layer_changed()
        import numpy as np
        self.assertFalse(np.allclose(dlg._image(), before, equal_nan=True))
        spec = dlg.view_spec()
        self.assertEqual(spec["layers"][0]["window"], [296.0, 298.0])

    def test_a_mode_is_saved_in_the_figure_and_comes_back(self):
        dlg = self.dialog()
        dlg.set_mode("composite")
        dlg.layer_rows[2][0].set("—")                 # only red and green
        dlg._layer_changed()
        fig = dlg.save_as_figure("rg")
        st = fig["state"]
        self.assertEqual(st["mode"], "composite")
        self.assertEqual(len(st["layers"]), 2)
        self.assertEqual(len(self.ws._figure_pngs(fig, dpi=50)), 1)
        for mode in ("grid", "ratio"):
            dlg.set_mode(mode)
            f = dlg.save_as_figure(mode)
            self.assertEqual(len(self.ws._deck_images(1, f)), 1, mode)
        again = self.dialog(view=st)
        self.assertEqual(again.effective_mode(), "composite")
        self.assertEqual(again.layer_rows[1][0].get(), "O 1s")
        path = os.path.join(tempfile.gettempdir(), "snapmap_modes.png")
        self.assertTrue(again.save_picture(path, spectrum=False))
        self.addCleanup(os.remove, path)

    def test_a_mosaic_and_a_mode_exclude_each_other(self):
        dlg = self.dialog()
        dlg.set_mode("grid")
        dlg.choose_tiles([(self.parser, self.regions[1])])
        self.assertEqual(dlg.mode.get(), "single")
        self.assertEqual(dlg.view_spec()["layers"], [])
        dlg.choose_tiles([])
        dlg.set_mode("grid")
        self.assertEqual(dlg.effective_mode(), "grid")
        dlg.choose_tiles([(self.parser, self.regions[1])])
        self.assertEqual(dlg.effective_mode(), "single")


@unittest.skipUnless(HAVE_PLOT and os.environ.get("XPS_SNAPMAP_CORPUS"),
                     "set XPS_SNAPMAP_CORPUS to the Vishay experiment folder")
class TestRealOverlappingMaps(unittest.TestCase):
    """Points 002a and 002b of the Vishay experiment are 1.446 mm apart on the
    stage and their 2 mm maps overlap by about 0.55 mm: the stage-position
    placement must line the overlap up (checked by correlating it)."""

    def test_the_stage_placement_lines_the_overlap_up(self):
        import numpy as np
        import snapmap
        from readers import load_file
        parser = load_file(os.environ["XPS_SNAPMAP_CORPUS"])

        def region(sample, name):
            return next(r for r in parser.regions
                        if r.sample.startswith(sample) and r.name == name
                        and r.extra.get("cube") is not None)

        def ncc(x, y):
            x, y = x - x.mean(), y - y.mean()
            return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()))

        for name in ("Al 2p", "O 1s", "F 1s"):
            ca = region("Point 002a", name).extra["cube"]
            cb = region("Point 002b", name).extra["cube"]
            win = snapmap.default_window(ca)
            a, b = ca.image(*win), cb.image(*win)
            _ox, oy = snapmapfig.tile_offset(cb, ca)
            s = int(round(oy / abs(ca.dy)))                  # 72 px
            here = ncc(a[s:], b[:ca.ny - s])
            flipped = ncc(a[:ca.ny - s], b[s:])              # the opposite sign
            self.assertGreater(here, 0.5, name)
            self.assertLess(abs(flipped), 0.2, name)
            best = max((ncc(a[k:], b[:ca.ny - k]), k) for k in range(60, 85))
            self.assertLessEqual(abs(best[1] - oy / abs(ca.dy)), 1.5, name)


if __name__ == "__main__":
    unittest.main()
