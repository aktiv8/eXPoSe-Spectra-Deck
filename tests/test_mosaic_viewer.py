"""The mosaic in the app's Images tab and in the offline HTML browser: the
sets found, the view in the app (with its points, zoom toolbar and the
refine switch), selecting a sample that only a mosaic shows, and the
stitched picture in the page's payload. Synthetic pictures; the app part
needs a display (skipped otherwise).

Run:  python -m unittest discover tests
"""

import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

try:
    import matplotlib
    matplotlib.use("Agg")
    import numpy  # noqa: F401
    from PIL import Image  # noqa: F401
    HAVE = True
except Exception:                                    # pragma: no cover
    HAVE = False

if HAVE:
    import htmlbrowser as hb
    import imagepages as ip
    import mosaic
    from test_mosaic import TestInThePages
    from test_htmlbrowser import doc, region


def make_doc(xs=(10.0, 11.0, 12.0)):
    helper = TestInThePages()
    return helper.doc(xs)


@unittest.skipUnless(HAVE, "numpy, Pillow and matplotlib needed")
class TestMosaicSets(unittest.TestCase):
    def setUp(self):
        mosaic._memo.clear()

    def test_one_set_per_cluster_with_points_and_sites(self):
        d = make_doc()
        sets = ip.mosaic_sets([d])
        self.assertEqual(len(sets), 1)
        self.assertIs(sets[0].parser, d)
        self.assertEqual(len(sets[0].cluster), 3)
        self.assertEqual(sets[0].title, "Mosaic of 3 pictures: S1")

    def test_far_apart_and_single_pictures_make_none(self):
        self.assertEqual(ip.mosaic_sets([make_doc(xs=(10.0, 40.0))]), [])
        self.assertEqual(ip.mosaic_sets([make_doc(xs=(10.0,))]), [])

    def test_skipped_pictures_are_not_in_a_set(self):
        d = make_doc()
        keys = [k for k, _l in ip.items([d])]
        sets = ip.mosaic_sets([d], skip=keys[:1])
        self.assertEqual(len(sets[0].cluster), 2)

    def test_refine_off_builds_a_picture_by_stage_position_alone(self):
        ms = ip.mosaic_sets([make_doc()])[0]
        m, pic = ip.mosaic_picture(ms.cluster, ms.positions, ms.sites,
                                   refine=False)
        self.assertEqual(m.matched, 0)
        self.assertIn("S1", pic.points)


@unittest.skipUnless(HAVE, "numpy, Pillow and matplotlib needed")
class TestNaming(unittest.TestCase):
    def test_names_that_differ_in_a_short_ending_are_folded(self):
        self.assertEqual(ip.sample_span(["Stained Pt #003a", "Stained Pt #003b",
                                         "Stained Pt #003c"]),
                         "Stained Pt #003a/b/c")
        self.assertEqual(ip.sample_span(["S1", "S1", ""]), "S1")
        self.assertEqual(ip.sample_span([]), "")

    def test_unrelated_names_are_listed_and_cut(self):
        self.assertEqual(ip.sample_span(["Alpha", "Beta"]), "Alpha, Beta")
        self.assertEqual(ip.sample_span(["A one", "B two", "C three", "D four"]),
                         "A one, B two, C three +1 more")


@unittest.skipUnless(HAVE, "numpy, Pillow and matplotlib needed")
class TestBeforeAnalysis(unittest.TestCase):
    def setUp(self):
        mosaic._memo.clear()

    def two_shots(self):
        d = make_doc()
        # photograph the first position again, later (after the analysis)
        again = TestInThePages().blob("S1 Pt #000a  12:00", 10.0)
        d.images = d.images + [again]
        return d

    def test_the_set_uses_the_latest_picture_unless_asked_for_the_earliest(self):
        d = self.two_shots()
        late = ip.mosaic_sets([d])[0]
        early = ip.mosaic_sets([d], keep="earliest")[0]
        self.assertIn("S1 Pt #000a  12:00", [p.name for p in late.cluster])
        self.assertIn("S1 Pt #000a  10:00", [p.name for p in early.cluster])
        self.assertEqual(late.key, early.key)

    def test_the_page_offers_a_before_analysis_mosaic_only_when_it_differs(self):
        names = [c["name"] for c in hb.build_payload([make_doc()])["cameras"]
                 if c.get("mosaic")]
        self.assertEqual(len(names), 1)
        cams = [c for c in hb.build_payload([self.two_shots()])["cameras"]
                if c.get("mosaic")]
        self.assertEqual([c["before"] for c in cams], [False, True])
        self.assertTrue(cams[1]["name"].endswith("before analysis"))
        self.assertIn("S1", cams[0]["name"])


@unittest.skipUnless(HAVE, "numpy, Pillow and matplotlib needed")
class TestMosaicInThePage(unittest.TestCase):
    def setUp(self):
        mosaic._memo.clear()

    def payload(self, **kw):
        return hb.build_payload([make_doc(**kw)])

    def test_the_page_lists_the_mosaic_after_the_pictures(self):
        cams = self.payload()["cameras"]
        self.assertEqual([bool(c.get("mosaic")) for c in cams],
                         [False, False, False, True])
        m = cams[-1]
        self.assertEqual(m["n"], 3)
        self.assertTrue(m["img"].startswith("data:image/jpeg;base64,"))
        self.assertIn("S1", m["points"])
        self.assertIn("pictures placed by stage position", m["desc"])

    def test_points_and_bar_are_in_the_shrunk_pictures_pixels(self):
        m = self.payload()["cameras"][-1]
        for x, y in m["points"].values():
            self.assertTrue(0 <= x <= m["w"] and 0 <= y <= m["h"])
        self.assertGreater(m["bar_px"], 0)
        self.assertLessEqual(max(m["w"], m["h"]), hb.MOSAIC_MAX_PX)
        json.dumps(m)

    def test_no_mosaic_when_the_pictures_do_not_overlap(self):
        cams = self.payload(xs=(10.0, 40.0))["cameras"]
        self.assertFalse(any(c.get("mosaic") for c in cams))

    def test_a_mosaic_over_budget_is_left_out_with_a_note(self):
        old = hb.CAMERA_BUDGET
        hb.CAMERA_BUDGET = 1
        try:
            p = self.payload()
        finally:
            hb.CAMERA_BUDGET = old
        self.assertFalse(any(c.get("mosaic") for c in p["cameras"]))
        self.assertTrue(any("mosaic" in n for n in p["build_notes"]))


try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_APP = ee.HAVE_MPL and HAVE
except Exception:                                    # pragma: no cover
    tk, ee, HAVE_APP = None, None, False


@unittest.skipUnless(HAVE_APP, "matplotlib / Tk not available")
class TestMosaicInTheApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.rc = dict(ee.matplotlib.rcParams)
        cls.ws = ee.Workspace(cls.root)

    @classmethod
    def tearDownClass(cls):
        ee.matplotlib.rcParams.update(cls.rc)
        cls.root.destroy()

    def setUp(self):
        mosaic._memo.clear()
        self.d = make_doc()
        self.ws.docs = [self.d]
        self.ws._cur_image = None
        self.ws._refresh_images()

    def tearDown(self):
        self.ws.docs = []
        self.ws._cur_image = None
        self.ws._refresh_images()

    def drawn_figures(self, show):
        """The figures the viewer puts on a canvas while ``show()`` runs."""
        made, real = [], ee.FigureCanvasTkAgg

        def spy(fig, *a, **k):
            made.append(fig)
            return real(fig, *a, **k)
        ee.FigureCanvasTkAgg = spy
        try:
            show()
        finally:
            ee.FigureCanvasTkAgg = real
        return made

    def test_showing_the_mosaic_draws_its_points_and_enables_refine(self):
        ms = ip.mosaic_sets([self.d])[0]
        figs = self.drawn_figures(lambda: self.ws._show_image(self.d, ms))
        self.assertEqual(str(self.ws.refine_cb.cget("state")), "normal")
        self.assertEqual(len(figs), 1)
        ax = figs[0].axes[0]
        self.assertEqual(len(ax.images), 1)
        self.assertIn("S1", [t.get_text() for t in ax.texts])
        titles = [ax.get_title(loc=loc) for loc in ("left", "center", "right")]
        self.assertTrue(any("Mosaic of 3 pictures" in t for t in titles))

    def test_the_brightness_and_before_choices_are_only_for_mosaics(self):
        self.ws._show_image(self.d, self.d.images[0])
        self.assertEqual(str(self.ws.gain_cb.cget("state")), "disabled")
        self.assertEqual(str(self.ws.keep_cb.cget("state")), "disabled")
        self.ws._show_image(self.d, ip.mosaic_sets([self.d])[0])
        self.assertEqual(str(self.ws.gain_cb.cget("state")), "normal")
        self.assertEqual(str(self.ws.keep_cb.cget("state")), "readonly")

    def test_choosing_before_analysis_redraws_with_the_earlier_pictures(self):
        again = TestInThePages().blob("S1 Pt #000a  12:00", 10.0)
        self.d.images = self.d.images + [again]
        ms = ip.mosaic_sets([self.d])[0]
        self.ws.keep_var.set(self.ws.KEEP_CHOICES[1])
        try:
            self.drawn_figures(lambda: self.ws._show_image(self.d, ms))
        finally:
            self.ws.keep_var.set(self.ws.KEEP_CHOICES[0])
        built = [m for m in mosaic._memo.values()]
        self.assertTrue(any("S1 Pt #000a  10:00" in m.names for m in built))

    def test_the_refine_switch_is_only_for_mosaics(self):
        self.ws._show_image(self.d, self.d.images[0])
        self.assertEqual(str(self.ws.refine_cb.cget("state")), "disabled")
        self.ws._show_image(self.d, ip.mosaic_sets([self.d])[0])
        self.assertEqual(str(self.ws.refine_cb.cget("state")), "normal")

    def test_selecting_a_sample_with_no_picture_of_its_own_shows_the_mosaic(self):
        # S2 has a position inside the pictures but no picture of its own
        self.ws.nb.select(self.ws.tab_images)

        class R:
            sample = "S2"
        self.ws.sel_regions = [R()]
        self.ws._follow_camera_image()
        cur = self.ws._cur_image
        self.assertTrue(cur and isinstance(cur[1], ip.MosaicSet))
        self.ws.sel_regions = []


if __name__ == "__main__":
    unittest.main()
