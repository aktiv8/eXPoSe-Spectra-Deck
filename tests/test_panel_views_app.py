"""The per-panel views in the real window (needs a display and matplotlib;
skipped otherwise): panels draw in their own view, old figures carry none,
saved looks round-trip.

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
class TestPanelViewsInApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # a Workspace applies its theme to matplotlib's global settings;
        # give them back so other tests see the defaults
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {
            "depth.vms": doc("depth.vms", [reg("C 1s", "Depth", 286 - i / 5,
                                               level=i + 1)
                                           for i in range(5)]),
            "bulk.vms": doc("bulk.vms", [reg("O 1s", "Bulk", 288 - i / 3)
                                         for i in range(3)])}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        # the tests assume the program's own defaults, not the settings the
        # person running them last saved (a saved fit_merge or colour scale
        # made three of them fail)
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
        ws.z_var.set("Auto")
        ws.scale_var.set("Binding")
        ws.panels_var.set("Auto")
        ws.traces_var.set("All")
        ws.panel_views = {}
        ws.trace_start = 0
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        self.keys = [k for k, _ in ws._groups()]

    def kinds(self):
        self.ws._render()
        return {self.ws._axmap[a]: v[1] for a, v in self.ws._axinfo.items()}

    def test_default_is_the_page_view_everywhere(self):
        self.assertEqual(set(self.kinds().values()), {"stack"})

    def test_panels_can_differ_on_one_page(self):
        depth, bulk = self.keys
        self.ws.panel_views = {depth: {"view": "Waterfall 3D"},
                               bulk: {"view": "Heatmap"}}
        k = self.kinds()
        self.assertEqual((k[depth], k[bulk]), ("waterfall 3d", "heatmap"))

    def test_page_default_still_applies_to_the_rest(self):
        depth, bulk = self.keys
        self.ws.view_var.set("Heatmap")
        self.ws.panel_views = {depth: {"view": "Stack"}}
        k = self.kinds()
        self.assertEqual((k[depth], k[bulk]), ("stack", "heatmap"))

    def test_series_panel_is_z_sorted_and_stack_can_reverse(self):
        depth, bulk = self.keys
        self.ws.panel_views = {depth: {"view": "Heatmap",
                                       "z_axis": "Etch level"},
                               bulk: {"reverse": True}}
        groups = dict(self.ws._groups())
        self.assertEqual([r.etch_level for r in groups[depth]],
                         [1, 2, 3, 4, 5])
        self.ws.panel_views = {}
        base = [id(r) for r in dict(self.ws._groups())[bulk]]
        self.ws.panel_views = {bulk: {"reverse": True}}
        self.assertEqual([id(r) for r in dict(self.ws._groups())[bulk]],
                         base[::-1])

    def test_right_click_hint_shows_when_something_is_plotted(self):
        self.ws._render()
        self.assertIn("Right-click a panel", self.ws.hint.cget("text"))
        self.ws.norm_var.set("At cursor")
        self.ws._render()
        self.assertIn("Click a panel to set the energy",
                      self.ws.hint.cget("text"))
        self.ws.norm_var.set("None")
        self.ws.checked = set()
        self.ws._render()
        self.assertEqual(self.ws.hint.cget("text"), "")
        self.ws.checked = {id(r) for d in self.ws.docs for r in d.regions}

    def test_fit_panel_shows_one_trace_and_moves_the_slider(self):
        depth, _bulk = self.keys
        self.ws.panel_views = {depth: {"view": "Fit"}}
        self.ws._render()
        self.assertEqual(self.ws._trace_limit_now, 1)
        self.assertEqual(int(float(self.ws.trace_scale.cget("to"))), 4)
        self.assertTrue(any("no fit stored" in n
                            for n in self.ws._view_notes))

    def test_zoom_is_kept_when_nothing_about_a_panel_changed(self):
        depth, _bulk = self.keys
        self.ws._render()
        ax = next(a for a, k in self.ws._axmap.items() if k == depth)
        ax.set_xlim(280.0, 284.0)
        ax.set_ylim(200.0, 800.0)
        self.ws._render()          # nothing changed: same ticks, same look
        ax2 = next(a for a, k in self.ws._axmap.items() if k == depth)
        self.assertEqual(ax2.get_xlim(), (280.0, 284.0))
        self.assertEqual(ax2.get_ylim(), (200.0, 800.0))

    def test_zoom_is_dropped_when_the_panels_own_data_changes(self):
        depth, _bulk = self.keys
        self.ws._render()
        ax = next(a for a, k in self.ws._axmap.items() if k == depth)
        ax.set_xlim(280.0, 284.0)
        ax.set_ylim(200.0, 800.0)
        self.ws.norm_var.set("Max = 1")    # the panel's own values change
        self.ws._render()
        ax2 = next(a for a, k in self.ws._axmap.items() if k == depth)
        self.assertNotEqual(ax2.get_ylim(), (200.0, 800.0))

    def test_old_look_without_panel_views_clears_them(self):
        depth, _bulk = self.keys
        st = self.ws.capture_state()
        old ={k: v for k, v in st.items() if k != "panel_views"}
        self.ws.panel_views = {depth: {"view": "Fit"}}
        with self.ws._temp_state(old):
            self.assertEqual(self.ws.panel_views, {})
        self.assertEqual(self.ws.panel_views, {depth: {"view": "Fit"}})

    def test_look_round_trips_and_is_sanitised(self):
        depth, bulk = self.keys
        self.ws.panel_views = {depth: {"view": "Waterfall 3D"}}
        st = self.ws.capture_state()
        self.assertEqual(st["panel_views"], self.ws.panel_views)
        self.ws.panel_views = {}
        self.ws.apply_state(st, render=False)
        self.assertEqual(self.ws.panel_views, {depth: {"view": "Waterfall 3D"}})
        self.ws.apply_state({"panel_views": {bulk: {"view": "nope"},
                                             depth: {"norm": "Max = 1"}}},
                            render=False)
        self.assertEqual(self.ws.panel_views, {depth: {"norm": "Max = 1"}})

    def test_report_figure_pages_follow_the_reports_page_size(self):
        import pdfstyle

        class FakePdf:
            def __init__(self):
                self.pages = []

            def savefig(self, page):
                self.pages.append(page)
        fig = {"name": "F", "caption": "", "state": self.ws.capture_state()}
        self.ws._pdf_figure_size = pdfstyle.page_size("a4")[1]
        pdf = FakePdf()
        self.ws._report_figure_pages(pdf, 1, fig)
        w, h = pdf.pages[0].get_size_inches()
        self.assertAlmostEqual(w, 11.69, places=1)
        self.assertAlmostEqual(h, 8.27, places=1)
        self.ws._pdf_figure_size = pdfstyle.page_size("letter")[1]
        pdf = FakePdf()
        self.ws._report_figure_pages(pdf, 1, fig)
        self.assertEqual(pdf.pages[0].get_size_inches().tolist(), [11.0, 8.5])

    def test_figure_pages_are_drawn_per_panel(self):
        depth, bulk = self.keys
        self.ws.panel_views = {depth: {"view": "Waterfall 3D"}}
        pages = []
        self.ws._render_figure_pages(
            {"name": "F", "caption": "", "state": self.ws.capture_state()},
            pages.append)
        self.assertTrue(pages)
        three_d = [ax for pg in pages for ax in pg.axes
                   if getattr(ax, "name", "") == "3d"]
        self.assertEqual(len(three_d), 1)
        flat = [ax for pg in pages for ax in pg.axes
                if getattr(ax, "name", "") != "3d"]
        self.assertTrue(flat)

    def test_menu_edits_one_panel_only(self):
        depth, bulk = self.keys
        self.ws._set_panel_view(depth, "view", "Heatmap")
        self.ws._set_panel_view(bulk, "norm", "Max = 1")
        self.ws._set_panel_view(bulk, "fit_show", ("envelope", False))
        self.assertEqual(self.ws.panel_views,
                         {depth: {"view": "Heatmap"},
                          bulk: {"norm": "Max = 1",
                                 "fit_show": {"envelope": False}}})
        self.ws.reset_panel_views(depth)
        self.assertEqual(list(self.ws.panel_views), [bulk])
        self.ws.reset_panel_views()
        self.assertEqual(self.ws.panel_views, {})

    def test_fit_merge_is_a_per_panel_choice(self):
        depth, bulk = self.keys
        self.ws._set_panel_view(depth, "fit_merge", True)
        self.assertEqual(self.ws.panel_views, {depth: {"fit_merge": True}})
        self.assertTrue(self.ws._look_for(depth)["fit_merge"])
        self.assertFalse(self.ws._look_for(bulk)["fit_merge"])

    def test_fit_merge_round_trips_and_a_missing_key_is_off(self):
        """CLAUDE.md: a saved figure applies its own state, so one saved
        before the choice existed must not inherit the live window's."""
        self.addCleanup(self.ws.fit_merge_var.set, False)
        self.ws.fit_merge_var.set(True)
        st = self.ws.capture_state()
        self.assertIs(st["fit_merge"], True)
        self.ws.fit_merge_var.set(False)
        self.ws.apply_state(st, render=False)
        self.assertTrue(self.ws.fit_merge_var.get())
        old = {k: v for k, v in st.items() if k != "fit_merge"}
        self.ws.apply_state(old, render=False)
        self.assertFalse(self.ws.fit_merge_var.get())
        # drawing a saved figure leaves the live choice as it was
        self.ws.fit_merge_var.set(True)
        with self.ws._temp_state(old):
            self.assertFalse(self.ws.fit_merge_var.get())
        self.assertTrue(self.ws.fit_merge_var.get())

    def test_residual_sd_round_trips_and_a_missing_key_is_off(self):
        """The Fit row's "Residual StDev" is part of the saved look: a figure
        saved before it existed must not inherit the live window's choice."""
        self.addCleanup(self.ws.fit_sd_var.set, False)
        self.ws.fit_sd_var.set(True)
        st = self.ws.capture_state()
        self.assertIs(st["fit_sd"], True)
        self.ws.fit_sd_var.set(False)
        self.ws.apply_state(st, render=False)
        self.assertTrue(self.ws.fit_sd_var.get())
        old = {k: v for k, v in st.items() if k != "fit_sd"}
        self.ws.apply_state(old, render=False)
        self.assertFalse(self.ws.fit_sd_var.get())
        self.ws.fit_sd_var.set(True)
        with self.ws._temp_state(old):
            self.assertFalse(self.ws.fit_sd_var.get())
        self.assertTrue(self.ws.fit_sd_var.get())

    def test_residual_sd_is_part_of_the_signature_and_reaches_the_results(self):
        self.addCleanup(self.ws.fit_sd_var.set, False)
        before = self.ws._signature()
        self.ws.fit_sd_var.set(True)
        self.assertNotEqual(before, self.ws._signature())
        self.assertIn("sd", self.ws._results().stats)     # reports follow it
        self.ws.fit_sd_var.set(False)
        self.assertNotIn("sd", self.ws._results().stats)

    def test_rms_and_chi2_round_trip_and_a_missing_key_is_on(self):
        """The RMS and Reduced chi-square boxes are on unless switched off,
        and a figure saved before they existed (which drew both) applies them
        on, whatever the live window has."""
        self.addCleanup(self.ws.fit_rms_var.set, True)
        self.addCleanup(self.ws.fit_chi2_var.set, True)
        self.ws.fit_rms_var.set(False)
        self.ws.fit_chi2_var.set(False)
        st = self.ws.capture_state()
        self.assertIs(st["fit_rms"], False)
        self.assertIs(st["fit_chi2"], False)
        self.ws.fit_rms_var.set(True)
        self.ws.fit_chi2_var.set(True)
        self.ws.apply_state(st, render=False)
        self.assertFalse(self.ws.fit_rms_var.get())
        self.assertFalse(self.ws.fit_chi2_var.get())
        old = {k: v for k, v in st.items() if k not in ("fit_rms", "fit_chi2")}
        self.ws.apply_state(old, render=False)
        self.assertTrue(self.ws.fit_rms_var.get())
        self.assertTrue(self.ws.fit_chi2_var.get())
        self.ws.fit_rms_var.set(False)
        with self.ws._temp_state(old):
            self.assertTrue(self.ws.fit_rms_var.get())
        self.assertFalse(self.ws.fit_rms_var.get())

    def test_the_three_statistics_reach_the_signature_and_the_results(self):
        for var in (self.ws.fit_rms_var, self.ws.fit_chi2_var):
            self.addCleanup(var.set, True)
        self.assertEqual(self.ws.fit_stats(), ("rms", "chi2"))
        before = self.ws._signature()
        self.ws.fit_rms_var.set(False)
        self.assertNotEqual(before, self.ws._signature())
        self.assertEqual(self.ws._results().stats, ("chi2",))
        self.ws.fit_chi2_var.set(False)
        self.assertEqual(self.ws._results().stats, ())

    def test_fit_merge_is_part_of_the_page_look_and_the_signature(self):
        self.addCleanup(self.ws.fit_merge_var.set, False)
        before = self.ws._signature()
        self.ws.fit_merge_var.set(True)
        self.assertTrue(self.ws._page_defaults()["fit_merge"])
        self.assertNotEqual(before, self.ws._signature())

    def test_menu_edits_ident_show_per_panel(self):
        depth, _bulk = self.keys
        self.ws._set_panel_view(depth, "ident_show", ("secondary", True))
        self.assertEqual(self.ws.panel_views,
                         {depth: {"ident_show": {"secondary": True}}})
        self.ws.reset_panel_views()
        self.assertEqual(self.ws.panel_views, {})

    def test_ident_show_round_trips_through_state(self):
        self.addCleanup(self.ws.ident_vars["secondary"].set, False)
        self.ws.ident_vars["secondary"].set(True)
        st = self.ws.capture_state()
        self.assertEqual(st["ident_show"],
                         {"secondary": True, "auger": False})
        self.ws.ident_vars["secondary"].set(False)
        self.ws.apply_state(st, render=False)
        self.assertTrue(self.ws.ident_vars["secondary"].get())

    def test_secondary_and_auger_markers_drawn_in_a_different_colour(self):
        # Isolate this from the real (and changeable) bundled table: a
        # small synthetic one with a known secondary line and an Auger
        # line landing on the same displayed BE as the primary marker.
        depth, _bulk = self.keys
        synthetic = [
            {"el": "C", "line": "1s", "be": 286.0, "rank": 1},
            {"el": "C", "line": "2s", "be": 287.0, "rank": 2},
            {"el": "Z", "line": "KLL", "ke": 1200.6, "rank": 3},
        ]
        orig_lines = self.ws._xps_lines
        self.ws._xps_lines = synthetic
        self.addCleanup(setattr, self.ws, "_xps_lines", orig_lines)
        self.addCleanup(self.ws.ann.markers.clear)
        self.addCleanup(self.ws.ident_vars["secondary"].set, False)
        self.addCleanup(self.ws.ident_vars["auger"].set, False)

        r = dict(self.ws._groups())[depth][0]
        self.ws.identify_add(r, 286.0, "C 1s")
        self.ws.ident_vars["secondary"].set(True)
        self.ws.ident_vars["auger"].set(True)
        self.ws._render()
        ax = next(a for a, k in self.ws._axmap.items() if k == depth)
        colours = {t.get_text(): t.get_color() for t in ax.texts}
        pal = self.ws.palette
        self.assertEqual(colours.get("C 1s"), pal["muted"])
        self.assertEqual(colours.get("C 2s"), pal["accent"])
        self.assertEqual(colours.get("Z KLL"), pal["cycle"][1])


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestFitStackInApp(unittest.TestCase):
    """The Fit stack view in the real window: ten fitted levels of one core
    level are shown eight at a time, each with its own fit."""

    @classmethod
    def setUpClass(cls):
        import copy
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        from test_casafit import fitted_region
        regs = []
        for i in range(10):
            r = copy.copy(fitted_region(n=61))
            r.counts = [c * (1 + 0.1 * i) for c in r.counts]
            r.sample, r.etch_level, r.etch_time = "Depth", i + 1, 30.0 * i
            regs.append(r)
        cls.docs = {"fits.vms": doc("fits.vms", regs)}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        # the tests assume the program's own defaults, not the settings the
        # person running them last saved (a saved fit_merge or colour scale
        # made three of them fail)
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
        ws.view_var.set("Fit stack")
        ws.norm_var.set("None")
        ws.scale_var.set("Binding")
        ws.panels_var.set("Auto")
        ws.traces_var.set("All")
        ws.colscale_var.set("Theme default")
        ws.panel_views = {}
        ws.trace_start = 0
        ws.fit_merge_var.set(False)
        for v in ws.fit_vars.values():
            v.set(True)
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        self.key = [k for k, _ in ws._groups()][0]

    def axis(self):
        self.ws._render()
        return next(a for a, k in self.ws._axmap.items() if k == self.key)

    def test_the_window_is_at_most_eight_traces(self):
        lim = self.ws._panel_limit
        look = {"view": "Fit stack"}
        self.assertEqual(lim(None, look), 8)
        self.assertEqual(lim(20, look), 8)
        self.assertEqual(lim(5, look), 5)
        self.assertEqual(lim(None, {"view": "Fit"}), 1)
        self.assertIsNone(lim(None, {"view": "Stack"}))
        self.assertEqual(lim(5, {"view": "Stack"}), 5)
        ax = self.axis()
        self.assertEqual(self.ws._axinfo[ax][2], 8)
        self.assertEqual(self.ws._axinfo[ax][1], "stack")
        self.assertIn("1–8 of 10", ax.get_title(loc="right"))

    def test_a_report_steps_the_window_until_every_trace_is_shown(self):
        ws = self.ws
        windows = ws._fit_stack_windows(ws._groups(), None, 0)
        self.assertEqual([first for _p, first in windows], [0, 8])
        self.assertEqual(len(windows[1][0]), 1)        # only the panel that goes on
        shown = []
        pages = ws._render_figure_pages(
            {"state": ws.capture_state(), "name": "Depth"},
            lambda page: (shown.append(page), page)[1])
        self.assertEqual(len(pages), 2)
        titles = [t.get_text() for p in shown for t in p.texts]
        head = ee.reportspec.figure_title(1, "Depth")
        self.assertIn(head, titles)
        self.assertIn(head + " (continued)", titles)
        notes = [a.get_title(loc="right") for p in shown for a in p.axes]
        self.assertTrue(any("1–8 of 10" in n for n in notes))
        self.assertTrue(any("3–10 of 10" in n for n in notes))   # a full last window

    def test_a_group_that_fits_one_window_gets_one_page(self):
        ws = self.ws
        ws.traces_var.set("All")
        ws.checked = {id(r) for d in ws.docs for r in d.regions[:5]}
        windows = ws._fit_stack_windows(ws._groups(), None, 0)
        self.assertEqual(len(windows), 1)

    def test_the_trace_slider_walks_the_window(self):
        self.ws.trace_start = 2
        ax = self.axis()
        self.assertEqual(self.ws._axinfo[ax][2], 8)
        self.assertIn("3–10 of 10", ax.get_title(loc="right"))

    def test_every_trace_is_drawn_with_its_fit_and_one_legend(self):
        ax = self.axis()
        from matplotlib.legend import Legend
        # 8 data + 8 x (background, components, envelope) + the scale bar
        self.assertGreater(len(ax.lines), 8 * 3)
        self.assertEqual(len([c for c in ax.get_children()
                              if isinstance(c, Legend)]), 1)
        self.assertGreater(len(ax.collections), 8)      # shaded components

    def test_the_fit_layers_and_merge_apply_to_the_stack(self):
        n_all = len(self.axis().lines)
        self.ws.fit_vars["background"].set(False)
        self.assertEqual(len(self.axis().lines), n_all - 8)
        self.ws.fit_vars["background"].set(True)
        self.ws.fit_vars["envelope"].set(False)
        self.ws.fit_vars["components"].set(False)
        self.assertEqual(len(self.axis().lines), n_all - 8 - 2 * 8)

    def test_the_plain_stack_view_still_draws_no_fit(self):
        self.ws.view_var.set("Stack")
        ax = self.axis()
        self.assertEqual(len(ax.collections), 0)
        self.assertEqual(len(ax.lines), 10 + 1)          # data and scale bar

    def test_the_data_keep_one_colour_so_they_do_not_clash_with_states(self):
        ax = self.axis()
        cols = {ax.lines[i].get_color() for i in range(8)}
        self.assertEqual(cols, {self.ws.palette["cycle"][0]})
        self.ws.colscale_var.set("Viridis")                 # the user's choice
        cols = {self.axis().lines[i].get_color() for i in range(8)}
        self.assertGreater(len(cols), 1)

    def test_a_trace_without_a_fit_is_noted_and_drawn_plain(self):
        reg = next(iter(self.docs.values())).regions[3]
        saved, reg.fit = reg.fit, None
        self.addCleanup(setattr, reg, "fit", saved)
        self.axis()
        self.assertIn(f"{self.key}: 1 of 8 spectra have no fit",
                      self.ws._view_notes)

    def test_a_group_with_no_fits_says_so(self):
        regs = next(iter(self.docs.values())).regions
        saved = [r.fit for r in regs]
        for r in regs:
            r.fit = None
        try:
            self.axis()
            self.assertIn(f"{self.key}: no fit stored, shown as a plain stack",
                          self.ws._view_notes)
        finally:
            for r, f in zip(regs, saved):
                r.fit = f

    def test_fit_overlay_gives_one_curve_set_per_trace_and_says_things_once(self):
        regs = [self.ws._display(r) for r in
                next(iter(self.docs.values())).regions[:4]]
        notes = []
        out = self.ws._fit_overlay(regs, regs, self.ws.palette, notes,
                                   {"components": True}, False, stack=True)
        self.assertEqual(len(out["per_trace"]), 4)
        self.assertTrue(all(out["per_trace"]))
        self.assertEqual(len(notes), len(set(notes)))
        single = self.ws._fit_overlay(regs[:1], regs[:1], self.ws.palette, [],
                                      {"components": True}, False)
        self.assertNotIn("per_trace", single)
        self.assertIsNone(self.ws._fit_overlay(regs, regs, self.ws.palette,
                                               [], {"components": True}))

    def test_a_panel_can_have_its_own_fit_stack(self):
        self.ws.view_var.set("Stack")
        self.ws._set_panel_view(self.key, "view", "Fit stack")
        ax = self.axis()
        self.assertEqual(self.ws._axinfo[ax][2], 8)
        self.assertGreater(len(ax.collections), 8)

    def test_the_view_survives_a_saved_state_and_draws_in_a_report(self):
        st = self.ws.capture_state()
        self.assertEqual(st["view_mode"], "Fit stack")
        self.ws.view_var.set("Stack")
        pages = []
        self.ws._render_figure_pages({"name": "F", "caption": "", "state": st},
                                     pages.append)
        self.assertTrue(pages)
        ax = next(a for pg in pages for a in pg.axes if a.collections)
        self.assertGreater(len(ax.collections), 8)
        self.assertEqual(self.ws.view_var.get(), "Stack")   # live view kept


if __name__ == "__main__":
    unittest.main()
