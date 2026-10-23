"""Manual panel order in the real window (needs a display and matplotlib;
skipped otherwise).

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
class TestPanelOrderInApp(unittest.TestCase):
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
        ws.panel_order = []
        ws.trace_start = 0
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        self.keys = [k for k, _ in ws._groups()]

    def kinds(self):
        self.ws._render()
        return {self.ws._axmap[a]: v[1] for a, v in self.ws._axinfo.items()}

    def order(self):
        return [k for k, _ in self.ws._groups()]

    def test_move_and_reset(self):
        a, b = self.keys
        self.ws.move_panel(b, "first")
        self.assertEqual(self.order(), [b, a])
        self.ws.move_panel(b, "later")
        self.assertEqual(self.order(), [a, b])
        self.ws.move_panel(a, "last")
        self.assertEqual(self.order(), [b, a])
        self.ws.reset_panel_order()
        self.assertEqual(self.order(), [a, b])

    def test_state_round_trip_and_old_state(self):
        a, b = self.keys
        self.ws.move_panel(b, "first")
        st = self.ws.capture_state()
        self.assertEqual(st["panel_order"], [b, a])
        self.ws.panel_order = []
        self.ws.apply_state(st, render=False)
        self.assertEqual(self.order(), [b, a])
        st.pop("panel_order")
        self.ws.apply_state(st, render=False)
        self.assertEqual(self.order(), [a, b])

    def test_drag_title_onto_another_panel(self):
        a, b = self.keys
        ws = self.ws
        ws._render()
        axes = {k: ax for ax, k in ws._axmap.items()}
        bb_a = axes[a].get_window_extent()
        bb_b = axes[b].get_window_extent()

        class Ev:
            inaxes = None
            button = 1
        press, drop = Ev(), Ev()
        press.x, press.y = (bb_b.x0 + bb_b.x1) / 2, bb_b.y1 + 8
        drop.x, drop.y = (bb_a.x0 + bb_a.x1) / 2, (bb_a.y0 + bb_a.y1) / 2
        ws._on_plot_click(press)
        self.assertIsNotNone(ws._drag)
        ws._on_plot_release(drop)
        self.assertEqual(self.order(), [b, a])


class TestOrderPure(unittest.TestCase):
    def test_order_groups(self):
        import panelview
        g = [("a", 1), ("b", 2), ("c", 3)]
        self.assertEqual([k for k, _ in panelview.order_groups(g, ["c", "x"])],
                         ["c", "a", "b"])
        self.assertEqual(panelview.moved(["a", "b", "c"], "c", "earlier"),
                         ["a", "c", "b"])
        self.assertEqual(panelview.moved_to(["a", "b", "c"], "c", "a"),
                         ["c", "a", "b"])
        self.assertEqual(panelview.sanitise_order(["a", "a", 3]), ["a"])


if __name__ == "__main__":
    unittest.main()
