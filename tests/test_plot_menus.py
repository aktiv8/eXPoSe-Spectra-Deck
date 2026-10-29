"""The plot right-click menus: sections, nothing lost, the page-wide entries
(real window; needs a display and matplotlib, skipped otherwise).

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


def entries(menu):
    """[(type, label, state)] of one menu's entries (separators have no
    label)."""
    out = []
    for i in range((menu.index("end") or -1) + 1 if menu.index("end")
                   is not None else 0):
        kind = menu.type(i)
        if kind == "separator":
            out.append((kind, "", ""))
        else:
            out.append((kind, menu.entrycget(i, "label").strip("✓ ").strip(),
                        menu.entrycget(i, "state")))
    return out


def walk(menu, path=()):
    """Every labelled entry of a menu tree as (path, label, state, kind)."""
    found = []
    for i in range((menu.index("end") if menu.index("end") is not None
                    else -1) + 1):
        kind = menu.type(i)
        if kind == "separator":
            continue
        label = menu.entrycget(i, "label").replace("✓", "").strip()
        state = menu.entrycget(i, "state")
        found.append((path, label, state, kind))
        if kind == "cascade":
            sub = menu.nametowidget(menu.entrycget(i, "menu"))
            found += walk(sub, path + (label,))
    return found


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestPlotMenus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {"a.vms": doc("a.vms", [reg("C 1s", "A", 286.0),
                                           reg("O 1s", "A", 288.0)])}
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
        ws.panels_var.set("Auto")
        ws.cols_var.set("Auto")
        ws.shape_var.set("Auto")
        ws.panel_views = {}
        ws.panel_order = []
        ws.panel_slots = []
        ws.panel_start = 0
        for p in self.docs:
            ws._add_file(p, file_id=p[:-4], refresh=False)
        ws._finish_adding()
        ws.checked = {id(r) for d in ws.docs for r in d.regions}
        self.keys = [k for k, _ in ws._groups()]

    def test_top_level_is_four_sections(self):
        menu = self.ws._build_panel_menu(self.keys[0])
        top = [(k, label) for k, label, _s in entries(menu) if k != "separator"]
        self.assertEqual(top, [
            ("command", f"Panel: {self.keys[0]}"), ("cascade", "Spectrum"),
            ("cascade", "Fit and peaks"), ("cascade", "Layout"),
            ("cascade", "Save")])

    def test_nothing_was_lost(self):
        labels = {label for _p, label, _s, _k in
                  walk(self.ws._build_panel_menu(self.keys[0]))}
        for want in ("View", "Normalise", "Reset Range end points", "Offset",
                     "Z axis", "Reverse order", "Use the page's view",
                     "Fit layers", "Other possible lines", "Identify peaks…",
                     "Move panel to first", "Move panel earlier",
                     "Move panel later", "Move panel to last",
                     "Reset panel order", "Shape of this position",
                     "Size of this position", "Reset page layout",
                     "Columns", "Page shape", "Spacing…",
                     "Save this panel as image…",
                     "Save plot image…", "Save this view as a figure…",
                     "Stack", "Heatmap", "Components", "Page setting",
                     "Tall (2 rows)", "Landscape 4:3"):
            self.assertIn(want, labels, want)

    def test_entries_sit_in_their_sections(self):
        by = {label: path for path, label, _s, _k in
              walk(self.ws._build_panel_menu(self.keys[0]))}
        self.assertEqual(by["View"], ("Spectrum",))
        self.assertEqual(by["Identify peaks…"], ("Fit and peaks",))
        self.assertEqual(by["Move panel to last"], ("Layout",))
        self.assertEqual(by["Reset page layout"], ("Layout",))
        self.assertEqual(by["Save plot image…"], ("Save",))
        self.assertEqual(by["Columns"], ("Layout",))

    def test_series_view_disables_fit_and_peaks(self):
        ws = self.ws
        top = {label: state for _p, label, state, _k in
               walk(ws._build_panel_menu(self.keys[0])) if _p == ()}
        self.assertEqual(top["Fit and peaks"], "normal")
        ws._set_panel_view(self.keys[0], "view", "Heatmap")
        got = {(p, label): state for p, label, state, _k in
               walk(ws._build_panel_menu(self.keys[0]))}
        self.assertEqual(got[((), "Fit and peaks")], "disabled")
        self.assertEqual(got[(("Spectrum",), "Z axis")], "normal")
        self.assertEqual(got[(("Spectrum",), "Offset")], "disabled")

    def test_layout_commands_work_from_the_menu(self):
        ws = self.ws
        menu = ws._build_panel_menu(self.keys[1])
        found = {}

        def find(m, want):
            for i in range((m.index("end") or -1) + 1):
                if m.type(i) == "separator":
                    continue
                label = m.entrycget(i, "label").replace("✓", "").strip()
                if m.type(i) == "cascade":
                    r = find(m.nametowidget(m.entrycget(i, "menu")), want)
                    if r:
                        return r
                if label == want:
                    return m, i
            return None
        m, i = find(menu, "Tall (2 rows)")                 # Size of position 2
        m.invoke(i)
        self.assertEqual(ws.panel_slots, [{}, {"span": "Tall (2 rows)"}])
        menu = ws._build_panel_menu(self.keys[1])
        m, i = find(menu, "Reset page layout")
        self.assertEqual(m.entrycget(i, "state"), "normal")
        m.invoke(i)
        self.assertEqual(ws.panel_slots, [])
        menu = ws._build_panel_menu(self.keys[0])
        cols_menu, i = find(menu, "Columns")
        sub = cols_menu.nametowidget(cols_menu.entrycget(i, "menu"))
        m, i = find(sub, "2")                              # Columns > 2
        m.invoke(i)
        self.assertEqual(ws.cols_var.get(), "2")
        shape_menu, i = find(menu, "Page shape")
        sub = shape_menu.nametowidget(shape_menu.entrycget(i, "menu"))
        m, i = find(sub, "Square")                         # Page shape > Square
        m.invoke(i)
        self.assertEqual(ws.shape_var.get(), "Square")
        self.assertEqual(ws.panel_slots, [])               # not the position's

    def test_position_header_follows_the_panel(self):
        ws = self.ws
        first = [label for _p, label, _s, _k in
                 walk(ws._build_panel_menu(self.keys[0]))
                 if label.startswith("Page layout")]
        second = [label for _p, label, _s, _k in
                  walk(ws._build_panel_menu(self.keys[1]))
                  if label.startswith("Page layout")]
        self.assertEqual((first, second), (["Page layout, position 1"],
                                            ["Page layout, position 2"]))

    def test_background_menu_has_layout_and_save_only(self):
        menu = self.ws._build_figure_menu()
        top = [label for k, label, _s in entries(menu) if k != "separator"]
        self.assertEqual(top, ["Layout", "Save"])
        labels = {label for _p, label, _s, _k in walk(menu)}
        self.assertIn("Reset panel order", labels)
        self.assertIn("Columns", labels)
        self.assertIn("Save plot image…", labels)
        self.assertIn("Save this view as a figure…", labels)
        self.assertNotIn("Save this panel as image…", labels)
        self.assertNotIn("Shape of this position", labels)
        self.assertNotIn("Move panel to first", labels)

    def test_menus_are_rebuilt_not_stacked(self):
        ws = self.ws
        a = len(walk(ws._build_panel_menu(self.keys[0])))
        b = len(walk(ws._build_panel_menu(self.keys[0])))
        self.assertEqual(a, b)
        c = len(walk(ws._build_figure_menu()))
        d = len(walk(ws._build_figure_menu()))
        self.assertEqual(c, d)


if __name__ == "__main__":
    unittest.main()
