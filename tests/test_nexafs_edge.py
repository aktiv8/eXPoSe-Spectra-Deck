"""NEXAFS edge normalisation: the maths (``nexafs``), the saved settings and
metadata (``Annotations``), the drawn copy (``Workspace._display``), the dialog
and the methods text. The "max" method is checked against the cached cells of
the author's own *NEXAFS Processing Template* (a spreadsheet), which gives
real numbers for it; "step" is checked on constructed edges whose step is
known.

Real data: ``XPS_NEXAFS_CORPUS=<folder>`` where ``<folder>/B07 XPS/NEXAFS
Processing Template.xlsx`` is the template.

Run:  python -m unittest discover tests
"""

import math
import os
import re
import shutil
import sys
import tempfile
import unittest
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import annotations  # noqa: E402
import methods  # noqa: E402
import nexafs  # noqa: E402
import readers  # noqa: E402

try:
    import h5py  # noqa: F401
    HAVE_H5 = True
except ImportError:                                 # pragma: no cover
    HAVE_H5 = False

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

from test_nexus_nexafs import write_scan  # noqa: E402


def edge(sign=1.0, e0=290.0, n=701, e_lo=270.0, step_ev=0.1, jump=3.2):
    """A scan with straight pre- and post-edge lines and a smooth rise
    between 289 and 291 eV: the lines are exact outside it, so the edge step
    at ``e0`` is known (``jump``)."""
    energy = [e_lo + step_ev * i for i in range(n)]

    def pre(e):
        return 2.0 + 0.01 * (e - 270.0)

    def post(e):
        return 2.0 + 0.01 * (e - 270.0) + jump + 0.02 * (e - e0)

    y = []
    for e in energy:
        t = min(1.0, max(0.0, (e - 289.0) / 2.0))
        w = t * t * (3 - 2 * t)
        y.append(sign * (pre(e) * (1 - w) + post(e) * w))
    return energy, y


class TestPolynomial(unittest.TestCase):
    def test_it_recovers_a_line_and_a_parabola_exactly(self):
        xs = [1500.0 + i for i in range(20)]
        line = nexafs.fit_poly(xs, [3.0 - 0.5 * (x - 1500) for x in xs], 1)
        self.assertAlmostEqual(nexafs.poly_at(line, 1510.0), -2.0, places=9)
        par = nexafs.fit_poly(xs, [1 + 2 * (x - 1505) ** 2 for x in xs], 2)
        self.assertAlmostEqual(nexafs.poly_at(par, 1505.0), 1.0, places=6)
        self.assertAlmostEqual(nexafs.poly_at(par, 1510.0), 51.0, places=5)

    def test_a_constant_is_the_mean(self):
        fit = nexafs.fit_poly([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], 0)
        self.assertAlmostEqual(nexafs.poly_at(fit, 7.0), 4.0)

    def test_too_few_points_or_a_singular_system_is_none(self):
        self.assertIsNone(nexafs.fit_poly([1.0], [1.0], 1))
        self.assertIsNone(nexafs.fit_poly([2.0, 2.0, 2.0], [1.0, 2.0, 3.0], 1))

    def test_a_not_a_number_point_is_ignored(self):
        fit = nexafs.fit_poly([1.0, 2.0, 3.0, 4.0],
                              [1.0, float("nan"), 3.0, 4.0], 1)
        self.assertAlmostEqual(nexafs.poly_at(fit, 2.0), 2.0)


class TestSettings(unittest.TestCase):
    def test_sanitise_keeps_only_complete_valid_settings(self):
        ok = {"mode": "step", "pre": [274, 270], "post": [330, 340],
              "e0": 290, "pre_order": 2}
        s = nexafs.sanitise_edge(ok)
        self.assertEqual(s["pre"], [270.0, 274.0])        # ordered
        self.assertEqual((s["pre_order"], s["post_order"]), (2, 1))
        for bad in ({}, None, {"mode": "other", "pre": [1, 2]},
                    {"mode": "max"}, {"mode": "max", "pre": [1]},
                    {"mode": "step", "pre": [1, 2]},
                    {"mode": "step", "pre": [1, 2], "post": [3, 4],
                     "e0": "x"},
                    {"mode": "max", "pre": [1, float("nan")]}):
            self.assertEqual(nexafs.sanitise_edge(bad), {})
        self.assertEqual(nexafs.sanitise_edge(
            {"mode": "max", "pre": [270, 274], "pre_order": 7})["pre_order"], 1)

    def test_each_reason_a_choice_cannot_be_applied(self):
        e, y = edge()
        good = {"mode": "step", "pre": [270, 274], "post": [330, 340],
                "e0": 290}
        self.assertEqual(nexafs.check_edge(good, e, y), [])
        self.assertTrue(nexafs.check_edge({}, e, y))
        self.assertTrue(any("fewer than" in m for m in nexafs.check_edge(
            dict(good, pre=[270, 270.05]), e, y)))
        self.assertTrue(any("outside the data" in m for m in nexafs.check_edge(
            dict(good, post=[500, 510]), e, y)))
        self.assertTrue(any("E0" in m for m in nexafs.check_edge(
            dict(good, e0=500), e, y)))
        self.assertTrue(any("pre-edge window must end" in m
                            for m in nexafs.check_edge(
                                dict(good, pre=[270, 300]), e, y)))
        self.assertTrue(any("post-edge window must start" in m
                            for m in nexafs.check_edge(
                                dict(good, post=[280, 340]), e, y)))

    def test_suggestions_follow_the_templates_rule(self):
        e, y = edge()
        s = nexafs.suggest_edge(e, y, "max")
        self.assertEqual(s, {"mode": "max", "pre": [270.0, 274.0],
                             "pre_order": 1})
        s = nexafs.suggest_edge(e, y, "step")
        self.assertEqual(s["pre"], [270.0, 274.0])
        self.assertAlmostEqual(s["e0"], 290.0, delta=1.0)
        self.assertLess(s["post"][0], 340.0 - 5 * 0.1)
        self.assertEqual(nexafs.check_edge(s, e, y), [])
        self.assertEqual(nexafs.suggest_edge(e[:10], y[:10]), {})

    def test_a_falling_edge_is_found_too(self):
        e, y = edge(sign=-1.0)
        self.assertAlmostEqual(nexafs.suggest_edge(e, y, "step")["e0"], 290.0,
                               delta=1.0)

    def test_the_text_states_the_windows(self):
        t = nexafs.edge_text({"mode": "step", "pre": [270, 274],
                              "post": [330, 340], "e0": 290.5,
                              "post_order": 2})
        self.assertEqual(t, "edge step: pre-edge 270-274 eV (linear), "
                            "post-edge 330-340 eV (quadratic), E0 290.5 eV; "
                            "step set to 1")
        self.assertIn("maximum of 1", nexafs.edge_text(
            {"mode": "max", "pre": [270, 274]}))


class TestEdgeStep(unittest.TestCase):
    P = {"mode": "step", "pre": [270, 280], "post": [300, 340], "e0": 290.0}

    def test_the_step_is_recovered_and_set_to_one(self):
        e, y = edge()
        out, info = nexafs.edge_normalise(e, y, self.P)
        self.assertAlmostEqual(info["step"], 3.2, places=9)
        i0 = e.index(270.0)
        self.assertAlmostEqual(out[i0], 0.0, places=9)        # pre-edge: 0
        # the post-edge line sits one step above the pre-edge line at E0
        i290 = min(range(len(e)), key=lambda i: abs(e[i] - 290.0))
        self.assertAlmostEqual(out[i290], 0.5, places=1)      # mid-rise
        self.assertAlmostEqual(out[-1], 1.0 + 0.02 * 50 / 3.2, places=9)

    def test_a_spectrum_of_negative_values_normalises_to_a_positive_edge(self):
        e, y = edge(sign=-1.0)
        out, info = nexafs.edge_normalise(e, y, self.P)
        self.assertAlmostEqual(info["step"], -3.2, places=9)
        self.assertAlmostEqual(out[0], 0.0, places=9)
        self.assertGreater(out[-1], 1.0)
        pos, _ = nexafs.edge_normalise(*edge(), self.P)
        for a, b in zip(out, pos):
            self.assertAlmostEqual(a, b, places=9)             # same result

    def test_a_high_energy_scale_does_not_cost_precision(self):
        e = [x + 1230.0 for x in edge()[0]]            # the same edge at 1520 eV
        out, info = nexafs.edge_normalise(
            e, edge()[1], {"mode": "step", "pre": [1500, 1510],
                           "post": [1530, 1570], "e0": 1520.0})
        self.assertAlmostEqual(info["step"], 3.2, places=9)

    def test_the_result_is_independent_of_the_scale_of_the_data(self):
        e, y = edge()
        a, _ = nexafs.edge_normalise(e, y, self.P)
        b, _ = nexafs.edge_normalise(e, [v * 1e-10 for v in y], self.P)
        for u, v in zip(a, b):
            self.assertAlmostEqual(u, v, places=9)

    def test_an_unusable_choice_returns_the_reasons_not_numbers(self):
        e, y = edge()
        out, why = nexafs.edge_normalise(e, y, dict(self.P, e0=900.0))
        self.assertIsNone(out)
        self.assertTrue(why)

    def test_a_flat_spectrum_has_no_step(self):
        e = [270 + 0.1 * i for i in range(701)]
        out, why = nexafs.edge_normalise(e, [5.0] * 701, self.P)
        self.assertIsNone(out)
        self.assertEqual(why, ["the edge step is zero"])

    def test_points_that_are_not_numbers_stay_not_numbers(self):
        e, y = edge()
        y[100] = float("nan")
        out, _ = nexafs.edge_normalise(e, y, self.P)
        self.assertTrue(math.isnan(out[100]))
        self.assertFalse(math.isnan(out[101]))


class TestMaxMethod(unittest.TestCase):
    P = {"mode": "max", "pre": [270, 274]}

    def test_it_spans_zero_to_one(self):
        e, y = edge()
        out, info = nexafs.edge_normalise(e, y, self.P)
        self.assertAlmostEqual(min(out), 0.0)
        self.assertAlmostEqual(max(out), 1.0)
        self.assertEqual(info["mode"], "max")

    def test_a_flat_result_is_refused(self):
        e = [270 + 0.1 * i for i in range(701)]
        out, why = nexafs.edge_normalise(e, [2.0 + 0.01 * (x - 270) for x in e],
                                         self.P)
        self.assertIsNone(out)
        self.assertEqual(why, ["the pre-edge-subtracted spectrum is flat"])


def _cached_column(xml, col, rows):
    out = {}
    for m in re.finditer(r'<c r="%s(\d+)"[^>]*?>(.*?)</c>' % col, xml, re.S):
        v = re.search(r"<v>(.*?)</v>", m.group(2))
        if v:
            out[int(m.group(1))] = v.group(1)
    return out


@unittest.skipUnless(
    os.environ.get("XPS_NEXAFS_CORPUS") and os.path.exists(os.path.join(
        os.environ.get("XPS_NEXAFS_CORPUS", ""), "B07 XPS",
        "NEXAFS Processing Template.xlsx")),
    "set XPS_NEXAFS_CORPUS to the folder holding 'B07 XPS' (the template)")
class TestAgainstTheTemplate(unittest.TestCase):
    """The spreadsheet's cached cells are real numbers: TEY / I0 (column AF,
    ca35b / ca18b), a line fitted over the automatic pre-edge window (the
    first 4 eV), subtracted (AM) and then ``(AM - min) x max`` (AN)."""

    @classmethod
    def setUpClass(cls):
        path = os.path.join(os.environ["XPS_NEXAFS_CORPUS"], "B07 XPS",
                            "NEXAFS Processing Template.xlsx")
        with zipfile.ZipFile(path) as z:
            xml = z.read("xl/worksheets/sheet1.xml").decode()
        cols = {c: _cached_column(xml, c, None) for c in
                ("V", "W", "X", "AF", "AM", "AN")}
        rows = sorted(r for r in cols["V"] if r >= 2 and r in cols["AN"])
        cls.energy = [float(cols["V"][r]) for r in rows]
        cls.ratio = [float(cols["X"][r]) / float(cols["W"][r]) for r in rows]
        cls.af = [float(cols["AF"][r]) for r in rows]
        cls.am = [float(cols["AM"][r]) for r in rows]
        cls.an = [float(cols["AN"][r]) for r in rows]

    def test_the_sheets_own_ratio_is_what_is_normalised(self):
        for a, b in zip(self.ratio, self.af):
            self.assertAlmostEqual(a / b, 1.0, places=9)

    def test_the_pre_edge_subtracted_curve_is_the_sheets_column_am(self):
        lo = min(self.energy)
        out, info = nexafs.edge_normalise(
            self.energy, self.ratio,
            {"mode": "max", "pre": [lo, lo + 4.0], "pre_order": 1})
        rng = max(self.ratio) - min(self.ratio)          # the sheet's AG6
        z_over_range = [(v * (info["maximum"] - info["minimum"])
                         + info["minimum"]) / rng for v in out]
        for mine, theirs in zip(z_over_range, self.am):
            self.assertAlmostEqual(mine, theirs, places=9)

    def test_the_final_column_differs_only_by_the_templates_last_factor(self):
        lo = min(self.energy)
        out, info = nexafs.edge_normalise(
            self.energy, self.ratio,
            {"mode": "max", "pre": [lo, lo + 4.0], "pre_order": 1})
        rng = max(self.ratio) - min(self.ratio)
        z_range = info["maximum"] - info["minimum"]
        factor = z_range / rng * info["maximum"] / rng   # (max-min) x max, as AM
        for mine, theirs in zip(out, self.an):
            self.assertAlmostEqual(mine * factor, theirs, places=9)
        self.assertAlmostEqual(factor, 1.0, delta=0.005)  # the 0.2 % in the doc


class TestProcess(unittest.TestCase):
    def region(self):
        from readers import Region
        e, y = edge()
        r = Region("ca35b", 0, 0, technique="NEXAFS", energy=e, counts=y,
                   energy_label="Photon Energy", decodable=True)
        ring = [300.0 - 0.02 * i for i in range(len(e))]
        r.extra["ring_current_points"] = ring
        return r

    def test_nothing_chosen_is_nothing_done(self):
        self.assertIsNone(nexafs.process(self.region(), False, {}))

    def test_the_ring_scaling_comes_first_and_the_edge_second(self):
        r = self.region()
        p = dict(TestEdgeStep.P)
        both = nexafs.process(r, True, p)
        scaled, mean = nexafs.scale_to_mean(r.counts, r.extra[
            "ring_current_points"])
        out, info = nexafs.edge_normalise(r.energy, scaled, p)
        self.assertEqual(both.mean, mean)
        self.assertEqual(both.counts, out)
        self.assertEqual(both.edge["step"], info["step"])

    def test_an_unusable_edge_leaves_the_counts_as_they_were(self):
        r = self.region()
        self.assertIsNone(nexafs.process(
            r, False, dict(TestEdgeStep.P, e0=900.0)))
        only_ring = nexafs.process(r, True, dict(TestEdgeStep.P, e0=900.0))
        self.assertIsNotNone(only_ring.mean)
        self.assertIsNone(only_ring.edge)


class TestAnnotations(unittest.TestCase):
    P = {"mode": "step", "pre": [270, 280], "post": [300, 340], "e0": 290.0}

    def test_saved_per_spectrum_and_read_back(self):
        a = annotations.Annotations()
        self.assertTrue(a.is_empty())
        a.set_edge("f1", "S", "ca35b", self.P)
        self.assertFalse(a.is_empty())
        self.assertEqual(a.edge_for("f1", "S", "ca35b")["e0"], 290.0)
        self.assertEqual(a.edge_for("f1", "S", "ca18b"), {})
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.edge_for("f1", "S", "ca35b"),
                         a.edge_for("f1", "S", "ca35b"))
        a.set_edge("f1", "S", "ca35b", {})
        self.assertTrue(a.is_empty())

    def test_bad_saved_settings_are_dropped_not_trusted(self):
        a = annotations.Annotations.from_json(
            {"nexafs_edge": {"f|S|x": {"mode": "step", "pre": [1, 2]},
                             "f|S|y": "text", "f|S|z": dict(self.P)}})
        self.assertEqual(list(a.nexafs_edge), ["f|S|z"])
        self.assertNotIn("nexafs_edge", annotations.Annotations().to_json())


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestWithAFile(unittest.TestCase):
    P = {"mode": "step", "pre": [270, 272], "post": [278, 280], "e0": 275.0}

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)                    # 270-280 eV, a current ramp
        cls.f = readers.load_file(cls.path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def annotated(self, ring=False, edge=None):
        a = annotations.Annotations()
        a.nexafs_ring = ring
        if edge:
            for r in self.f.regions:
                a.set_edge("f1", r.sample, r.name, edge)
        self.f.annotations, self.f.file_id = a, "f1"
        return a

    def tearDown(self):
        self.f.annotations = None

    def test_the_metadata_says_what_was_done(self):
        self.annotated(edge=self.P)
        md = self.f.region_metadata(self.f.regions[0])
        self.assertIn("edge step: pre-edge 270-272 eV (linear)",
                      md["Edge normalisation"])
        self.assertIn("E0 275 eV", md["Edge normalisation"])
        self.assertTrue(md["Edge step"].endswith(" A"))
        self.assertNotIn("Normalisation", md)            # no ring scaling

    def test_nothing_is_stated_without_a_saved_choice(self):
        self.annotated()
        md = self.f.region_metadata(self.f.regions[0])
        self.assertNotIn("Edge normalisation", md)
        self.assertIn("the spectra are not normalised",
                      methods.generate(self.f.metadata_rows()))

    def test_the_methods_text_states_it_once_and_drops_not_normalised(self):
        self.annotated(ring=True, edge=self.P)
        text = methods.generate(self.f.metadata_rows())
        self.assertIn("The spectra were normalised (edge step: pre-edge "
                      "270-272 eV (linear), post-edge 278-280 eV (linear), "
                      "E0 275 eV; step set to 1)", text)
        self.assertIn("scaled to the mean ring current", text)
        self.assertNotIn("not normalised", text)

    def test_different_windows_are_counted_not_listed(self):
        a = self.annotated(edge=self.P)
        a.set_edge("f1", self.f.regions[1].sample, self.f.regions[1].name,
                   dict(self.P, e0=276.0))
        text = methods.generate(self.f.metadata_rows())
        self.assertIn("with the windows listed in the metadata (2 different "
                      "settings)", text)


@unittest.skipUnless(HAVE_H5 and HAVE_MPL, "h5py / matplotlib / Tk missing")
class TestInTheWindow(unittest.TestCase):
    P = {"mode": "step", "pre": [270, 272], "post": [278, 280], "e0": 275.0}

    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror, ee.messagebox.askyesno)
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None
        ee.messagebox.askyesno = lambda *a, **k: True
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path, channels=("ca18b", "ca35b", "ca36b"))

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror, ee.messagebox.askyesno) = cls._boxes
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def setUp(self):
        self.ws = ee.Workspace(self.root)
        self.ws._add_files([self.path])
        self.regions = self.ws.docs[0].regions
        self.ws.checked = {id(r) for r in self.regions}

    def dialog(self):
        from nexafs_ui import NexafsEdgeDialog
        dlg = NexafsEdgeDialog(self.root, self.ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        return dlg

    def test_nothing_is_normalised_until_a_choice_is_saved(self):
        self.assertIs(self.ws._display(self.regions[0]), self.regions[0])

    def test_a_saved_choice_changes_the_drawn_copy_only(self):
        r = self.regions[0]
        original = list(r.counts)
        self.ws.nexafs_edge_set([r], self.P)
        shown = self.ws._display(r)
        self.assertEqual(r.counts, original)
        out, info = nexafs.edge_normalise(r.energy, original, self.P)
        self.assertEqual(shown.counts, out)
        self.assertEqual((shown.count_label, shown.count_units),
                         ("Normalised intensity", "edge step = 1"))
        self.assertEqual(shown.extra["edge_norm"]["step"], info["step"])
        self.ws.nexafs_edge_set([r], {})
        self.assertIs(self.ws._display(r), r)

    def test_it_draws_with_the_new_label(self):
        self.ws.nexafs_edge_set(self.regions, dict(self.P, mode="max",
                                                   post=None, e0=None))
        self.ws._render()
        self.root.update_idletasks()
        labels = [ax.get_ylabel() for ax in self.ws.fig.axes]
        self.assertTrue(any("Normalised intensity" in t and "0 to 1" in t
                            for t in labels), labels)

    def test_the_dialog_suggests_but_saves_nothing_by_itself(self):
        dlg = self.dialog()
        self.assertEqual(len(dlg.regions), 3)
        self.assertEqual(dlg.params()["mode"], "step")     # a suggestion
        self.assertEqual(self.ws.nexafs_edge_get(self.regions[0]), {})
        self.assertIn("Not saved yet", dlg.status.cget("text"))
        self.assertIs(self.ws._display(self.regions[0]), self.regions[0])

    def test_use_for_this_spectrum_saves_and_says_so(self):
        dlg = self.dialog()
        dlg.set_values(pre=(270, 272), post=(278, 280), e0=275.0)
        self.assertTrue(dlg.use_here())
        saved = self.ws.nexafs_edge_get(self.regions[0])
        self.assertEqual((saved["pre"], saved["post"], saved["e0"]),
                         ([270.0, 272.0], [278.0, 280.0], 275.0))
        self.assertEqual(self.ws.nexafs_edge_get(self.regions[1]), {})
        self.assertIn("Saved with this file", dlg.status.cget("text"))

    def test_an_unusable_choice_is_not_saved_and_says_why(self):
        dlg = self.dialog()
        dlg.set_values(pre=(270, 272), post=(278, 280), e0=500.0)
        self.assertFalse(dlg.use_here())
        self.assertEqual(self.ws.nexafs_edge_get(self.regions[0]), {})
        self.assertIn("Cannot apply", dlg.status.cget("text"))
        self.assertIn("E0", dlg.status.cget("text"))

    def test_use_for_all_listed_saves_where_the_windows_fit(self):
        dlg = self.dialog()
        dlg.set_values(pre=(270, 272), post=(278, 280), e0=275.0)
        fits = dlg.use_all()
        self.assertEqual(len(fits), 3)
        for r in self.regions:
            self.assertEqual(self.ws.nexafs_edge_get(r)["e0"], 275.0)

    def test_the_max_method_needs_only_the_pre_edge(self):
        dlg = self.dialog()
        dlg.set_mode("max")
        dlg.set_values(pre=(270, 272))
        out, info = dlg.preview()
        self.assertAlmostEqual(min(out), 0.0)
        self.assertAlmostEqual(max(out), 1.0)
        self.assertTrue(dlg.use_here())
        self.assertEqual(self.ws.nexafs_edge_get(self.regions[0])["mode"],
                         "max")

    def test_remove_forgets_the_choice(self):
        self.ws.nexafs_edge_set(self.regions, self.P)
        dlg = self.dialog()
        self.assertIn("Saved with this file", dlg.status.cget("text"))
        dlg.remove()
        self.assertEqual(self.ws.nexafs_edge_get(self.regions[0]), {})
        self.assertIs(self.ws._display(self.regions[0]), self.regions[0])

    def test_dragging_on_the_plot_sets_the_window_chosen(self):
        dlg = self.dialog()
        dlg.drag.set("post")
        dlg.set_window("post", 277.0, 279.5)
        self.assertEqual((dlg.vars["post_lo"].get(), dlg.vars["post_hi"].get()),
                         ("277", "279.5"))

    def test_it_travels_with_the_workbook_annotations(self):
        self.ws.nexafs_edge_set(self.regions, self.P)
        saved = self.ws.ann.to_json()
        self.ws.ann = annotations.Annotations.from_json(saved)
        self.ws._ann_changed(relabel=False)
        self.assertIsNot(self.ws._display(self.regions[0]), self.regions[0])

    def test_the_menu_has_the_entry_and_it_dims_without_data(self):
        self.assertIn("NEXAFS edge step…", ee.TOOLS_NEEDS_DATA)


if __name__ == "__main__":
    unittest.main()
