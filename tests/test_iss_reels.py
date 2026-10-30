"""ISS kinematics and element identification, the REELS band gap, their
storage beside the data, and how they are drawn.

Run:  python -m unittest discover tests
"""

import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import annotations  # noqa: E402
import elements as el  # noqa: E402
import reels  # noqa: E402
from readers import Region  # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

HE, AU, CU = 4.0026, 196.9666, 62.9296


class TestKinematics(unittest.TestCase):
    def test_head_on_backscatter_is_the_textbook_formula(self):
        for m2 in (CU, AU, 12.0):
            a = m2 / HE
            self.assertAlmostEqual(el.kinematic_factor(HE, m2, 180.0),
                                   ((a - 1) / (a + 1)) ** 2, 12)

    def test_ninety_degrees_is_the_other_textbook_formula(self):
        for m2 in (CU, AU, 12.0):
            a = m2 / HE
            self.assertAlmostEqual(el.kinematic_factor(HE, m2, 90.0),
                                   (a - 1) / (a + 1), 12)

    def test_heavier_targets_scatter_more_energy(self):
        ks = [el.kinematic_factor(HE, m, 123.03)
              for m in (12.0, 28.0, 63.0, 197.0)]
        self.assertEqual(ks, sorted(ks))
        self.assertTrue(all(0 < k < 1 for k in ks))

    def test_a_lighter_target_cannot_backscatter(self):
        self.assertIsNone(el.kinematic_factor(20.0, HE, 123.0))
        self.assertIsNone(el.kinematic_factor(20.0, 19.0, 150.0))
        self.assertIsNone(el.kinematic_factor(20.0, HE, 30.0))   # beyond arcsin(A)
        self.assertIsNotNone(el.kinematic_factor(20.0, HE, 5.0))  # small angles ok
        self.assertIsNone(el.kinematic_factor(0.0, 10.0, 100.0))

    def test_energy_scales_with_the_beam(self):
        self.assertAlmostEqual(el.iss_energy(2000.0, HE, CU, 123.03)
                               / el.iss_energy(1000.0, HE, CU, 123.03), 2.0)

    def test_known_peak_positions(self):
        # He+ 1 keV, 123.03 deg (worked out from the formula)
        for sym, mass in (("Cu", CU), ("Au", AU)):
            a = mass / HE
            th = math.radians(123.03)
            want = 1000 * ((math.cos(th) + math.sqrt(a * a - math.sin(th) ** 2))
                           / (1 + a)) ** 2
            e = el.iss_energy(1000.0, el.ion_mass("He+"), el.ELEMENTS[sym][1],
                              123.03)
            self.assertAlmostEqual(e, want, 6)
        self.assertGreater(el.iss_energy(1000, HE, AU, 123.03),
                           el.iss_energy(1000, HE, CU, 123.03))

    def test_table_is_sorted_and_excludes_light_targets(self):
        t = el.iss_table(1000.0, "He+", 123.03)
        e = [d["energy"] for d in t]
        self.assertEqual(e, sorted(e))
        syms = {d["symbol"] for d in t}
        self.assertNotIn("H", syms)
        self.assertNotIn("He", syms)                 # not heavier than He
        self.assertIn("Cu", syms)
        self.assertNotIn("He", {d["symbol"] for d in el.iss_table(
            1000, "Ne+", 123.03)})
        self.assertNotIn("C", {d["symbol"] for d in el.iss_table(
            1000, "Ar+", 123.03)})                   # lighter than the ion

    def test_candidates_near_a_measured_peak(self):
        cu = el.iss_energy(1000.0, HE, CU, 123.03)
        hits = el.candidates(cu + 1.2, 1000.0, "He+", 123.03)
        self.assertEqual(hits[0]["symbol"], "Cu")
        self.assertAlmostEqual(hits[0]["delta"], -1.2, 6)
        self.assertEqual([abs(h["delta"]) for h in hits],
                         sorted(abs(h["delta"]) for h in hits))
        self.assertEqual(el.candidates(10.0, 1000.0), [])
        wide = el.candidates(cu, 1000.0, window=60)
        self.assertGreater(len(wide), len(hits))

    def test_mass_from_energy_inverts_the_kinematics(self):
        for m in (12.0, CU, AU):
            e = el.iss_energy(1000.0, HE, m, 123.03)
            self.assertAlmostEqual(el.mass_from_energy(e, 1000.0, "He+",
                                                       123.03), m, 4)
        self.assertIsNone(el.mass_from_energy(-5.0, 1000.0))       # not physical
        self.assertIsNone(el.mass_from_energy(1000.0, 1000.0))     # above E0 x k_max
        self.assertIsNone(el.mass_from_energy(500.0, 100.0))

    def test_ion_names(self):
        self.assertEqual(el.ion_mass("He"), HE)
        self.assertEqual(el.ion_mass("Ne+"), el.ION_MASS["Ne+"])
        with self.assertRaises(ValueError):
            el.ion_mass("Zz+")

    def test_symbols_and_marks(self):
        self.assertEqual(el.parse_symbols("cu, AU  xx ni Cu"),
                         ["Cu", "Au", "Ni"])
        self.assertEqual(el.parse_symbols(""), [])
        self.assertEqual(el.parse_symbols(None), [])
        m = el.marks_for(["Cu", "Au", "C"], 1000.0, "He+", 123.03)
        self.assertEqual([s for s, _e in m], ["Cu", "Au", "C"])
        m = el.marks_for(["Cu", "Au"], 1000.0, "He+", 123.03, lo=900, hi=960)
        self.assertEqual([s for s, _e in m], ["Au"])   # only the ones in range
        self.assertEqual(el.marks_for(["He"], 1000.0, "Ne+", 123.03), [])

    def test_element_table_is_sane(self):
        zs = [v[0] for v in el.ELEMENTS.values()]
        self.assertEqual(zs, sorted(zs))
        self.assertEqual(zs[0], 1)
        for sym, (z, mass, name) in el.ELEMENTS.items():
            self.assertGreaterEqual(mass, z * 0.99, sym)     # mass >= ~Z
            self.assertTrue(name)


class TestReels(unittest.TestCase):
    def spectrum(self, gap=3.2, slope=400.0, base=100.0, el_ke=1000.0):
        """Loss spectrum: elastic peak, flat baseline, then a linear onset."""
        ke, y = [], []
        for i in range(0, 1201):
            loss = i * 0.01                          # 0..12 eV
            v = base + (slope * (loss - gap) if loss > gap else 0.0)
            v += 20000 * math.exp(-(loss / 0.25) ** 2)
            ke.append(el_ke - loss)
            y.append(v)
        return ke, y

    def test_the_gap_is_recovered_from_two_points_on_the_edge(self):
        ke, y = self.spectrum(gap=3.2)
        loss = reels.loss_axis(ke, 1000.0)
        p1 = (4.0, 100 + 400 * 0.8)
        p2 = (5.0, 100 + 400 * 1.8)
        res = reels.band_gap(loss, y, p1, p2)
        self.assertAlmostEqual(res["gap"], 3.2, 6)
        self.assertAlmostEqual(res["slope"], 400.0, 6)
        self.assertAlmostEqual(res["base"], 100.0, 6)

    def test_order_of_the_two_points_does_not_matter(self):
        ke, y = self.spectrum()
        loss = reels.loss_axis(ke, 1000.0)
        a = reels.band_gap(loss, y, (4.0, 420.0), (5.0, 820.0))
        b = reels.band_gap(loss, y, (5.0, 820.0), (4.0, 420.0))
        self.assertAlmostEqual(a["gap"], b["gap"], 9)

    def test_an_explicit_baseline_wins(self):
        ke, y = self.spectrum()
        loss = reels.loss_axis(ke, 1000.0)
        res = reels.band_gap(loss, y, (4.0, 420.0), (5.0, 820.0), base=0.0)
        self.assertAlmostEqual(res["gap"], 4.0 - 420.0 / 400.0, 9)

    def test_bad_picks_give_none(self):
        ke, y = self.spectrum()
        loss = reels.loss_axis(ke, 1000.0)
        self.assertIsNone(reels.band_gap(loss, y, (4.0, 420.0), (4.0, 820.0)))
        self.assertIsNone(reels.band_gap(loss, y, (4.0, 820.0), (5.0, 420.0)))
        self.assertIsNone(reels.band_gap(loss, y, (4.0, 100.0), (5.0, 100.0)))
        self.assertIsNone(reels.band_gap(loss, y, (float("nan"), 1.0),
                                         (2.0, 2.0)))
        # a tangent that would cross the baseline at a negative loss
        self.assertIsNone(reels.band_gap(loss, y, (0.1, 900.0), (0.2, 1300.0),
                                         base=100.0))

    def test_baseline_ignores_the_elastic_tail_and_falls_back(self):
        ke, y = self.spectrum(base=100.0)
        loss = reels.loss_axis(ke, 1000.0)
        self.assertAlmostEqual(reels.baseline(loss, y, (4.0, 5.0)), 100.0, 6)
        # onset right after the elastic tail: too few points between
        near = reels.baseline(loss, y, (0.6, 5.0))
        self.assertIsNotNone(near)
        self.assertIsNone(reels.baseline([], [], (1.0, 2.0)))

    def test_elastic_peak_and_axes(self):
        ke, y = self.spectrum(el_ke=1486.6)
        self.assertAlmostEqual(reels.elastic_peak(ke, y), 1486.6, 6)
        self.assertIsNone(reels.elastic_peak([], []))
        self.assertIsNone(reels.elastic_peak([1, 2], [1]))
        self.assertEqual(reels.loss_axis([10, 9], 10), [0, 1])

    def test_tangent_points_run_from_the_gap_up_the_edge(self):
        res = {"gap": 3.2, "slope": 400.0, "base": 100.0}
        a, b = reels.tangent_points(res, (4.0, 420.0), (5.0, 820.0))
        self.assertEqual(a, (3.2, 100.0))
        self.assertGreater(b[0], 5.0)
        self.assertAlmostEqual(b[1], 100 + 400 * (b[0] - 3.2), 6)


class TestStorage(unittest.TestCase):
    def test_kinetic_markers_round_trip_and_stay_apart(self):
        a = annotations.Annotations()
        a.add_marker("f1", "S", "ISS", 850.5, "Cu", kin=True)
        a.add_marker("f1", "S", "ISS", 850.5, "Cu", kin=True)      # duplicate
        a.add_marker("f1", "S", "ISS", 285.0, "C 1s")               # binding
        ms = a.markers_for("f1", "S", "ISS")
        self.assertEqual(len(ms), 2)
        self.assertTrue(ms[0]["kin"])
        self.assertNotIn("kin", ms[1])
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.markers_for("f1", "S", "ISS"), ms)

    def test_reels_round_trip_metadata_and_clearing(self):
        a = annotations.Annotations()
        self.assertIsNone(a.reels_for("f1", "S", "R"))
        d = {"elastic": 1486.6, "p1": [4.0, 420.0], "p2": [5.0, 820.0],
             "gap": 3.2, "base": 100.0, "slope": 400.0}
        a.set_reels("f1", "S", "R", d)
        self.assertFalse(a.is_empty())
        self.assertEqual(a.reels_for("f1", "S", "R")["gap"], 3.2)
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.reels_for("f1", "S", "R"), a.reels_for("f1", "S", "R"))
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S", "Region": "R"})
        self.assertEqual(md["REELS band gap (eV)"], "3.20")
        self.assertEqual(md["REELS elastic peak (eV KE)"], "1486.60")
        a.set_reels("f1", "S", "R", None)
        self.assertTrue(a.is_empty())

    def test_a_construction_without_a_gap_is_kept_but_not_reported(self):
        a = annotations.Annotations()
        a.set_reels("f1", "S", "R", {"elastic": 1.0, "p1": [1, 2],
                                     "p2": [3, 4]})
        self.assertIsNone(a.reels_for("f1", "S", "R")["gap"])
        r = Region(name="R", index=0, offset=0, sample="S")
        self.assertNotIn("REELS band gap (eV)",
                         a.apply_metadata("f1", 0, r, {"Sample": "S"}))

    def test_loading_is_tolerant(self):
        b = annotations.Annotations.from_json({"reels": {
            "a": {"elastic": "x"}, "b": "junk", "c": {"elastic": 1,
                                                       "p1": [1], "p2": [1, 2]},
            "d": {"elastic": float("nan"), "p1": [1, 2], "p2": [1, 2]},
            "ok": {"elastic": 5, "p1": [1, 2], "p2": [2, 3], "gap": "bad"}}})
        self.assertEqual(list(b.reels), ["ok"])
        self.assertIsNone(b.reels["ok"]["gap"])
        self.assertEqual(annotations.Annotations.from_json(
            {"reels": [1]}).reels, {})
        b = annotations.Annotations.from_json({"markers": {"k": [
            {"be": 1, "label": "x", "kin": "yes"}, {"be": "z", "label": "y"}]}})
        self.assertEqual(b.markers["k"], [{"be": 1.0, "label": "x"}])


def region(kinetic=True, hv=1000.0, n=501):
    """An ISS-like spectrum: KE axis, peaks at Cu and Au energies."""
    cu = el.iss_energy(hv, HE, CU, 123.03)
    au = el.iss_energy(hv, HE, AU, 123.03)
    ke = [600.0 + i * 0.75 for i in range(n)]
    y = [50 + 900 * math.exp(-((k - cu) / 6) ** 2)
         + 1500 * math.exp(-((k - au) / 6) ** 2) for k in ke]
    if kinetic:
        return Region(name="ISS", index=0, offset=0, energy=ke, counts=y,
                      decodable=True, sample="S", photon_energy=hv,
                      energy_label="Kinetic Energy", technique="ISS")
    return Region(name="ISS", index=0, offset=0,
                  energy=[hv - k for k in ke], counts=y, decodable=True,
                  sample="S", photon_energy=hv, energy_label="Binding Energy")


@unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
class TestDrawing(unittest.TestCase):
    def draw(self, r, markers=(), reels_arg=None, norm="None"):
        import plots
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        plots.draw_stack(ax, [r], 0.6, norm, None, ["#1f77b4"], "ISS", "",
                         (), False, markers=markers, reels=reels_arg,
                         scale="Kinetic" if "kinetic" in r.energy_label.lower()
                         else "Binding")
        fig.canvas.draw()
        return fig, ax

    def marker_xs(self, ax):
        return sorted(round(float(t.get_position()[0]), 1) for t in ax.texts
                      if t.get_text() in ("Cu", "Au"))

    def test_kinetic_markers_land_on_the_peak_on_a_kinetic_axis(self):
        r = region(kinetic=True)
        cu = el.iss_energy(1000.0, HE, CU, 123.03)
        au = el.iss_energy(1000.0, HE, AU, 123.03)
        _fig, ax = self.draw(r, [(cu, "Cu", True), (au, "Au", True)])
        self.assertEqual(self.marker_xs(ax), sorted([round(cu, 1),
                                                     round(au, 1)]))

    def test_the_same_markers_follow_the_axis_when_it_is_binding(self):
        r = region(kinetic=False)
        cu = el.iss_energy(1000.0, HE, CU, 123.03)
        _fig, ax = self.draw(r, [(cu, "Cu", True)])
        self.assertEqual(self.marker_xs(ax), [round(1000.0 - cu, 1)])

    def test_binding_markers_are_unchanged(self):
        r = region(kinetic=False)
        _fig, ax = self.draw(r, [(400.0, "Cu")])            # 2-tuple as before
        self.assertEqual(self.marker_xs(ax), [400.0])

    def reels_dict(self):
        return {"elastic": 1000.0, "p1": [4.0, 420.0], "p2": [5.0, 820.0],
                "gap": 3.2, "base": 100.0, "slope": 400.0}

    def reels_region(self):
        ke, y = TestReels().spectrum(el_ke=1000.0)
        return Region(name="REELS", index=0, offset=0, energy=ke, counts=y,
                      decodable=True, sample="S", photon_energy=1000.0,
                      energy_label="Kinetic Energy")

    def test_reels_overlay_draws_the_tangent_and_names_the_gap(self):
        r = self.reels_region()
        _fig, ax = self.draw(r, reels_arg=self.reels_dict())
        self.assertTrue(any("Eg = 3.20 eV" in t.get_text()
                            for t in ax.texts))
        tangent = [ln for ln in ax.lines if len(ln.get_xdata()) == 2
                   and ln.get_linestyle() == "-"]
        self.assertTrue(tangent)
        xs = list(tangent[-1].get_xdata())
        self.assertAlmostEqual(max(xs), 1000.0 - 3.2, 5)        # KE = el - gap

    def test_reels_overlay_follows_a_binding_axis(self):
        ke, y = TestReels().spectrum(el_ke=1000.0)
        r = Region(name="REELS", index=0, offset=0,
                   energy=[1000.0 - k for k in ke], counts=y, decodable=True,
                   sample="S", photon_energy=1000.0,
                   energy_label="Binding Energy")
        _fig, ax = self.draw(r, reels_arg=self.reels_dict())
        tangent = [ln for ln in ax.lines if len(ln.get_xdata()) == 2
                   and ln.get_linestyle() == "-"][-1]
        # loss L sits at BE = hv - (el - L) = L on this axis (el = hv)
        self.assertAlmostEqual(min(tangent.get_xdata()), 3.2, 5)

    def test_reels_overlay_scales_with_normalisation(self):
        r = self.reels_region()
        _fig, ax = self.draw(r, reels_arg=self.reels_dict(), norm="Max = 1")
        tangent = [ln for ln in ax.lines if len(ln.get_xdata()) == 2
                   and ln.get_linestyle() == "-"][-1]
        self.assertLess(max(tangent.get_ydata()), 5.0)          # not raw counts

    def test_incomplete_construction_draws_only_the_points(self):
        r = self.reels_region()
        d = {"elastic": 1000.0, "p1": [4.0, 420.0], "p2": [5.0, 820.0],
             "gap": None}
        _fig, ax = self.draw(r, reels_arg=d)
        self.assertFalse(any("Eg" in t.get_text() for t in ax.texts))

    def test_no_overlay_on_a_stack(self):
        import plots
        fig = Figure()
        ax = fig.add_subplot(111)
        r1, r2 = self.reels_region(), self.reels_region()
        plots.draw_stack(ax, [r1, r2], 0.6, "None", None, ["#000", "#111"],
                         "R", "", (), False, reels=self.reels_dict(),
                         scale="Kinetic")
        self.assertFalse(any("Eg" in t.get_text() for t in ax.texts))


def plasmon_spectrum(centre=27.0, sigma=3.0, ramp=0.0, level=0.0, n=1001,
                     el_ke=1000.0):
    """Loss spectrum 0..50 eV: elastic peak, a π+σ-like Gaussian and a linear
    ramp (descending-loss, i.e. ascending-KE order, as the app holds it)."""
    loss = [i * 0.05 for i in range(n)]
    y = [level + ramp * x + 1000 * math.exp(-0.5 * ((x - centre) / sigma) ** 2)
         + 20000 * math.exp(-(x / 0.25) ** 2) for x in loss]
    return [el_ke - x for x in loss], y, loss


class TestPlasmon(unittest.TestCase):
    def test_gaussian_on_a_linear_background_is_exact(self):
        ke, y, loss = plasmon_spectrum(ramp=5.0, level=50.0)
        res = reels.plasmon(loss, y, (10, 50), "Linear")
        self.assertAlmostEqual(res["ep"], 27.0, delta=0.01)
        self.assertEqual(res["method"], "gauss")
        self.assertEqual(res["background"], "Linear")
        self.assertAlmostEqual(res["fwhm"], 3.0 * reels.FWHM_PER_SIGMA,
                               delta=0.1)

    def test_a_flat_background_is_recovered_exactly(self):
        ke, y, loss = plasmon_spectrum(level=100.0)
        res = reels.plasmon(loss, y, (10, 50), "Linear")
        self.assertAlmostEqual(res["ep"], 27.0, delta=0.05)

    def test_input_order_does_not_matter(self):
        ke, y, loss = plasmon_spectrum()
        a = reels.plasmon(loss, y, (15, 45), "Linear")
        b = reels.plasmon(loss[::-1], y[::-1], (45, 15), "Linear")
        self.assertAlmostEqual(a["ep"], b["ep"], 9)

    def test_second_derivative_ignores_a_linear_background(self):
        _ke, y0, loss = plasmon_spectrum(ramp=0.0)
        _ke, y1, _l = plasmon_spectrum(ramp=40.0, level=300.0)
        a = reels.plasmon(loss, y0, (15, 45), "None", "d2", smooth=1.5)
        b = reels.plasmon(loss, y1, (15, 45), "None", "d2", smooth=1.5)
        self.assertAlmostEqual(a["ep"], 27.0, delta=0.1)
        self.assertAlmostEqual(a["ep"], b["ep"], delta=0.02)
        self.assertIsNone(a["fwhm"])
        self.assertEqual(a["background"], "None")
        self.assertEqual(a["bg"], [])

    def test_shirley_moves_a_skewed_peak_the_way_the_paper_finds(self):
        # a peak with a long high-loss tail on a rising background: the
        # Shirley background takes more from the high-loss side than a
        # straight line, so Ep and the width come out smaller (Filippi &
        # Calliari 2006: 26.2 eV / 14.6 eV against 27.1 eV / 15.6 eV)
        loss = [i * 0.05 for i in range(1001)]
        y = []
        for x in loss:
            g = 1000 * math.exp(-0.5 * ((x - 27.0) / (5.0 if x < 27 else 8.0)) ** 2)
            y.append(g + 20000 * math.exp(-(x / 0.25) ** 2))
        lin = reels.plasmon(loss, y, (12, 48), "Linear")
        shi = reels.plasmon(loss, y, (12, 48), "Shirley")
        self.assertIsNotNone(shi)
        self.assertLess(shi["ep"], lin["ep"])
        self.assertLess(shi["fwhm"], lin["fwhm"])

    def test_the_background_curve_is_kept_for_drawing(self):
        _ke, y, loss = plasmon_spectrum(level=100.0)
        res = reels.plasmon(loss, y, (10, 50), "Linear")
        self.assertTrue(10 < len(res["bg"]) <= 2 * reels.BG_POINTS)
        self.assertTrue(all(10 <= a <= 50 for a, _b in res["bg"]))
        self.assertAlmostEqual(res["bg"][0][1], 100.0, 3)

    def test_refusals(self):
        _ke, y, loss = plasmon_spectrum()
        self.assertIsNone(reels.plasmon(loss, y, (15, 15), "Linear"))
        self.assertIsNone(reels.plasmon(loss, y, (15, float("nan")), "Linear"))
        self.assertIsNone(reels.plasmon(loss, y, (15, 45), "Cubic"))
        self.assertIsNone(reels.plasmon(loss, y, (15, 45), "Linear", "x"))
        self.assertIsNone(reels.plasmon(loss, y, (15, 45), "Linear",
                                        frac=1.5))
        self.assertIsNone(reels.plasmon(loss, y, (60, 70), "Linear"))   # no data
        self.assertIsNone(reels.plasmon([1, 2, 3], [1, 2, 3], (0, 4),
                                        "Linear"))                      # few
        self.assertIsNone(reels.plasmon(loss, y[:-1], (15, 45), "Linear"))
        flat = [5.0] * len(loss)
        self.assertIsNone(reels.plasmon(loss, flat, (15, 45), "Linear"))

    def test_shirley_with_the_step_the_wrong_way_round_is_refused(self):
        loss = [i * 0.05 for i in range(1001)]
        y = [(60 - x) * 10 + 1000 * math.exp(-0.5 * ((x - 27) / 4) ** 2)
             for x in loss]          # falls with loss: high-loss end is lower
        self.assertIsNone(reels.plasmon(loss, y, (15, 45), "Shirley"))
        self.assertIsNotNone(reels.plasmon(loss, y, (15, 45), "Linear"))

    def test_non_finite_points_are_dropped(self):
        _ke, y, loss = plasmon_spectrum(level=100.0)
        y = list(y)
        y[400] = float("nan")
        res = reels.plasmon(loss, y, (10, 50), "Linear")
        self.assertAlmostEqual(res["ep"], 27.0, delta=0.05)

    def test_method_text_names_everything_that_moves_the_number(self):
        _ke, y, loss = plasmon_spectrum()
        g = reels.method_text(reels.plasmon(
            loss, [v + 5 * x for v, x in zip(y, loss)], (10, 50), "Shirley"))
        self.assertIn("90 %", g)
        self.assertIn("shirley", g)
        self.assertIn("10–50", g)
        d = reels.method_text(reels.plasmon(loss, y, (15, 45), "None", "d2",
                                            smooth=2.0))
        self.assertIn("second derivative", d)
        self.assertIn("2", d)

    def test_density_reproduces_the_published_pairs(self):
        # Yang et al., PCCP 23, 25335 (2021): 25.5 eV -> 2.046 g/cm3 and the
        # raw REELS peak 27.0 eV -> 2.294 g/cm3 (m* = 0.87 m, 4 electrons)
        self.assertAlmostEqual(reels.density_from_plasmon(25.5), 2.046,
                               delta=0.001)
        self.assertAlmostEqual(reels.density_from_plasmon(27.0), 2.294,
                               delta=0.001)
        self.assertIsNone(reels.density_from_plasmon(0))
        self.assertIsNone(reels.density_from_plasmon(25.5, m_eff=0))


class TestPlasmonStorage(unittest.TestCase):
    def result(self):
        _ke, y, loss = plasmon_spectrum(level=100.0, ramp=5.0)
        return {**reels.plasmon(loss, y, (10, 50), "Shirley"),
                "elastic": 1000.0}

    def test_round_trip_metadata_and_clearing(self):
        a = annotations.Annotations()
        self.assertIsNone(a.plasmon_for("f1", "S", "R"))
        a.set_plasmon("f1", "S", "R", self.result())
        self.assertFalse(a.is_empty())
        got = a.plasmon_for("f1", "S", "R")
        self.assertAlmostEqual(got["ep"], self.result()["ep"], 9)
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.plasmon_for("f1", "S", "R"), got)
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S", "Region": "R"})
        self.assertEqual(md["REELS plasmon energy (eV)"], f"{got['ep']:.2f}")
        self.assertIn("shirley", md["REELS plasmon method"])
        self.assertEqual(md["REELS elastic peak (eV KE)"], "1000.00")
        a.set_plasmon("f1", "S", "R", None)
        self.assertTrue(a.is_empty())

    def test_it_is_independent_of_the_band_gap(self):
        a = annotations.Annotations()
        a.set_plasmon("f1", "S", "R", self.result())
        a.set_reels("f1", "S", "R", None)                 # clearing the gap
        self.assertIsNotNone(a.plasmon_for("f1", "S", "R"))
        a.set_reels("f1", "S", "R", {"elastic": 1000.0, "p1": [4, 420],
                                     "p2": [5, 820], "gap": 3.2})
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S"})
        self.assertIn("REELS band gap (eV)", md)
        self.assertIn("REELS plasmon energy (eV)", md)

    def test_loading_is_tolerant(self):
        good = self.result()
        b = annotations.Annotations.from_json({"plasmon": {
            "ok": good, "junk": "x", "nowin": {**good, "window": [1]},
            "badm": {**good, "method": "fft"},
            "badbg": {**good, "background": "Spline"},
            "nan": {**good, "ep": float("nan")},
            "noel": {k: v for k, v in good.items() if k != "elastic"}}})
        self.assertEqual(list(b.plasmon), ["ok"])
        self.assertEqual(annotations.Annotations.from_json(
            {"plasmon": [1]}).plasmon, {})
        c = annotations.Annotations.from_json({"plasmon": {"k": {
            **good, "bg": [[1, 2], ["x", 3], [4]], "fwhm": "bad"}}})
        self.assertEqual(c.plasmon["k"]["bg"], [[1.0, 2.0]])
        self.assertIsNone(c.plasmon["k"]["fwhm"])


@unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
class TestPlasmonDrawing(unittest.TestCase):
    def setUp(self):
        ke, y, _loss = plasmon_spectrum(level=100.0)
        self.region = Region(name="REELS", index=0, offset=0, energy=ke,
                             counts=y, decodable=True, sample="S",
                             photon_energy=1000.0,
                             energy_label="Kinetic Energy")
        _k, _y, loss = plasmon_spectrum(level=100.0)
        self.pl = {**reels.plasmon(loss, y, (10, 50), "Linear"),
                   "elastic": 1000.0}

    def draw(self, regions, norm="None", scale="Kinetic"):
        import plots
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        plots.draw_stack(ax, regions, 0.6, norm, None,
                         ["#1f77b4"] * len(regions), "R", "", (), False,
                         plasmon=self.pl, scale=scale)
        fig.canvas.draw()
        return ax

    def lines(self, ax, style):
        return [ln for ln in ax.lines if ln.get_linestyle() == style]

    def test_the_energy_the_window_and_the_background_are_drawn(self):
        ax = self.draw([self.region])
        self.assertTrue(any("Ep = 27.0 eV" in t.get_text() for t in ax.texts))
        xs = [ln.get_xdata()[0] for ln in ax.lines
              if len(ln.get_xdata()) == 2 and ln.get_linestyle() == "-"]
        self.assertTrue(any(abs(x - (1000.0 - 27.0)) < 0.2 for x in xs))
        self.assertTrue(self.lines(ax, "--"))              # the background
        self.assertTrue(ax.patches or ax.collections)      # the shaded window

    def test_it_follows_a_binding_axis(self):
        r = Region(name="REELS", index=0, offset=0,
                   energy=[1000.0 - k for k in self.region.energy],
                   counts=self.region.counts, decodable=True, sample="S",
                   photon_energy=1000.0, energy_label="Binding Energy")
        ax = self.draw([r], scale="Binding")
        xs = [ln.get_xdata()[0] for ln in ax.lines
              if len(ln.get_xdata()) == 2 and ln.get_linestyle() == "-"]
        self.assertTrue(any(abs(x - 27.0) < 0.2 for x in xs))

    def test_the_background_scales_with_normalisation(self):
        ax = self.draw([self.region], norm="Max = 1")
        bg = self.lines(ax, "--")[-1]
        self.assertLess(max(bg.get_ydata()), 5.0)

    def test_no_overlay_on_a_stack(self):
        ax = self.draw([self.region, self.region])
        self.assertFalse(any("Ep =" in t.get_text() for t in ax.texts))


LOSS_TRUTH = {                       # key: (loss, FWHM, area) of a synthetic spectrum
    "pi": (5.8, 3.0, 900.0), "inter": (12.5, 7.0, 300.0),
    "sp2S": (19.5, 7.0, 500.0), "sp3S": (23.0, 7.0, 600.0),
    "sp2B": (27.0, 7.0, 1500.0), "sp3B": (33.8, 7.0, 500.0),
    "ss35": (35.0, 7.0, 100.0), "ss42": (42.0, 7.0, 150.0)}
LAM_6_5 = 6.5 / (math.cos(math.radians(35)) / (1 + math.cos(math.radians(35)))
                 * math.log(100))     # the mean free path that gives their 6.5 nm


def loss_spectrum(truth=LOSS_TRUTH, level=20.0, ramp=0.2, el_ke=2000.0):
    """Elastic peak plus the loss-region Gaussians on a linear background, as
    plain lists in the kinetic-energy order the app holds."""
    loss = [i * 0.05 for i in range(1201)]
    y = []
    for x in loss:
        v = level + ramp * x + 50000 * math.exp(-(x / 0.3) ** 2)
        for c, w, a in truth.values():
            v += a * (2 * math.sqrt(math.log(2) / math.pi) / w) * math.exp(
                -4 * math.log(2) * ((x - c) / w) ** 2)
        y.append(v)
    return [el_ke - x for x in loss], y, loss


class TestLossFit(unittest.TestCase):
    def test_the_components_and_the_ratio_are_recovered(self):
        _ke, y, loss = loss_spectrum()
        res = reels.decompose(loss, y, (3, 50), "Linear")
        by = {c["key"]: c for c in res["components"]}
        self.assertEqual(len(by), 8)
        self.assertLess(res["rms"], 0.02)
        self.assertAlmostEqual(by["pi"]["loss"], 5.8, delta=0.2)
        self.assertAlmostEqual(by["pi"]["fwhm"], 3.0, delta=0.5)
        self.assertAlmostEqual(by["sp2B"]["loss"], 27.0, delta=0.5)
        self.assertAlmostEqual(by["sp2S"]["loss"], 19.5, delta=1.0)
        # neighbouring Gaussians trade area, so R is good to about 10 %
        self.assertAlmostEqual(res["r"], 1500 / 500, delta=0.5)
        self.assertEqual(res["background"], "Linear")

    def test_positions_held_when_the_tolerance_is_zero(self):
        _ke, y, loss = loss_spectrum()
        res = reels.decompose(loss, y, (3, 50), "Linear", tol=0)
        for (key, _n, c0, _g), c in zip(reels.LOSS_COMPONENTS,
                                        res["components"]):
            self.assertEqual(c["key"], key)
            self.assertEqual(c["loss"], c0)

    def test_positions_stay_within_the_tolerance(self):
        _ke, y, loss = loss_spectrum()
        res = reels.decompose(loss, y, (3, 50), "Linear", tol=0.4)
        for (_k, _n, c0, _g), c in zip(reels.LOSS_COMPONENTS,
                                       res["components"]):
            self.assertLessEqual(abs(c["loss"] - c0), 0.4 + 1e-9)

    def test_areas_are_never_negative_and_a_missing_component_is_zero(self):
        truth = {k: v for k, v in LOSS_TRUTH.items() if k != "sp3B"}
        _ke, y, loss = loss_spectrum(truth)
        res = reels.decompose(loss, y, (3, 50), "Linear")
        self.assertTrue(all(c["area"] >= 0 for c in res["components"]))

    def test_input_order_does_not_matter(self):
        _ke, y, loss = loss_spectrum()
        a = reels.decompose(loss, y, (3, 50), "Linear")
        b = reels.decompose(loss[::-1], y[::-1], (50, 3), "Linear")
        self.assertAlmostEqual(a["r"], b["r"], 6)

    def test_refusals(self):
        _ke, y, loss = loss_spectrum()
        self.assertIsNone(reels.decompose(loss, y, (3, 3), "Linear"))
        self.assertIsNone(reels.decompose(loss, y, (3, 50), "Cubic"))
        self.assertIsNone(reels.decompose(loss, y, (3, float("nan")),
                                          "Linear"))
        self.assertIsNone(reels.decompose(loss, y, (70, 80), "Linear"))
        self.assertIsNone(reels.decompose(loss[:10], y[:10], (0, 0.4),
                                          "Linear"))
        self.assertIsNone(reels.decompose(loss, y[:-1], (3, 50), "Linear"))
        self.assertIsNone(reels.decompose(loss, [5.0] * len(loss), (3, 50),
                                          "Linear"))
        self.assertIsNone(reels.decompose(loss, y, (3, 50), "Linear",
                                          components=[]))

    def test_no_sp2_surface_gives_no_ratio(self):
        _ke, y, loss = loss_spectrum()
        comps = [c for c in reels.LOSS_COMPONENTS if c[0] != "sp2S"]
        res = reels.decompose(loss, y, (3, 50), "Linear", components=comps)
        self.assertIsNone(res["r"])

    def test_the_curves_add_up(self):
        _ke, y, loss = loss_spectrum()
        res = reels.decompose(loss, y, (3, 50), "Linear")
        import numpy as np
        x = np.linspace(3, 50, 200)
        total, parts = reels.loss_curves(res, x)
        bg = np.interp(x, [p[0] for p in res["bg"]], [p[1] for p in res["bg"]])
        self.assertTrue(np.allclose(total, bg + sum(p for _c, p in parts)))
        self.assertEqual(len(parts), 8)

    def test_gaussian_unit_has_unit_area_and_the_stated_width(self):
        import numpy as np
        x = np.linspace(-30, 30, 60001)
        g = reels.gaussian_unit(x, 2.0, 5.0)
        self.assertAlmostEqual(float(g.sum() * (x[1] - x[0])), 1.0, 6)
        half = x[g >= g.max() / 2]
        self.assertAlmostEqual(float(half[-1] - half[0]), 5.0, 3)


class TestLayers(unittest.TestCase):
    def test_the_information_depth_formula(self):
        self.assertAlmostEqual(reels.information_depth(LAM_6_5, 0, 35), 6.5, 9)
        # normal incidence and emission: lam/2 ln 100
        self.assertAlmostEqual(reels.information_depth(3.0, 0, 0),
                               1.5 * math.log(100), 9)
        for bad in ((0, 0, 0), (3, 90, 0), (3, 0, 90), (-1, 0, 0),
                    ("x", 0, 0)):
            self.assertIsNone(reels.information_depth(*bad))
        self.assertIsNone(reels.information_depth(3, 0, 0, p=100))

    def test_the_published_ratios_give_their_layer_counts(self):
        # Stobinski et al. 2014: R = 2.72 -> 8 layers (commercial GO) and
        # R = 1.55 -> 4 layers (commercial rGO), 2 keV, 0 and 35 degrees,
        # graphite information depth 6.5 nm
        go = reels.layers(2.72, LAM_6_5, 0, 35)
        rgo = reels.layers(1.55, LAM_6_5, 0, 35)
        self.assertEqual(round(go["n"]), 8)
        self.assertEqual(round(rgo["n"]), 4)
        self.assertAlmostEqual(go["id"], 6.5, 9)
        self.assertAlmostEqual(go["n_graphite"], 6.5 / 0.335, 9)

    def test_one_layer_at_zero_and_graphite_at_the_graphite_ratio(self):
        self.assertAlmostEqual(reels.layers(0, LAM_6_5, 0, 35)["n"], 1.0, 9)
        g = reels.layers(4.2, LAM_6_5, 0, 35)
        self.assertIsNone(g["n"])
        self.assertAlmostEqual(g["at_least"], 6.5 / 0.335, 9)
        self.assertIsNone(reels.layers(9.0, LAM_6_5, 0, 35)["n"])

    def test_more_layers_for_a_larger_ratio(self):
        ns = [reels.layers(r, LAM_6_5, 0, 35)["n"] for r in (0.5, 1, 2, 3, 4)]
        self.assertEqual(ns, sorted(ns))

    def test_a_longer_path_reads_more_layers_for_the_same_ratio(self):
        # a deeper probe reaches the bulk sooner only if the film is thicker
        short = reels.layers(2.0, 2.0, 0, 35)["n"]
        long_ = reels.layers(2.0, 4.0, 0, 35)["n"]
        self.assertLess(short, long_)

    def test_undefined_inputs(self):
        self.assertIsNone(reels.layers(-1, LAM_6_5, 0, 35))
        self.assertIsNone(reels.layers(float("nan"), LAM_6_5, 0, 35))
        self.assertIsNone(reels.layers(2, LAM_6_5, 90, 35))
        self.assertIsNone(reels.layers(2, 0.1, 0, 35))     # under one layer deep
        self.assertIsNone(reels.layers(2, LAM_6_5, 0, 35, d=0))
        self.assertIsNone(reels.layers(2, LAM_6_5, 0, 35, r_graphite=0))


class TestLossFitStorage(unittest.TestCase):
    def fit(self, geometry=True):
        _ke, y, loss = loss_spectrum()
        d = {**reels.decompose(loss, y, (3, 50), "Linear"),
             "elastic": 2000.0}
        if geometry:
            d["geometry"] = {"lam": LAM_6_5, "a_in": 0.0, "a_out": 35.0,
                             "d": 0.335, "r_graphite": 4.2, "p": 99.0}
            d["layers"] = reels.layers(d["r"], LAM_6_5, 0, 35)
        return d

    def test_round_trip_metadata_and_clearing(self):
        a = annotations.Annotations()
        self.assertIsNone(a.loss_fit_for("f1", "S", "R"))
        d = self.fit()
        a.set_loss_fit("f1", "S", "R", d)
        self.assertFalse(a.is_empty())
        got = a.loss_fit_for("f1", "S", "R")
        self.assertAlmostEqual(got["r"], d["r"], 9)
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.loss_fit_for("f1", "S", "R"), got)
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S", "Region": "R"})
        self.assertEqual(md["REELS sp2 bulk / sp2 surface (R)"],
                         f"{d['r']:.2f}")
        self.assertIn("REELS π loss (eV)", md)
        self.assertIn("REELS C sp2 surface loss (eV)", md)
        self.assertIn("linear background", md["REELS loss fit"])
        self.assertIn("graphite R 4.2", md["REELS layers (estimate)"])
        a.set_loss_fit("f1", "S", "R", None)
        self.assertTrue(a.is_empty())

    def test_no_layers_row_without_a_geometry(self):
        a = annotations.Annotations()
        a.set_loss_fit("f1", "S", "R", self.fit(geometry=False))
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S"})
        self.assertNotIn("REELS layers (estimate)", md)
        self.assertIn("REELS sp2 bulk / sp2 surface (R)", md)

    def test_a_graphite_like_ratio_says_or_more(self):
        d = self.fit()
        d["layers"] = reels.layers(9.0, LAM_6_5, 0, 35)
        a = annotations.Annotations()
        a.set_loss_fit("f1", "S", "R", d)
        r = Region(name="R", index=0, offset=0, sample="S")
        md = a.apply_metadata("f1", 0, r, {"Sample": "S"})
        self.assertIn("or more", md["REELS layers (estimate)"])

    def test_it_is_independent_of_the_other_reels_results(self):
        a = annotations.Annotations()
        a.set_loss_fit("f1", "S", "R", self.fit())
        a.set_reels("f1", "S", "R", None)
        a.set_plasmon("f1", "S", "R", None)
        self.assertIsNotNone(a.loss_fit_for("f1", "S", "R"))

    def test_loading_is_tolerant(self):
        good = self.fit()
        b = annotations.Annotations.from_json({"loss_fit": {
            "ok": good, "junk": "x", "nowin": {**good, "window": [1]},
            "badbg": {**good, "background": "Spline"},
            "nocomp": {**good, "components": []},
            "noel": {k: v for k, v in good.items() if k != "elastic"}}})
        self.assertEqual(list(b.loss_fit), ["ok"])
        self.assertEqual(annotations.Annotations.from_json(
            {"loss_fit": [1]}).loss_fit, {})
        c = annotations.Annotations.from_json({"loss_fit": {"k": {
            **good, "components": good["components"] + [
                {"key": "x", "name": "x", "loss": "bad", "fwhm": 1,
                 "area": 1},
                {"key": "y", "name": "y", "loss": 1, "fwhm": -1, "area": 1}],
            "geometry": {"lam": 3}}}})
        self.assertEqual(len(c.loss_fit["k"]["components"]), 8)
        self.assertNotIn("geometry", c.loss_fit["k"])     # incomplete: dropped


@unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
class TestLossFitDrawing(unittest.TestCase):
    def setUp(self):
        ke, y, loss = loss_spectrum()
        self.region = Region(name="REELS", index=0, offset=0, energy=ke,
                             counts=y, decodable=True, sample="S",
                             photon_energy=2000.0,
                             energy_label="Kinetic Energy")
        self.fit = {**reels.decompose(loss, y, (3, 50), "Linear"),
                    "elastic": 2000.0}

    def draw(self, regions, norm="None"):
        import plots
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        plots.draw_stack(ax, regions, 0.6, norm, None,
                         ["#1f77b4"] * len(regions), "R", "", (), False,
                         lossfit=self.fit, scale="Kinetic")
        fig.canvas.draw()
        return ax

    def test_the_ratio_the_components_and_the_sum_are_drawn(self):
        ax = self.draw([self.region])
        self.assertTrue(any(t.get_text().startswith("R = ")
                            for t in ax.texts))
        self.assertTrue(any(t.get_text() == "C sp2 bulk" for t in ax.texts))
        # the sum and every component with an area (the 35 eV one has none)
        n_drawn = sum(1 for c in self.fit["components"] if c["area"] > 0)
        self.assertEqual(len([ln for ln in ax.lines
                              if len(ln.get_xdata()) == 400]), n_drawn + 1)
        self.assertGreaterEqual(len(ax.collections), 3)        # the fills

    def test_the_sum_follows_the_normalisation(self):
        ax = self.draw([self.region], norm="Max = 1")
        sums = [ln for ln in ax.lines if len(ln.get_xdata()) == 400]
        self.assertLess(max(float(max(ln.get_ydata())) for ln in sums), 5.0)

    def test_no_overlay_on_a_stack(self):
        ax = self.draw([self.region, self.region])
        self.assertFalse(any(t.get_text().startswith("R = ")
                             for t in ax.texts))


class TestBrowserPayload(unittest.TestCase):
    def test_kinetic_markers_are_flagged_and_not_shifted(self):
        import htmlbrowser
        from readers import SpectrumFile
        r = region(kinetic=True)
        f = SpectrumFile()
        f.path = "a.vms"
        f.format_name = "Test"
        f.regions = [r]
        f._finish()
        ann = annotations.Annotations()
        f.annotations, f.file_id = ann, "f1"
        ann.add_marker("f1", "S", "ISS", 821.3, "Cu", kin=True)
        ann.add_marker("f1", "S", "ISS", 285.0, "C 1s")
        ann.set_shift(annotations.sample_key("f1", "S"), 0.8)
        marks = htmlbrowser.build_payload([f])["samples"][0]["regions"][0][
            "markers"]
        self.assertEqual(marks[0], {"be": 821.3, "label": "Cu", "kin": True})
        self.assertEqual(marks[1], {"be": 285.8, "label": "C 1s"})


if __name__ == "__main__":
    unittest.main()
