"""Fitting a SnapMap: the solvers, the model, the refinement on the mean
spectrum and the per-pixel partition into chemical states (snapfit.py), on
synthetic maps whose answer is known.

Run:  python -m unittest discover tests
"""

import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

try:
    import numpy as np
    HAVE = True
except Exception:                                   # pragma: no cover
    HAVE = False

import snapfit  # noqa: E402

if HAVE:
    import lineshapes
    import snapmap


class TestSpec(unittest.TestCase):
    def test_a_bad_spec_becomes_a_valid_one(self):
        s = snapfit.sanitise({"components": [
            {"name": "Si", "be": 99.4, "fwhm": 99, "shape": "nonsense"},
            {"be": "x"}, {"name": "Si", "be": 103.4,
                          "doublet": {"split": 0.6, "ratio": 0.5}}],
            "bin": 7, "background": "weird", "window": [110, 95]})
        self.assertEqual([c["name"] for c in s["components"]], ["Si", "Si 2"])
        self.assertEqual(s["components"][0]["fwhm"], snapfit.MAX_FWHM)
        self.assertEqual(s["components"][0]["shape"], snapfit.DEFAULT_SHAPE)
        self.assertEqual(s["components"][1]["doublet"],
                         {"split": 0.6, "ratio": 0.5})
        self.assertEqual((s["bin"], s["background"]), (1, "shirley"))
        self.assertEqual(s["window"], [95.0, 110.0])
        self.assertEqual(snapfit.sanitise(json.loads(json.dumps(s))), s)
        self.assertEqual(snapfit.sanitise(None)["components"], [])

    def test_doublets_follow_the_degeneracies(self):
        self.assertAlmostEqual(snapfit.doublet_for("2p", 0.61)["ratio"], 0.5)
        self.assertAlmostEqual(snapfit.doublet_for("3d", 0.7)["ratio"], 2 / 3)
        self.assertAlmostEqual(snapfit.doublet_for("4f", 1.0)["ratio"], 0.75)
        self.assertIsNone(snapfit.doublet_for("1s", 0.5))
        self.assertIsNone(snapfit.doublet_for("2p", 0))


class TestSeeds(unittest.TestCase):
    def test_si_2p_is_seeded_from_the_bundled_states(self):
        comps = snapfit.seed_components("Si 2p")
        self.assertEqual([c["name"] for c in comps], ["Si", "SiO2"])
        self.assertEqual([c["be"] for c in comps], [99.4, 103.4])
        d = comps[0]["doublet"]
        self.assertAlmostEqual(d["split"], 100.03 - 99.4, places=6)  # the line table's
        self.assertAlmostEqual(d["ratio"], 0.5)
        self.assertEqual(snapfit.seed_components("Si2p")[0]["name"], "Si")

    def test_other_levels_are_seeded_by_their_main_line(self):
        c = snapfit.seed_components("C 1s")            # many table states: one start
        self.assertEqual([(x["name"], x["be"], x["doublet"]) for x in c],
                         [("C", 284.8, None)])
        al = snapfit.seed_components("Al 2p")           # no table states at all
        self.assertEqual((al[0]["name"], al[0]["be"]), ("Al", 72.7))
        self.assertAlmostEqual(al[0]["doublet"]["split"], 0.44, places=6)
        f = snapfit.seed_components("F 1s")
        self.assertEqual(f[0]["be"], 686.0)
        self.assertEqual(snapfit.seed_components("Xx 9z"), [])

    def test_the_seed_moves_to_the_measured_peak_but_only_nearby(self):
        be = np.linspace(292.0, 280.0, 241)
        y = 20 + 500 * np.exp(-0.5 * ((be - 285.5) / 0.7) ** 2)
        c = snapfit.seed_components("C 1s", be=be, y=y)
        self.assertAlmostEqual(c[0]["be"], 285.5, delta=0.1)
        far = 20 + 500 * np.exp(-0.5 * ((be - 290.9) / 0.4) ** 2)  # 6 eV away
        self.assertEqual(snapfit.seed_components("C 1s", be=be, y=far)[0]["be"],
                         284.8)
        # a window with no peak in it (only the slope of one) is not a position
        self.assertEqual(snapfit.peak_position(be, y, 284.8, span=0.5), 284.8)

    def test_the_library_lists_each_state_once_by_binding_energy(self):
        lib = snapfit.library_states("C 1s")
        names = [e["name"] for e in lib]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual([e["be"] for e in lib], sorted(e["be"] for e in lib))
        self.assertIn("C-C", names)
        self.assertEqual([e["name"] for e in snapfit.library_states("Si 2p")],
                         ["Si", "SiO2"])

    def test_a_window_surrounds_the_components_and_stays_in_the_data(self):
        be = np.linspace(110.0, 95.0, 301)
        comps = snapfit.seed_components("Si 2p")
        lo, hi = snapfit.suggest_window(be, comps)
        self.assertLess(lo, 99.4 - 0.5)
        self.assertGreater(hi, 103.4 + 0.5 + 0.63)
        self.assertLessEqual(hi, 110.0)
        narrow = snapfit.suggest_window(np.linspace(105, 98, 50), comps)
        self.assertEqual((narrow[0], narrow[1]), (98.0, 105.0))

    def test_positions_may_move_half_an_electronvolt_by_default(self):
        for c in snapfit.seed_components("Si 2p"):
            self.assertEqual(c["be_range"], 0.5)
        self.assertEqual(snapfit.new_component("X", 100.0)["be_range"], 0.5)
        self.assertEqual(snapfit.sanitise({"components": [
            {"be": 100.0, "be_range": -1}]})["components"][0]["be_range"], 0.5)

    def test_a_level_without_states_or_lines_is_just_empty_or_singlet(self):
        self.assertEqual(snapfit.seed_components("Xx 9z"), [])
        s_only = snapfit.seed_components("Si 2s", states=[
            {"core_level": "Si 2s", "state": "Si2s Si", "be": 150.0,
             "fwhm": 1.0}])
        self.assertEqual(len(s_only), 1)
        self.assertIsNone(s_only[0]["doublet"])

    def test_the_doublet_comes_from_the_given_line_table(self):
        lines = [{"el": "Pt", "line": "4f7/2", "be": 71.0},
                 {"el": "Pt", "line": "4f5/2", "be": 74.33}]
        d = snapfit.doublet_of_level("Pt 4f", lines)
        self.assertAlmostEqual(d["split"], 3.33)
        self.assertAlmostEqual(d["ratio"], 0.75)
        self.assertIsNone(snapfit.doublet_of_level("Pt 4f", []))


@unittest.skipUnless(HAVE, "numpy not available")
class TestSolvers(unittest.TestCase):
    def test_nnls_matches_least_squares_when_the_answer_is_positive(self):
        rng = np.random.RandomState(1)
        A = rng.rand(30, 4)
        x = np.array([1.0, 2.0, 0.5, 3.0])
        got = snapfit.nnls(A, A @ x)
        self.assertTrue(np.allclose(got, x, atol=1e-8))

    def test_nnls_clamps_what_would_go_negative(self):
        rng = np.random.RandomState(2)
        A = rng.rand(30, 3)
        b = A @ np.array([2.0, -1.5, 1.0])
        x = snapfit.nnls(A, b)
        self.assertTrue((x >= 0).all())
        free = np.linalg.lstsq(A, b, rcond=None)[0]
        self.assertGreater(np.linalg.norm(A @ x - b) + 1e-12,
                           np.linalg.norm(A @ free - b) - 1e-9)
        # KKT: the gradient is non-positive where x is zero, zero where positive
        g = A.T @ (b - A @ x)
        self.assertTrue((g[x > 1e-9] < 1e-7).all())
        self.assertTrue((g[x <= 1e-9] <= 1e-7).all())

    def test_nnls_many_agrees_with_one_at_a_time(self):
        rng = np.random.RandomState(3)
        A = rng.rand(25, 3)
        Y = A @ rng.randn(3, 40) + 0.05 * rng.randn(25, 40)
        X = snapfit.nnls_many(A, Y)
        self.assertTrue((X >= -1e-12).all())
        for j in range(Y.shape[1]):
            self.assertTrue(np.allclose(X[:, j], snapfit.nnls(A, Y[:, j]),
                                        atol=1e-7))


@unittest.skipUnless(HAVE, "numpy not available")
class TestModel(unittest.TestCase):
    BE = np.linspace(110.0, 95.0, 301) if HAVE else None

    def test_a_state_sums_to_one_count_over_its_channels(self):
        c = snapfit.new_component("Si", 99.4, 1.0, doublet={"split": 0.61,
                                                             "ratio": 0.5})
        col = snapfit.design(self.BE, 1486.6, snapfit.sanitise(
            {"components": [c]}))[:, 0]
        self.assertAlmostEqual(float(col.sum()), 1.0, places=2)

    def test_the_second_line_of_a_doublet_sits_higher_and_weaker(self):
        c = snapfit.new_component("Si", 99.4, 0.6, doublet={"split": 0.61,
                                                             "ratio": 0.5})
        prof = snapfit.component_profile(self.BE, 1486.6, c)
        i1 = int(np.argmin(abs(self.BE - 99.4)))
        i2 = int(np.argmin(abs(self.BE - 100.01)))
        self.assertGreater(prof[i1], 1.5 * prof[i2])         # 2 : 1 areas
        self.assertGreater(prof[i2], 3 * prof[int(np.argmin(
            abs(self.BE - 100.7)))])                          # a real second peak

    def test_the_shirley_step_gives_the_codes_shirley_background(self):
        be = np.linspace(110.0, 95.0, 301)
        y = 30 + 400 * np.exp(-0.5 * ((be - 102.0) / 0.8) ** 2)
        y += np.where(be < 102.0, 60 * (1 - (be - 95) / 7.0), 0)
        cols = snapfit.background_columns("shirley", be, y)
        hi_be = float(y[:3].mean())                          # high-BE end
        lo_be = float(y[-3:].mean())
        mine = cols @ np.array([lo_be, hi_be - lo_be])     # base + step
        order = np.argsort(hi_be * 0 + (1486.6 - be))        # ascending KE
        ref = lineshapes.shirley(y[order], avg=3)
        self.assertTrue(np.allclose(mine[order], ref, atol=1e-6 * y.max()))

    def test_the_end_averaging_is_a_validated_setting_that_shapes_the_step(self):
        self.assertEqual(snapfit.sanitise({})["avg"], 3)
        for bad in (0, 11, "x", True, None, -2):
            self.assertEqual(snapfit.sanitise({"avg": bad})["avg"], 3, bad)
        self.assertEqual(snapfit.sanitise({"avg": 7})["avg"], 7)
        s = snapfit.sanitise({"avg": 5, "components": [{"be": 100.0}]})
        self.assertEqual(snapfit.sanitise(json.loads(json.dumps(s)))["avg"], 5)
        # a noisy end moves the Shirley step with the number of channels averaged
        rng = np.random.RandomState(4)
        be = np.linspace(110.0, 95.0, 301)
        y = 30 + 400 * np.exp(-0.5 * ((be - 102.0) / 0.8) ** 2)
        y = y + rng.normal(0, 25, y.size)
        a1 = snapfit.background_columns("shirley", be, y, avg=1)[:, 1]
        a3 = snapfit.background_columns("shirley", be, y, avg=3)[:, 1]
        a9 = snapfit.background_columns("shirley", be, y, avg=9)[:, 1]
        self.assertFalse(np.allclose(a1, a9))
        self.assertTrue(np.allclose(
            a3, snapfit.background_columns("shirley", be, y)[:, 1]))   # 3 is the default
        # and it reaches the design through the spec
        comps = [snapfit.new_component("S", 102.0, 1.0)]
        d3 = snapfit.design(be, 1486.6, snapfit.sanitise({"components": comps}),
                            None, y)
        d9 = snapfit.design(be, 1486.6, snapfit.sanitise(
            {"components": comps, "avg": 9}), None, y)
        self.assertTrue(np.allclose(d3[:, :2], d9[:, :2]))
        self.assertFalse(np.allclose(d3[:, 2:], d9[:, 2:]))

    def test_a_linear_background_joins_the_end_levels(self):
        be = np.linspace(110.0, 95.0, 11)
        cols = snapfit.background_columns("linear", be, np.zeros(11))
        b = cols @ np.array([10.0, 30.0])         # low-BE level, high-BE level
        self.assertAlmostEqual(float(b[-1]), 10.0)
        self.assertAlmostEqual(float(b[0]), 30.0)
        self.assertAlmostEqual(float(b[5]), 20.0)


def two_state_map(nx=40, ny=30, seed=5, scale=1.0, hv=1486.6):
    """A Si 2p map: elemental Si strong on the left, oxide on the right (with a
    gradient), both spin-orbit doublets, on a Shirley-like step, Poisson noise.
    Returns ``(cube, truth)``; ``truth`` has the areas and the be / fwhm used."""
    be = np.linspace(110.0, 95.0, 121)
    spec = snapfit.sanitise({"components": [
        snapfit.new_component("Si", 99.4, 0.85, doublet={"split": 0.61,
                                                          "ratio": 0.5}),
        snapfit.new_component("SiO2", 103.4, 1.5, doublet={"split": 0.61,
                                                            "ratio": 0.5})],
        "window": [95.0, 110.0]})
    step = abs(be[1] - be[0])
    prof = [snapfit.component_profile(be, hv, c) * step
            for c in spec["components"]]
    xs = np.linspace(0, 1, nx)[None, :] * np.ones((ny, 1))
    a_si = 900.0 * scale * (1.0 - xs) ** 1.2
    a_ox = 700.0 * scale * (0.2 + 0.8 * xs)
    nominal = prof[0] * 450 + prof[1] * 450 + 10.0   # the mean map, as the fit sees it
    bg_step = snapfit.background_columns("shirley", be, nominal)
    rng = np.random.RandomState(seed)
    blocks = []
    for iy in range(ny):
        for ix in range(nx):
            m = a_si[iy, ix] * prof[0] + a_ox[iy, ix] * prof[1]
            m = m + bg_step @ np.array([6.0, 8.0]) * scale
            blocks.append(((ix, iy), rng.poisson(m).astype(float)))
    cube = snapmap.build(list(be), nx, ny, 0.0, 10.0, 0.0, 10.0, blocks)
    return cube, {"a_si": a_si, "a_ox": a_ox, "spec": spec, "be": be}


@unittest.skipUnless(HAVE, "numpy not available")
class TestNoiseGain(unittest.TestCase):
    def test_integer_counts_have_a_gain_of_one_and_scaled_ones_of_the_scale(self):
        cube, truth = two_state_map()
        self.assertAlmostEqual(snapfit.noise_gain(cube), 1.0, delta=0.15)
        scaled = snapmap.MapCube(
            list(cube.energy), cube.nx, cube.ny, 0.0, 10.0, 0.0, 10.0,
            snapmap.array("f", (cube.array3d() * 100.0).astype("float32").ravel()))
        self.assertAlmostEqual(snapfit.noise_gain(scaled), 100.0, delta=15.0)
        rate = snapmap.MapCube(
            list(cube.energy), cube.nx, cube.ny, 0.0, 10.0, 0.0, 10.0,
            snapmap.array("f", (cube.array3d() / 4.0).astype("float32").ravel()))
        self.assertAlmostEqual(snapfit.noise_gain(rate), 0.25, delta=0.05)

    def test_the_fit_does_not_depend_on_the_scale_and_chi2_is_about_one(self):
        cube, truth = two_state_map()
        scaled = snapmap.MapCube(
            list(cube.energy), cube.nx, cube.ny, 0.0, 10.0, 0.0, 10.0,
            snapmap.array("f", (cube.array3d() * 100.0).astype("float32").ravel()))
        be = truth["be"]
        mean = cube.array3d().reshape(-1, len(be)).mean(axis=0)
        start = {"components": [dict(c, be=c["be"] + 0.2)
                                for c in truth["spec"]["components"]],
                 "window": [95.0, 110.0]}
        a = snapfit.fit_mean(be, mean, 1486.6, start, n_pix=1200)
        g = snapfit.noise_gain(scaled)
        b = snapfit.fit_mean(be, mean * 100.0, 1486.6, start, n_pix=1200, gain=g)
        for ca, cb in zip(a.spec["components"], b.spec["components"]):
            self.assertAlmostEqual(ca["be"], cb["be"], delta=0.03)
        self.assertAlmostEqual(b.chi2red, a.chi2red, delta=0.5)
        self.assertLess(b.chi2red, 2.0)
        # without the gain the scaled data would look 100 times too noisy-free
        wrong = snapfit.fit_mean(be, mean * 100.0, 1486.6, start, n_pix=1200)
        self.assertGreater(wrong.chi2red, 20 * b.chi2red)
        maps = snapfit.fit_maps(scaled, be, 1486.6, b.spec)
        self.assertAlmostEqual(float(maps.chi2.mean()), 1.0, delta=0.3)

    def test_the_width_limit_is_reported_and_clipped_parameters_do_not_spoil_errors(self):
        cube, truth = two_state_map()
        be = truth["be"]
        mean = cube.array3d().reshape(-1, len(be)).mean(axis=0)
        s = {"components": [dict(c, fwhm_range=[0.3, 0.5])
                            for c in truth["spec"]["components"]],
             "window": [95.0, 110.0]}
        res = snapfit.fit_mean(be, mean, 1486.6, s, n_pix=1200)
        self.assertTrue(any("width at its limit" in n for n in res.notes))
        self.assertNotIn("could not be estimated", " ".join(res.notes))


@unittest.skipUnless(HAVE, "numpy not available")
class TestMeanFit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cube, cls.truth = two_state_map()
        cls.be = cls.truth["be"]
        cls.mean = cls.cube.array3d().reshape(-1, len(cls.be)).mean(axis=0)

    def start(self, **kw):
        """The true spec, started off by 0.3 eV and 0.3 eV of width."""
        comps = []
        for c in self.truth["spec"]["components"]:
            comps.append(dict(c, be=c["be"] + 0.3, fwhm=c["fwhm"] * 1.3, **kw))
        return {"components": comps, "window": [95.0, 110.0],
                "background": "shirley"}

    def test_the_positions_and_widths_come_back(self):
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, self.start(),
                               n_pix=1200)
        got = {c["name"]: c for c in res.spec["components"]}
        for c in self.truth["spec"]["components"]:
            self.assertAlmostEqual(got[c["name"]]["be"], c["be"], delta=0.06)
            self.assertAlmostEqual(got[c["name"]]["fwhm"], c["fwhm"], delta=0.12)
        self.assertLess(res.chi2red, 2.0)
        self.assertGreater(res.chi2red, 0.4)
        self.assertEqual(res.notes, [])
        self.assertGreater(res.errors[("Si", "be")], 0)

    def test_the_fitted_areas_sum_to_the_counts_in_the_map(self):
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, self.start(),
                               n_pix=1200)
        truth = (self.truth["a_si"].mean(), self.truth["a_ox"].mean())
        self.assertAlmostEqual(res.amplitudes[0] / truth[0], 1.0, delta=0.05)
        self.assertAlmostEqual(res.amplitudes[1] / truth[1], 1.0, delta=0.05)
        # the pieces add up to the envelope the figure draws
        self.assertTrue(np.allclose(res.components.sum(axis=0) + res.bg,
                                    res.envelope))

    def test_fixed_parameters_do_not_move(self):
        s = self.start(free_be=False, free_fwhm=False)
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, s, n_pix=1200)
        for c, c0 in zip(res.spec["components"], s["components"]):
            self.assertEqual((c["be"], c["fwhm"]), (c0["be"], c0["fwhm"]))

    def test_a_shared_width_is_one_parameter(self):
        s = dict(self.start(), shared_fwhm=True)
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, s, n_pix=1200)
        w = {c["fwhm"] for c in res.spec["components"]}
        self.assertEqual(len(w), 1)

    def test_a_charging_shift_within_half_an_ev_is_followed(self):
        # the peaks are 0.45 eV from the library values: inside the room
        s = self.start(free_fwhm=True)
        for c, shift in zip(s["components"], (0.0, 0.0)):
            c["be"] = c["be"] - 0.3 - 0.15                 # seed 0.45 eV low
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, s, n_pix=1200)
        for c, c0 in zip(res.spec["components"], self.truth["spec"]["components"]):
            self.assertAlmostEqual(c["be"], c0["be"], delta=0.06)
        self.assertFalse(any("limit" in n for n in res.notes))

    def test_a_shift_beyond_half_an_ev_stops_at_the_limit_and_says_so(self):
        s = self.start()
        for c in s["components"]:
            c["be"] = c["be"] - 0.3 - 0.45                 # seed 0.75 eV low
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, s, n_pix=1200)
        self.assertTrue(any("limit" in n for n in res.notes))

    def test_a_position_that_hits_its_limit_is_reported(self):
        s = self.start()
        s["components"][1]["be"] = 105.5        # 2 eV off, limit is 0.6
        res = snapfit.fit_mean(self.be, self.mean, 1486.6, s, n_pix=1200)
        self.assertTrue(any("limit" in n for n in res.notes))

    def test_too_few_channels_and_no_components_are_refused(self):
        with self.assertRaises(ValueError):
            snapfit.fit_mean(self.be, self.mean, 1486.6, self.start(),
                             window=[100.0, 100.5])
        with self.assertRaises(ValueError):
            snapfit.fit_mean(self.be, self.mean, 1486.6, {"components": []})


def stripe_map(nx=40, ny=30, seed=3, hv=1486.6, scale=100.0, amp=150.0):
    """Oxide everywhere, elemental Si (0.3 eV from the library and broader
    than its seed) in ONE row only (a thin stripe), counts scaled by 100 as
    the Avantage maps are."""
    be = np.linspace(110.0, 95.0, 121)
    truth = snapfit.sanitise({"components": [
        snapfit.new_component("Si", 99.7, 1.3, doublet={"split": 0.61,
                                                         "ratio": 0.5}),
        snapfit.new_component("SiO2", 103.3, 1.6, doublet={"split": 0.61,
                                                            "ratio": 0.5})],
        "window": [95.0, 110.0]})
    step = abs(be[1] - be[0])
    prof = [snapfit.component_profile(be, hv, c) * step
            for c in truth["components"]]
    rng = np.random.RandomState(seed)
    blocks = []
    for iy in range(ny):
        for ix in range(nx):
            m = (amp if iy == 13 else 0.0) * prof[0] + 600.0 * prof[1] + 9.0
            blocks.append(((ix, iy), rng.poisson(m).astype(float) * scale))
    return snapmap.build(list(be), nx, ny, 0.0, 10.0, 0.0, 10.0, blocks), truth


@unittest.skipUnless(HAVE, "numpy not available")
class TestStrongestPixels(unittest.TestCase):
    """A state present in only a few pixels barely shows in the mean of the
    map: its width is biased there, but the pixels where it is strongest
    pin it."""

    @classmethod
    def setUpClass(cls):
        cls.cube, cls.truth = stripe_map()
        cls.be = np.array(cls.cube.energy)
        cls.gain = snapfit.noise_gain(cls.cube)
        cls.mean = cls.cube.array3d().reshape(-1, 121).mean(axis=0)
        cls.seed = {"components": snapfit.seed_components("Si 2p"),
                    "window": [95.0, 110.0]}

    def test_enriched_spectra_pick_the_stripe(self):
        got = snapfit.enriched_spectra(self.cube, self.be, 1486.6, self.seed,
                                       gain=self.gain)
        names = [g[2] for g in got]
        self.assertEqual(names, ["Si"])              # the dominant oxide is left alone
        y, n, _ = got[names.index("Si")]
        self.assertGreaterEqual(n, 30)
        band = (self.be > 99) & (self.be < 100.5)
        self.assertGreater(y[band].mean(), 2 * self.mean[band].mean())

    def test_the_enriched_fit_is_closer_to_the_truth_than_the_mean_alone(self):
        extra = snapfit.enriched_spectra(self.cube, self.be, 1486.6, self.seed,
                                         gain=self.gain)
        mean_only = snapfit.fit_mean(self.be, self.mean, 1486.6, self.seed,
                                     n_pix=1200, gain=self.gain)
        both = snapfit.fit_multi(
            self.be, [(self.mean, 1200)] + [(y, n) for y, n, _ in extra],
            1486.6, self.seed, gain=self.gain)
        w = lambda r: r.spec["components"][0]["fwhm"]
        self.assertLess(abs(w(both) - 1.3), abs(w(mean_only) - 1.3))
        got = {c["name"]: c for c in both.spec["components"]}
        self.assertAlmostEqual(got["Si"]["be"], 99.7, delta=0.08)
        self.assertAlmostEqual(got["Si"]["fwhm"], 1.3, delta=0.1)
        self.assertAlmostEqual(got["SiO2"]["be"], 103.3, delta=0.05)
        self.assertEqual(both.notes, [])
        maps = snapfit.fit_maps(self.cube, self.be, 1486.6, both.spec,
                                gain=self.gain)
        stripe = maps.areas[0][13].mean()
        other = np.concatenate([maps.areas[0][:13].ravel(),
                                maps.areas[0][14:].ravel()]).mean()
        self.assertGreater(stripe, 8 * max(other, 1.0))

    def test_a_state_that_is_nowhere_prominent_adds_no_spectrum(self):
        cube, truth = two_state_map()                 # both states everywhere
        be = np.array(cube.energy)
        spec = {"components": [snapfit.new_component("Ghost", 105.5, 1.0)] +
                list(truth["spec"]["components"]), "window": [95.0, 110.0]}
        got = snapfit.enriched_spectra(cube, be, 1486.6, spec)
        self.assertNotIn("Ghost", [g[2] for g in got])

    def test_one_spectrum_through_fit_multi_is_fit_mean(self):
        cube, truth = two_state_map()
        be = truth["be"]
        mean = cube.array3d().reshape(-1, len(be)).mean(axis=0)
        s = {"components": [dict(c, be=c["be"] + 0.2)
                            for c in truth["spec"]["components"]],
             "window": [95.0, 110.0]}
        a = snapfit.fit_mean(be, mean, 1486.6, s, n_pix=1200)
        b = snapfit.fit_multi(be, [(mean, 1200)], 1486.6, s)
        for ca, cb in zip(a.spec["components"], b.spec["components"]):
            self.assertEqual((ca["be"], ca["fwhm"]), (cb["be"], cb["fwhm"]))


@unittest.skipUnless(HAVE, "numpy not available")
class TestMaps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cube, cls.truth = two_state_map()
        cls.be = cls.truth["be"]
        cls.spec = cls.truth["spec"]               # the true shapes, frozen

    def test_the_partition_recovers_the_amplitude_maps(self):
        maps = snapfit.fit_maps(self.cube, self.be, 1486.6, self.spec)
        self.assertEqual(maps.names, ["Si", "SiO2"])
        self.assertEqual(maps.areas.shape, (2, 30, 40))
        for k, key in enumerate(("a_si", "a_ox")):
            truth = self.truth[key]
            got = maps.areas[k]
            r = np.corrcoef(truth.ravel(), got.ravel())[0, 1]
            self.assertGreater(r, 0.95, key)
            self.assertAlmostEqual(float(got.mean() / truth.mean()), 1.0,
                                   delta=0.05)
        self.assertTrue((maps.areas >= 0).all())
        self.assertAlmostEqual(float(maps.chi2.mean()), 1.0, delta=0.25)

    def test_the_oxide_fraction_follows_the_true_fraction(self):
        maps = snapfit.fit_maps(self.cube, self.be, 1486.6, self.spec)
        f = snapfit.fraction(maps, ["SiO2"])
        truth = self.truth["a_ox"] / (self.truth["a_ox"] + self.truth["a_si"])
        ok = np.isfinite(f)
        self.assertGreater(ok.mean(), 0.9)
        self.assertLess(float(np.abs(f[ok] - truth[ok]).mean()), 0.05)
        self.assertLess(float(f[:, :5].mean()), 0.3)           # elemental side
        self.assertGreater(float(f[:, -5:].mean()), 0.8)       # oxide side

    def test_binning_reduces_the_noise_without_moving_the_mean(self):
        raw = snapfit.fit_maps(self.cube, self.be, 1486.6, self.spec)
        binned = snapfit.fit_maps(self.cube, self.be, 1486.6,
                                  dict(self.spec, bin=3))
        err = lambda m: float(np.abs(m.areas[0] - self.truth["a_si"]).mean())
        self.assertLess(err(binned), err(raw))
        self.assertAlmostEqual(float(binned.areas[0].mean() / raw.areas[0].mean()),
                               1.0, delta=0.03)
        self.assertTrue(any("binned" in n for n in binned.notes))

    def test_a_window_with_no_data_is_refused_not_a_crash(self):
        spec = dict(self.spec, window=[60.0, 70.0])           # the cube is 95-110
        with self.assertRaises(ValueError):
            snapfit.fit_maps(self.cube, self.be, 1486.6, spec)
        spec = dict(self.spec, window=[100.0, 100.4])         # a few channels
        with self.assertRaises(ValueError):
            snapfit.fit_maps(self.cube, self.be, 1486.6, spec)

    def test_the_cube_is_not_touched(self):
        before = bytes(self.cube.data)
        snapfit.fit_maps(self.cube, self.be, 1486.6, dict(self.spec, bin=2))
        self.assertEqual(bytes(self.cube.data), before)

    def test_a_blank_total_is_not_a_fraction(self):
        maps = snapfit.fit_maps(self.cube, self.be, 1486.6, self.spec)
        maps.areas[:, 3:6, 3:6] = 0.0
        maps.total = maps.areas.sum(axis=0)
        f = snapfit.fraction(maps, ["SiO2"])
        self.assertTrue(np.isnan(f[3:6, 3:6]).all())

    def test_the_energy_axis_may_be_ascending(self):
        asc = snapmap.MapCube(list(self.be[::-1]), self.cube.nx, self.cube.ny,
                              0.0, 10.0, 0.0, 10.0,
                              snapmap.array("f", self.cube.array3d()[:, :, ::-1]
                                            .astype("float32").ravel()))
        a = snapfit.fit_maps(self.cube, self.be, 1486.6, self.spec)
        b = snapfit.fit_maps(asc, self.be[::-1], 1486.6, self.spec)
        self.assertTrue(np.allclose(a.areas, b.areas, rtol=1e-3, atol=1e-3))

    def test_a_noisy_map_still_separates_with_the_fitted_shapes(self):
        cube, truth = two_state_map(scale=0.25, seed=9)
        mean = cube.array3d().reshape(-1, len(truth["be"])).mean(axis=0)
        start = {"components": [dict(c, be=c["be"] + 0.25)
                                for c in truth["spec"]["components"]],
                 "window": [95.0, 110.0]}
        res = snapfit.fit_mean(truth["be"], mean, 1486.6, start,
                               n_pix=1200)
        maps = snapfit.fit_maps(cube, truth["be"], 1486.6, res.spec,
                                y_mean=mean)
        f = snapfit.fraction(maps, ["SiO2"])
        self.assertGreater(float(np.nanmean(f[:, -5:])),
                           float(np.nanmean(f[:, :5])) + 0.4)


def level_map(centres, widths, patterns, lo, hi, n=161, nx=30, ny=24, seed=3,
              bg=8.0, hv=1486.6):
    """A map of one core level with singlet states at ``centres`` whose
    spatial amplitude maps are ``patterns`` (functions of x in 0..1)."""
    be = np.linspace(hi, lo, n)
    spec = snapfit.sanitise({"components": [
        snapfit.new_component(f"S{i}", c, w) for i, (c, w) in
        enumerate(zip(centres, widths))], "window": [lo, hi]})
    step = abs(be[1] - be[0])
    prof = [snapfit.component_profile(be, hv, c) * step
            for c in spec["components"]]
    xs = np.linspace(0, 1, nx)
    rng = np.random.RandomState(seed)
    base = snapfit.background_columns(
        "shirley", be, sum(p * 300 for p in prof) + 5.0) @ np.array([bg, bg * 0.8])
    blocks, truth = [], np.zeros((len(centres), ny, nx))
    for iy in range(ny):
        for ix in range(nx):
            amps = [f(xs[ix]) for f in patterns]
            truth[:, iy, ix] = amps
            m = sum(a * p for a, p in zip(amps, prof)) + base
            blocks.append(((ix, iy), rng.poisson(m).astype(float)))
    return snapmap.build(list(be), nx, ny, 0.0, 10.0, 0.0, 10.0, blocks), truth, spec


@unittest.skipUnless(HAVE, "numpy not available")
class TestOtherLevels(unittest.TestCase):
    """Nothing in the method is Si-specific: any core level, one to several
    states, seeded from the line table or by hand."""

    def test_c_1s_three_states_from_one_seed_and_two_added(self):
        cube, truth, spec = level_map(
            (284.8, 286.4, 288.5), (1.1, 1.2, 1.2),
            (lambda x: 500 * (1 - 0.7 * x), lambda x: 150 + 200 * x,
             lambda x: 60 + 120 * x ** 2), 281.0, 292.0)
        be = np.array(cube.energy)
        mean = cube.array3d().reshape(-1, len(be)).mean(axis=0)
        seed = snapfit.seed_components("C 1s", be=be, y=mean)      # C-C only
        self.assertEqual(len(seed), 1)
        comps = seed + [snapfit.new_component("C-O", 286.2, 1.2),
                        snapfit.new_component("C=O", 288.3, 1.2)]
        start = {"components": comps, "window": [281.0, 292.0]}
        res = snapfit.fit_mean(be, mean, 1486.6, start, n_pix=720)
        got = [c["be"] for c in res.spec["components"]]
        for g, want in zip(got, (284.8, 286.4, 288.5)):
            self.assertAlmostEqual(g, want, delta=0.12)
        maps = snapfit.fit_maps(cube, be, 1486.6, res.spec, y_mean=mean)
        for k, floor in enumerate((0.9, 0.9, 0.7)):     # the weak state is noisier
            r = np.corrcoef(truth[k].ravel(), maps.areas[k].ravel())[0, 1]
            self.assertGreater(r, floor, k)
        f = snapfit.fraction(maps, [res.spec["components"][2]["name"]])
        self.assertGreater(float(np.nanmean(f[:, -5:])),
                           float(np.nanmean(f[:, :5])))

    def test_a_single_state_level_gives_one_amplitude_map(self):
        cube, truth, spec = level_map((686.0,), (1.8,),
                                      (lambda x: 100 + 400 * x,), 680.0, 692.0)
        be = np.array(cube.energy)
        mean = cube.array3d().reshape(-1, len(be)).mean(axis=0)
        seed = snapfit.seed_components("F 1s", be=be, y=mean)
        self.assertAlmostEqual(seed[0]["be"], 686.0, delta=0.2)  # the measured peak
        res = snapfit.fit_mean(be, mean, 1486.6,
                               {"components": seed, "window": [680.0, 692.0]},
                               n_pix=720)
        maps = snapfit.fit_maps(cube, be, 1486.6, res.spec, y_mean=mean)
        self.assertEqual(maps.areas.shape[0], 1)
        self.assertGreater(np.corrcoef(truth[0].ravel(),
                                       maps.areas[0].ravel())[0, 1], 0.95)
        # one state is its own fraction wherever there is signal
        f = snapfit.fraction(maps, [maps.names[0]])
        self.assertTrue(np.allclose(f[np.isfinite(f)], 1.0))

    def test_an_al_2p_metal_and_oxide_pair_with_the_line_tables_doublet(self):
        seed = snapfit.seed_components("Al 2p")
        self.assertEqual(seed[0]["be"], 72.7)
        oxide = snapfit.new_component("Al2O3", 74.5, 1.3,
                                      doublet=seed[0]["doublet"])
        self.assertAlmostEqual(oxide["doublet"]["split"], 0.44, places=6)
        spec = snapfit.sanitise({"components": seed + [oxide],
                                 "window": [68.0, 80.0]})
        self.assertEqual([c["name"] for c in spec["components"]],
                         ["Al", "Al2O3"])


if __name__ == "__main__":
    unittest.main()
