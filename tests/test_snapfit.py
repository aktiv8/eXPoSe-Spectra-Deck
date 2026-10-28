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
        mine = cols @ np.array([lo_be, hi_be])
        order = np.argsort(hi_be * 0 + (1486.6 - be))        # ascending KE
        ref = lineshapes.shirley(y[order], avg=3)
        self.assertTrue(np.allclose(mine[order], ref, atol=1e-6 * y.max()))

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
            m = m + bg_step @ np.array([6.0, 14.0]) * scale
            blocks.append(((ix, iy), rng.poisson(m).astype(float)))
    cube = snapmap.build(list(be), nx, ny, 0.0, 10.0, 0.0, 10.0, blocks)
    return cube, {"a_si": a_si, "a_ox": a_ox, "spec": spec, "be": be}


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


if __name__ == "__main__":
    unittest.main()
