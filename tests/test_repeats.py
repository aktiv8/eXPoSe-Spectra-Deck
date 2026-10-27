"""Repeat groups: mean, sd (n - 1) and n of the atomic percent over samples
that repeat one another (repeats.py), their storage in the Annotations and the
suggestions from sample names.

The statistic is XPSView / HarwellXPS's (its repeat-group table: mean, sd with
n - 1, n). Its numbers are not a reference for the composition itself: it sums
the survey's regions and the high-resolution scans into one total (O 1s twice).

Run:  python -m unittest discover tests
"""

import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import annotations  # noqa: E402
import quant  # noqa: E402
import repeats  # noqa: E402
import resultspages as rp  # noqa: E402
from test_results import row  # noqa: E402


def sample(label, regions, key=None, levels=1, kind="regions"):
    """A quantified sample; ``regions`` = [(name, rsf, area, states)]."""
    s = rp.Sample(key=key or f"f/{label}", label=label, kind=kind)
    for li in range(levels):
        entries = [{"spectrum": "sp", "row": row(n, rsf, a, st)}
                   for n, rsf, a, st in regions]
        lv = rp.Level(level=li if levels > 1 else None, entries=entries)
        lv.include = [True] * len(entries)
        lv.why = [""] * len(entries)
        lv.res = quant.normalise([e["row"] for e in entries])
        s.levels.append(lv)
    return s


class TestStats(unittest.TestCase):
    def test_xpsview_numbers(self):
        # ATH PR034 C 1s in XPSView's repeat-group table: 10.566 +- 0.337, n 3
        m, sd, n = repeats.stats([10.921, 10.251, 10.525])
        self.assertAlmostEqual(m, 10.5657, 3)
        self.assertAlmostEqual(sd, 0.337, 3)
        self.assertEqual(n, 3)

    def test_sample_sd_not_population(self):
        self.assertAlmostEqual(repeats.stats([1.0, 3.0])[1], math.sqrt(2.0))

    def test_one_value_has_no_sd(self):
        self.assertEqual(repeats.stats([5.0]), (5.0, None, 1))
        self.assertEqual(repeats.stats([]), (None, None, 0))

    def test_non_finite_values_are_ignored(self):
        self.assertEqual(repeats.stats([1.0, float("nan"), None, 3.0])[2], 2)


class TestCompute(unittest.TestCase):
    def setUp(self):
        self.a = sample("S", [("C 1s", 1.0, 100.0, ()),
                              ("O 1s", 1.0, 100.0, ())])
        self.b = sample("S - loc 1", [("C 1s", 1.0, 120.0, ()),
                                      ("O 1s", 1.0, 80.0, ())])
        self.c = sample("S - loc 2", [("C 1s", 1.0, 110.0, ()),
                                      ("O 1s", 1.0, 90.0, ())])
        self.all = [self.a, self.b, self.c]
        self.groups = {"S": {"kind": "locations",
                             "samples": [s.key for s in self.all]}}

    def test_mean_sd_n_per_region(self):
        (g,) = repeats.compute(self.all, self.groups)
        c = next(r for r in g.rows if r.region == "C 1s")
        self.assertEqual(c.n, 3)
        self.assertAlmostEqual(c.mean, (50 + 60 + 55) / 3, 6)
        self.assertAlmostEqual(c.sd, repeats.stats([50, 60, 55])[1], 6)
        self.assertEqual([lab for lab, _ in c.values],
                         ["S", "S - loc 1", "S - loc 2"])
        self.assertEqual(g.members, ["S", "S - loc 1", "S - loc 2"])

    def test_the_means_of_a_total_add_to_100(self):
        (g,) = repeats.compute(self.all, self.groups)
        self.assertAlmostEqual(sum(r.mean for r in g.rows if not r.state), 100)

    def test_a_state_only_some_members_have(self):
        a = sample("S", [("Al 2p", 1.0, 100.0,
                          [("Al2O3", 60.0), ("AlOH", 40.0)])])
        b = sample("S2", [("Al 2p", 1.0, 100.0, [("Al2O3", 100.0)])])
        (g,) = repeats.compute([a, b], {"g": {"kind": "locations",
                                              "samples": [a.key, b.key]}})
        states = {r.state: r for r in g.rows if r.state}
        self.assertEqual(states["Al2O3"].n, 2)
        self.assertEqual(states["AlOH"].n, 1)
        self.assertIsNone(states["AlOH"].sd)
        self.assertEqual(states["AlOH"].of, 2)

    def test_one_state_everywhere_adds_no_rows(self):
        a = sample("S", [("Al 2p", 1.0, 100.0, [("Al2O3", 100.0)])])
        b = sample("S2", [("Al 2p", 1.0, 100.0, [("Al2O3", 100.0)])])
        (g,) = repeats.compute([a, b], {"g": {"kind": "locations",
                                              "samples": [a.key, b.key]}})
        self.assertEqual([(r.region, r.state) for r in g.rows],
                         [("Al 2p", "")])

    def test_depth_profile_and_missing_keys_are_left_out_with_notes(self):
        prof = sample("P", [("C 1s", 1.0, 100.0, ())], levels=3)
        (g,) = repeats.compute(self.all + [prof], {"g": {
            "kind": "repeats",
            "samples": [self.a.key, self.b.key, prof.key, "f/gone"]}})
        self.assertEqual(g.members, ["S", "S - loc 1"])
        self.assertTrue(any("depth profile" in n for n in g.notes))
        self.assertTrue(any("f/gone" in n for n in g.notes))

    def test_one_member_says_so(self):
        (g,) = repeats.compute(self.all, {"g": {"kind": "locations",
                                                "samples": [self.a.key]}})
        self.assertTrue(any("Only one member" in n for n in g.notes))
        self.assertTrue(all(r.sd is None for r in g.rows))

    def test_a_region_left_out_of_a_member_is_not_counted(self):
        s = sample("X", [("C 1s", 1.0, 100.0, ()), ("O 1s", 0.0, 50.0, ())])
        (g,) = repeats.compute([s], {"g": {"kind": "locations",
                                           "samples": [s.key]}})
        self.assertEqual([r.region for r in g.rows], ["C 1s"])

    def test_csv_layout(self):
        rows = repeats.csv_rows(repeats.compute(self.all, self.groups))
        self.assertEqual(rows[0][:7], ["Group", "Kind", "Region", "State",
                                       "at % mean", "at % sd (n-1)", "n"])
        self.assertEqual(rows[1][:3], ["S", "different locations on the "
                                       "sample", "C 1s"])
        self.assertEqual(rows[1][6], "3")
        self.assertIn("S - loc 1=60.000", rows[1][7])


class TestAnnotations(unittest.TestCase):
    def test_round_trip_and_empty(self):
        a = annotations.Annotations()
        self.assertTrue(a.is_empty())
        a.set_repeat_group("S", "repeats", ["f/a", "f/b", "f/a"])
        self.assertFalse(a.is_empty())
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.repeat_groups,
                         {"S": {"kind": "repeats", "samples": ["f/a", "f/b"]}})

    def test_empty_list_removes_and_bad_input_is_ignored(self):
        a = annotations.Annotations()
        a.set_repeat_group("S", "locations", ["f/a"])
        a.set_repeat_group("S", "locations", [])
        self.assertEqual(a.repeat_groups, {})
        b = annotations.Annotations.from_json(
            {"repeat_groups": {"x": {"samples": "nope"}, "y": 3,
                               "z": {"kind": "weird", "samples": ["k"]}}})
        self.assertEqual(b.repeat_groups,
                         {"z": {"kind": "locations", "samples": ["k"]}})

    def test_older_files_load(self):
        self.assertEqual(
            annotations.Annotations.from_json({"version": 1}).repeat_groups,
            {})


class TestShortLabels(unittest.TestCase):
    def test_words_every_label_has_are_dropped(self):
        got = repeats.short_labels([
            "ATH PR035 - Spot sequence - 13 Sep 2026 03:34",
            "ATH PR035 - Spot sequence - location 1 - 13 Sep 2026 03:54",
            "ATH PR035 - Spot sequence - location 2 - 13 Sep 2026 04:17"])
        self.assertEqual(got, ["03:34", "location 1 03:54",
                               "location 2 04:17"])

    def test_edges(self):
        self.assertEqual(repeats.short_labels(["A", "B"]), ["A", "B"])
        self.assertEqual(repeats.short_labels(["Same"]), ["Same"])
        self.assertEqual(repeats.short_labels(["S x", "S"]), ["x", "S"])
        self.assertEqual(repeats.short_labels([]), [])
        self.assertEqual(len(repeats.short_labels(
            ["a" + "x" * 60, "b" + "y" * 60])[0]), 26)


class TestSuggest(unittest.TestCase):
    def test_base_label(self):
        for label, want in (
                ("ATH PR034 - Spot sequence - location 1", "ATH PR034"),
                ("Stained Pt <Pt #001a>", "Stained Pt"),
                ("Wafer 2 location 3", "Wafer 2"),
                ("Plain", "Plain")):
            self.assertEqual(repeats.base_label(label), want)

    def test_suggestions_skip_profiles_and_keep_surveys_apart(self):
        s1 = sample("A - one", [("C 1s", 1.0, 1.0, ())])
        s2 = sample("A - two", [("C 1s", 1.0, 1.0, ())])
        prof = sample("A - dp", [("C 1s", 1.0, 1.0, ())], levels=3)
        sv = sample("A - one (survey)", [("C 1s", 1.0, 1.0, ())],
                    kind="survey")
        sv2 = sample("A - two (survey)", [("C 1s", 1.0, 1.0, ())],
                     kind="survey")
        lone = sample("B", [("C 1s", 1.0, 1.0, ())])
        got = repeats.suggest([s1, s2, prof, sv, sv2, lone])
        self.assertEqual(list(got), ["A", "A (survey)"])
        self.assertEqual(got["A"]["samples"], [s1.key, s2.key])
        self.assertEqual(got["A (survey)"]["samples"], [sv.key, sv2.key])


if __name__ == "__main__":
    unittest.main()
