"""Element identification: the line table, candidate lookup and auto-labelling.

Run:  python -m unittest discover tests
"""

import json
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import xpslines as xl  # noqa: E402

LINES = xl.load_lines()


def label(be, window=2.0, hv=None, split=False):
    c = xl.candidates(be, window, LINES, hv, split=split)
    return xl.label_of(c[0][1]) if c else None


class TestTable(unittest.TestCase):
    def test_bundled_table_is_sane(self):
        self.assertGreater(len(LINES), 80)
        els = {e["el"] for e in LINES}
        for el in ("C", "O", "N", "Si", "Al", "Mo", "S", "Au", "Ag", "Cu"):
            self.assertIn(el, els)
        for e in LINES:
            if "be" in e:
                # up to the deepest level any anode of the Kratos reader can
                # reach (Cr K-alpha 5414.7 eV); xpslines.reachable hides the
                # levels a given photon energy cannot excite
                self.assertTrue(9.0 < e["be"] < 5500, e)
            else:
                # Auger kinetic energy is set by the atom's own level
                # spacing, not by the exciting photon energy, so it can
                # legitimately exceed a survey's BE range (up to ~2.4 keV
                # in the bundled table).
                self.assertTrue(0 < e["ke"] < 2500, e)

    def test_missing_or_broken_file_gives_an_empty_table(self):
        self.assertEqual(xl.load_lines("/no/such/file.json"), [])
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.json")
            with open(p, "w") as fh:
                fh.write("{not json")
            self.assertEqual(xl.load_lines(p), [])
            with open(p, "w") as fh:
                json.dump({"lines": [{"el": "X", "line": "1s", "be": 5},
                                     {"el": "Y"}, "junk"]}, fh)
            self.assertEqual(len(xl.load_lines(p)), 1)


class TestCandidates(unittest.TestCase):
    def test_common_peaks(self):
        self.assertEqual(label(285.2), "C 1s")
        self.assertEqual(label(532.4), "O 1s")
        self.assertEqual(label(399.6), "N 1s")
        # a doublet is named once by default, by its components on request
        self.assertEqual(label(228.9), "Mo 3d")
        self.assertEqual(label(161.4), "S 2p")
        self.assertEqual(label(83.9), "Au 4f")
        self.assertEqual(label(228.9, split=True), "Mo 3d5/2")
        self.assertEqual(label(161.4, split=True), "S 2p3/2")
        self.assertEqual(label(83.9, split=True), "Au 4f7/2")

    def test_nothing_in_the_window(self):
        self.assertIsNone(label(700.0))
        self.assertEqual(xl.candidates(287.5, 0.1, LINES), [])

    def test_sorted_nearest_first_with_secondary_lines_penalised(self):
        c = xl.candidates(119.0, 3.0, LINES)          # Al 2s 118 / Tl 4f 118
        self.assertTrue(c)
        self.assertLessEqual(abs(c[0][0]), abs(c[-1][0]) + 1.6)

    def test_o1s_at_the_edge_of_the_window_beats_a_rare_secondary_line(self):
        # O 1s at 531.0 clicked at 533.0 (2 eV away) used to lose to a
        # rank-2 line of a rare element (At 4d3/2) sitting right on the click
        for be in (532.5, 533.0):
            self.assertEqual(label(be), "O 1s", be)
        e, y = survey([(532.9, 9000)])
        self.assertEqual([lbl for _be, lbl in xl.auto_label(e, y, LINES)],
                         ["O 1s"])

    def test_a_main_line_only_loses_at_its_own_energy_to_a_true_overlap(self):
        # the four genuine coincidences (a rare element's main line sitting
        # on a common element's): nothing else may beat a main line on its
        # own energy
        allowed = {"Pr 3d5/2", "Ho 4d", "Hg 4f7/2", "Th 4f7/2"}
        lost = set()
        for e in LINES:
            if e.get("rank", 1) == 1 and "be" in e:
                c = xl.candidates(e["be"], 2.0, LINES, 1486.6, split=True)
                if c and c[0][1] is not e:
                    lost.add(xl.label_of(e))
        self.assertLessEqual(lost, allowed)

    def test_rare_secondary_lines_are_still_offered(self):
        # the penalty reorders; it never hides a line from the list
        labels = [xl.label_of(e) for _d, e in
                  xl.candidates(533.0, 2.0, LINES, 1486.6)]
        self.assertIn("O 1s", labels)
        self.assertIn("At 4d", labels)
        split = [xl.label_of(e) for _d, e in
                 xl.candidates(533.0, 2.0, LINES, 1486.6, split=True)]
        self.assertIn("At 4d3/2", split)

    def test_auger_lines_follow_the_photon_energy(self):
        # O KL1 (ke=506): an Auger line's binding-energy position is
        # hv-dependent (BE = hv - KE), unlike every photoelectron line. The
        # richer table now has hundreds of Auger entries, so rather than
        # pin whichever label happens to rank first globally at one BE
        # (fragile as the table grows), check the specific O KL1 entry's
        # own computed BE shifts exactly with hv, and is found in the
        # candidate list near the BE it implies for a given hv.
        o_kl1 = next(e for e in LINES if e["el"] == "O" and e["line"] == "KL1")
        be_at_al_ka = xl.line_be(o_kl1, hv=1486.6)
        be_at_mg_ka = xl.line_be(o_kl1, hv=1253.6)
        self.assertAlmostEqual(be_at_al_ka - be_at_mg_ka, 1486.6 - 1253.6)
        cands = xl.candidates(be_at_al_ka, 0.5, LINES, hv=1486.6)
        self.assertIn("O KL1", [xl.label_of(e) for _d, e in cands])


class TestNearbyLines(unittest.TestCase):
    HV = 1486.6
    TABLE = [
        {"el": "X", "line": "1s", "be": 100.0, "rank": 1},
        {"el": "X", "line": "2s", "be": 101.0, "rank": 2},
        {"el": "X", "line": "3s", "be": 102.0, "rank": 3},
        {"el": "Y", "line": "KLL", "ke": 1386.6, "rank": 3},   # be = hv-ke = 100.0
    ]

    def test_excludes_the_primary_and_returns_secondary_nearest_first(self):
        got = xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                              exclude="X 1s", secondary=True)
        self.assertEqual([g[1] for g in got], ["X 2s", "X 3s"])
        self.assertTrue(all(g[2] == "secondary" for g in got))

    def test_auger_only_when_asked(self):
        without = xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                                  exclude="X 1s", secondary=True)
        self.assertNotIn("Y KLL", [g[1] for g in without])
        got = xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                              exclude="X 1s", auger=True)
        self.assertEqual([(g[1], g[2]) for g in got], [("Y KLL", "auger")])

    def test_max_extra_caps_secondary_not_auger(self):
        got = xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                              exclude="X 1s", secondary=True, auger=True,
                              max_extra=1)
        secondary = [g for g in got if g[2] == "secondary"]
        self.assertEqual([g[1] for g in secondary], ["X 2s"])
        self.assertIn("Y KLL", [g[1] for g in got])

    def test_neither_flag_gives_nothing(self):
        self.assertEqual(
            xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                            exclude="X 1s"), [])

    def test_be_of_each_candidate(self):
        got = xl.nearby_lines(100.0, 3.0, self.TABLE, hv=self.HV,
                              exclude="X 1s", secondary=True, auger=True)
        by_label = {g[1]: g[0] for g in got}
        self.assertAlmostEqual(by_label["X 2s"], 101.0)
        self.assertAlmostEqual(by_label["Y KLL"], 100.0)


def survey(peaks, lo=0.0, hi=1100.0, step=0.5, width=1.2, noise=0.0):
    n = int((hi - lo) / step) + 1
    e = [hi - i * step for i in range(n)]
    y = [200.0 for _ in e]
    for be, amp in peaks:
        y = [v + amp * math.exp(-((x - be) / width) ** 2)
             for v, x in zip(y, e)]
    return e, y


class TestPeaks(unittest.TestCase):
    def test_finds_the_peaks_strongest_first(self):
        e, y = survey([(285.0, 3000), (532.0, 9000), (99.0, 1500)])
        peaks = xl.find_peaks(e, y)
        self.assertEqual([round(p[0]) for p in peaks][:3], [532, 285, 99])

    def test_close_peaks_are_merged_and_noise_is_ignored(self):
        e, y = survey([(285.0, 3000), (286.5, 2900)])
        near = [p for p in xl.find_peaks(e, y, min_sep=3.0)
                if 280 < p[0] < 292]
        self.assertEqual(len(near), 1)
        import random
        rnd = random.Random(2)
        e, y = survey([(532.0, 9000)])
        y = [v + rnd.gauss(0, 20) for v in y]
        self.assertEqual(len(xl.find_peaks(e, y)), 1)

    def test_flat_or_tiny_input(self):
        self.assertEqual(xl.find_peaks([1, 2], [1, 1]), [])
        self.assertEqual(xl.find_peaks(list(range(50)), [5.0] * 50), [])

    def test_auto_label(self):
        e, y = survey([(285.3, 3000), (532.6, 9000), (74.0, 800),
                       (228.5, 2500)])
        got = dict((lbl, be) for be, lbl in xl.auto_label(e, y, LINES))
        for lbl, be in (("O 1s", 532.6), ("C 1s", 285.3), ("Al 2p", 74.0),
                        ("Mo 3d", 228.5)):
            self.assertIn(lbl, got)
            self.assertAlmostEqual(got[lbl], be, delta=0.5)
        got = dict((lbl, be) for be, lbl in
                   xl.auto_label(e, y, LINES, split=True))
        for lbl, be in (("Al 2p3/2", 74.0), ("Mo 3d5/2", 228.5)):
            self.assertIn(lbl, got)
            self.assertAlmostEqual(got[lbl], be, delta=0.5)


class TestBundledTable(unittest.TestCase):
    """What the KherveFitting-derived additions must keep true."""

    def test_keys_are_unique_and_every_row_is_valid(self):
        keys = [(e["el"], e["line"]) for e in LINES]
        self.assertEqual(len(keys), len(set(keys)))
        for e in LINES:
            self.assertIn(e.get("rank", 1), (1, 2, 3), e)
            self.assertIn(e.get("src", "avantage"), ("avantage", "orange"))

    def test_added_rows_never_override_and_gases_keep_their_solid_state_value(self):
        self.assertEqual([e["be"] for e in LINES
                          if (e["el"], e["line"]) == ("O", "1s")], [531.0])
        self.assertEqual([e["be"] for e in LINES
                          if (e["el"], e["line"]) == ("N", "1s")], [400.0])
        self.assertFalse([e for e in LINES
                          if e.get("src") == "orange" and e["el"] in ("H", "He")])

    def test_each_element_has_one_main_line_reachable_with_al_k_alpha(self):
        from collections import Counter
        main = Counter(e["el"] for e in LINES if "be" in e
                       and e["be"] < 1481.6 and e.get("rank", 1) == 1)
        have = {e["el"] for e in LINES if "be" in e and e["be"] < 1481.6}
        self.assertEqual(have - set(main), set())
        self.assertEqual([el for el, n in main.items() if n > 1], [])

    def test_elements_the_avantage_library_lacks_are_present(self):
        els = {e["el"] for e in LINES}
        for el in ("Ac", "Fr", "Pa", "Po"):
            self.assertIn(el, els)

    def test_the_lanthanide_levels_that_were_missing_are_offered(self):
        got = {(e["el"], e["line"]) for e in LINES}
        for k in (("Gd", "4d5/2"), ("Er", "4p3/2"), ("Ce", "5p1/2")):
            self.assertIn(k, got)


class TestPhotonReach(unittest.TestCase):
    def test_a_level_beyond_the_photon_is_not_offered(self):
        # a deep Pt level (3d5/2, 2122 eV) exists only for harder sources
        deep = next(e for e in LINES if (e["el"], e["line"]) == ("Pt", "3d5/2"))
        self.assertGreater(deep["be"], 2000)
        for hv in (1253.6, 1486.6):
            got = [xl.label_of(e) for _d, e in
                   xl.candidates(deep["be"], 5.0, LINES, hv, split=True)]
            self.assertNotIn("Pt 3d5/2", got, hv)
        got = [xl.label_of(e) for _d, e in
               xl.candidates(deep["be"], 5.0, LINES, 2984.2, split=True)]
        self.assertIn("Pt 3d5/2", got)

    def test_mg_k_alpha_loses_what_al_k_alpha_just_reaches(self):
        entry = {"el": "X", "line": "1s", "be": 1300.0, "rank": 1}
        self.assertTrue(xl.reachable(entry, 1486.6))
        self.assertFalse(xl.reachable(entry, 1253.6))
        self.assertTrue(xl.reachable(dict(entry, be=1248.0), 1253.6))
        self.assertFalse(xl.reachable(dict(entry, be=1249.0), 1253.6))   # 5 eV margin

    def test_auger_lines_are_not_filtered_by_reach(self):
        self.assertTrue(xl.reachable({"el": "O", "line": "KL1", "ke": 506}, 1253.6))
        self.assertTrue(xl.reachable({"el": "O", "line": "KL1", "ke": 506}, None))


class TestDoublets(unittest.TestCase):
    def test_split_line_and_base_label(self):
        self.assertEqual(xl.split_line("2p3/2"), ("2p", "3/2"))
        self.assertEqual(xl.split_line("4f7/2"), ("4f", "7/2"))
        self.assertEqual(xl.split_line("3p"), ("3p", ""))
        self.assertEqual(xl.split_line("KL1"), ("KL1", ""))
        self.assertEqual(xl.base_label("Ti 2p3/2"), "Ti 2p")
        self.assertEqual(xl.base_label("Ti 2p"), "Ti 2p")
        self.assertEqual(xl.base_label("C 1s"), "C 1s")

    def test_a_pair_is_one_candidate_named_by_the_pair(self):
        c = xl.candidates(454.3, 3.0, LINES, 1486.6)
        ti = [e for _d, e in c if e["el"] == "Ti"]
        self.assertEqual([xl.label_of(e) for e in ti], ["Ti 2p"])
        self.assertEqual(ti[0]["component"], "2p3/2")
        both = {xl.label_of(e) for _d, e in
                xl.candidates(457.0, 8.0, LINES, 1486.6, split=True)
                if e["el"] == "Ti"}
        self.assertTrue({"Ti 2p3/2", "Ti 2p1/2"} <= both)

    def test_the_weak_component_still_finds_the_pair_ranked_by_the_strong_one(self):
        # 2p1/2 is a rank-2 line; the pair is rank 1, so a click on the
        # 460.2 eV peak finds Ti 2p first, as a click on the 2p3/2 does
        self.assertEqual(label(460.2), "Ti 2p")
        self.assertEqual(label(460.2, split=True), "Ti 2p1/2")
        e = next(e for _d, e in xl.candidates(460.2, 1.0, LINES, 1486.6)
                 if e["el"] == "Ti")
        self.assertEqual(e["rank"], 1)

    def test_unsplit_entries_and_other_lines_are_untouched(self):
        got = [xl.label_of(e) for _d, e in
               xl.candidates(33.4, 0.3, LINES, 1486.6)]
        self.assertIn("Ti 3p", got)                    # no j in the table
        self.assertEqual(label(285.0), "C 1s")
        by_id = {id(e) for e in LINES}
        for _d, e in xl.candidates(285.0, 2.0, LINES, 1486.6):
            if e["line"] in ("1s",):
                self.assertIn(id(e), by_id)               # not a copy

    def test_auto_label_names_a_doublet_once(self):
        e, y = survey([(454.0, 9000), (460.2, 4500)])
        names = [lbl for _be, lbl in xl.auto_label(e, y, LINES)]
        self.assertEqual(names.count("Ti 2p"), 1)
        split = [lbl for _be, lbl in xl.auto_label(e, y, LINES, split=True)]
        self.assertEqual(sorted(split), ["Ti 2p1/2", "Ti 2p3/2"])

    def test_nearby_lines_leave_out_the_markers_own_doublet(self):
        for exclude in ("Ti 2p", "Ti 2p3/2", "Ti 2p1/2"):
            for split in (False, True):
                got = xl.nearby_lines(457.0, 8.0, LINES, 1486.6, exclude,
                                      secondary=True, auger=False,
                                      max_extra=50, split=split)
                self.assertFalse([g for g in got if g[1].startswith("Ti 2p")],
                                 (exclude, split))


if __name__ == "__main__":
    unittest.main()
