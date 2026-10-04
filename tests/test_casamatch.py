"""CasaXPS's own quantification tied to the ticked fit regions (casamatch)."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import casamatch as cm  # noqa: E402
import quant  # noqa: E402

REAL = r"D:\Temp\for claude files\PtCl2"


def fit_row(region, comps=(), lo=None, hi=None):
    """A ``quant.fit_rows``-shaped row: ``comps`` are component positions."""
    return {"region": region, "rsf": 1.0, "area": 100.0, "area_t": None,
            "be_lo": lo, "be_hi": hi, "source": "high-res",
            "components": [{"name": region, "be": b, "area": 10.0, "rsf": 1.0}
                           for b in comps]}


def entries(*rows):
    return [{"spectrum": r["region"], "row": r} for r in rows]


def file_rows(*rows):
    return [{"name": n, "position": p, "at_pct": a} for n, p, a in rows]


class TestNames(unittest.TestCase):
    def test_a_component_label_belongs_to_its_region(self):
        names = ["Cl 2p", "Pt 4f", "C 1s"]
        self.assertEqual(cm.region_of("Cl 2p x", names), "Cl 2p")
        self.assertEqual(cm.region_of("Cl 2p #", names), "Cl 2p")
        self.assertEqual(cm.region_of("Cl 2p ", names), "Cl 2p")
        self.assertEqual(cm.region_of("Pt 4f loss", names), "Pt 4f")
        self.assertIsNone(cm.region_of("Pt 4d", names))
        self.assertIsNone(cm.region_of("Cl 2pp", names))

    def test_the_longest_region_name_wins(self):
        self.assertEqual(cm.region_of("Pt 4f7", ["Pt 4f", "Pt 4f7"]), "Pt 4f7")

    def test_a_region_listed_again_starts_a_new_block(self):
        rows = file_rows(("C 1s", 285.0, 1), ("C 1s", 286.0, 2),
                         ("O 1s", 532.0, 3), ("C 1s", 285.1, 4))
        got = cm.blocks(rows, ["C 1s", "O 1s"])
        self.assertEqual([(r, len(b)) for r, b in got],
                         [("C 1s", 2), ("O 1s", 1), ("C 1s", 1)])


class TestTagRegions(unittest.TestCase):
    def test_rows_are_summed_per_region_and_shared_out(self):
        es = entries(fit_row("C 1s", [285.0], 280, 290),
                     fit_row("Pt 4f", [73.2, 76.5], 60, 90))
        notes = cm.tag_regions(es, file_rows(
            ("C 1s", 284.9, 20.0), ("Pt 4f", 73.2, 30.0),
            ("Pt 4f", 76.5, 20.0), ("Pt 4f loss", 78.4, 30.0)))
        self.assertEqual([e["row"]["casa_pct"] for e in es], [20.0, 80.0])
        res = quant.normalise([e["row"] for e in es])
        self.assertAlmostEqual(res[0]["at_pct"], 20.0)
        self.assertAlmostEqual(res[1]["at_pct"], 80.0)
        self.assertEqual([n for n in notes if "component" in n], [
            n for n in notes])           # 3 file rows vs 2 components: said

    def test_only_the_ticked_regions_count_and_the_rest_renormalise(self):
        es = entries(fit_row("C 1s", [285.0]))          # Pt 4d not ticked
        cm.tag_regions(es, file_rows(("C 1s", 285.0, 10.0),
                                     ("Pt 4d", 316.0, 30.0)))
        res = quant.normalise([e["row"] for e in es])
        self.assertAlmostEqual(res[0]["at_pct"], 100.0)

    def test_a_region_without_a_file_row_is_left_out(self):
        es = entries(fit_row("C 1s", [285.0]), fit_row("N 1s", [400.0]))
        cm.tag_regions(es, file_rows(("C 1s", 285.0, 10.0)))
        self.assertEqual(es[1]["row"]["casa_why"], cm.WHY_NONE)
        res = quant.normalise([e["row"] for e in es])
        self.assertIsNone(res[1]["at_pct"])
        self.assertEqual(res[1]["why"], cm.WHY_NONE)
        self.assertAlmostEqual(res[0]["at_pct"], 100.0)

    def test_no_file_rows_leaves_the_rows_alone(self):
        es = entries(fit_row("C 1s", [285.0]))
        self.assertEqual(cm.tag_regions(es, []), [])
        self.assertNotIn("casa_pct", es[0]["row"])
        self.assertNotIn("casa_why", es[0]["row"])

    def test_a_region_with_no_components_matches_by_name(self):
        es = entries(fit_row("C 1s", [], 278, 298))
        cm.tag_regions(es, file_rows(("C 1s", 284.7, 16.8)))
        self.assertAlmostEqual(es[0]["row"]["casa_pct"], 16.8)

    def test_rows_far_outside_the_fitted_region_do_not_match(self):
        es = entries(fit_row("C 1s", [], 278, 298))
        cm.tag_regions(es, file_rows(("C 1s", 330.0, 16.8)))
        self.assertIn("outside", es[0]["row"]["casa_why"])

    def test_two_acquisitions_pair_in_file_order(self):
        es = entries(fit_row("C 1s", [285.0]), fit_row("O 1s", [532.0]),
                     fit_row("C 1s", [285.1]), fit_row("O 1s", [532.1]))
        cm.tag_regions(es, file_rows(
            ("C 1s", 285.0, 1.0), ("O 1s", 532.0, 9.0),
            ("C 1s", 285.1, 2.0), ("O 1s", 532.1, 8.0)))
        self.assertEqual([e["row"]["casa_pct"] for e in es],
                         [1.0, 9.0, 2.0, 8.0])

    def test_one_block_for_two_spectra_goes_to_the_closer_fit(self):
        es = entries(fit_row("C 1s", [290.0]), fit_row("C 1s", [285.0]))
        cm.tag_regions(es, file_rows(("C 1s", 285.1, 5.0)))
        self.assertEqual(es[0]["row"]["casa_why"], cm.WHY_NONE)
        self.assertAlmostEqual(es[1]["row"]["casa_pct"], 5.0)


class TestTagSurvey(unittest.TestCase):
    def test_elements_take_the_survey_values(self):
        es = entries(fit_row("O 1s"), fit_row("C 1s"), fit_row("Cl 2p"))
        cm.tag_survey(es, [{"element": "O 1s", "pct": 1.82},
                           {"element": "C 1s", "pct": 27.46}])
        self.assertEqual(es[0]["row"]["casa_pct"], 1.82)
        self.assertEqual(es[1]["row"]["casa_pct"], 27.46)
        self.assertEqual(es[2]["row"]["casa_why"], cm.WHY_NO_SURVEY)

    def test_no_rows_leaves_the_survey_recomputed(self):
        es = entries(fit_row("O 1s"))
        cm.tag_survey(es, [])
        self.assertNotIn("casa_why", es[0]["row"])


@unittest.skipUnless(os.path.isfile(os.path.join(REAL, "quant_regions.txt")),
                     "the PtCl2 reference files are not on this machine")
class TestRealPtCl2(unittest.TestCase):
    """The real pair: the numbers must come back as CasaXPS wrote them."""

    @classmethod
    def setUpClass(cls):
        import casaquant
        import readers
        cls.doc = readers.load_file(os.path.join(REAL, "PtCl2_quantified.vms"))
        cls.cq = casaquant.load(REAL)

    def collect(self, **kw):
        import resultspages as rp
        return rp.collect([self.doc], casa_quant=self.cq, **kw)

    def test_regions_reproduce_the_regions_file(self):
        res = self.collect()
        s = next(x for x in res.samples
                 if x.label == "PtCl2" and x.kind == "regions")
        got = {e["row"]["region"]: x["at_pct"]
               for e, x in zip(s.levels[0].entries, s.levels[0].res)}
        self.assertAlmostEqual(got["C 1s"], 16.84, places=1)
        self.assertAlmostEqual(got["Cl 2p"], 33.23, places=1)
        self.assertAlmostEqual(got["Pt 4d"], 17.41, places=1)
        self.assertAlmostEqual(got["Pt 4f"], 17.54, places=1)
        self.assertAlmostEqual(got["O 1s"], 2.20, places=1)
        self.assertAlmostEqual(got["Pt 4s"], 12.77, places=1)

    def test_unticking_a_region_renormalises_the_rest(self):
        import resultspages as rp
        res = self.collect()
        s = next(x for x in res.samples
                 if x.label == "PtCl2" and x.kind == "regions")
        lv = s.levels[0]
        i = next(i for i, e in enumerate(lv.entries)
                 if e["row"]["region"] == "Pt 4d")
        key = rp.entry_key(s.key, lv.level, lv.entries, i)
        res2 = self.collect(overrides={key: False})
        lv2 = next(x for x in res2.samples
                   if x.label == "PtCl2" and x.kind == "regions").levels[0]
        total = sum(x["at_pct"] for x in lv2.res if x["at_pct"] is not None)
        self.assertAlmostEqual(total, 100.0)
        c1s = next(x for e, x in zip(lv2.entries, lv2.res)
                   if e["row"]["region"] == "C 1s")
        self.assertAlmostEqual(c1s["at_pct"], 100 * 16.84 / (100 - 17.41), places=1)

    def test_the_survey_is_its_own_total_from_the_survey_file(self):
        res = self.collect()
        s = next(x for x in res.samples
                 if x.label == "PtCl2 (survey)" and x.kind == "survey")
        got = {e["row"]["region"]: x["at_pct"]
               for e, x in zip(s.levels[0].entries, s.levels[0].res)}
        self.assertEqual(set(got), {"O 1s", "C 1s", "Cl 2p", "Pt 4f"})
        self.assertAlmostEqual(got["O 1s"], 1.82, places=1)
        self.assertAlmostEqual(got["Cl 2p"], 46.80, places=1)
        regions = next(x for x in res.samples
                       if x.label == "PtCl2" and x.kind == "regions")
        self.assertNotIn("survey", {e["row"].get("source")
                                    for e in regions.levels[0].entries})


if __name__ == "__main__":
    unittest.main()
