"""CasaXPS's "Residual StDev" (Eq. 3 of Major et al., Surf. Interface Anal.
2021): ``sqrt((1/n) * sum(((d - e) / sqrt(d))^2))`` over the region's points,
``d`` the counts, ``e`` the fitted envelope -- computed beside
``casafit.Curves.chi2_red`` and shown only when the user's toggle (the Fit
row's "Residual StDev") is on: the plot's RMS note, the Quantification table,
the PDF / slides / Word tables and the HTML page's fit header.

Run:  python -m unittest discover tests
"""
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import casafit  # noqa: E402
import resultspages as rp  # noqa: E402

try:
    import numpy as np
    HAVE_NP = True
except ImportError:                                  # pragma: no cover
    HAVE_NP = False

from test_results import level, row  # noqa: E402


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class TestTheStatistic(unittest.TestCase):
    HV = 1486.71
    DWELL, SCANS = 0.27, 25

    def setUp(self):
        from test_casafit import CASA, model_data
        self.fit = casafit.parse(CASA)
        self.be, self.counts = model_data(self.fit, self.HV, self.DWELL,
                                          self.SCANS)

    def _curves(self, counts, dwell=None, scans=None):
        return casafit.curves(self.fit, self.be, list(counts), self.HV,
                              self.DWELL if dwell is None else dwell,
                              self.SCANS if scans is None else scans)[0]

    def test_counting_noise_alone_gives_about_one(self):
        # the statistic's point: a fit that explains the data to counting
        # statistics has a Residual StDev of 1
        rng = np.random.default_rng(7)
        noisy = rng.poisson(np.array(self.counts)).astype(float)
        cv = self._curves(noisy)
        self.assertGreater(cv.residual_sd, 0.85)
        self.assertLess(cv.residual_sd, 1.15)

    def test_it_is_chi_square_over_n_instead_of_n_minus_one(self):
        rng = np.random.default_rng(3)
        noisy = rng.poisson(np.array(self.counts) * 2.0).astype(float) / 2.0
        cv = self._curves(noisy)
        n = int(np.sum(~np.isnan(np.array(cv.envelope, dtype=float))))
        self.assertAlmostEqual(cv.residual_sd ** 2,
                               cv.chi2_red * (n - 1) / n, places=9)

    def test_it_needs_true_counts_like_the_chi_square(self):
        cv = casafit.curves(self.fit, self.be, self.counts, self.HV, None, 1)[0]
        self.assertIsNone(cv.chi2_red)
        self.assertIsNone(cv.residual_sd)


class TestTheTable(unittest.TestCase):
    W = [34, 24, 16, 28, 24, 18, 18, 18]
    BOTH = ("rms", "chi2")
    ALL = ("rms", "chi2", "sd")

    def _level(self, sd=1.07, **kw):
        r = row("C 1s", 1.0, 200.0, states=[("a", 100.0), ("b", 100.0)], **kw)
        r["residual_sd"] = sd
        return level(None, [r])

    def test_the_default_table_is_unchanged(self):
        lv = self._level()
        self.assertEqual(rp.composition_header(lv), rp.COMPOSITION_HEADER)
        for _kind, cells in rp.composition_cells(lv):
            self.assertEqual(len(cells), len(rp.COMPOSITION_HEADER))

    def test_the_toggle_adds_one_last_column(self):
        lv = self._level()
        header = rp.composition_header(lv, self.ALL)
        self.assertEqual(header, rp.COMPOSITION_HEADER + (rp.SD_HEADER,))
        rows = rp.composition_cells(lv, self.ALL)
        self.assertEqual(rows[0][0], "region")
        self.assertEqual(rows[0][1][-1], "1.07")
        for kind, cells in rows:
            self.assertEqual(len(cells), len(header))
            if kind == "state":
                self.assertEqual(cells[-1], "")        # a state has none

    def test_a_report_table_pads_its_weights_and_keeps_the_width(self):
        lv = self._level()
        header, rows, weights = rp.composition_table(lv, self.W, self.ALL)
        self.assertEqual(header[-1], rp.SD_HEADER)
        self.assertEqual(len(weights), len(header))
        self.assertAlmostEqual(sum(weights), sum(self.W) + self.W[-1])

    def test_an_empty_column_is_left_out_as_the_others_are(self):
        lv = self._level(sd=None)
        header, _rows, weights = rp.composition_table(lv, self.W, self.ALL)
        self.assertNotIn(rp.SD_HEADER, header)
        self.assertEqual(len(weights), len(header))

    def test_each_statistic_can_be_dropped_on_its_own(self):
        lv = self._level()
        fixed = rp.COMPOSITION_HEADER[:rp.FIXED_COLUMNS]
        for stats in [(), ("rms",), ("chi2",), ("sd",), ("rms", "sd"),
                      ("chi2", "sd"), self.BOTH, self.ALL]:
            header, rows, weights = rp.composition_table(lv, self.W, stats)
            want = fixed + tuple(rp._STAT_HEADERS[k] for k in rp.STATS
                                 if k in stats)
            self.assertEqual(tuple(header), want, stats)
            self.assertEqual(len(weights), len(header), stats)
            for _kind, cells in rows:
                self.assertEqual(len(cells), len(header), stats)

    def test_none_means_the_default_pair(self):
        lv = self._level()
        self.assertEqual(rp.composition_header(lv, None),
                         rp.COMPOSITION_HEADER)
        self.assertEqual(rp.composition_header(lv, self.BOTH),
                         rp.COMPOSITION_HEADER)

    def test_casaxps_own_numbers_show_no_fit_quality(self):
        r = row("C 1s", 1.0, 200.0)
        r["casa_pct"] = 50.0
        r["residual_sd"] = 1.0
        lv = level(None, [r])
        self.assertIs(rp.composition_header(lv, self.ALL), rp.CASA_HEADER)
        for _kind, cells in rp.composition_cells(lv, self.ALL):
            self.assertEqual(len(cells), len(rp.CASA_HEADER))


class TestThePlotNote(unittest.TestCase):
    class _Cv:
        def __init__(self, region, rms, sd):
            self.region, self.residual_rms, self.residual_sd = region, rms, sd

    def _text(self, stats, curves):
        try:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib.figure import Figure
        except ImportError:
            self.skipTest("matplotlib not installed")
        import plots
        ax = Figure().subplots()
        plots._fit_rms_text(ax, curves, "grey", stats)
        return [t.get_text() for t in ax.texts]

    def test_off_it_is_the_old_note(self):
        self.assertEqual(self._text(("rms",), [self._Cv("C 1s", 0.012, 1.07)]),
                         ["RMS 1.2%"])

    def test_the_rms_line_can_be_switched_off(self):
        cv = [self._Cv("C 1s", 0.012, 1.07)]
        self.assertEqual(self._text((), cv), [])
        self.assertEqual(self._text(("chi2",), cv), [])
        self.assertEqual(self._text(("sd",), cv), ["Residual StDev 1.07"])

    def test_on_a_second_line_gives_the_statistic(self):
        self.assertEqual(self._text(("rms", "sd"),
                                    [self._Cv("C 1s", 0.012, 1.07)]),
                         ["RMS 1.2%\nResidual StDev 1.07"])

    def test_a_curve_without_it_is_left_out(self):
        got = self._text(("rms", "sd"), [self._Cv("C 1s", 0.012, None),
                                         self._Cv("O 1s", 0.02, 1.5)])
        self.assertEqual(got, ["RMS C 1s: 1.2%, O 1s: 2.0%\n"
                               "Residual StDev O 1s: 1.50"])


class TestTheHtmlPage(unittest.TestCase):
    def test_the_payload_flag_follows_the_toggle(self):
        import htmlbrowser as hb
        from test_htmlbrowser import doc, region
        d = [doc("a.vms", [region()])]
        self.assertNotIn("fit_stats", hb.build_payload(d))
        self.assertNotIn("fit_stats",
                         hb.build_payload(d, fit_stats=("rms", "chi2")))
        self.assertEqual(
            hb.build_payload(d, fit_stats=("rms", "sd"))["fit_stats"],
            ["rms", "sd"])
        self.assertEqual(hb.build_payload(d, fit_stats=())["fit_stats"], [])

    def test_the_page_shows_it_only_when_the_flag_is_there(self):
        with open(os.path.join(ROOT, "viewer", "viewer.js"),
                  encoding="utf-8") as fh:
            js = fh.read()
        self.assertIn("S.data.fit_stats || ['rms', 'chi2']", js)
        self.assertIn("stats.indexOf('sd') >= 0 && row.residual_sd", js)
        self.assertIn("'Residual StDev '", js)


if __name__ == "__main__":
    unittest.main()
