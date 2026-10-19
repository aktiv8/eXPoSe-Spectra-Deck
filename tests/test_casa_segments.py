"""Several same-named CasaXPS regions on one spectrum (a segmented Shirley
background) are one fit, not fragments."""

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import casafit  # noqa: E402

try:
    import numpy  # noqa: F401
    HAVE_NP = True
except ImportError:
    HAVE_NP = False

HV = 1486.6
CORPUS = os.environ.get("XPS_CE_CORPUS", "")


def reg(a, b, rsf=8.81, name="Ce 3d"):
    return casafit.FitRegion(name=name, background="Shirley", start_ke=a,
                             end_ke=b, rsf=rsf, avg=5)


def spectrum(lo, hi, step=0.1):
    """Binding energies (descending) covering KE lo..hi, with one peak."""
    n = int(round((hi - lo) / step)) + 1
    ke = [lo + i * step for i in range(n)]
    be = [HV - k for k in ke]
    counts = [100.0 + 1000.0 * (1.0 / (1.0 + ((k - 595.0) / 1.5) ** 2))
              for k in ke]
    return be, counts


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class TestSegments(unittest.TestCase):
    def test_sliver_overlap_is_not_a_duplicate(self):
        a, b = reg(595.17, 614.70), reg(576.79, 595.20)
        self.assertEqual(casafit.distinct_regions([a, b]), [a, b])

    def test_nested_regions_still_collapse(self):
        narrow, wide = reg(10, 20), reg(8, 24)
        self.assertEqual(casafit.distinct_regions([narrow, wide]), [wide])

    def test_groups_chain_touching_regions_only(self):
        a, b, c = reg(595.17, 614.7), reg(576.79, 595.2), reg(563.6, 576.79)
        far = reg(520, 530)
        groups = casafit.segment_groups([a, b, c, far])
        self.assertEqual([[g.start_ke for g in grp] for grp in groups],
                         [[563.6, 576.79, 595.17], [520]])

    def test_three_segments_with_no_components_make_one_background(self):
        fit = casafit.Fit(regions=[reg(595.17, 614.7), reg(576.79, 595.2),
                                   reg(563.6, 576.79)])
        be, counts = spectrum(563.0, 615.0)
        cvs = casafit.curves(fit, be, counts, HV, dwell=1.0, scans=1)
        self.assertEqual(len(cvs), 1)
        bg = [v for v in cvs[0].background if v == v]
        inside = sum(1 for k in (HV - e for e in be) if 563.6 <= k <= 614.7)
        self.assertEqual(len(bg), inside)           # no gap in the middle
        self.assertIsNone(cvs[0].envelope)

    def test_components_are_drawn_across_the_whole_group(self):
        fit = casafit.Fit(
            regions=[reg(594.44, 613.95, 8.81), reg(570.9, 594.58, 51.62)],
            components=[
                casafit.FitComponent(name="Ce 3d", shape="GL(30)", area=50.0,
                                     fwhm=2.0, pos_ke=605.0, region="Ce 3d"),
                casafit.FitComponent(name="Ce 3d", shape="GL(30)", area=80.0,
                                     fwhm=2.0, pos_ke=586.0, region="Ce 3d")])
        be, counts = spectrum(570.0, 615.0)
        cvs = casafit.curves(fit, be, counts, HV, dwell=1.0, scans=1)
        self.assertEqual(len(cvs), 1)
        cv = cvs[0]
        self.assertEqual(len(cv.components), 2)
        for _c, vals in cv.components:          # full length over the union
            self.assertEqual(sum(1 for v in vals if v == v),
                             sum(1 for v in cv.background if v == v))
        self.assertIsNotNone(cv.envelope)
        self.assertEqual(cv.fit_region.rsf, 0.0)    # the two RSFs differ
        self.assertAlmostEqual(cv.fit_region.start_ke, 570.9)
        self.assertAlmostEqual(cv.fit_region.end_ke, 613.95)

    def test_equal_rsf_is_kept(self):
        m = casafit._merged_region([reg(1, 5), reg(5, 9)])
        self.assertEqual(m.rsf, 8.81)


@unittest.skipUnless(HAVE_NP and CORPUS, "set XPS_CE_CORPUS to the Ceria folder")
class TestRealFile(unittest.TestCase):
    def test_every_ce3d_block_is_one_complete_fit(self):
        import readers
        d = readers.load_file(os.path.join(
            CORPUS, "Ce4 and Ce3 Suitable fits.vms"))
        for r in d.regions:
            dw, sc = r.dwell_and_scans()
            cvs = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                 dw, sc)
            self.assertEqual(len(cvs), 1, r.sample)
            lo, hi = (cvs[0].fit_region.start_ke, cvs[0].fit_region.end_ke)
            inside = sum(1 for e in r.energy
                         if lo <= r.photon_energy - e <= hi)
            self.assertAlmostEqual(
                sum(1 for v in cvs[0].background if v == v), inside,
                delta=2, msg=r.sample)           # edge points, calibration


if __name__ == "__main__":
    unittest.main()
