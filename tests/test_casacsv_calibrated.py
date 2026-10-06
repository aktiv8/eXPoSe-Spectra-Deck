"""CasaXPS ASCII exports of a charge-corrected file (casacsv.py).

The export's B.E. column is the CALIBRATED axis (what CasaXPS shows); a
region's own axis is the file's raw one, ``calibration_shift`` lower (the real
PET example: 3.008 eV). Matching must allow for it and the curves must land on
the raw kinetic-energy axis ``casafit.curves`` works in. Also the
``is_export`` test the Open path uses, and the real PET files when
``XPS_CASACSV_CORPUS`` names their folder.

Run:  python -m unittest discover tests
"""

import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import casacsv  # noqa: E402
import casafit  # noqa: E402
import readers  # noqa: E402
from test_casacsv import _mk_region  # noqa: E402

try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class TestChargeCorrectedFile(unittest.TestCase):
    HV, SHIFT, DWELL, SCANS = 1486.6, 3.0, 0.3, 2

    def make(self, flagged=False, csv_shift=None, cps_scale=1.0,
             region_shift=None, scramble=False):
        """``(region, block)``: a Gaussian peak at 531.0 eV calibrated
        (528.0 raw) on a flat background. ``csv_shift`` = how far the CSV's
        axis sits above the raw one (default: the file's own correction)."""
        hv, sh = self.HV, self.SHIFT
        raw = np.linspace(540.0, 520.0, 201)                  # descending
        csv_shift = sh if csv_shift is None else csv_shift
        cal = raw + csv_shift

        def comp(be_axis):
            return 100.0 * np.exp(
                -((be_axis - (528.0 + csv_shift)) / 0.8) ** 2)
        bg = np.full(len(raw), 50.0)
        cps = bg + comp(cal)
        counts = cps * self.DWELL * self.SCANS * cps_scale
        if flagged:        # limits and positions stored raw ("Regions Comps")
            fr = casafit.FitRegion(name="R", background="Linear",
                                   start_ke=hv - 542.0, end_ke=hv - 512.0)
            fc = casafit.FitComponent(name="A", pos_ke=hv - 528.0, region="R")
            fit = casafit.Fit(regions=[fr], components=[fc], calib_shift=sh,
                              region_shift=0.0, comp_shift=0.0)
        else:              # stored in the calibrated frame
            fr = casafit.FitRegion(name="R", background="Linear",
                                   start_ke=hv - 545.0, end_ke=hv - 515.0)
            fc = casafit.FitComponent(name="A", pos_ke=hv - 531.0, region="R")
            fit = casafit.Fit(regions=[fr], components=[fc], calib_shift=sh)
        region = _mk_region("R", "S", raw, counts, hv, dwell=self.DWELL,
                            scans=self.SCANS, fit=fit)
        region.calibration_shift = sh if region_shift is None else region_shift
        if scramble:
            cps = cps[np.random.RandomState(1).permutation(len(cps))]
        block = casacsv.CsvBlock(
            source="x.csv", block_index=0, layout="rows", cycle=0,
            sample="S", scan_name="R Scan", hv=None, dwell_total=None,
            ke=None, counts=None, be=tuple(cal), cps=tuple(cps),
            background_cps=tuple(bg), envelope_cps=tuple(bg + comp(cal)),
            components=[casacsv.CsvComponent(
                name="A", position_be=528.0 + csv_shift, fwhm=1.0, area=1.0,
                lineshape="GL(30)", curve_cps=tuple(comp(cal)))])
        return region, block

    def curves(self, region):
        dw, sc = region.dwell_and_scans()
        return casafit.curves(region.fit, region.energy, region.counts,
                              region.photon_energy, dw, sc,
                              prefer_csv=True)[0]

    def test_rows_layout_matches_across_the_correction_and_lands_on_the_peak(
            self):
        for flagged in (False, True):
            region, block = self.make(flagged=flagged)
            report = casacsv.match_to_regions([block], [region])
            res = report.results[0]
            self.assertEqual(res.reason, "", flagged)
            self.assertIs(res.region, region)
            self.assertEqual(
                (res.n_components_total, res.n_components_aligned), (1, 1))
            self.assertIn("charge correction", res.frame_note)
            casacsv.apply_matches(report)
            cv = self.curves(region)
            self.assertFalse(cv.approximate)
            self.assertLess(cv.residual_rms, 1e-6, flagged)
            # the curves' own axis is the raw kinetic energy
            ke = np.asarray(res.csv_curves.ke)
            self.assertAlmostEqual(float(ke.min()), self.HV - 540.0, places=6)

    def test_the_axis_the_old_code_built_would_misplace_every_curve(self):
        region, block = self.make()
        report = casacsv.match_to_regions([block], [region])
        casacsv.apply_matches(report)
        good = self.curves(region).residual_rms
        cc = region.fit.regions[0].csv_curves
        region.fit.regions[0].csv_curves = casafit.CsvCurves(
            ke=tuple(self.HV - b for b in block.be),         # calibrated KE
            background=cc.background, components=cc.components,
            envelope=cc.envelope, source="old")
        self.assertGreater(self.curves(region).residual_rms,
                           100 * good + 0.05)

    def test_a_csv_on_the_raw_axis_still_matches_without_a_note(self):
        region, block = self.make(csv_shift=0.0)
        report = casacsv.match_to_regions([block], [region])
        res = report.results[0]
        self.assertEqual(res.reason, "")
        self.assertEqual(res.frame_note, "")
        casacsv.apply_matches(report)
        self.assertLess(self.curves(region).residual_rms, 1e-6)

    def test_a_different_correction_is_inferred_only_when_the_cps_agrees(self):
        # the file carries no correction, the CSV was exported at +2.5 eV
        region, block = self.make(csv_shift=2.5, region_shift=0.0)
        report = casacsv.match_to_regions([block], [region])
        res = report.results[0]
        self.assertEqual(res.reason, "")
        self.assertIn("inferred", res.frame_note)
        casacsv.apply_matches(report)
        self.assertLess(self.curves(region).residual_rms, 1e-6)
        # the same spans, but a CPS column that is not this spectrum
        region, block = self.make(csv_shift=2.5, region_shift=0.0,
                                  scramble=True)
        res = casacsv.match_to_regions([block], [region]).results[0]
        self.assertIsNone(res.csv_curves)
        self.assertTrue(res.reason.startswith("no matching region"))

    def test_the_cps_scale_does_not_matter_for_an_inferred_match(self):
        # an unknown dwell: the CSV CPS is the counts times one constant
        region, block = self.make(csv_shift=2.5, region_shift=0.0,
                                  cps_scale=3.7)
        self.assertEqual(
            casacsv.match_to_regions([block], [region]).results[0].reason, "")

    def test_a_failed_match_names_both_ranges(self):
        region, block = self.make(csv_shift=9.0, region_shift=0.0,
                                  scramble=True)
        reason = casacsv.match_to_regions([block], [region]).results[0].reason
        self.assertIn("CSV 529.0-549.0 eV, 201 points", reason)
        self.assertIn("'R' 520.0-540.0 eV, 201 points", reason)

    def test_no_fitted_region_in_the_sample_is_said_so(self):
        region, block = self.make()
        region.fit = None
        reason = casacsv.match_to_regions([block], [region]).results[0].reason
        self.assertIn("has a CasaXPS fit", reason)

    def test_columns_layout_takes_the_offset_from_its_own_data(self):
        region, block = self.make()
        raw_be = np.asarray(region.energy)
        block = casacsv.CsvBlock(
            source="x.csv", block_index=0, layout="columns", cycle=0,
            sample="S", scan_name="R Scan", hv=self.HV,
            dwell_total=self.DWELL * self.SCANS,
            ke=tuple(self.HV - raw_be), counts=tuple(region.counts),
            be=block.be, cps=block.cps, background_cps=block.background_cps,
            envelope_cps=block.envelope_cps, components=block.components)
        report = casacsv.match_to_regions([block], [region])
        res = report.results[0]
        self.assertEqual(res.reason, "")
        self.assertIn("K.E. column", res.frame_note)
        casacsv.apply_matches(report)
        self.assertLess(self.curves(region).residual_rms, 1e-6)

    def test_fit_windows_are_compared_in_the_calibrated_frame(self):
        """Two named windows on one scan, stored raw (flagged) or calibrated:
        both must be found inside the calibrated block."""
        hv, sh = self.HV, self.SHIFT
        raw = np.linspace(560.0, 520.0, 81)
        cal = raw + sh
        for flagged in (False, True):
            lift = sh if flagged else 0.0      # a raw KE is `sh` higher

            def window(name, hi, lo):
                return casafit.FitRegion(
                    name=name, background="Linear",
                    start_ke=hv - hi + lift, end_ke=hv - lo + lift)
            fit = casafit.Fit(
                regions=[window("A", 545.0, 541.0),
                         window("B", 530.0, 526.0)],
                components=[
                    casafit.FitComponent(name="a1", pos_ke=hv - 543.0 + lift,
                                         region="A"),
                    casafit.FitComponent(name="b1", pos_ke=hv - 528.0 + lift,
                                         region="B")],
                calib_shift=sh, region_shift=0.0 if flagged else None,
                comp_shift=0.0 if flagged else None)
            region = _mk_region("R", "S", raw, np.ones(len(raw)), hv, fit=fit)
            region.calibration_shift = sh
            block = casacsv.CsvBlock(
                source="x.csv", block_index=0, layout="rows", cycle=0,
                sample="S", scan_name="A B Scan", hv=None, dwell_total=None,
                ke=None, counts=None, be=tuple(cal),
                cps=tuple(np.ones(len(raw))), background_cps=None,
                envelope_cps=None, components=[
                    casacsv.CsvComponent("a1", 543.0, 1.0, 1.0, "GL(30)",
                                         tuple([0.1] * len(raw))),
                    casacsv.CsvComponent("b1", 528.0, 1.0, 1.0, "GL(30)",
                                         tuple([0.2] * len(raw)))])
            report = casacsv.match_to_regions([block], [region])
            got = {r.fit_region.name: r for r in report.results}
            self.assertEqual(set(got), {"A", "B"}, flagged)
            for r in got.values():
                self.assertEqual(r.reason, "", flagged)

    def test_summary_mentions_the_correction(self):
        region, block = self.make()
        text = casacsv.summarise(casacsv.match_to_regions([block], [region]))
        self.assertIn("1 of 1 CSV block(s) matched.", text)
        self.assertIn("charge correction (+3.000 eV)", text)


class TestIsExport(unittest.TestCase):
    ROWS = ('Name,,"A","B",,Name,,"C"\nPosition,,"1","2",,Position,,"3"\n'
            'FWHM,,"1","1",,FWHM,,"1"\nArea,,"1","1",,Area,,"1"\n'
            'Lineshape,,"GL(30)","GL(30)",,Lineshape,,"GL(30)"\n\n'
            'B.E.,"Cycle 0:a.vms:R:CPS","R:A","R:B"\n1.0,2.0,3.0,4.0\n')
    COLUMNS = ('Cycle 0:s.vms:R Scan\n,Characteristic Energy eV,1486.6,'
               'Acquisition Time s,0.6\nName,,A\nPosition,,1\n'
               'K.E.,Counts,B.E.,CPS,A\n1,2,3,4,5\n')

    def test_both_layouts_are_recognised(self):
        self.assertTrue(casacsv.is_export(self.ROWS.encode()))
        self.assertTrue(casacsv.is_export(self.COLUMNS.encode()))
        self.assertTrue(casacsv.is_export(b"\xef\xbb\xbf" + self.ROWS.encode()))

    def test_other_text_is_not(self):
        self.assertFalse(casacsv.is_export(b""))
        self.assertFalse(casacsv.is_export(b"BE,Counts/s\n1,2\n3,4\n"))
        self.assertFalse(casacsv.is_export(               # a table of names
            b"Name,Area\n" + b"".join(b"C 1s,%d\n" % i for i in range(30))))
        self.assertFalse(casacsv.is_export(b"\x00\x01\x02"))

    def test_reader_for_says_what_it_is_and_what_to_do(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "fit.csv")
        with open(p, "w", encoding="utf-8", newline="") as fh:
            fh.write(self.ROWS)
        try:
            with self.assertRaises(readers.UnsupportedFormat) as cm:
                readers.reader_for(p)
            self.assertIn("CasaXPS ASCII export", str(cm.exception))
            self.assertIn("Import CasaXPS CSV export", str(cm.exception))
        finally:
            os.remove(p)
            os.rmdir(d)


@unittest.skipUnless(os.environ.get("XPS_CASACSV_CORPUS") and HAVE_NP,
                     "set XPS_CASACSV_CORPUS to the 'PET EXAMPLE' folder")
class TestRealPetExample(unittest.TestCase):
    """``PET 10 eV Fitted.vms`` (3.008 eV charge correction) and its rows-layout
    ``PET_CasaFit.csv``: the case that used to match nothing."""

    def setUp(self):
        d = os.environ["XPS_CASACSV_CORPUS"]
        self.doc = readers.load_file(os.path.join(d, "PET 10 eV Fitted.vms"))
        self.imp = casacsv.parse(os.path.join(d, "PET_CasaFit.csv"))

    def test_both_blocks_match_with_every_component(self):
        report = casacsv.match_to_regions(self.imp.blocks, self.doc.regions)
        self.assertEqual([r.reason for r in report.results], ["", ""])
        self.assertEqual([(r.n_components_aligned, r.n_components_total)
                          for r in report.results], [(5, 5), (6, 6)])
        for r in report.results:
            self.assertIn("charge correction (+3.008 eV)", r.frame_note)

    def test_the_csv_curves_fit_the_data_better_than_the_reconstruction(self):
        report = casacsv.match_to_regions(self.imp.blocks, self.doc.regions)
        casacsv.apply_matches(report)
        for r in self.doc.regions:
            dw, sc = r.dwell_and_scans()
            lit = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                 dw, sc, prefer_csv=True)[0]
            rec = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                                 dw, sc)[0]
            self.assertLess(lit.residual_rms, 0.03, r.name)
            self.assertLess(lit.residual_rms, rec.residual_rms, r.name)
            self.assertFalse(lit.approximate)


if __name__ == "__main__":
    unittest.main()
