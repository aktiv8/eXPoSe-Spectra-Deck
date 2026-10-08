"""Ion scattering (ISS) spectra from Thermo Avantage files: the ion beam the
file records, how it is described, and how the ISS dialog uses it.

Synthetic files run everywhere; the real-data tests need the Nexsa ISS runs:

    set XPS_THERMO_ISS_CORPUS=<...\\Experiment_XPS ISS>
    python -m unittest tests.test_iss_thermo

Thermo records the ion (He+), the set beam energy (1000 eV), the scattering
angle (123.028 degrees) and an "ISS calibration" of the beam energy (966.313
eV). The calibrated energy is the one the peaks fit: In falls at 867.7 eV for
it and 897.9 eV for 1000 eV, and the spectra peak at 866-870 eV.
"""

import glob
import os
import shutil
import sys
import tempfile
import types
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import annotations  # noqa: E402
import elements as el  # noqa: E402
import methods  # noqa: E402
from readers import load_file  # noqa: E402
from test_readers import AVG_TEMPLATE  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

CORPUS = os.environ.get("XPS_THERMO_ISS_CORPUS", "")
HE = el.ION_MASS["He+"]
IN = el.ELEMENTS["In"][1]


def iss_avg(title="ISS Survey", ion=True):
    """An ISS ``.avg`` as Avantage writes it: technique 6, no photon energy,
    energy-scale flag 0, the ion gun block."""
    rows = "\n".join("LIST@ %3d=  %s" % (k, ",  ".join(
        "%.6f" % (100 + i) for i in range(k, k + 4))) for k in (0, 4))
    text = AVG_TEMPLATE.format(title=title, rows=rows)
    text = (text.replace("DS_SOPROPID_ENERGY                          "
                         ": VT_R4   = 1486.680054\n", "")
            .replace("DS_ACPROPID_EV_SCALE                        : VT_I2   = 1",
                     "DS_ACPROPID_EV_SCALE                        : VT_I2   = 0")
            .replace("DS_ANPROPID_LENS_MODE_NAME                  "
                     ": VT_BSTR = 'Standard'",
                     "DS_ANPROPID_LENS_MODE_NAME                  "
                     ": VT_BSTR = 'ISS'")
            .replace("DS_GEPROPID_INSTRUMENT ",
                     "DS_GEPROPID_TECHNIQUE                       : VT_I4   = 6\n"
                     "DS_GEPROPID_INSTRUMENT "))
    if ion:
        text = text.replace(
            "DS_ANPROPID_PASS ",
            "DS_SOURCE_IONGUNPROPID_ENERGY               : VT_R4   = 1000.000000\n"
            "DS_SOURCE_IONGUNPROPID_IONTYPE              : VT_BSTR = 'He+'\n"
            "DS_SOURCE_IONGUNPROPID_DESCRIPTION          : VT_BSTR = '1000, Low (iss)'\n"
            "DS_SOURCE_IONGUNPROPID_ISS_ANGLE            : VT_R4   = 123.028000\n"
            "DS_SOURCE_IONGUNPROPID_ISS_CALIBRATION      : VT_R4   = 966.312988\n"
            "DS_ANPROPID_PASS ")
    return text


class Tmp(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, name, text):
        p = os.path.join(self.dir.name, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p


class TestThermoIssFile(Tmp):
    def test_it_is_an_iss_spectrum_without_noise(self):
        f = load_file(self.write("ISS Survey.avg", iss_avg()))
        r = f.regions[0]
        self.assertEqual(f.warnings, [])        # kinetic by nature
        self.assertEqual((r.technique, r.name), ("ISS", "ISS Survey"))
        self.assertTrue(r.is_iss)
        self.assertFalse(r.is_survey)
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertIsNone(r.photon_energy)
        self.assertEqual(r.anode, "")

    def test_the_recorded_beam_is_kept_and_listed(self):
        f = load_file(self.write("ISS Survey.avg", iss_avg()))
        r = f.regions[0]
        self.assertEqual(r.extra["iss"], {
            "ion": "He+", "e0": 1000.0, "e0_cal": 966.313, "theta": 123.028,
            "description": "1000, Low (iss)"})
        md = f.region_metadata(r)
        self.assertEqual(
            (md["Technique"], md["ISS ion"], md["ISS beam energy (eV)"],
             md["ISS calibrated beam energy (eV)"],
             md["ISS scattering angle (°)"]),
            ("ISS", "He+", "1000", "966.313", "123.028"))
        self.assertIn("KE start (eV)", md)

    def test_a_file_without_the_gun_block_still_reads_as_iss(self):
        f = load_file(self.write("ISS Survey.avg", iss_avg(ion=False)))
        self.assertEqual(f.regions[0].technique, "ISS")
        self.assertNotIn("iss", f.regions[0].extra)

    def test_an_xps_file_is_untouched(self):
        rows = "LIST@   0=  1.0,  2.0,  3.0,  4.0\nLIST@   4=  5.0,  6.0,  7.0,  8.0"
        f = load_file(self.write("C1s Scan.avg", AVG_TEMPLATE.format(
            title="C1s Scan", rows=rows)))
        self.assertEqual((f.regions[0].technique, f.regions[0].name),
                         ("XPS", "C 1s"))

    def test_methods_state_the_recorded_beam_with_its_calibration(self):
        f = load_file(self.write("ISS Survey.avg", iss_avg()))
        text = methods.generate([f.region_metadata(r) for r in f.regions])
        self.assertIn("He+ ions, a beam energy of 1000 eV (calibrated "
                      "966.313 eV) and a scattering angle of 123.028°", text)
        self.assertNotIn("X-ray", text)

    def test_a_saved_energy_replaces_the_recorded_pair(self):
        f = load_file(self.write("ISS Survey.avg", iss_avg()))
        r = f.regions[0]
        a = annotations.Annotations()
        a.set_iss("f1", {"ion": "He+", "e0": 966.313, "theta": 123.028})
        md = a.apply_metadata("f1", 0, r, f.region_metadata(r))
        self.assertEqual(md["ISS beam energy (eV)"], "966.313")
        self.assertNotIn("ISS calibrated beam energy (eV)", md)


class TestEnergyFit(unittest.TestCase):
    def test_in_sits_where_the_calibration_says(self):
        cal = el.iss_energy(966.313, HE, IN, 123.028)
        nominal = el.iss_energy(1000.0, HE, IN, 123.028)
        self.assertAlmostEqual(cal, 867.7, delta=0.1)
        self.assertAlmostEqual(nominal, 897.9, delta=0.1)


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestIssDialogOnThermo(Tmp):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def open(self):
        ws = ee.Workspace(self.root)
        ws._add_files([self.write("ISS Survey.avg", iss_avg())])
        region = ws.docs[0].regions[0]
        ws.checked = {id(region)}
        return ws, region

    def dialog(self, ws):
        from iss_ui import IssReelsDialog
        dlg = IssReelsDialog(self.root, ws)
        self.addCleanup(lambda: dlg.winfo_exists() and dlg.destroy())
        return dlg

    def test_the_file_beats_the_defaults(self):
        ws, region = self.open()
        self.assertEqual(ws.iss_recorded(region),
                         {"ion": "He+", "e0": 966.313, "e0_set": 1000.0,
                          "theta": 123.028})
        self.assertEqual(ws.iss_offer(region), {})       # not a Kratos file
        dlg = self.dialog(ws)
        self.assertEqual((dlg.ion.get(), dlg.e0.get(), dlg.theta.get()),
                         ("He+", "966.313", "123.028"))
        self.assertIn("Recorded in the file", dlg.iss_note.cget("text"))
        self.assertIn("set 1000 eV", dlg.iss_note.cget("text"))

    def test_a_saved_value_beats_the_file(self):
        ws, region = self.open()
        dlg = self.dialog(ws)
        dlg.theta.set("124")
        dlg._remember()
        dlg.destroy()
        self.assertEqual(self.dialog(ws).theta.get(), "124")

    def test_in_is_a_candidate_with_the_calibrated_energy_only(self):
        ws, _region = self.open()
        dlg = self.dialog(ws)
        dlg._on_iss_click(867.7)
        self.assertIn("In", [c["symbol"] for c in dlg.cands])
        self.assertNotIn("In", [c["symbol"] for c in el.candidates(
            867.7, 1000.0, "He+", 123.028)])


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_THERMO_ISS_CORPUS to the Experiment_XPS ISS folder")
class TestRealThermoIss(unittest.TestCase):
    def pairs(self):
        return [(a, a[:-4] + ".VGD") for a in sorted(glob.glob(
            os.path.join(CORPUS, "ISS Source", "*", "*.avg")))]

    def test_every_pair_is_iss_and_the_two_readers_agree(self):
        n = 0
        for avg, vgd in self.pairs():
            a, v = load_file(avg), load_file(vgd)
            ra, rv = a.regions[0], v.regions[0]
            with self.subTest(os.path.relpath(avg, CORPUS)):
                self.assertEqual((a.warnings, v.warnings), ([], []))
                self.assertEqual((ra.technique, rv.technique), ("ISS", "ISS"))
                self.assertEqual(ra.name, rv.name)
                self.assertEqual(ra.extra["iss"], rv.extra["iss"])
                self.assertEqual(ra.extra["iss"]["e0_cal"], 966.313)
                self.assertEqual(ra.energy, rv.energy)
                for x, y in zip(ra.counts, rv.counts):
                    self.assertAlmostEqual(x, y, delta=1e-3 * max(1.0, abs(x)))
            n += 1
        self.assertEqual(n, 40)

    def test_the_xps_files_of_the_experiment_stay_xps(self):
        for avg in sorted(glob.glob(os.path.join(
                CORPUS, "X-Ray001*", "*", "*.avg")))[:12]:
            with self.subTest(os.path.basename(avg)):
                r = load_file(avg).regions[0]
                self.assertEqual(r.technique, "XPS")
                self.assertFalse(r.is_iss)

    def test_whole_experiment_loads_without_warnings(self):
        e = load_file(CORPUS)
        self.assertEqual(e.warnings, [])
        self.assertEqual(sum(1 for r in e.regions if r.technique == "ISS"), 40)

    def test_indium_peak_fits_the_calibrated_energy(self):
        want = el.iss_energy(966.313, HE, IN, 123.028)
        found = 0
        for avg, _vgd in self.pairs():
            r = load_file(avg).regions[0]
            hi = [i for i, e in enumerate(r.energy) if 800 <= e <= 940]
            i = max(hi, key=r.counts.__getitem__)
            others = sorted(r.counts)[len(r.counts) // 2]
            if r.counts[i] < 3 * others:
                continue                              # no clear peak here
            if abs(r.energy[i] - want) < 6:
                found += 1
                self.assertIn("In", [c["symbol"] for c in el.candidates(
                    r.energy[i], 966.313, "He+", 123.028)])
        self.assertGreaterEqual(found, 3)


if __name__ == "__main__":
    unittest.main()
