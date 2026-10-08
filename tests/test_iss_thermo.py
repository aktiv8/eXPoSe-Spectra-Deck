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
CORPUS2 = os.environ.get("XPS_THERMO_ISS_CORPUS2", "")
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


class TestUnacquiredAndSeries(Tmp):
    def test_a_header_with_only_the_iss_lens_is_still_iss_and_quiet(self):
        # an unacquired scan: no technique code, no ion gun block, no values
        text = iss_avg(ion=False).replace(
            "DS_GEPROPID_TECHNIQUE                       : VT_I4   = 6\n", "")
        text = text.split("$DATA=*")[0] + "$DATA=*\n"
        f = load_file(self.write("ISS Survey.avg", text))
        r = f.regions[0]
        self.assertEqual((r.technique, r.decodable), ("ISS", False))
        self.assertEqual(f.warnings, [])
        self.assertIn("Header only", r.note)

    def test_an_iteration_series_is_not_called_depth_profiling(self):
        rows = [{"Sample": "s", "Region": "ISS Survey", "Technique": "ISS",
                 "Instrument": "Nexsa", "Lens mode": "ISS",
                 "Etch level": str(k), "Pass energy (eV)": "200",
                 "KE start (eV)": "300", "KE end (eV)": "1000"}
                for k in range(50)]
        text = methods.generate(rows)
        self.assertIn("repeated in succession: 50 iterations per series", text)
        for word in ("sputter", "Depth profiling", "etching"):
            self.assertNotIn(word, text)

    def test_real_depth_profiles_keep_their_wording(self):
        rows = [{"Sample": "s", "Region": "C 1s", "Technique": "XPS",
                 "Etch level": str(k), "Etch time (s)": str(10 * k),
                 "BE start (eV)": "295", "BE end (eV)": "280"}
                for k in range(5)]
        self.assertIn("Depth profiling was performed", methods.generate(rows))


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




@unittest.skipUnless(CORPUS2 and os.path.isdir(CORPUS2),
                     "set XPS_THERMO_ISS_CORPUS2 to the 'YX AB' folder")
class TestSecondThermoSet(unittest.TestCase):
    """Four Nexsa ISS iteration series (8, 50, an unacquired one, 50 levels):
    both readers agree, the unacquired one is header-only in both, the
    recorded beam is read, and the peaks sit where the calibrated beam energy
    puts them, with the scatter these (probably charging) samples show."""

    def path(self, d, ext):
        return os.path.join(CORPUS2, "ISS Source", d, "Iteration",
                            "ISS Survey." + ext)

    def test_the_two_readers_agree_on_every_level(self):
        for d, n in (("1", 8), ("1b", 50), ("2b", 50)):
            with self.subTest(d):
                a, v = load_file(self.path(d, "avg")), load_file(self.path(d, "VGD"))
                self.assertEqual((len(a.regions), len(v.regions)), (n, n))
                self.assertEqual((a.warnings, v.warnings), ([], []))
                for x, y in zip(a.regions, v.regions):
                    self.assertEqual((x.technique, x.etch_level, x.energy),
                                     ("ISS", y.etch_level, y.energy))
                    self.assertEqual(x.extra["iss"], y.extra["iss"])
                    self.assertEqual(x.extra["iss"]["e0_cal"], 966.313)
                    for p, q in zip(x.counts, y.counts):
                        self.assertAlmostEqual(p, q, delta=1e-3 * max(1.0, abs(p)))

    def test_the_unacquired_file_is_header_only_in_both(self):
        for ext in ("avg", "VGD"):
            with self.subTest(ext):
                f = load_file(self.path("2", ext))
                r = f.regions[0]
                self.assertEqual((r.decodable, r.technique), (False, "ISS"))
                self.assertEqual(f.warnings, [])

    def test_the_methods_describe_iterations_not_sputtering(self):
        f = load_file(self.path("1b", "avg"))
        text = methods.generate([f.region_metadata(r) for r in f.regions])
        self.assertIn("50 iterations per series", text)
        self.assertNotIn("sputter", text)

    def test_peaks_sit_where_the_calibrated_energy_puts_them(self):
        cal, nom = 966.313, 1000.0

        def peak(r, lo, hi):
            """(energy, clear?) of the maximum of a 5-point (10 eV) average in
            a window. Clear = the bump over the mean of the window's two ends
            is at least 4 times the noise (rms of the counts around that
            average), so a few counts of noise in a weak spectrum are not
            taken for a peak."""
            c = r.counts
            sm = [sum(c[max(0, i - 2):i + 3]) / len(c[max(0, i - 2):i + 3])
                  for i in range(len(c))]
            noise = (sum((x - y) ** 2 for x, y in zip(c, sm)) / len(c)) ** 0.5
            idx = [i for i, e in enumerate(r.energy) if lo <= e <= hi]
            i = max(idx, key=sm.__getitem__)
            edge = (sm[idx[0]] + sm[idx[-1]]) / 2
            return r.energy[i], (sm[i] - edge) >= 4 * noise
        tested = skipped = 0
        for d in ("1", "1b", "2b"):
            f = load_file(self.path(d, "avg"))
            for r in (f.regions[0], f.regions[len(f.regions) // 2],
                      f.regions[-1]):
                for sym, lo, hi in (("Al", 560, 640), ("Ti", 700, 760)):
                    m = el.ELEMENTS[sym][1]
                    pc = el.iss_energy(cal, HE, m, 123.028)
                    pn = el.iss_energy(nom, HE, m, 123.028)
                    e, clear = peak(r, lo, hi)
                    if not clear:
                        skipped += 1
                        continue
                    tested += 1
                    with self.subTest(f"{d} level {r.etch_level} {sym}"):
                        self.assertLess(abs(e - pc), 30)      # charging scatter
                        self.assertLess(abs(e - pc), abs(e - pn))
                o, clear = peak(r, 400, 500)
                if clear:
                    with self.subTest(f"{d} level {r.etch_level} O"):
                        self.assertLess(abs(o - el.iss_energy(
                            cal, HE, 15.9949, 123.028)), 20)
        self.assertGreaterEqual(tested, 12, f"only {tested} windows tested "
                                            f"({skipped} had no clear peak)")


if __name__ == "__main__":
    unittest.main()
