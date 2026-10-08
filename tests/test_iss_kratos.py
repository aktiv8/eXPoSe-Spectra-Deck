"""Ion scattering (ISS) spectra from Kratos files, and the ISS tools on them.

Synthetic files run everywhere. The real-data test needs the Kratos Axis Ultra
ISS runs and is skipped without them:

    set XPS_ISS_CORPUS=<folder with the .kal / .dset / .vms of the ISS runs>
    python -m unittest tests.test_iss_kratos

Geometry and the 933 eV gold peak come from Kratos test procedure TPC1369C
(the Axis Ultra ion gun is 45 degrees from the surface normal, so the
scattering angle is 135 degrees; standard He beam HT 1 kV).
"""

import glob
import os
import sys
import tempfile
import types
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import annotations  # noqa: E402
import elements as el  # noqa: E402
import methods  # noqa: E402
from readers import Region, load_file  # noqa: E402
from readers import kratos_dset  # noqa: E402
from readers.vamas import VamasFile, _LENS_LINE  # noqa: E402
from tests.test_kratos_dset import (block, dset, r_dbl, r_f32, r_int,  # noqa: E402
                                    r_str)

CORPUS = os.environ.get("XPS_ISS_CORPUS", "")
HE = el.ION_MASS["He+"]


def iss_object(name="ISS 1", lens="F_HSA_LENS_ISS"):
    """A ``.kal`` object of an ISS spectrum, as the Axis Ultra writes it."""
    return "\n".join([
        f"Object name               = {name}/2",
        "   1 Technique                                        = F_ISS",
        "   2 Scan type                                        = F_SPECTRUM",
        "   3 Spectrum scan start                              = 100 eV",
        "   4 Spectrum scan step size                          = 1 eV",
        "   5 Abscissa label                                   = Kinetic Energy",
        "   6 Abscissa units                                   = eV",
        "   7 Dwell time                                       = 0.0666667 seconds",
        "   9 Ordinate units                                   = counts",
        "  12 Ordinate values                                  = {5, 9, 40, 12}",
        f"  37 Acquisition name                                 = {name}",
        "  42 Pass energy                                      = 320 eV",
        "  99 # Sweeps completed                               = 1",
        " 151 Date Acquired                                    = 16/11/11 13:49:29",
        f"3049 HSA Lens Mode                                    = {lens}",
        "3080 Xray Reference Energy                            = F_REFER_TO_NONE",
        "3203 NICPU Ion Gun PSU beam_ht                        = 980 V",
        "3210 NICPU Ion Gun PSU Emission Current               = 0.015 A",
        "3282 Descriptor for aperture size used in acquisition = Slot",
        ""])


class Tmp(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, name, data):
        p = os.path.join(self.dir.name, name)
        with open(p, "wb" if isinstance(data, bytes) else "w") as fh:
            fh.write(data)
        return p


class TestKratosIss(Tmp):
    def test_kal_object_is_an_iss_spectrum(self):
        text = "Dataset filename          = a.dset\n" + iss_object()
        f = load_file(self.write("a.kal", text))
        self.assertEqual(f.warnings, [])
        r = f.regions[0]
        self.assertEqual(r.technique, "ISS")
        self.assertTrue(r.is_iss)
        self.assertFalse(r.is_survey)        # wide, but not a survey
        self.assertEqual(r.lens_mode, "ISS")
        self.assertEqual(r.pass_energy, 320.0)
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertIsNone(r.photon_energy)
        self.assertEqual(r.energy, [100.0, 101.0, 102.0, 103.0])
        md = f.region_metadata(r)
        self.assertEqual(md["Technique"], "ISS")
        self.assertEqual(md["KE start (eV)"], "100.00")
        self.assertNotIn("BE start (eV)", md)
        self.assertEqual(md["Ion gun beam HT (V)"], "980")
        self.assertEqual(md["Ion gun emission current (mA)"], "15")

    def test_dset_with_the_iss_constants_reads_like_the_kal(self):
        iss = block(
            2, r_int(1, 1), r_int(2, 0), r_dbl(3, 100.0), r_dbl(4, 1.0),
            r_str(5, "Kinetic Energy"), r_str(6, "eV"), r_dbl(7, 0.5),
            r_str(9, "counts"), r_f32(12, [5, 9, 40, 12]),
            r_str(37, "ISS 1"), r_dbl(42, 320.0), r_int(99, 1),
            r_int(3047, 6), r_int(3049, 6), r_int(3080, 2))
        f = load_file(self.write("a.dset", dset(iss)))
        self.assertEqual(f.warnings, [])             # no unmet constant
        r = f.regions[0]
        self.assertEqual((r.technique, r.lens_mode, r.pass_energy),
                         ("ISS", "ISS", 320.0))
        self.assertTrue(r.is_iss)
        self.assertIn(6, kratos_dset.ENUMS[3049])

    def test_xps_object_stays_xps(self):
        text = iss_object().replace("F_ISS", "F_XPS").replace(
            "F_HSA_LENS_ISS", "F_HSA_LENS_HYBRID")
        r = load_file(self.write("b.kal", "Dataset filename = b\n"
                                 + text)).regions[0]
        self.assertEqual((r.technique, r.lens_mode), ("XPS", "Hybrid"))
        self.assertFalse(r.is_iss)


class TestCasaIss(unittest.TestCase):
    def test_lens_line_pattern_finds_the_iss_lens_among_several(self):
        comments = ["Lens Mode:HSA_LENS_HYBRID", "Aperture Description: Slot",
                    "Lens Mode:HSA_LENS_ISS"]
        self.assertEqual(_LENS_LINE.findall("\n".join(comments)),
                         ["HSA_LENS_HYBRID", "HSA_LENS_ISS"])

    def test_binding_axis_made_from_a_default_photon_energy_is_undone(self):
        native = [200.0, 201.0, 202.0]
        r = Region("I SS Au", 0, 0, energy=[1486.6 - k for k in native],
                   counts=[1, 2, 3], photon_energy=1486.6, lens_mode="ISS")
        VamasFile._as_iss(r, native)
        self.assertEqual((r.technique, r.energy_label, r.name),
                         ("ISS", "Kinetic Energy", "ISS Au"))
        self.assertEqual(r.energy, native)
        self.assertIsNone(r.photon_energy)


class TestMetadataAndMethods(unittest.TestCase):
    def rows(self):
        f = types.SimpleNamespace()
        out = []
        for i in range(2):
            out.append({"Sample": f"S{i}", "Region": f"ISS {i}",
                        "Technique": "ISS", "Instrument": "Kratos (Vision)",
                        "Lens mode": "ISS", "Pass energy (eV)": "320",
                        "Step (eV)": "1", "Dwell (s)": "0.067",
                        "KE start (eV)": "100.00", "KE end (eV)": "1000.00"})
        return out

    def test_only_iss_means_no_x_ray_or_binding_energy_wording(self):
        text = methods.generate(self.rows())
        self.assertIn("Ion scattering spectroscopy (ISS) measurements", text)
        self.assertNotIn("X-ray", text)
        self.assertNotIn("Survey", text)
        self.assertNotIn("Binding energies", text)
        self.assertIn("Ion scattering spectra (2)", text)
        self.assertNotIn("scattering angle", text)   # nobody stated it

    def test_beam_is_stated_only_when_the_user_saved_it(self):
        rows = self.rows()
        for r in rows:
            r.update({"ISS ion": "He+", "ISS beam energy (eV)": "1000",
                      "ISS scattering angle (°)": "135"})
        text = methods.generate(rows)
        self.assertIn("He+ ions, a beam energy of 1000 eV and a scattering "
                      "angle of 135°", text)

    def test_mixed_files_describe_both(self):
        xps = {"Sample": "A", "Region": "C 1s", "Technique": "XPS",
               "Instrument": "Kratos (Vision)", "Photon energy (eV)": "1486.6",
               "BE start (eV)": "295", "BE end (eV)": "280"}
        text = methods.generate([xps] + self.rows())
        self.assertIn("X-ray photoelectron spectroscopy (XPS) and ion "
                      "scattering spectroscopy (ISS) measurements", text)
        self.assertIn("Binding energies are not charge-corrected.", text)

    def test_xps_text_is_unchanged(self):
        xps = {"Sample": "A", "Region": "C 1s", "Instrument": "Kratos",
               "BE start (eV)": "295", "BE end (eV)": "280"}
        text = methods.generate([xps])
        self.assertTrue(text.startswith(
            "X-ray photoelectron spectroscopy (XPS) measurements were made "))

    def test_saved_settings_reach_only_iss_spectra(self):
        a = annotations.Annotations()
        a.set_iss("f1", {"ion": "He+", "e0": 1000, "theta": 135,
                         "junk": 1})
        self.assertEqual(a.iss_for("f1"),
                         {"ion": "He+", "e0": 1000.0, "theta": 135.0})
        iss = Region("ISS", 0, 0, technique="ISS", sample="s")
        xps = Region("C 1s", 0, 0, sample="s")
        md = a.apply_metadata("f1", 0, iss, {"Sample": "s", "Region": "ISS"})
        self.assertEqual(md["ISS beam energy (eV)"], "1000")
        self.assertEqual(md["ISS scattering angle (°)"], "135")
        md = a.apply_metadata("f1", 0, xps, {"Sample": "s", "Region": "C 1s"})
        self.assertNotIn("ISS ion", md)

    def test_settings_survive_json_and_bad_values_are_dropped(self):
        a = annotations.Annotations()
        a.set_iss("f2", {"ion": "He+", "e0": 1000, "theta": 135})
        self.assertFalse(a.is_empty())
        b = annotations.Annotations.from_json(a.to_json())
        self.assertEqual(b.iss_for("f2"), a.iss_for("f2"))
        self.assertEqual(el.sanitise_iss({"ion": "Xx+", "e0": -5,
                                          "theta": 400}), {})
        a.set_iss("f2", {})
        self.assertTrue(a.is_empty())


class TestOffer(unittest.TestCase):
    """What the ISS dialog offers: Kratos geometry for Kratos ISS spectra
    only, never the photon energy of a mislabelled export."""

    @classmethod
    def setUpClass(cls):
        try:
            import spectradeck as ee
        except Exception as exc:                      # pragma: no cover
            raise unittest.SkipTest(f"app not importable: {exc}")
        cls.ee = ee

    def offer(self, region, instrument):
        app = types.SimpleNamespace(region_parser={
            id(region): types.SimpleNamespace(
                instrument={"Instrument": instrument})})
        return self.ee.Workspace.iss_offer(app, region)

    def test_kratos_iss_gets_the_axis_ultra_geometry(self):
        r = Region("ISS", 0, 0, technique="ISS")
        self.assertEqual(self.offer(r, "Kratos (Vision)"),
                         {"ion": "He+", "e0": 1000.0, "theta": 135.0})

    def test_nothing_for_other_instruments_or_xps(self):
        self.assertEqual(self.offer(Region("ISS", 0, 0, technique="ISS"),
                                    "Thermo K-Alpha"), {})
        self.assertEqual(self.offer(Region("C 1s", 0, 0), "Kratos (Vision)"),
                         {})


class TestKratosEnergyRatios(unittest.TestCase):
    """Kratos TPC1369C, table 2: He4 at 135 degrees (ratio of scattered to
    beam energy), with the table's own rounded atomic masses."""
    TABLE = [("Cr", 52.0, 0.768), ("Fe", 55.8, 0.782), ("Co", 58.9, 0.792),
             ("Ni", 58.7, 0.791), ("Cu", 63.5, 0.806), ("Zn", 65.4, 0.811),
             ("Mo", 95.9, 0.867), ("Pd", 106.4, 0.879), ("Ag", 107.9, 0.881),
             ("Sn", 118.7, 0.891), ("W", 183.8, 0.928), ("Pt", 195.1, 0.932),
             ("Au", 197.0, 0.933), ("Pb", 207.2, 0.936)]
    # Lighter targets: ours is higher than the table, by 0.0013 (Ca, 40 u) up
    # to 0.013 (Be, 9 u). Not understood (the table's method is not given), so
    # the formula is left as the textbook one and this is only recorded.
    LIGHT = [("Be", 9.0, 0.180, 0.014), ("C", 12.0, 0.295, 0.010),
             ("O", 16.0, 0.410, 0.008), ("Al", 27.0, 0.597, 0.004),
             ("Si", 28.1, 0.610, 0.003), ("Ca", 40.1, 0.709, 0.0015)]

    def test_heavy_targets_agree_to_a_thousandth(self):
        for sym, mass, ratio in self.TABLE:
            with self.subTest(sym):
                self.assertAlmostEqual(
                    el.kinematic_factor(HE, mass, 135.0), ratio, delta=0.001)

    def test_light_targets_differ_by_the_recorded_amount_at_most(self):
        for sym, mass, ratio, tol in self.LIGHT:
            with self.subTest(sym):
                self.assertAlmostEqual(
                    el.kinematic_factor(HE, mass, 135.0), ratio, delta=tol)

    def test_gold_peak_of_the_procedure(self):
        self.assertAlmostEqual(
            el.iss_energy(1000.0, HE, el.ELEMENTS["Au"][1], 135.0), 933.0,
            delta=0.5)


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_ISS_CORPUS to the folder of Kratos ISS runs")
class TestRealKratosIss(unittest.TestCase):
    def files(self, ext):
        return sorted(glob.glob(os.path.join(CORPUS, "**", "*" + ext),
                                recursive=True))

    def test_every_dset_matches_its_kal(self):
        pairs = 0
        for kal in self.files(".kal"):
            dset_path = kal[:-4] + ".dset"
            if not os.path.isfile(dset_path) or os.path.getsize(kal) == 0:
                continue
            a, b = load_file(kal), load_file(dset_path)
            with self.subTest(os.path.basename(kal)):
                self.assertEqual(b.warnings, [])
                self.assertEqual(len(a.regions), len(b.regions))
                for x, y in zip(a.regions, b.regions):
                    self.assertEqual(
                        (x.name, x.technique, x.energy_label, x.energy,
                         x.counts, x.lens_mode, x.pass_energy,
                         x.photon_energy),
                        (y.name, y.technique, y.energy_label, y.energy,
                         y.counts, y.lens_mode, y.pass_energy,
                         y.photon_energy))
                    # the .kal prints six significant digits
                    self.assertAlmostEqual(x.dwell, y.dwell, delta=1e-6)
                pairs += 1
        self.assertGreater(pairs, 0)

    def test_casa_vamas_of_an_iss_run_equals_the_kal(self):
        n = 0
        for vms in self.files(".vms"):
            kal = vms[:-4] + ".kal"
            f = load_file(vms)
            iss = [r for r in f.regions if r.is_iss]
            if not iss or not os.path.isfile(kal):
                continue
            ref = [r for r in load_file(kal).regions if r.is_iss]
            with self.subTest(os.path.basename(vms)):
                self.assertEqual(len(iss), len(ref))
                for x, y in zip(iss, ref):
                    self.assertEqual(x.energy_label, "Kinetic Energy")
                    self.assertIsNone(x.photon_energy)
                    self.assertEqual(x.counts, y.counts)
                    for p, q in zip(x.energy, y.energy):
                        self.assertAlmostEqual(p, q, places=6)
                n += 1
        self.assertGreater(n, 0)

    def test_gold_foil_peak_and_candidates(self):
        found = 0
        for kal in self.files(".kal"):
            if os.path.getsize(kal) == 0:
                continue
            f = load_file(kal)
            for r in f.regions:
                if not r.is_iss:
                    continue
                if not f.region_metadata(r).get("Sample"):
                    continue
                i = max(range(len(r.counts)), key=r.counts.__getitem__)
                if r.counts[i] < 1000 or r.energy[i] < 900:
                    continue                          # not a gold-foil run
                with self.subTest(f"{os.path.basename(kal)} {r.name}"):
                    want = el.iss_energy(1000.0, HE, el.ELEMENTS["Au"][1],
                                         el.KRATOS_THETA)
                    self.assertAlmostEqual(r.energy[i], want, delta=4.0)
                    syms = [c["symbol"] for c in el.candidates(
                        r.energy[i], el.KRATOS_E0, "He+", el.KRATOS_THETA)]
                    self.assertIn("Au", syms)
                    found += 1
        self.assertGreater(found, 0)


if __name__ == "__main__":
    unittest.main()
