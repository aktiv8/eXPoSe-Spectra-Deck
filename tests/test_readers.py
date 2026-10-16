"""Reader regression tests. Run:  python -m unittest discover tests

Each format is exercised with a small synthetic file built here, so the tests
need no instrument data. (``test_corpus.py`` checks real files if you point
XPS_CORPUS at a folder of them.)
"""

import math
import os
import struct
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from readers import load_file, reader_for, UnsupportedFormat, Region  # noqa: E402
from readers.base import (canon_region_name, guess_region_name,  # noqa: E402
                          kv_from_lines, clean_text, is_survey_span,
                          is_survey_name, is_survey_region)
from exporters import export_vamas  # noqa: E402


def peak(n, lo, hi, centre, amp=1000.0):
    e = [hi - i * (hi - lo) / (n - 1) for i in range(n)]      # descending BE
    y = [50 + amp * math.exp(-((x - centre) ** 2) / 2) for x in e]
    return e, y


class Tmp(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def path(self, name):
        return os.path.join(self.dir.name, name)

    def write(self, name, data):
        mode = "wb" if isinstance(data, bytes) else "w"
        with open(self.path(name), mode, **({} if mode == "wb" else
                                            {"newline": "", "encoding": "latin-1"})) as fh:
            fh.write(data)
        return self.path(name)


class TestNames(unittest.TestCase):
    def test_canonical_region_names(self):
        self.assertEqual(canon_region_name("C1s"), "C 1s")
        self.assertEqual(canon_region_name("Cu2p3/2"), "Cu 2p3/2")
        self.assertEqual(canon_region_name("Fe 2p"), "Fe 2p")
        for s in ("XPS Survey", "wide/1", "SUR", "Survey"):
            self.assertEqual(canon_region_name(s), "Survey")
        self.assertEqual(canon_region_name("Not Specified"), "")

    def test_wide_core_level_is_a_survey(self):
        wide = [0.0, 1300.0]
        self.assertEqual(guess_region_name("C 1s", wide), "Survey")
        self.assertEqual(guess_region_name("C 1s", [280.0, 295.0]), "C 1s")

    def test_is_survey_span(self):
        self.assertTrue(is_survey_span([0.0, 1300.0]))
        self.assertFalse(is_survey_span([280.0, 295.0]))
        self.assertFalse(is_survey_span([280.0]))
        self.assertFalse(is_survey_span(None))
        self.assertFalse(is_survey_span([]))

    def test_is_survey_name(self):
        for s in ("XPS Survey", "wide/1", "SUR", "Survey"):
            self.assertTrue(is_survey_name(s))
        self.assertFalse(is_survey_name("C 1s"))
        self.assertFalse(is_survey_name(""))

    def test_is_survey_region(self):
        # span only
        r = Region("C 1s", 0, 0, energy=[0.0, 1300.0])
        self.assertTrue(is_survey_region(r))
        # name only
        r = Region("Survey", 0, 0, energy=[280.0, 295.0])
        self.assertTrue(is_survey_region(r))
        # both
        r = Region("Survey", 0, 0, energy=[0.0, 1300.0])
        self.assertTrue(is_survey_region(r))
        # neither
        r = Region("C 1s", 0, 0, energy=[280.0, 295.0])
        self.assertFalse(is_survey_region(r))
        self.assertFalse(r.is_survey)
        self.assertTrue(Region("Survey", 0, 0, energy=[0.0, 1300.0]).is_survey)

    def test_kv_from_lines(self):
        kv = kv_from_lines(["Step(meV): 700.0    Dwell(ms): 76   Sweeps: 2",
                            "Lens Mode:Hybrid", "no colon here"])
        self.assertEqual(kv["dwell(ms)"], "76")
        self.assertEqual(kv["lens mode"], "Hybrid")
        self.assertEqual(clean_text(" Not Specified "), "")


class TestVamas(Tmp):
    def _regions(self):
        out = []
        for s, c in (("S1", 285.0), ("S2", 286.0)):
            e, y = peak(101, 280, 295, c)
            out.append(Region(name="C 1s", index=len(out), offset=0, energy=e,
                              counts=y, decodable=True, sample=s,
                              photon_energy=1486.6, pass_energy=20.0,
                              dwell=0.1, step=0.15, anode="Al Ka"))
        e, y = peak(61, 525, 540, 532)
        out.append(Region(name="O 1s", index=2, offset=0, energy=e, counts=y,
                          decodable=True, sample="S1", photon_energy=1486.6,
                          pass_energy=20.0, dwell=0.1, step=0.25))
        return out

    def test_round_trip_through_our_exporter(self):
        regs = self._regions()
        p = self.path("out.vms")
        export_vamas(regs, p, include_transmission=False)
        back = load_file(p)
        self.assertEqual(back.format_name, "VAMAS (ISO 14976)")
        self.assertEqual(len(back.regions), 3)
        for a, b in zip(regs, back.regions):
            self.assertEqual(a.n_points, b.n_points)
            self.assertEqual(b.name, a.name)
            self.assertEqual(b.energy_label, "Binding Energy")
            self.assertAlmostEqual(b.photon_energy, 1486.6, places=3)
            self.assertAlmostEqual(b.pass_energy, 20.0)
            self.assertLess(max(abs(x - y) for x, y in zip(a.energy, b.energy)), 1e-3)
            self.assertLess(max(abs(x - y) for x, y in zip(a.counts, b.counts)), 1e-3)

    def test_transmission_is_read_as_second_variable(self):
        regs = self._regions()[:1]
        r = regs[0]
        ke = [r.photon_energy - x for x in r.energy]
        r.tf_ke = [min(ke) - 5, max(ke) + 5]
        r.tf_values = [1.0, 0.5]
        p = self.path("tf.vms")
        export_vamas(regs, p, include_transmission=True)
        back = load_file(p).regions[0]
        self.assertIsNotNone(back.transmission())
        self.assertEqual(len(back.transmission()), back.n_points)

    def test_truncated_file_keeps_the_good_blocks(self):
        regs = self._regions()
        p = self.path("cut.vms")
        export_vamas(regs, p, include_transmission=False)
        with open(p, encoding="latin-1", newline="") as fh:
            text = fh.read()
        cut = text[:int(len(text) * 0.75)]
        self.write("cut2.vms", cut)
        back = load_file(self.path("cut2.vms"))
        self.assertGreaterEqual(len(back.regions), 1)
        self.assertTrue(back.warnings)

    def test_not_vamas(self):
        p = self.write("junk.vms", "hello\nworld\n")
        with self.assertRaises(Exception):
            load_file(p)


AVG_TEMPLATE = """;================
;Dump of DataSpace 'F:\\x\\{title}.VGD'
; on 12/5/2008 at 14:26:48
;================
$FORMAT=3
$PROPERTIES=SUM
DS_EXT_SUPROPID_TITLE       : VT_BSTR = '{title}'
DS_EXT_SUPROPID_CREATED     : VT_DATE = 12/5/2008   17:34:00
$PROPERTIES=STD
DS_GEPROPID_INSTRUMENT                      : VT_BSTR = 'Alpha'
DS_GEPROPID_VALUE_LABEL                     : VT_BSTR = 'Counts'
DS_SOPROPID_MONO                            : VT_BOOL = True
DS_SOPROPID_ENERGY                          : VT_R4   = 1486.680054
DS_SOPROPID_VOLTAGE                         : VT_R4   = 12000.000000
DS_SOPROPID_CURRENT                         : VT_R4   = 0.002100
DS_ACPROPID_ACQ_TIME                        : VT_R4   = 0.050000
DS_ACPROPID_EV_SCALE                        : VT_I2   = 1
DS_ANPROPID_PASS                            : VT_R4   = 50.000000
DS_ANPROPID_LENS_MODE_NAME                  : VT_BSTR = 'Standard'
$DATAAXES=1,#empty#
    0=        0,        7,          1
$SPACEAXES=1
    0=    1188.600000,       0.100000,        8,  ENERGY,   LINEAR,  'E',   'eV',    'Energy'
$DATA=*
{rows}
"""


class TestThermoAvg(Tmp):
    def test_energy_axis_is_kinetic_and_converted(self):
        vals = [100 + i for i in range(8)]
        rows = "\n".join(
            "LIST@ %3d=  %s" % (k, ",  ".join("%.6f" % v for v in vals[k:k + 4]))
            for k in (0, 4))
        p = self.write("C1s Scan.avg", AVG_TEMPLATE.format(title="C1s Scan",
                                                          rows=rows))
        f = load_file(p)
        r = f.regions[0]
        self.assertEqual(r.name, "C 1s")
        self.assertEqual(r.energy_label, "Binding Energy")
        self.assertAlmostEqual(r.photon_energy, 1486.68, places=3)
        self.assertAlmostEqual(r.energy[0], 1486.68 - 1188.6, places=6)
        self.assertAlmostEqual(r.energy[-1], 1486.68 - 1189.3, places=6)
        self.assertEqual(r.counts, [float(v) for v in vals])
        self.assertEqual(r.pass_energy, 50.0)
        self.assertEqual(r.date, "2008-05-12 17:34:00")            # D/M/Y
        self.assertAlmostEqual(float(r.conditions["X-ray Power"].split()[0]), 25.2, 1)

    def test_ion_gun_properties_fill_the_sputter_hint(self):
        vals = [100 + i for i in range(8)]
        rows = "\n".join(
            "LIST@ %3d=  %s" % (k, ",  ".join("%.6f" % v for v in vals[k:k + 4]))
            for k in (0, 4))
        extra = (
            "DS_DEPTHPROFILE_IONGUNPROPID_CURRENT        : VT_R4   = 10.000000\n"
            "DS_DEPTHPROFILE_IONGUNPROPID_ENERGY         : VT_R4   = 2000.000000\n"
            "DS_DEPTHPROFILE_IONGUNPROPID_RASTER_WIDTH   : VT_R4   = 1.000000\n"
            "DS_DEPTHPROFILE_IONGUNPROPID_RASTER_HEIGHT  : VT_R4   = 0.000000\n"
            "DS_DEPTHPROFILE_IONGUNPROPID_SPUTTERRATE    : VT_R4   = 2.204846\n"
            "DS_DEPTHPROFILE_IONGUNPROPID_IONTYPE        : VT_BSTR = 'Ar+'\n")
        text = AVG_TEMPLATE.format(title="C1s Scan", rows=rows)
        text = text.replace("$DATAAXES=", extra + "$DATAAXES=")
        f = load_file(self.write("C1s Scan.avg", text))
        hint = f.sputter_hint
        self.assertEqual(hint["ion"], "Ar+")
        self.assertEqual(hint["energy_ev"], 2000.0)
        self.assertEqual(hint["current"], 10.0)
        self.assertEqual(hint["current_unit"], "µA")
        self.assertEqual((hint["raster_x"], hint["raster_y"]), (1.0, 1.0))
        self.assertAlmostEqual(hint["etch_rate"], 2.204846)
        self.assertEqual(hint["rate_unit"], "nm/min")

    def test_header_only_dump_is_no_data_not_an_error(self):
        rows = "\n".join("LIST@ %3d=  #empty#,  #empty#,  #empty#,  #empty#"
                         % k for k in (0, 4))
        f = load_file(self.write("Empty.avg", AVG_TEMPLATE.format(
            title="O1s Scan", rows=rows)))
        self.assertFalse(f.regions[0].decodable)
        self.assertIn("Header only", f.regions[0].note)


class TestPhiSpe(Tmp):
    def test_multi_region_layout(self):
        defs = [("C1s", 6, 5, -0.5, 290.0, 1.5, 23.5),
                ("Si2p", 14, 4, -0.5, 104.0, 0.5, 23.5)]
        head = ["SOFH", "Technique: XPS", "FileDesc: test",
                "XraySource: Al 1486.6 mono", "XrayPower: 25.3W",
                "AnalyserWorkFcn: 4.2 eV", "NoSpectralReg: %d" % len(defs)]
        for i, (n, z, npts, step, start, dwell, pe) in enumerate(defs, 1):
            head.append("SpectralRegDef: %d 1 %s %d %d %.3f %.3f %.3f %.3f "
                        "%.3f %.3f %.2f none" % (i, n, z, npts, step, start,
                                                 start + step * (npts - 1), start,
                                                 start + step * (npts - 1), dwell, pe))
        head += ["NoSpatialArea: 1", "EOFH"]
        blob = "\r\n".join(head).encode("latin-1") + b"\r\n"
        data = [[10.0 + i + 100 * k for i in range(d[2])]
                for k, d in enumerate(defs)]
        blob += struct.pack("<4i", 1, len(defs), 96 * len(defs), 16)
        blob += b"\0" * (96 * len(defs))
        for col in data:
            blob += struct.pack("<%dd" % len(col), *col)
        f = load_file(self.write("t.SPE", blob))
        self.assertEqual([r.name for r in f.regions], ["C 1s", "Si 2p"])
        self.assertEqual(f.regions[1].counts, data[1])
        self.assertAlmostEqual(f.regions[0].energy[0], 290.0)
        self.assertAlmostEqual(f.regions[0].energy[-1], 288.0)
        self.assertEqual(f.regions[0].photon_energy, 1486.6)
        self.assertEqual(f.regions[1].dwell, 0.5)
        self.assertEqual(f.regions[0].count_units, "counts/s")

    def test_wrong_length_is_rejected(self):
        blob = b"SOFH\r\nNoSpectralReg: 1\r\nSpectralRegDef: 1 1 C1s 6 50 -0.1 " \
               b"290 285 290 285 1.0 20 none\r\nEOFH\r\n" + b"\0" * 40
        with self.assertRaises(ValueError):
            load_file(self.write("bad.spe", blob))

    # -- block headers that say how the data are stored (Quantes, K-Alpha+) --
    @staticmethod
    def phi_blob(regions, order="<", tag=b"f8", source="Al 1486.6 mono",
                 head_extra=(), date="2022 8 8", pad=0, desc="test"):
        """A PHI .spe whose 96-byte block headers carry the data type (+72),
        length (+76) and offset (+80), as the real files do. ``regions`` is
        [(name, [counts], step, start)]; ``pad`` bytes sit before the data."""
        code, size = {b"f8": ("d", 8), b"f4": ("f", 4)}.get(tag, ("d", 8))
        head = ["SOFH", "Technique: XPS", "FileDesc: " + desc,
                "FileDate: " + date, "AcqFileDate: " + date,
                "XraySource: " + source, "XrayPower: 25.0 W",
                "XrayBeamDiameter: 100.0 um", "AnalyserMode: FAT",
                "SourceAnalyserAngle: 45.0 d",
                "IntensityCalCoeff: 107.072 0.591",
                "AnalyserWorkFcn: 4.45 eV", "NoSpectralRegFull: 1",
                "SpectralRegDefFull: 1 1 Su1s 111 9 -1.0 8.0 0.0 8.0 0.0 "
                "0.0 280.00 TOTAL", "NoSpectralReg: %d" % len(regions)]
        for i, (n, c, step, start) in enumerate(regions, 1):
            head += ["SpectralRegDef: %d %d %s 40 %d %.4f %.4f %.4f %.4f "
                     "%.4f %.6f %.2f TOTAL" % (
                         i, i, n, len(c), step, start, start + step * (len(c) - 1),
                         start, start + step * (len(c) - 1), 0.5, 23.5),
                     "SpectralRegDef2: %d 40.0 3 0 6 1" % i]
        head += ["NoSpatialArea: 1", *head_extra, "EOFH"]
        blob = "\r\n".join(head).encode("latin-1") + b"\r\n"
        n = len(regions)
        blob += struct.pack(order + "4i", 1, n, 96 * n, 16)
        offset = 16 + 96 * n + pad
        for i, (name, c, _s, _st) in enumerate(regions, 1):
            h = bytearray(96)
            struct.pack_into(order + "i", h, 0, i)
            struct.pack_into(order + "i", h, 20, len(c))
            h[32:35], h[40:43], h[56:59] = b"chn", b"sar", b"c/s"
            h[72:72 + len(tag)] = tag
            struct.pack_into(order + "i", h, 76, len(c) * size)
            struct.pack_into(order + "i", h, 80, offset)
            blob += bytes(h)
            offset += len(c) * size
        blob += b"\x7f" * pad
        for _n, c, _s, _st in regions:
            blob += struct.pack("%s%d%s" % (order, len(c), code), *c)
        return blob

    REGIONS = [("Su1s", [10.5 + i for i in range(9)], -1.0, 8.0),
               ("Zr2p3", [100.25 + i for i in range(5)], -0.125, 2255.0)]

    def test_float32_data_type_from_the_block_header(self):
        f = load_file(self.write("q.spe", self.phi_blob(self.REGIONS,
                                                       tag=b"f4")))
        self.assertEqual(f.warnings, [])
        self.assertEqual(f.regions[1].counts, self.REGIONS[1][1])
        self.assertEqual(f.regions[0].counts, self.REGIONS[0][1])

    def test_big_endian_file(self):
        f = load_file(self.write("be.spe", self.phi_blob(self.REGIONS,
                                                        order=">")))
        self.assertEqual(f.regions[1].counts, self.REGIONS[1][1])
        f = load_file(self.write("be4.spe", self.phi_blob(
            self.REGIONS, order=">", tag=b"f4")))
        self.assertEqual(f.regions[0].counts, self.REGIONS[0][1])

    def test_data_offset_comes_from_the_block_header(self):
        f = load_file(self.write("pad.spe", self.phi_blob(
            self.REGIONS, tag=b"f4", pad=24)))
        self.assertEqual(f.warnings, [])
        self.assertEqual(f.regions[1].counts, self.REGIONS[1][1])

    def test_full_and_def2_lines_are_not_regions(self):
        f = load_file(self.write("full.spe", self.phi_blob(self.REGIONS)))
        self.assertEqual(len(f.regions), 2)
        self.assertEqual([r.n_points for r in f.regions], [9, 5])

    def test_phi_region_names(self):
        f = load_file(self.write("n.spe", self.phi_blob(self.REGIONS)))
        self.assertEqual([r.name for r in f.regions], ["Survey", "Zr 2p3"])

    def test_unknown_data_type_is_refused(self):
        blob = self.phi_blob(self.REGIONS, tag=b"f8")
        i = blob.find(b"EOFH") + 6 + 16 + 72
        blob = blob[:i] + b"i2" + blob[i + 2:]
        with self.assertRaises(ValueError) as cm:
            load_file(self.write("t.spe", blob))
        self.assertIn("unsupported data type", str(cm.exception))

    def test_implausible_values_are_refused_not_returned(self):
        bad = [("C1s", [1e300, 2e300, 3e300], -1.0, 3.0)]
        with self.assertRaises(ValueError) as cm:
            load_file(self.write("bad.spe", self.phi_blob(bad)))
        self.assertIn("plausible", str(cm.exception))

    def test_source_with_an_element_names_the_anode(self):
        for text, anode, hv in (("Al 1486.6 mono", "Al mono", 1486.6),
                                ("Cr 5414.8 mono", "Cr mono", 5414.8)):
            f = load_file(self.write("s.spe", self.phi_blob(
                self.REGIONS, source=text)))
            r = f.regions[0]
            self.assertEqual((r.anode, r.photon_energy), (anode, hv))
            self.assertEqual(f.instrument["X-ray source"],
                             "%s (%g eV)" % (anode, hv))

    def test_energy_only_source_names_no_anode(self):
        """Thermo's ``KA1486 1486.68``: the energy is the whole number, no
        anode is made up, and nothing reads "KA1486"."""
        import methods
        f = load_file(self.write("ka.spe", self.phi_blob(
            self.REGIONS, source="KA1486 1486.68", desc="CasaXPS MULTIPAK output")))
        r = f.regions[0]
        self.assertEqual((r.anode, r.photon_energy), ("", 1486.68))
        self.assertNotIn("X-ray source", f.instrument)
        md = f.region_metadata(r)
        text = methods._source_sentence([md])
        self.assertNotIn("KA1486", text)
        self.assertIn("1486.68", text)
        self.assertEqual(f.instrument["Acquisition software"],
                         "CasaXPS (MultiPak-format output)")
    def test_dates_with_any_spacing(self):
        for text, want in (("2022 8 8", "2022-08-08"),
                           ("2025 10  6", "2025-10-06"),
                           ("2009 03 19", "2009-03-19")):
            f = load_file(self.write("d.spe", self.phi_blob(self.REGIONS,
                                                           date=text)))
            self.assertEqual(f.regions[0].date, want)

    def test_recorded_settings_become_metadata(self):
        f = load_file(self.write("m.spe", self.phi_blob(self.REGIONS)))
        r = f.regions[0]
        md = f.region_metadata(r)
        self.assertEqual(md["Analyser mode"], "Constant analyser energy (CAE)")
        self.assertEqual(md["X-ray spot (µm)"], "100")
        self.assertIn("107.072, 0.591", md["Intensity calibration (a, b)"])
        self.assertEqual(md["Source–analyser angle (°)"], "45")
        self.assertNotIn("Take-off angle (°)", md)  # a different quantity
        self.assertIsNone(r.tf_ke)              # recorded, never applied
        self.assertEqual(r.counts, self.REGIONS[0][1])

    def test_npl_is_offered_beside_vms(self):
        import readers
        vamas = dict(readers.supported_patterns())["VAMAS (ISO 14976)"]
        self.assertIn("*.npl", vamas)
        try:
            import spectradeck as ee
        except Exception as exc:                # noqa: BLE001 - no Tk/display
            self.skipTest(f"app not importable here: {exc}")
        types = ee.Workspace._open_filetypes()
        self.assertIn("*.npl", types[0][1])              # all supported
        self.assertTrue(any(n.startswith("VAMAS") and "*.npl" in p
                            for n, p in types))

    def test_nul_in_a_vamas_species_label(self):
        """PHI MultiPak's ``.npl`` export ends each species label with a NUL."""
        e, y = peak(21, 280, 290, 285)
        p = self.path("nul.vms")
        export_vamas([Region(name="C 1s", index=1, offset=0, energy=e,
                             counts=y, decodable=True, sample="S1",
                             photon_energy=1486.6, pass_energy=20.0,
                             dwell=0.1, step=0.5)], p,
                     include_transmission=False)
        with open(p, "rb") as fh:
            blob = fh.read()
        blob = blob.replace(b"C 1s", b"C 1s\x00")
        self.assertIn(b"C 1s\x00", blob)
        with open(self.path("nul.npl"), "wb") as fh:
            fh.write(blob)
        self.assertEqual(load_file(self.path("nul.npl")).regions[0].name,
                         "C 1s")


class TestScienta(Tmp):
    def test_regions_and_units(self):
        e = [545.0 - 0.5 * i for i in range(6)]
        y = [100.0 + i for i in range(6)]
        txt = ["[Info]", "Number of Regions=1", "Version=1.3.1", "",
               "[Region 1]", "Region Name=O1s",
               "Dimension 1 name=Binding Energy [eV]", "Dimension 1 size=6",
               "Dimension 1 scale=" + " ".join("%.5f" % x for x in e), "",
               "[Info 1]", "Excitation Energy=700.0", "Pass Energy=100",
               "Step Time=200", "Number of Sweeps=20", "Lens Mode=Transmission",
               "Date=2011-08-23", "Time=22:22:44", "", "[Data 1]"]
        txt += [" %.6E  %.6E" % (a, b) for a, b in zip(e, y)]
        f = load_file(self.write("s.txt", "\n".join(txt) + "\n"))
        r = f.regions[0]
        self.assertEqual(r.name, "O 1s")
        self.assertEqual(r.energy, e)
        self.assertEqual(r.counts, y)
        self.assertEqual(r.photon_energy, 700.0)
        self.assertAlmostEqual(r.dwell, 0.2)
        self.assertEqual(r.date, "2011-08-23 22:22:44")

    def test_detector_columns_are_summed(self):
        txt = ["[Info]", "Number of Regions=1", "", "[Region 1]",
               "Region Name=Au4f", "Dimension 1 name=Binding Energy [eV]",
               "Dimension 1 size=3", "Dimension 1 scale=90 89 88", "",
               "[Info 1]", "Excitation Energy=700", "", "[Data 1]",
               "90 1 2 3", "89 4 5 6", "88 7 8 9"]
        r = load_file(self.write("d.txt", "\n".join(txt))).regions[0]
        self.assertEqual(r.counts, [6.0, 15.0, 24.0])
        self.assertIn("summed", r.note)


KAL = """Dataset filename          = t.dset
Object name               = hy160/1
   1 Technique                 = F_XPS
   2 Scan type                 = F_SPECTRUM
   3 Spectrum scan start       = 100 eV
   4 Spectrum scan step size   = 2 eV
   5 Abscissa label            = Kinetic Energy
   7 Dwell time                = 1 seconds
   9 Ordinate units            = Counts
  12 Ordinate values           = {10, 11,
12, 13}
  42 Pass energy               = 160 eV
 151 Date Acquired             = 98/07/07 19:11:35
3053 Xray Gun Anode            = F_MONO_ANODE
3080 Xray Reference Energy     = F_REFER_TO_XRAY_AG
3113 Chemical symbol or formula = Mo
3114 Transition or charge state = 3d
5587 Transmission Function Object (ke,t) = {
    5615 Transmission Function Kinetic Energy = {100, 110}
    5616 Transmission Function Value = {1.0, 0.5}
}
"""


class TestKratosKal(Tmp):
    def test_object_fields(self):
        f = load_file(self.write("t.kal", KAL))
        r = f.regions[0]
        self.assertEqual(r.name, "Mo 3d")
        self.assertEqual(r.counts, [10.0, 11.0, 12.0, 13.0])          # multi-line list
        self.assertAlmostEqual(r.photon_energy, 2984.2)
        self.assertAlmostEqual(r.energy[0], 2984.2 - 100)
        self.assertEqual(r.date, "1998-07-07 19:11:35")
        self.assertEqual(r.tf_ke, [100.0, 110.0])
        self.assertEqual(r.pass_energy, 160.0)

    @staticmethod
    def _scans(n, reference):
        """``n`` kinetic-energy scans whose X-ray reference line is
        ``reference`` ('' leaves the field out)."""
        out = ["Dataset filename = t.dset"]
        for i in range(n):
            out += [f"Object name = pe{i}/{i + 1}",
                    "   3 Spectrum scan start = 100 eV",
                    "   4 Spectrum scan step size = 1 eV",
                    "   5 Abscissa label = Kinetic Energy",
                    "   7 Dwell time = 1 seconds",
                    "  12 Ordinate values = {1, 2, 3}"]
            if reference:
                out.append(f"3080 Xray Reference Energy = {reference}")
        return "\n".join(out) + "\n"

    def test_reference_none_is_a_kinetic_scan_not_a_missing_source(self):
        # Kratos transmission test: the file itself says there is no X-ray
        # reference, so the kinetic-energy axis is as recorded
        f = load_file(self.write("n.kal", self._scans(3, "F_REFER_TO_NONE")))
        self.assertEqual(len(f.regions), 3)
        self.assertEqual({r.energy_label for r in f.regions}, {"Kinetic Energy"})
        self.assertEqual(f.warnings, [])

    def test_missing_source_warns_once_for_the_file(self):
        for ref in ("", "F_REFER_TO_SOMETHING_ELSE"):
            f = load_file(self.write("m.kal", self._scans(30, ref)))
            self.assertEqual(len(f.warnings), 1, ref)
            self.assertIn("X-ray energy unknown for 30 spectra", f.warnings[0])
            self.assertIn("pe0/1", f.warnings[0])         # names a few
            self.assertIn("(30 in all)", f.warnings[0])
        one = load_file(self.write("o.kal", self._scans(1, "")))
        self.assertIn("for 1 spectrum (pe0/1)", one.warnings[0])

    def test_known_anode_still_gives_binding_energy_and_no_warning(self):
        f = load_file(self.write("a.kal",
                                 self._scans(2, "F_REFER_TO_XRAY_MONO_AL")))
        self.assertEqual({r.energy_label for r in f.regions}, {"Binding Energy"})
        self.assertEqual(f.warnings, [])


class TestRegistry(Tmp):
    def test_content_sniffing_beats_extension(self):
        e, y = peak(21, 280, 290, 285)
        p = self.path("real.vms")
        export_vamas([Region(name="C 1s", index=0, offset=0, energy=e, counts=y,
                             decodable=True, photon_energy=1486.6, step=0.5,
                             pass_energy=20, dwell=0.1)], p,
                     include_transmission=False)
        q = self.path("renamed.dat")
        os.replace(p, q)
        self.assertEqual(reader_for(q).format_name, "VAMAS (ISO 14976)")

    def test_unknown_file(self):
        with self.assertRaises(UnsupportedFormat):
            reader_for(self.write("x.bin", b"\x00\x01\x02"))


if __name__ == "__main__":
    unittest.main()
