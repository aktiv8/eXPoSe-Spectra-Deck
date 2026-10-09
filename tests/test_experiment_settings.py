"""What the Kratos ESCApe ``.experiment`` reader takes from the acquisition
records of each spectrum: the charge neutraliser (filament current, charge
balance, bias), the X-ray emission current, and the stage position x / y / z,
plus the operator.

The records are read from the last one before the spectrum. Synthetic bytes
test the rules everywhere; the real-file test needs HarwellXPS's own export of
the same experiment and compares every region:

    set XPS_EXPERIMENT_CORPUS=<folder with MI-LD-20264066-26-21.experiment
                               and its "... data.vms">
    python -m unittest tests.test_experiment_settings

(checked on all 144 regions: neutraliser, emission current, stage x / y / z,
power and operator all agree).

Run:  python -m unittest discover tests
"""

import os
import re
import struct
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import methods  # noqa: E402
import readers  # noqa: E402
from readers import Region  # noqa: E402
from readers.kratos_experiment import EscapeParser  # noqa: E402

CORPUS = os.environ.get("XPS_EXPERIMENT_CORPUS", "")
NEU, XRY, AN2 = (EscapeParser.NEUTRALISER_MARKER, EscapeParser.XRAY_MARKER,
                 EscapeParser.ANALYSIS_PAIR)


def neutraliser(fil=0.45, bal=5.0, bias=1.4, flag=1):
    return (NEU + struct.pack("<i", 1) + bytes([flag])
            + struct.pack("<ddd", fil, bal, bias) + b"\x00" * 8)


def xray(emission=0.015):
    return XRY + b"\x02\x00\x00\x00\x01\x0cSpectroscopy" + struct.pack(
        "<d", emission) + b"\x00" * 16


def stage(x=0.0422982, y=0.0182745, z=-0.0002697):
    return AN2 + struct.pack("<ddd", x, y, z) + b"\x00" * 8


def parser(raw):
    p = EscapeParser()
    p.raw = raw
    p.corruption = {"corrupted": False}
    return p


class TestRecords(unittest.TestCase):
    def read(self, raw, off=None):
        p = parser(raw)
        reg = Region("C 1s", 0, 0)
        p._acquisition_records(reg, len(raw) if off is None else off)
        return reg

    def test_every_field_is_read(self):
        reg = self.read(neutraliser() + xray() + stage())
        self.assertEqual(reg.extra["neutraliser"],
                         "on (filament current 0.45 A, charge balance 5 V, "
                         "bias 1.4 V)")
        self.assertEqual(reg.conditions["Emission current (mA)"], "15")
        self.assertAlmostEqual(reg.pos_x, 42.2982)
        self.assertAlmostEqual(reg.pos_y, 18.2745)
        self.assertAlmostEqual(reg.extra["pos_z"], -0.2697)
        self.assertTrue(reg.extra["stage_own"])

    def test_the_last_record_before_the_spectrum_wins(self):
        raw = (neutraliser(0.45) + stage(0.01) + b"XX" + neutraliser(0.30, 3.0, 1.1)
               + xray(0.015) + stage(0.02) + b"SPECTRUM" + neutraliser(0.55))
        off = raw.index(b"SPECTRUM")
        reg = self.read(raw, off)
        self.assertIn("filament current 0.3 A", reg.extra["neutraliser"])
        self.assertAlmostEqual(reg.pos_x, 20.0)          # the second stage record
        # and nothing later than the spectrum is used
        self.assertNotIn("0.55", reg.extra["neutraliser"])

    def test_a_missing_or_implausible_record_leaves_the_field_out(self):
        self.assertEqual(self.read(b"nothing here").extra, {})
        reg = self.read(neutraliser(flag=0))             # flag not "on"
        self.assertNotIn("neutraliser", reg.extra)
        reg = self.read(neutraliser(fil=float("nan")))
        self.assertNotIn("neutraliser", reg.extra)
        reg = self.read(neutraliser(fil=-1.0))
        self.assertNotIn("neutraliser", reg.extra)
        self.assertNotIn("Emission current (mA)",
                         self.read(xray(5.0)).conditions)     # 5 A: not an emission
        self.assertNotIn("Emission current (mA)",
                         self.read(xray(0.0)).conditions)
        reg = self.read(stage(0.0, 0.0, 0.0))
        self.assertNotIn("pos_z", reg.extra)
        reg = self.read(stage(9.0, 0.0, 0.0))             # 9 m: not a stage
        self.assertIsNone(reg.pos_x)

    def test_a_record_cut_short_by_the_end_of_the_file_is_not_read(self):
        raw = neutraliser()[:len(NEU) + 12]
        self.assertNotIn("neutraliser", self.read(raw).extra)

    def test_a_damaged_file_gives_nothing(self):
        p = parser(neutraliser() + xray() + stage())
        p.corruption = {"corrupted": True}
        reg = Region("C 1s", 0, 0)
        p._acquisition_records(reg, len(p.raw))
        self.assertEqual((reg.extra, reg.conditions), ({}, {}))


class TestWording(unittest.TestCase):
    def rows(self, n):
        return [{"Sample": "s", "Region": "C 1s", "Instrument": "x",
                 "Charge neutraliser": f"on (filament current {0.1 * i:g} A)"}
                for i in range(1, n + 1)]

    def test_a_few_settings_are_listed(self):
        text = methods.generate(self.rows(2))
        self.assertIn("filament current 0.1 A", text)
        self.assertIn("filament current 0.2 A", text)

    def test_many_settings_are_counted_not_listed(self):
        text = methods.generate(self.rows(9))
        self.assertIn("9 different settings were used (listed in the "
                      "metadata)", text)
        self.assertNotIn("filament current", text)


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_EXPERIMENT_CORPUS to the folder of the experiment "
                     "and its Harwell export")
class TestAgainstHarwell(unittest.TestCase):
    def test_every_region_agrees_with_harwells_export(self):
        f = readers.load_file(os.path.join(
            CORPUS, "MI-LD-20264066-26-21.experiment"))
        v = readers.load_file(os.path.join(
            CORPUS, "MI-LD-20264066-26-21 data.vms"))
        self.assertEqual(len(f.regions), len(v.regions))
        self.assertGreater(len(f.regions), 100)
        self.assertEqual(f.instrument["Operator"],
                         v.regions[0].extra["comments"]["operator"])
        for r, h in zip(f.regions, v.regions):
            c, md = h.extra["comments"], f.region_metadata(r)
            with self.subTest(r.index):
                self.assertEqual(
                    md["Charge neutraliser"],
                    "on (filament current "
                    f"{float(c['neutraliser filament (a)']):g} A, charge "
                    "balance "
                    f"{float(c['neutraliser charge balance (v)']):g} V, bias "
                    f"{float(c['neutraliser bias (v)']):g} V)")
                self.assertAlmostEqual(
                    float(md["Emission current (mA)"]),
                    1000 * float(c["emission current (a)"]), places=6)
                self.assertEqual(md["Source power (W)"] + "W",
                                 c["x-ray power"])
                x, y, z = (float(g) for g in re.match(
                    r"x=(-?[\d.]+) y=(-?[\d.]+) z=(-?[\d.]+)",
                    c["stage position (mm)"]).groups())
                self.assertAlmostEqual(float(md["Position X (mm)"]), x, delta=6e-4)
                self.assertAlmostEqual(float(md["Position Y (mm)"]), y, delta=6e-4)
                self.assertAlmostEqual(float(md["Position Z (mm)"]), z, delta=6e-4)


if __name__ == "__main__":
    unittest.main()
