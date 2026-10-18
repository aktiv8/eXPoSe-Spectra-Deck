"""Thermo Avantage stage-tilt (angle-resolved) series: the second space axis of
type code 2 / 'ANGLE' (0, 15 ... 75 degrees in the K-Alpha+ files it was
found in) is a sample tilt, not an etch level -- read the same from the
binary ``.VGD`` and the ``.avg`` dump, shown as "Sample tilt (°)", a z axis of
its own, a methods sentence and kept through a VAMAS round trip.

``XPS_THERMO_TILT_CORPUS=<folder of .VGD/.avg pairs>`` checks real files.

Run:  python -m unittest tests.test_thermo_tilt
"""
import glob
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import exporters  # noqa: E402
import methods  # noqa: E402
import viewdata  # noqa: E402
from readers import load_file  # noqa: E402
from readers.thermo_avg import AXIS_CODES  # noqa: E402
from test_thermo_kinds import Tmp, dump  # noqa: E402

TILTS = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0]


def tilt_text(ne=4):
    axes = [(1188.0, 0.5, ne, "ENERGY", "LINEAR", "E", "eV", "Energy"),
            (0.0, 15.0, len(TILTS), "ANGLE", "LINEAR", "\xb0", "\xb0",
             "Angle")]
    blocks = []
    for k, t in enumerate(TILTS):
        labels = ["$AXISVALUE= DATAXIS=1 SPACEAXIS=1 LABEL='Angle' "
                  "POINT=%d VALUE=%f;" % (k, t)]
        blocks.append(((k,), labels, [100.0 - k + c for c in range(ne)]))
    return dump("Si2p Scan", axes, [(0, ne - 1, 1), (0, len(TILTS) - 1, 1)],
                blocks)


class TestTiltSeries(Tmp):
    def load(self):
        return load_file(self.write("Si2p Scan.avg", tilt_text()))

    def test_code_two_is_the_angle_axis(self):
        self.assertEqual(AXIS_CODES[2], "ANGLE")

    def test_each_spectrum_carries_its_tilt(self):
        f = self.load()
        self.assertEqual([r.extra["tilt"] for r in f.regions], TILTS)
        md = f.region_metadata(f.regions[3])
        self.assertEqual(md["Sample tilt (\xb0)"], "45")

    def test_it_is_not_called_an_etch_level_or_a_depth_profile(self):
        f = self.load()
        for r in f.regions:
            self.assertNotIn("Etch level", f.region_metadata(r))
        self.assertIn("not sputtering", f.depth_profile["etch_source"])

    def test_z_axis_offers_tilt_and_auto_picks_it(self):
        f = self.load()
        self.assertIn("Tilt", viewdata.Z_MODES)
        zi = viewdata.resolve_z(f.regions)
        self.assertEqual((zi.mode, zi.label), ("Tilt", "Tilt (\xb0)"))
        self.assertEqual(zi.values, TILTS)
        self.assertEqual(viewdata.resolve_z(f.regions, "Etch level").mode,
                         "Etch level")

    def test_methods_say_tilt_not_sputtering(self):
        f = self.load()
        rows = [f.region_metadata(r) for r in f.regions]
        s = methods._tilt_sentence(rows)
        self.assertIn("6 sample tilts from 0\xb0 to 75\xb0 in 15\xb0 steps", s)
        self.assertNotIn("sputter", methods._depth_sentence(rows).lower())

    def test_a_single_tilt_says_nothing(self):
        self.assertEqual(methods._tilt_sentence(
            [{"Sample": "a", "Sample tilt (\xb0)": "45"}]), "")

    def test_vamas_round_trip_keeps_the_tilt(self):
        f = self.load()
        md = [f.region_metadata(r) for r in f.regions]
        out = os.path.join(self.dir.name, "t.vms")
        exporters.export_vamas(f.regions, out, metadata=md)
        back = load_file(out)
        self.assertEqual([r.extra.get("tilt") for r in back.regions], TILTS)
        self.assertEqual(back.regions[3].conditions["Sample tilt (\xb0)"], "45")


    def test_html_payload_gives_each_spectrum_its_tilt(self):
        import htmlbrowser
        f = self.load()
        payload = htmlbrowser.build_payload([f])
        spectra = payload["samples"][0]["regions"]
        self.assertEqual([s.get("tilt") for s in spectra], TILTS)
        js = open(os.path.join(ROOT, "viewer", "viewer.js"),
                  encoding="utf-8").read()
        self.assertIn("V.tiltTag = function", js)


class TestOtherTiltsAreUntouched(unittest.TestCase):
    def test_one_constant_tilt_is_not_a_z_axis(self):
        from readers.base import Region
        regs = [Region(name="C 1s", index=i, offset=i, energy=[1, 2],
                       counts=[1, 2]) for i in range(3)]
        for r in regs:
            r.extra["tilt"] = 30.0
        self.assertNotEqual(viewdata.resolve_z(regs).mode, "Tilt")


@unittest.skipUnless(os.environ.get("XPS_THERMO_TILT_CORPUS"),
                     "set XPS_THERMO_TILT_CORPUS to a folder of tilt-series "
                     ".VGD/.avg pairs")
class TestRealTiltFiles(unittest.TestCase):
    def test_every_file_reads_the_same_from_both_formats(self):
        root = os.environ["XPS_THERMO_TILT_CORPUS"]
        vgds = glob.glob(os.path.join(root, "**", "*.VGD"), recursive=True)
        self.assertTrue(vgds)
        for v in vgds:
            a = v[:-4] + ".avg"
            if not os.path.exists(a):
                continue
            fv, fa = load_file(v), load_file(a)
            with self.subTest(v):
                self.assertEqual([r.extra.get("tilt") for r in fv.regions],
                                 TILTS)
                self.assertEqual([r.extra.get("tilt") for r in fa.regions],
                                 TILTS)
                for x, y in zip(fv.regions, fa.regions):
                    self.assertEqual(len(x.energy), len(y.energy))
                    for p, q in zip(x.energy, y.energy):
                        self.assertAlmostEqual(p, q, delta=1e-9)
                    self.assertEqual((x.pos_x, x.pos_y), (y.pos_x, y.pos_y))
                    for p, q in zip(x.counts, y.counts):
                        self.assertAlmostEqual(p, q, delta=2e-6 * max(1, abs(p)))


if __name__ == "__main__":
    unittest.main()
