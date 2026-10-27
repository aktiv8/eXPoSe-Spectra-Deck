"""Block type: what kind of measurement a spectrum is, as the instrument file
says (blocktype.py). Kratos ESCApe run markers, and Avantage levels that are
iterations, a depth profile or a tilt series.

Set ``XPS_BLOCK_CORPUS=<D:\\Temp\\hxps>`` for the real Kratos sessions and
``XPS_THERMO_ITERATION=<an Avantage experiment folder of iteration series>``
for the Thermo case.

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import blocktype  # noqa: E402
import readers  # noqa: E402
from readers import Region, SpectrumFile, load_file  # noqa: E402
from test_thermo_kinds import Tmp, dump  # noqa: E402

RUN = b"Kernel.Core.Instrument.Acquisition."


class TestText(unittest.TestCase):
    def test_wording(self):
        b = blocktype.make("depth_profile", "the levels have etch times",
                           run=3, size=39)
        self.assertEqual(blocktype.text(b),
                         "Depth profile (the levels have etch times; run 3, "
                         "39 spectra)")
        self.assertEqual(blocktype.text(blocktype.make("etch", "")), "Etch")

    def test_nothing_for_nothing_or_nonsense(self):
        for b in (None, {}, {"kind": "mystery"}, "x"):
            self.assertEqual(blocktype.text(b), "")

    def test_levels(self):
        self.assertEqual(blocktype.levels_block(5, etch=True)["kind"],
                         "depth_profile")
        self.assertEqual(blocktype.levels_block(5, tilt=True)["kind"],
                         "angle_resolved")
        b = blocktype.levels_block(50)
        self.assertEqual(b["kind"], "iteration_series")
        self.assertIn("50 levels", b["reason"])
        # an etch time wins over a tilt, as in the reader
        self.assertEqual(blocktype.levels_block(2, True, True)["kind"],
                         "depth_profile")


class TestKratosRuns(unittest.TestCase):
    def regions(self, *offsets):
        return [Region(f"R{i}", i, o) for i, o in enumerate(offsets)]

    def test_each_spectrum_takes_the_run_before_it(self):
        raw = (b"x" * 100 + RUN + b"UserSequenceData\x00" + b"y" * 100
               + RUN + b"DepthProfileData\x00" + b"z" * 100)
        a, b, c, d = self.regions(10, 150, 300, 310)
        n = blocktype.kratos_runs(raw, [a, b, c, d])
        self.assertEqual(n, 2)
        self.assertNotIn("block", a.extra)           # before the first run
        self.assertEqual((b.extra["block"]["kind"], b.extra["block"]["run"],
                          b.extra["block"]["size"]), ("sequence", 1, 1))
        self.assertEqual((c.extra["block"]["kind"], c.extra["block"]["size"]),
                         ("depth_profile", 2))
        self.assertEqual(d.extra["block"]["run"], 2)
        self.assertIn("recorded by the instrument as a depth-profile run",
                      c.extra["block"]["reason"])

    def test_every_class_has_a_kind(self):
        for name, (kind, _why) in blocktype.RUN_CLASSES.items():
            raw = RUN + name.encode()
            (r,) = self.regions(len(raw) + 5)
            blocktype.kratos_runs(raw + b"....", [r])
            self.assertEqual(r.extra["block"]["kind"], kind, name)
            self.assertIn(kind, blocktype.LABELS)

    def test_an_unknown_class_or_no_marker_labels_nothing(self):
        (r,) = self.regions(500)
        self.assertEqual(blocktype.kratos_runs(
            RUN + b"SomethingNewData" + b"." * 600, [r]), 0)
        self.assertNotIn("block", r.extra)
        self.assertEqual(blocktype.kratos_runs(b"nothing", [r]), 0)

    def test_the_metadata_row(self):
        f = SpectrumFile()
        r = Region("C 1s", 0, 0)
        r.extra["block"] = blocktype.make("sequence", "why", run=1, size=4)
        self.assertEqual(f.region_metadata(r)["Block type"],
                         "Acquisition sequence (why; run 1, 4 spectra)")
        self.assertNotIn("Block type",
                         f.region_metadata(Region("C 1s", 0, 0)))


class TestAvantageLevels(Tmp):
    NE = 4

    def levels(self, name, second_axis, labels_for, n=3):
        axes = [(1188.0, 0.5, self.NE, "ENERGY", "LINEAR", "E", "eV",
                 "Energy")] + second_axis
        blocks = [((k,), labels_for(k),
                   [10.0 * k + c for c in range(self.NE)]) for k in range(n)]
        return load_file(self.write(name, dump(
            name, axes, [(0, self.NE - 1, 1), (0, n - 1, 1)], blocks)))

    def kinds(self, f):
        return {(r.extra.get("block") or {}).get("kind") for r in f.regions}

    def test_levels_with_neither_etch_time_nor_tilt_are_iterations(self):
        f = self.levels(
            "C1s Scan.avg",
            [(0.0, 1.0, 3, "ETCHLEVEL", "LINEAR", "EtchLevel", "",
              "Iteration")],
            lambda k: ["$AXISVALUE= DATAXIS=1 SPACEAXIS=1 LABEL='Iteration' "
                       "POINT=%d VALUE=%f;" % (k, k)])
        self.assertEqual(self.kinds(f), {"iteration_series"})
        self.assertIn("3 levels", f.regions[0].extra["block"]["reason"])
        self.assertIn("Iteration series",
                      f.region_metadata(f.regions[0])["Block type"])

    def test_etch_times_make_a_depth_profile(self):
        f = self.levels(
            "Al2p Snap.avg",
            [(0.0, 5.0, 3, "ETCHTIME", "NON-LINEAR", "EtchTime", "s",
              "Etch Time")],
            lambda k: ["$AXISVALUE= DATAXIS=1 SPACEAXIS=1 LABEL='Etch Time' "
                       "POINT=%d VALUE=%f;" % (k, 5.0 * k)])
        self.assertEqual(self.kinds(f), {"depth_profile"})

    def test_a_tilt_axis_makes_an_angle_resolved_series(self):
        f = self.levels(
            "Si2p Tilt.avg",
            [(0.0, 15.0, 3, "ANGLE", "LINEAR", "T", "\u00b0", "Angle")],
            lambda k: ["$AXISVALUE= DATAXIS=1 SPACEAXIS=1 LABEL='Angle' "
                       "POINT=%d VALUE=%f;" % (k, 15.0 * k)])
        self.assertEqual(self.kinds(f), {"angle_resolved"})

    def test_a_single_scan_has_no_block(self):
        axes = [(1188.0, 0.5, 4, "ENERGY", "LINEAR", "E", "eV", "Energy")]
        f = load_file(self.write("C1s Scan.avg", dump(
            "C1s Scan", axes, [(0, 3, 1)], [((), [], [1, 2, 3, 4])])))
        self.assertEqual(self.kinds(f), {None})
        self.assertNotIn("Block type", f.region_metadata(f.regions[0]))


BLOCK = os.environ.get("XPS_BLOCK_CORPUS", "")
ITER = os.environ.get("XPS_THERMO_ITERATION", "")


@unittest.skipUnless(BLOCK and os.path.isdir(BLOCK),
                     "set XPS_BLOCK_CORPUS to the XPSView projects folder")
class TestKratosCorpus(unittest.TestCase):
    """The run markers against XPSView's own run list for the same files."""

    def load(self, *parts):
        return readers.load_file(os.path.join(BLOCK, *parts))

    def runs(self, d):
        out = {}
        for r in d.regions:
            b = r.extra.get("block")
            if b:
                out[b["run"]] = (b["kind"], b["size"])
        return [out[k] for k in sorted(out)]

    def test_sr_has_six_sequences_and_eight_depth_profiles(self):
        d = self.load("2026-Q3 SR-20263976-26-14",
                      "SR-20263976-26-14.experiment")
        self.assertEqual(self.runs(d),
                         [("sequence", 4)] * 6 + [("depth_profile", 39)] * 8)
        self.assertTrue(all(r.extra.get("block") for r in d.regions))

    def test_ea_runs_match_xpsview_and_the_angle_resolved_run_is_named(self):
        d = self.load("2026-Q3-EA-20264028-26-83",
                      "EA-20264028-26-83 G1B6.experiment")
        runs = self.runs(d)
        self.assertEqual(runs.count(("sequence", 7)), 6)
        self.assertEqual(runs.count(("depth_profile", 112)), 1)
        # XPSView guessed "an exposure test?" for these 35 spectra
        self.assertEqual(runs.count(("angle_resolved", 35)), 1)


@unittest.skipUnless(ITER and os.path.isdir(ITER),
                     "set XPS_THERMO_ITERATION to an Avantage experiment")
class TestThermoCorpus(unittest.TestCase):
    def test_iteration_folders_are_iteration_series(self):
        d = readers.load_file(ITER)
        by = {}
        for r in d.regions:
            b = r.extra.get("block") or {}
            by.setdefault(r.sample, set()).add((b.get("kind"), b.get("size")))
        self.assertEqual(by["CrO3 Quick Scan"], {(None, None)})
        self.assertEqual(by["CrO3 Dual Neutralisation"],
                         {("iteration_series", 50)})
        self.assertEqual(by["CrO3 Neutralisation Via Electrons 100 uA"],
                         {("iteration_series", 10)})


if __name__ == "__main__":
    unittest.main()
