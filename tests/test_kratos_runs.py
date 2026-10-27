"""Kratos runs told apart (readers/kratos_experiment.py): each depth profile has
its own levels and etch times, a sequence has none, a sample measured in several
runs becomes one sample per run, and the annotations of the plain sample carry
over to the new ones.

Synthetic Kratos ESCApe bytes are built here (length-prefixed strings, the run, graph
and location markers); ``XPS_KRATOS_RUNS=<D:\\Temp\\hxps>`` checks the real SR
and EA files against XPSView's own sample, level and etch time of every spectrum.

Run:  python -m unittest discover tests
"""

import collections
import glob
import gzip
import json
import os
import struct
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import annotations  # noqa: E402
import blocktype  # noqa: E402
import readers  # noqa: E402
from readers import Region  # noqa: E402
from readers.kratos_experiment import EscapeParser  # noqa: E402

RUN = b"Kernel.Core.Instrument.Acquisition."


def ls(text):
    """A length-prefixed string, as Kratos ESCApe writes them."""
    return bytes([len(text)]) + text.encode("ascii")


def location(sample, loc=""):
    """An analysis-location record: GUID, location name, the two "Analysis"
    strings, a stage position and the sample's name."""
    return (b"InstrumentAnalysisLocation\x00\x00" + ls("0123-4567-89ab-cdef")
            + ls(loc or sample) + ls("Analysis") + ls("Analysis")
            + struct.pack("<3d", 0.03, 0.01, -0.0003) + ls(sample) + b"\x00" * 8)


def beam(text):
    return b"Mb6EtchSettings\x00\x00" + text.encode("ascii") + b"\x00" * 4


def graph(*durations):
    """A duration graph: (etch number, duration) pairs, 16 bytes apart."""
    body = b"".join(struct.pack("<2d", k + 1, d)
                    for k, d in enumerate(durations))
    return (b"DataTypes.NonUniformGraphData\x00" + ls("Etch number")
            + ls("Duration") + b"\x00" * 20 + body)


class Build:
    """A synthetic file: runs of (kind, sample, region names, beam, durations)."""

    def __init__(self, runs, legacy_marker=False):
        raw, self.regions, self.starts = b"HEAD" + b"\x00" * 40, [], []
        for kind, sample, names, settings, durs in runs:
            raw += location(sample)
            if settings:
                raw += beam(settings)
            self.starts.append(len(raw))
            raw += RUN + kind.encode() + b"\x00" * 30
            for name in names:
                self.regions.append(Region(name, len(self.regions),
                                           len(raw), sample=sample))
                raw += b"SPECTRUM" + b"\x00" * 40
            if durs:
                raw += graph(*durs)
            raw += b"\x00" * 40
        if legacy_marker:
            raw += b"Depth Profile"
        self.raw = raw

    def parser(self):
        p = EscapeParser()
        p.raw = self.raw
        p.corruption = {"corrupted": False}
        p._scan_strings()
        p.regions = self.regions
        p.path = "t.experiment"
        return p

    def run(self):
        p = self.parser()
        p._runs = blocktype.kratos_runs(p.raw, p.regions)
        p._assign_run_samples()
        p._parse_depth_profile()
        p._qualify_samples()
        return p


DP = "DepthProfileData"
SEQ = "UserSequenceData"
TWO_LEVELS = ["C 1s", "O 1s", "C 1s", "O 1s"]


def two_profiles():
    return Build([
        (SEQ, "Sample A", ["wide", "C 1s"], "", None),
        (DP, "Sample A", TWO_LEVELS, "10 keV Ar1000+", (30.0,)),
        (DP, "Sample B", TWO_LEVELS, "5 keV Ar+", (60.0,)),
        (SEQ, "Sample B", ["wide", "C 1s"], "", None),
    ])


class TestRuns(unittest.TestCase):
    def test_each_profile_has_its_own_levels_and_etch_times(self):
        p = two_profiles().run()
        dp = [r for r in p.regions if r.extra["block"]["kind"] == "depth_profile"]
        self.assertEqual([r.etch_level for r in dp], [0, 0, 1, 1] * 2)
        self.assertEqual([r.etch_time for r in dp[:4]], [0.0, 0.0, 30.0, 30.0])
        self.assertEqual([r.etch_time for r in dp[4:]], [0.0, 0.0, 60.0, 60.0])
        self.assertEqual(dp[0].extra["acq_metadata"]["Etch beam"],
                         "10 keV Ar1000+")
        self.assertEqual(dp[4].extra["acq_metadata"]["Etch beam"],
                         "5 keV Ar+")

    def test_a_sequence_has_no_levels_or_etch_time(self):
        p = two_profiles().run()
        seq = [r for r in p.regions if r.extra["block"]["kind"] == "sequence"]
        self.assertEqual(len(seq), 4)
        self.assertTrue(all(r.etch_level is None and r.etch_time is None
                            for r in seq))

    def test_the_summary_describes_the_profiles(self):
        dp = two_profiles().run().depth_profile
        self.assertTrue(dp["is_profile"])
        self.assertEqual((dp["n_levels"], dp["regions_per_level"],
                          dp["n_runs"]), (2, 2, 2))
        self.assertEqual(dp["total_etch_time"], 60.0)

    def test_a_file_without_a_profile_has_none(self):
        p = Build([(SEQ, "Sample A", ["wide"], "", None),
                   (SEQ, "Sample A", ["wide"], "", None)]).run()
        self.assertFalse(p.depth_profile["is_profile"])

    def test_a_file_without_run_markers_keeps_the_old_reading(self):
        b = Build([], legacy_marker=True)
        b.raw += location("A")
        for lv in range(3):
            for name in ("C 1s", "O 1s"):
                b.regions.append(Region(name, len(b.regions), len(b.raw),
                                        sample="A"))
                b.raw += b"SPECTRUM" + b"\x00" * 40
        p = b.run()
        self.assertEqual([r.etch_level for r in p.regions], [0, 0, 1, 1, 2, 2])
        self.assertEqual(p.sample_aliases, {})


class TestSamples(unittest.TestCase):
    def test_the_sample_comes_from_the_runs_own_location_record(self):
        b = two_profiles()
        for r in b.regions:           # the old reading: the last sample wins
            r.sample = "Sample B"
        p = b.run()
        plain = {r.extra.get("sample_plain") or r.sample for r in p.regions[:2]}
        self.assertEqual(plain, {"Sample A"})
        self.assertEqual({r.extra["sample_plain"] for r in p.regions[2:6]},
                         {"Sample A"})
        self.assertEqual({r.extra["sample_plain"] for r in p.regions[6:10]},
                         {"Sample B"})

    def test_a_sample_with_several_runs_becomes_one_sample_per_run(self):
        p = two_profiles().run()
        self.assertEqual(
            sorted({r.sample for r in p.regions}),
            ["Sample A - Depth profile 1", "Sample A - Sequence 1", "Sample B - Depth profile 2",
             "Sample B - Sequence 2"])
        self.assertEqual(p.sample_aliases["Sample A - Depth profile 1"], "Sample A")

    def test_a_sample_with_one_run_keeps_its_name(self):
        p = Build([(SEQ, "Sample A", ["wide"], "", None),
                   (DP, "Sample B", TWO_LEVELS, "", (30.0,))]).run()
        self.assertEqual({r.sample for r in p.regions}, {"Sample A", "Sample B"})
        self.assertEqual(p.sample_aliases, {})

    def test_numbering_is_by_kind_through_the_file(self):
        p = Build([(DP, "Sample A", TWO_LEVELS, "", (30.0,)),
                   (SEQ, "Sample A", ["wide"], "", None),
                   (DP, "Sample A", TWO_LEVELS, "", (30.0,)),
                   (SEQ, "Sample A", ["wide"], "", None)]).run()
        self.assertEqual(sorted(set(p.sample_aliases)),
                         ["Sample A - Depth profile 1", "Sample A - Depth profile 2",
                          "Sample A - Sequence 1", "Sample A - Sequence 2"])

    def test_the_tree_lays_out_a_profile_sample_with_levels_only(self):
        p = two_profiles().run()
        p._build_tree()
        kinds = {n.label: [c.type_name for c in n.children]
                 for n in p.tree.children}
        self.assertEqual(set(kinds["Sample A - Depth profile 1"]), {"regionfolder"})
        self.assertEqual(set(kinds["Sample A - Sequence 1"]), {"EscaSpectrum"})
        folder = next(c for n in p.tree.children
                      if n.label == "Sample A - Depth profile 1"
                      for c in n.children if c.label == "C 1s")
        self.assertEqual([c.label for c in folder.children],
                         ["Level 0", "Level 1"])


class TestAnnotations(unittest.TestCase):
    def test_the_new_samples_start_with_the_plain_ones_edits(self):
        a = annotations.Annotations()
        a.sample_notes[annotations.sample_key("f1", "A")] = "a note"
        a.set_shift(annotations.sample_key("f1", "A"), 0.4)
        a.shifts[annotations.region_key("f1", "A", "C 1s")] = 0.7
        a.sample_names[annotations.sample_key("f1", "A")] = "Renamed"
        a.add_marker("f1", "A", "O 1s", 531.0, "O")
        n = a.adopt_aliases("f1", {"Sample A - Depth profile 1": "A"})
        self.assertGreater(n, 3)
        new = annotations.sample_key("f1", "Sample A - Depth profile 1")
        self.assertEqual(a.sample_notes[new], "a note")
        self.assertEqual(a.shifts[new], 0.4)
        self.assertEqual(a.shifts[annotations.region_key(
            "f1", "Sample A - Depth profile 1", "C 1s")], 0.7)
        self.assertEqual(len(a.markers_for("f1", "Sample A - Depth profile 1",
                                           "O 1s")), 1)
        self.assertNotIn(new, a.sample_names)        # the run keeps its name
        self.assertIn(annotations.sample_key("f1", "A"), a.sample_notes)

    def test_nothing_is_overwritten_and_other_samples_are_left_alone(self):
        a = annotations.Annotations()
        new = annotations.sample_key("f1", "Sample A - Depth profile 1")
        a.shifts[new] = 1.5
        a.shifts[annotations.sample_key("f1", "A")] = 0.4
        a.shifts[annotations.sample_key("f1", "AB")] = 9.0   # a longer name
        a.shifts[annotations.sample_key("f2", "A")] = 8.0    # another file
        a.adopt_aliases("f1", {"Sample A - Depth profile 1": "A"})
        self.assertEqual(a.shifts[new], 1.5)
        self.assertNotIn(annotations.sample_key("f1", "AB - Depth profile 1"),
                         a.shifts)
        self.assertEqual(a.adopt_aliases("f1", {}), 0)
        self.assertEqual(a.adopt_aliases("f1", None), 0)


try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_APP = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_APP = None, None, False


@unittest.skipUnless(HAVE_APP, "Tk / matplotlib not available")
class TestInApp(unittest.TestCase):
    def test_registering_a_file_adopts_the_plain_samples_edits(self):
        import matplotlib
        saved = matplotlib.rcParams.copy()
        try:
            root = tk.Tk()
        except tk.TclError:
            self.skipTest("no display")
        root.withdraw()
        self.addCleanup(lambda: (root.destroy(),
                                 matplotlib.rcParams.update(saved)))
        p = two_profiles().run()
        p._build_summary()
        ws = ee.Workspace(root)
        ws.ann.sample_notes[annotations.sample_key("f1", "Sample A")] = "kept"
        ws._register_parser(p, "t.experiment", file_id="f1")
        self.assertEqual(ws.ann.sample_notes[annotations.sample_key(
            "f1", "Sample A - Depth profile 1")], "kept")


CORPUS = os.environ.get("XPS_KRATOS_RUNS", "")


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_KRATOS_RUNS to the XPSView projects folder")
class TestAgainstXpsview(unittest.TestCase):
    """Every spectrum's sample, level and etch time against XPSView's own."""

    def compare(self, xpsv, files):
        import numpy as np
        xp = json.loads(gzip.open(xpsv).read())["experiment"]
        sm = {s["id"]: s["name"] for s in xp["samples"]}
        pool = collections.defaultdict(list)
        for f in files:
            for r in readers.load_file(f).regions:
                pool[len(r.counts)].append(r)
        used, bad, seen = set(), [], 0
        for sp in xp["spectra"]:
            for li, row in enumerate(sp["y"]):
                c = np.array([np.nan if v is None else v for v in row], float)
                o = next((o for o in pool[len(c)] if id(o) not in used
                          and np.allclose(np.asarray(o.counts, float), c,
                                          equal_nan=True)), None)
                self.assertIsNotNone(o, sp["name"])
                used.add(id(o))
                seen += 1
                plain = o.extra.get("sample_plain") or o.sample
                if (not plain.startswith("Sample ")
                        and plain != sm[sp["sample_id"]]):
                    bad.append(("sample", plain, sm[sp["sample_id"]]))
                xe = sp["levels"][li].get("etch_time")
                if len(sp["levels"]) > 1 and xe is not None and (
                        o.etch_time is None or abs(o.etch_time - xe) > 1e-6):
                    bad.append(("etch", o.sample, sp["name"], o.etch_time, xe))
        return seen, bad

    def test_sr_agrees_for_every_spectrum(self):
        base = os.path.join(CORPUS, "2026-Q3 SR-20263976-26-14")
        seen, bad = self.compare(
            os.path.join(base, "20263976-26-14.xpsv"),
            [os.path.join(base, "SR-20263976-26-14.experiment")])
        self.assertEqual((seen, bad), (336, []))

    def test_ea_agrees_except_one_profile_with_an_extra_spectrum(self):
        base = os.path.join(CORPUS, "2026-Q3-EA-20264028-26-83")
        seen, bad = self.compare(
            os.path.join(base, "XPSView", "20264028-26-83.xpsv"),
            [f for f in glob.glob(os.path.join(base, "*.experiment"))
             if "-18" not in f])
        self.assertEqual(seen, 1549)
        # the one-etch profile of Cu whose run holds a third spectrum
        self.assertLessEqual(len(bad), 1, bad)


if __name__ == "__main__":
    unittest.main()
