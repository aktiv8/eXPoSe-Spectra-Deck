"""Kratos Vision2 ``.dset`` reader: synthetic files built here, plus (optional)
every real ``.dset`` / ``.kal`` pair in a folder.

    set XPS_DSET_CORPUS=C:\\path\\to\\folder\\with\\dset\\and\\kal      (any depth)
    python -m unittest tests.test_kratos_dset

The pair test checks each ``.dset`` against the ``.kal`` DumpDataset made of it:
the same regions, energies, counts, photon / pass energy, dwell, date,
metadata and transmission function (the ``.kal`` prints numbers to about six
significant digits, so dwell and transmission compare to that precision).
"""

import glob
import os
import struct
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from readers import load_file, reader_for  # noqa: E402
from readers import kratos_dset  # noqa: E402

CORPUS = os.environ.get("XPS_DSET_CORPUS", "")


# ----------------------------------------------------------------------
# a test-only writer (the app never writes .dset)
# ----------------------------------------------------------------------

def w32(*xs):
    return b"".join(struct.pack(">I", x) for x in xs)


def r_int(i, v):
    return w32(i, v & 0xFFFFFFFF)


def r_dbl(i, x):
    return w32(i) + struct.pack(">d", x)


def r_str(i, s):
    return w32(i, len(s) + 1, *[ord(c) for c in s], 0)


def r_f32(i, xs):
    return w32(i, len(xs)) + b"".join(struct.pack(">f", x) for x in xs)


def r_f64(i, xs):
    return w32(i, len(xs)) + b"".join(struct.pack(">d", x) for x in xs)


def r_box(i, *recs):
    """A container: three header words, its records, its own closing 0, 0."""
    return w32(i, 0, 0, 0) + b"".join(recs) + w32(0, 0)


def block(ordinal, *recs):
    body = w32(ordinal, 0, 0) + b"".join(recs) + w32(0, 0)
    return w32(len(body)) + body


def dset(*blocks, header=True, listed=None):
    """``listed`` is the object count the index states (word 5); the default
    is the number of blocks given."""
    head = bytearray(kratos_dset.MAGIC + b"\x00" * (kratos_dset.FIRST_BLOCK - 8))
    struct.pack_into(">I", head, kratos_dset.INDEX_COUNT,
                     len(blocks) if listed is None else listed)
    return bytes(head) + b"".join(blocks)


def spectrum(ordinal=2, name="Mo 3d", anode=8, extra=()):
    return block(
        ordinal,
        r_int(1, 3), r_int(2, 0),
        r_dbl(3, 100.0), r_dbl(4, 2.0),
        r_str(5, "Kinetic Energy"), r_str(6, "eV"),
        r_dbl(7, 0.5), r_str(9, "counts"),
        r_f32(12, [10, 11, 12, 13.5]),
        r_str(37, name), r_dbl(42, 40.0), r_int(99, 3),
        r_str(151, "19/04/01 15:40:59"),
        r_int(3049, 0), r_int(3080, anode),
        r_str(3113, "Mo"), r_str(3114, "3d"),
        r_dbl(3192, 0.01), r_dbl(3193, 14000.0),
        *extra,
        r_box(5587, r_f64(5615, [100.0, 110.0]), r_f64(5616, [1.0, 0.5])),
    )


class Tmp(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, name, data):
        p = os.path.join(self.dir.name, name)
        with open(p, "wb") as fh:
            fh.write(data)
        return p


class TestDset(Tmp):
    def test_decodes_a_spectrum(self):
        positions = block(1, r_str(37, "Positions"), r_int(38, 1))
        f = load_file(self.write("a.dset", dset(positions, spectrum())))
        self.assertEqual(type(f).__name__, "KratosDsetFile")
        self.assertEqual(len(f.regions), 1)           # Positions has no data
        r = f.regions[0]
        self.assertEqual(r.name, "Mo 3d")
        self.assertEqual(r.counts, [10.0, 11.0, 12.0, 13.5])   # float32 exact
        self.assertAlmostEqual(r.photon_energy, 1486.6)
        self.assertAlmostEqual(r.energy[0], 1486.6 - 100.0)
        self.assertEqual(r.pass_energy, 40.0)
        self.assertEqual(r.dwell, 0.5)
        self.assertEqual(r.date, "2019-04-01 15:40:59")
        self.assertEqual(r.tf_ke, [100.0, 110.0])         # nested container
        self.assertEqual(r.tf_values, [1.0, 0.5])
        self.assertEqual(r.extra["n_scans"], 3)
        self.assertEqual(r.conditions["Anode voltage (kV)"], "14")

    def test_registry_sniffs_content_not_extension(self):
        p = self.write("renamed.bin", dset(spectrum()))
        self.assertEqual(reader_for(p).format_name, "Kratos Vision (.dset)")

    def test_index_only_file_says_so(self):
        p = self.write("index.dset", dset())
        with self.assertRaises(ValueError) as cm:
            load_file(p)
        self.assertIn("only its index", str(cm.exception))

    def test_index_only_helper(self):
        p = self.write("index.dset", dset())
        self.assertTrue(kratos_dset.is_index_only(p))
        self.assertFalse(kratos_dset.is_index_only(
            self.write("full.dset", dset(spectrum()))))
        self.assertFalse(kratos_dset.is_index_only(
            self.write("other.bin", b"" * 3000)))
        self.assertFalse(kratos_dset.is_index_only(p + ".missing"))

    def test_truncated_file_is_refused(self):
        data = dset(spectrum())
        with self.assertRaises(ValueError) as cm:
            load_file(self.write("t.dset", data[:-40]))
        self.assertIn("chain", str(cm.exception))

    # what follows the last indexed block in a transmission-function file:
    # not block-shaped (``size, 1024, 25792, 0``), claims more than is there
    TRAILER = w32(0x43980, 1024, 25792, 0, 1, 3, 2, 0) + b"\x00" * 200

    def test_trailing_section_is_skipped_when_every_indexed_block_was_read(self):
        a, b = spectrum(2, "Hyb_pe05"), spectrum(3, "Hyb_pe10")
        plain = load_file(self.write("p.dset", dset(a, b)))
        f = load_file(self.write("x.dset", dset(a, b) + self.TRAILER))
        self.assertEqual([r.counts for r in f.regions],
                         [r.counts for r in plain.regions])
        self.assertEqual(len(f.regions), 2)
        self.assertEqual(len(f.warnings), 1)
        self.assertIn("were not read", f.warnings[0])
        self.assertIn("all of those objects were loaded", f.warnings[0])
        self.assertFalse(plain.warnings)

    def test_trailing_section_is_refused_when_indexed_blocks_are_missing(self):
        # the index lists 3 objects, 2 chain: not an extra section but a loss
        data = dset(spectrum(2), spectrum(3), listed=3) + self.TRAILER
        with self.assertRaises(ValueError) as cm:
            load_file(self.write("m.dset", data))
        msg = str(cm.exception)
        self.assertIn("do not chain", msg)
        self.assertIn("2 of the 3 objects", msg)

    def test_truncated_last_block_is_still_refused_with_its_count(self):
        # a normal block header with too little behind it is truncation, even
        # though the earlier blocks are all there
        data = dset(spectrum(2), spectrum(3))
        with self.assertRaises(ValueError) as cm:
            load_file(self.write("t2.dset", data[:-40]))
        self.assertIn("1 of the 2 objects", str(cm.exception))

    def test_ion_gun_block_with_the_high_numbered_pah_records_decodes(self):
        gun = block(3, r_str(37, "Etch"), r_int(38, 0),
                    r_dbl(3230, 4.0), r_int(3251, 0), r_dbl(3278, 0.0),
                    r_int(0x20200CE0, 1), r_int(0x20200CE7, 0),
                    r_dbl(0x20500CDB, 0.0), r_dbl(0x20500CDC, 0.0),
                    r_dbl(0x20500CE3, 0.0))
        f = load_file(self.write("g.dset", dset(spectrum(), gun)))
        self.assertEqual(len(f.regions), 1)
        self.assertFalse(f.warnings)                 # decoded, no unknown enum

    def test_unknown_enum_number_is_left_out_not_guessed(self):
        f = load_file(self.write("u.dset", dset(spectrum(anode=99))))
        r = f.regions[0]
        self.assertIsNone(r.photon_energy)                # no invented anode
        self.assertEqual(r.energy_label, "Kinetic Energy")
        text = " ".join(f.warnings)
        self.assertIn("Xray Reference Energy = 99", text)

    def test_unknown_fields_of_every_shape_are_skipped(self):
        extra = (r_int(4001, 7), r_dbl(4002, 2.5), r_str(4003, "hello"),
                 r_f32(4004, [1, 2, 3]))
        f = load_file(self.write("k.dset", dset(spectrum(extra=extra))))
        r = f.regions[0]
        self.assertEqual(r.counts, [10.0, 11.0, 12.0, 13.5])
        self.assertEqual(r.tf_ke, [100.0, 110.0])         # still in step after
        self.assertFalse(f.warnings)

    def test_a_block_that_cannot_be_decoded_costs_only_that_block(self):
        bad = block(3, r_str(37, "Bad"), w32(4005, 0xFFFFFFF0, 1, 2))
        f = load_file(self.write("b.dset", dset(spectrum(), bad)))
        self.assertEqual(len(f.regions), 1)
        self.assertTrue(any("not decoded" in w for w in f.warnings))


# ----------------------------------------------------------------------
# real pairs
# ----------------------------------------------------------------------

def _close(a, b, tol):
    if a is None or b is None:
        return a == b
    return abs(a - b) <= tol * max(abs(b), 1e-30)


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "XPS_DSET_CORPUS not set")
class TestRealPairs(unittest.TestCase):
    ATTRS = ("name", "sample", "photon_energy", "pass_energy", "step",
             "lens_mode", "date", "anode", "aperture", "energy_label",
             "count_units", "technique")

    def test_each_dset_matches_its_kal(self):
        pairs = [(d, d[:-5] + ".kal") for d in
                 glob.glob(os.path.join(CORPUS, "**", "*.dset"),
                           recursive=True)
                 if os.path.exists(d[:-5] + ".kal")]
        self.assertTrue(pairs, "no .dset / .kal pairs found")
        checked = 0
        for dpath, kpath in pairs:
            try:
                a = load_file(dpath)
            except ValueError as exc:
                # an index-only file, or one whose chain stops short of the
                # objects its index lists (07092017_NPL_Trans): refused
                self.assertTrue("only its index" in str(exc)
                                or "objects the index lists" in str(exc), dpath)
                continue
            b = load_file(kpath)
            self.assertEqual(len(a.regions), len(b.regions), dpath)
            for ra, rb in zip(a.regions, b.regions):
                tag = f"{os.path.basename(dpath)}:{rb.name}"
                for at in self.ATTRS:
                    self.assertEqual(getattr(ra, at), getattr(rb, at),
                                     f"{tag} {at}")
                self.assertEqual(ra.energy, rb.energy, tag)
                self.assertEqual(ra.counts, rb.counts, tag)
                self.assertTrue(_close(ra.dwell, rb.dwell, 1e-5), tag)
                self.assertEqual(len(ra.tf_ke), len(rb.tf_ke), tag)
                for x, y in zip(ra.tf_ke + ra.tf_values,
                                rb.tf_ke + rb.tf_values):
                    self.assertTrue(_close(x, y, 1e-4), f"{tag} transmission")
                self.assertEqual(ra.conditions, rb.conditions, tag)
                self.assertEqual(
                    {k: v for k, v in ra.extra.items() if k != "fields"},
                    {k: v for k, v in rb.extra.items() if k != "fields"}, tag)
                ma = a.region_metadata(ra)
                mb = b.region_metadata(rb)
                for k in ("File format", "Source file"):
                    ma.pop(k, None)
                    mb.pop(k, None)
                self.assertEqual(ma, mb, tag)
                checked += 1
        print(f"\n{len(pairs)} pairs, {checked} regions identical")


if __name__ == "__main__":
    unittest.main()
