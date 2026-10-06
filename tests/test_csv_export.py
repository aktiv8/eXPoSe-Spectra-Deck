"""The spectra CSV: UTF-8 with a byte-order mark, optional provenance header.

Run:  python -m unittest discover tests
"""

import csv
import datetime
import math
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import exporters  # noqa: E402
from readers import Region  # noqa: E402


def reg(name="C 1s", sample="Probe µ-1 α", n=21):
    e = [292.0 - i * 12.0 / (n - 1) for i in range(n)]
    c = [100 + 900 * math.exp(-((x - 286) / 1.0) ** 2) for x in e]
    return Region(name=name, index=0, offset=0, energy=e, counts=c,
                  decodable=True, sample=sample, photon_energy=1486.6,
                  pass_energy=20.0, dwell=0.1, step=0.2, source="a.vms",
                  count_units="counts/s")


class TestCsv(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.path = os.path.join(self.dir, "x.csv")

    def raw(self):
        with open(self.path, "rb") as fh:
            return fh.read()

    def rows(self):
        with open(self.path, encoding="utf-8-sig", newline="") as fh:
            return list(csv.reader(fh))

    def test_a_greek_sample_name_is_written_and_read_back(self):
        # used to raise UnicodeEncodeError under the Windows default (cp1252)
        exporters.export_csv([reg()], self.path)
        self.assertTrue(self.raw().startswith(b"\xef\xbb\xbf"))
        head = self.rows()[0]
        self.assertTrue(head[0].startswith("Probe µ-1 α C 1s"), head)

    def test_the_default_file_has_no_comment_lines(self):
        exporters.export_csv([reg()], self.path)
        self.assertFalse(any(r and r[0].startswith("#") for r in self.rows()))
        self.assertEqual(len(self.rows()), 1 + 21)

    def test_a_comment_comes_first_and_the_table_is_unchanged(self):
        exporters.export_csv([reg()], self.path, comment="one\ntwo")
        text = self.raw().decode("utf-8-sig")
        self.assertTrue(text.startswith("# one\r\n# two\r\n"))
        rows = self.rows()
        self.assertEqual(rows[0], ["# one"])
        self.assertEqual(rows[1], ["# two"])
        plain = os.path.join(self.dir, "p.csv")
        exporters.export_csv([reg()], plain)
        with open(plain, encoding="utf-8-sig", newline="") as fh:
            self.assertEqual(rows[2:], list(csv.reader(fh)))


class TestProvenance(unittest.TestCase):
    when = datetime.datetime(2026, 10, 6, 21, 5, tzinfo=datetime.timezone.utc)

    def test_it_names_the_app_the_time_the_files_and_the_corrections(self):
        text = exporters.provenance_comment("App", "1.2", ["a.vms", "b.kal"],
                                            True, self.when)
        lines = text.splitlines()
        self.assertEqual(lines[0], "App 1.2, written 2026-10-06 21:05 UTC")
        self.assertEqual(lines[1], "Source file(s): a.vms, b.kal")
        self.assertIn("are applied", lines[2])

    def test_uncorrected_and_unknown_are_said_plainly(self):
        text = exporters.provenance_comment("App", "1", [], False, self.when)
        self.assertIn("(unknown)", text)
        self.assertIn("not applied", text)


if __name__ == "__main__":
    unittest.main()
