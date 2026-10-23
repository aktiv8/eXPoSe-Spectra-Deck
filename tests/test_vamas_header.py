"""VAMAS block header fields written from what the source recorded: the
acquisition date (not the export time), work function and beam width."""
import datetime
import os
import tempfile
import unittest

import exporters
from readers.base import Region


def region(**kw):
    e = [290.0 - 0.5 * i for i in range(8)]
    return Region(name="C 1s", index=1, offset=0, energy=e,
                  counts=[10.0 + i for i in range(8)], decodable=True,
                  sample="S1", photon_energy=1486.68, pass_energy=50.0,
                  dwell=0.05, step=0.5, **kw)


class TestBlockHeader(unittest.TestCase):
    def export(self, r, md=None):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.vms")
            exporters.export_vamas([r], p, include_transmission=False,
                                   metadata=[md] if md else None)
            with open(p, encoding="latin-1") as fh:
                lines = fh.read().split("\n")
        # block starts after the experiment header: name, sample, date
        i = lines.index("C 1s")
        return lines, i

    def test_date_is_acquisition_not_export(self):
        r = region(date="2024-06-25 14:18:05")
        lines, i = self.export(r)
        self.assertEqual(lines[i + 2:i + 8],
                         ["2024", "6", "25", "14", "18", "5"])

    def test_run_start_wins_over_date(self):
        r = region(date="2020-01-01 00:00:00")
        r.extra["t_start"] = "2024-06-25 14:18:05"
        lines, i = self.export(r)
        self.assertEqual(lines[i + 2], "2024")

    def test_date_falls_back_to_metadata_then_export_time(self):
        lines, i = self.export(region(), {"Date acquired": "2023-02-03 04:05:06"})
        self.assertEqual(lines[i + 2:i + 8],
                         ["2023", "2", "3", "4", "5", "6"])
        lines, i = self.export(region())
        self.assertEqual(int(lines[i + 2]), datetime.datetime.now().year)

    def test_work_function_and_beam_width_from_metadata(self):
        md = {"Work function (eV)": "4.2", "X-ray spot (µm)": "200 × 100"}
        lines, i = self.export(region(), md)
        j = lines.index("XPS", i)
        # technique, source label, hv, power, beam x, beam y ...
        self.assertEqual(lines[j + 4:j + 6], ["200", "100"])
        k = lines.index("FAT", j)
        self.assertEqual(lines[k + 3], "4.2")        # work function

    def test_unknown_values_are_not_invented(self):
        lines, i = self.export(region())
        j = lines.index("XPS", i)
        self.assertEqual(lines[j + 4:j + 6], ["1E+37", "1E+37"])
        k = lines.index("FAT", j)
        self.assertEqual(lines[k + 3], "1E+37")


COEFFS = [3.82653, 0.749286, -0.612674, 0.057846]


def tx_region(removed=True):
    r = region()
    r.extra["txfn"] = {"removed": removed, "coeffs": COEFFS}
    return r


class TestCasaForm(unittest.TestCase):
    def test_transmission_formula(self):
        import txfn
        # Casa's own column for C 1s at pass 50: 3528.62 at KE 1186.6 eV
        self.assertAlmostEqual(txfn.transmission(1186.6, 50, COEFFS) / 3528.62,
                               1.0, places=5)
        self.assertAlmostEqual(txfn.transmission(136.6, 150, COEFFS) / 9357.38,
                               1.0, places=5)

    def test_counts_and_column(self):
        import txfn
        r = tx_region()
        counts, t = txfn.casa_form(r)
        ke = r.kinetic_energy
        ref = txfn.transmission(txfn.REF_KE, 50.0, COEFFS)
        for c0, c, tt, e in zip(r.counts, counts, t, ke):
            self.assertAlmostEqual(c, c0 * tt / ref)
            self.assertAlmostEqual(c / tt * ref, c0)   # T divides back out

    def test_only_when_the_source_removed_it(self):
        import txfn
        self.assertIsNone(txfn.casa_form(tx_region(removed=False)))
        self.assertIsNone(txfn.casa_form(region()))      # no coefficients

    def write(self, r, **kw):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.vms")
            exporters.export_vamas([r], p, include_transmission=False, **kw)
            with open(p, encoding="latin-1") as fh:
                return fh.read().split("\n")

    def test_export_writes_second_variable(self):
        lines = self.write(tx_region(), casa_transmission=True)
        self.assertIn("Transmission", lines)
        plain = self.write(tx_region())
        self.assertNotIn("Transmission", plain)
        # unflagged data is untouched even with the option on
        self.assertNotIn("Transmission",
                         self.write(region(), casa_transmission=True))

    def test_mfp_parameters(self):
        lines = self.write(region(), escape_depth={"exponent": -0.6, "type": 5})
        i = lines.index("MFP Exponent")
        self.assertEqual(lines[i - 1], "2")
        self.assertEqual(lines[i:i + 6],
                         ["MFP Exponent", "d", "-0.6",
                          "ESCAPE DEPTH TYPE", "d", "5"])
        self.assertNotIn("MFP Exponent", self.write(region()))

    def test_zero_exponent_is_written_as_casaxps_does(self):
        # a file saved by CasaXPS with a zero exponent holds "0", never "-0"
        for zero in (0.0, -0.0):
            lines = self.write(region(),
                               escape_depth={"exponent": zero, "type": 1})
            i = lines.index("MFP Exponent")
            self.assertEqual(lines[i:i + 6], ["MFP Exponent", "d", "0",
                                              "ESCAPE DEPTH TYPE", "d", "1"])

    def test_every_escape_type_code_is_written(self):
        for kind in range(7):
            lines = self.write(region(),
                               escape_depth={"exponent": 2.0, "type": kind})
            i = lines.index("ESCAPE DEPTH TYPE")
            self.assertEqual(lines[i + 2], str(kind))
            self.assertEqual(lines[lines.index("MFP Exponent") + 2], "2")


@unittest.skipUnless(os.environ.get("XPS_AVG_CASA"),
                     "set XPS_AVG_CASA=<folder with the Naresh files>")
class TestAgainstCasa(unittest.TestCase):
    """The export of the Avantage files equals CasaXPS's own import."""

    def test_counts_and_transmission_match(self):
        import readers
        folder = os.environ["XPS_AVG_CASA"]
        with open(os.path.join(folder, "casa from avg.vms"),
                  encoding="latin-1") as fh:
            casa = [l.rstrip("\r") for l in fh.read().split("\n")]
        f = readers.load_file(os.path.join(folder, "C1s Scan.avg"))
        r = f.regions[0]
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.vms")
            exporters.export_vamas([r], p, casa_transmission=True,
                                   metadata=[f.region_metadata(r)])
            with open(p, encoding="latin-1") as fh:
                ours = [l.rstrip("\r") for l in fh.read().split("\n")]

        # 502 values (251 points x 2 variables), then min/max of each
        ci = casa.index("502") + 5
        oi = ours.index("502") + 5
        got = [float(x) for x in ours[oi:oi + 2 * 251]]
        want = [float(x) for x in casa[ci:ci + 2 * 251]]
        for g, w in zip(got, want):
            self.assertLess(abs(g - w) / abs(w), 1e-5)


if __name__ == "__main__":
    unittest.main()
