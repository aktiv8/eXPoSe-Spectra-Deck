"""Plain numbers in columns (CSV / ASC / TXT / DAT): the table finder, the
column guesses, the options, the sniffer and the reader. Shapes taken from the
files of Surface Science Spectra submissions (tab or comma separated, header
``BE_Cl2p`` / ``CPS_Cl2p`` or ``Binding Energy (eV),c/s`` or none at all,
intensity before energy, stray lines under the table).

Set ``XPS_COLUMN_CORPUS=<folder>`` to check that every .csv / .asc / .txt file
in a folder loads.

Run:  python -m unittest discover tests
"""

import math
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import annotations  # noqa: E402
import columntext as ct  # noqa: E402
import readers  # noqa: E402
from readers import column_text  # noqa: E402


def peak(n=60, lo=280.0, hi=292.0, centre=285.0):
    """(be descending, counts) of one Gaussian on a flat background."""
    step = (hi - lo) / (n - 1)
    be = [round(hi - i * step, 4) for i in range(n)]
    return be, [round(100 + 900 * math.exp(-((x - centre) / 0.8) ** 2), 3)
                for x in be]


def text(rows, header=None, sep="\t", eol="\r\n", pre=(), post=()):
    lines = list(pre)
    if header:
        lines.append(sep.join(header))
    lines += [sep.join(str(v) for v in r) for r in rows]
    lines += list(post)
    return eol.join(lines) + eol


class TestTable(unittest.TestCase):
    def test_tab_separated_with_a_header(self):
        be, y = peak()
        t = ct.read_table(text(zip(be, y), ["BE_C1s", "CPS_C1s"]))
        self.assertEqual(t.delimiter, "\t")
        self.assertEqual(t.header, ["BE_C1s", "CPS_C1s"])
        self.assertEqual((t.n_rows, t.n_cols), (60, 2))
        self.assertEqual(t.columns[0][:2], be[:2])
        self.assertEqual(t.first_line, 2)
        self.assertEqual(t.trailing, 0)

    def test_comma_semicolon_and_white_space(self):
        be, y = peak()
        for sep, delim in ((",", ","), (";", ";"), (" ", None), ("  ", None)):
            t = ct.read_table(text(zip(be, y), sep=sep))
            self.assertEqual(t.delimiter, delim, repr(sep))
            self.assertEqual(t.n_rows, 60)
            self.assertEqual(t.header, [])

    def test_a_decimal_comma_with_semicolons_or_tabs(self):
        rows = [(f"{285 - i * 0.1:.1f}".replace(".", ","),
                 f"{100 + i}".replace(".", ",") + ",5") for i in range(20)]
        for sep in (";", "\t"):
            t = ct.read_table(text(rows, sep=sep))
            self.assertAlmostEqual(t.columns[0][1], 284.9)
            self.assertAlmostEqual(t.columns[1][0], 100.5)

    def test_no_header_and_intensity_first(self):
        be, y = peak(40)
        t = ct.read_table(text(zip(y, be)))
        self.assertEqual(t.header, [])
        self.assertEqual(ct.find_energy_col(t), (1, True))

    def test_lines_under_the_table_are_set_aside_and_counted(self):
        be, y = peak(30)
        t = ct.read_table(text(zip(be, y), ["BE", "Counts/s"],
                               sep=",", post=["65406.7,", "65298.3,", "7,"]))
        self.assertEqual(t.n_rows, 30)
        self.assertEqual(t.trailing, 3)
        self.assertEqual(t.trailing_at, 1 + 30 + 1)   # header, rows

    def test_text_above_the_header_is_kept_as_a_preamble(self):
        be, y = peak(30)
        t = ct.read_table(text(zip(be, y), ["BE", "CPS"],
                               pre=["# Sample A", "Al Ka, 20 eV pass"]))
        self.assertEqual(t.header, ["BE", "CPS"])
        self.assertEqual(t.preamble, ["# Sample A", "Al Ka, 20 eV pass"])

    def test_blank_lines_inside_do_not_split_the_table(self):
        be, y = peak(20)
        rows = [f"{a}\t{b}" for a, b in zip(be, y)]
        rows.insert(10, "")
        t = ct.read_table("\r\n".join(rows))
        self.assertEqual(t.n_rows, 20)

    def test_the_longest_run_wins_over_a_short_one(self):
        be, y = peak(30)
        short = [f"{i}\t{i * 2}\t{i * 3}" for i in range(9)]
        t = ct.read_table("\n".join(short + ["x"] + [
            f"{a}\t{b}" for a, b in zip(be, y)]))
        self.assertEqual((t.n_rows, t.n_cols), (30, 2))

    def test_too_short_or_not_numbers_is_none(self):
        self.assertIsNone(ct.read_table(""))
        self.assertIsNone(ct.read_table("a,b\n1,2\n3,4\n"))
        self.assertIsNone(ct.read_table("name,value\n" + "\n".join(
            f"k{i},v{i}" for i in range(30))))

    def test_nan_inf_and_underscored_cells_are_not_numbers(self):
        rows = ["1\t2"] * 10
        self.assertIsNone(ct.read_table("\n".join(
            rows[:5] + ["nan\t1"] + rows[5:])))      # two runs of five
        self.assertIsNone(ct._num("1_0"))
        self.assertIsNone(ct._num("inf"))
        self.assertEqual(ct._num(" 1e3 "), 1000.0)

    def test_utf8_bom_and_windows_1252(self):
        self.assertEqual(ct.decode(b"\xef\xbb\xbfBE"), "BE")
        self.assertEqual(ct.decode("Kα".encode("utf-8")), "Kα")
        self.assertEqual(ct.decode(b"\xb5m"), "µm")


class TestHeaders(unittest.TestCase):
    def test_axis_name_and_unit(self):
        h = ct.header_info("BE_Cl2p")
        self.assertEqual((h["axis"], h["energy"], h["name"]), ("BE", True,
                                                                "Cl2p"))
        self.assertEqual(ct.header_info("KE_C 1s")["axis"], "KE")
        self.assertEqual(ct.header_info("KE_C 1s")["name"], "C 1s")
        h = ct.header_info("CPS_Mo 3d")
        self.assertEqual((h["units"], h["name"]), ("counts/s", "Mo 3d"))
        self.assertEqual(ct.header_info("Binding Energy (eV)")["axis"], "BE")
        self.assertEqual(ct.header_info("Binding Energy (eV0")["axis"], "BE")
        self.assertEqual(ct.header_info("c/s")["units"], "counts/s")
        self.assertEqual(ct.header_info("Counts/s")["units"], "counts/s")
        self.assertEqual(ct.header_info("Counts")["units"], "counts")
        self.assertEqual(ct.header_info("Intensity (a.u.)")["units"], "a.u.")
        self.assertIsNone(ct.header_info("Column 3")["units"])
        self.assertIsNone(ct.header_info("beam")["axis"])     # not "be"

    def test_names(self):
        self.assertEqual(ct.region_name("C1s"), "C 1s")
        self.assertEqual(ct.region_name("Gen"), "Survey")
        self.assertEqual(ct.region_name("Survey"), "Survey")
        self.assertEqual(ct.region_name("VB"), "VB")
        self.assertEqual(ct.region_name("", [0, 1300]), "Survey")
        self.assertEqual(ct.region_name("", [280, 290]), "")


class TestColumns(unittest.TestCase):
    def test_the_steady_column_is_the_energy_wherever_it_stands(self):
        be, y = peak(50)
        a = ct.read_table(text(zip(be, y)))
        b = ct.read_table(text(zip(y, be)))
        self.assertEqual(ct.find_energy_col(a), (0, True))
        self.assertEqual(ct.find_energy_col(b), (1, True))

    def test_a_header_breaks_a_tie(self):
        n = 30
        rows = [(i * 0.1, 500 - i * 0.1) for i in range(n)]       # both steady
        t = ct.read_table(text(rows, ["CPS", "BE"]))
        self.assertEqual(ct.find_energy_col(t)[0], 1)

    def test_no_steady_column_is_reported(self):
        rows = [(i % 7, (i * 13) % 11) for i in range(40)]
        t = ct.read_table(text(rows))
        self.assertFalse(ct.find_energy_col(t)[1])

    def test_axis_quality(self):
        self.assertEqual(ct.axis_quality([1, 2, 3, 4]), (1.0, 1.0))
        mono, reg = ct.axis_quality([1, 2, 3, 3.5, 6, 8])
        self.assertEqual(mono, 1.0)
        self.assertLess(reg, ct.REGULAR)
        self.assertEqual(ct.axis_quality([3]), (0.0, 0.0))


class TestOptions(unittest.TestCase):
    def table(self, header, rows):
        return ct.read_table(text(rows, header))

    def test_guesses_from_the_header(self):
        be, y = peak()
        t = self.table(["BE_Cl2p", "CPS_Cl2p"], zip(be, y))
        o = ct.guess_options(t)
        self.assertEqual((o["energy_col"], o["intensity_cols"]), (0, [1]))
        self.assertEqual((o["axis"], o["units"], o["names"]),
                         ("BE", "counts/s", ["Cl2p"]))
        self.assertIsNone(o["photon_energy"])

    def test_kinetic_header_and_unknown_units(self):
        ke = [1180 + i * 0.1 for i in range(30)]
        o = ct.guess_options(self.table(["KE_C 1s", "Signal"],
                                        zip(ke, range(30))))
        self.assertEqual((o["axis"], o["units"]), ("KE", "a.u."))

    def test_several_intensity_columns_are_several_spectra(self):
        be, y = peak(30)
        t = self.table(["BE", "Sample A", "Sample B"],
                       [(a, b, b * 2) for a, b in zip(be, y)])
        o = ct.guess_options(t)
        self.assertEqual(o["intensity_cols"], [1, 2])
        self.assertEqual(o["names"], ["Sample A", "Sample B"])

    def test_sanitise_repairs_anything(self):
        d = ct.sanitise_options(None)
        self.assertEqual(d["energy_col"], 0)
        self.assertEqual(d["intensity_cols"], [1])
        bad = ct.sanitise_options(
            {"energy_col": "x", "intensity_cols": [0, 0, "q", -1, 9],
             "axis": "sideways", "units": 4, "names": [5, "C 1s"],
             "photon_energy": -3, "pass_energy": "nan", "anode": "Cu",
             "sample": ["no"]}, n_cols=3)
        self.assertEqual(bad["energy_col"], 0)
        self.assertEqual(bad["intensity_cols"], [1])   # none valid: first other
        self.assertEqual((bad["axis"], bad["units"]), ("BE", "a.u."))
        self.assertIsNone(bad["photon_energy"])
        self.assertIsNone(bad["pass_energy"])
        self.assertEqual((bad["anode"], bad["sample"]), ("", ""))

    def test_the_energy_column_is_never_an_intensity_column(self):
        o = ct.sanitise_options({"energy_col": 1, "intensity_cols": [1, 0]},
                                n_cols=2)
        self.assertEqual((o["energy_col"], o["intensity_cols"]), (1, [0]))
        o = ct.sanitise_options({"energy_col": 1, "intensity_cols": [1]},
                                n_cols=2)
        self.assertEqual(o["intensity_cols"], [0])           # falls back

    def test_photon_energy_text(self):
        self.assertEqual(ct.parse_photon("Al Kα (1486.6 eV)"), 1486.6)
        self.assertEqual(ct.parse_photon("Mg"), 1253.6)
        self.assertEqual(ct.parse_photon("1486,7"), 1486.7)
        self.assertEqual(ct.parse_photon("  21.2 eV"), 21.2)
        for none in ("", "Unknown", "(not set)", "abc", "0", "-5"):
            self.assertIsNone(ct.parse_photon(none), none)
        self.assertEqual(ct.photon_label(1486.6), "Al Kα (1486.6 eV)")
        self.assertEqual(ct.photon_label(21.2), "21.2 eV")
        self.assertEqual(ct.photon_label(None), "unknown")


class TestSuggestion(unittest.TestCase):
    def test_one_common_strongest_line_in_the_window(self):
        self.assertEqual(ct.suggest_name(279.0, 291.0, 1486.6), "C 1s")
        self.assertEqual(ct.suggest_name(521.0, 542.0, 1486.6), "O 1s")

    def test_ambiguous_wide_or_rare_windows_give_nothing(self):
        self.assertEqual(ct.suggest_name(276.0, 297.0, 1486.6), "")  # C 1s, K 2p
        self.assertEqual(ct.suggest_name(0, 1300), "")           # too wide
        self.assertEqual(ct.suggest_name(259.5, 277.5), "")      # only Fr 4f
        self.assertEqual(ct.suggest_name(1200, 1240), "")


class TestSniff(unittest.TestCase):
    def sniff(self, body, ext=".csv"):
        return column_text.sniff(body.encode("utf-8"), ext)

    def test_a_spectrum_in_any_of_the_extensions(self):
        be, y = peak()
        body = text(zip(be, y), ["BE_C1s", "CPS_C1s"])
        for ext in (".csv", ".asc", ".txt", ".dat", ".tsv"):
            self.assertTrue(self.sniff(body, ext), ext)

    def test_other_extensions_are_not_taken(self):
        be, y = peak()
        self.assertFalse(self.sniff(text(zip(be, y)), ".xlsx"))
        self.assertFalse(self.sniff(text(zip(be, y)), ".vms"))

    def test_tables_that_are_not_spectra_are_left_alone(self):
        rows = [(i % 7, (i * 13) % 11) for i in range(40)]
        self.assertFalse(self.sniff(text(rows)))                 # no axis
        self.assertFalse(self.sniff("a,b\n1,2\n3,4\n"))          # too short
        self.assertFalse(self.sniff("Region,Area\n" + "\n".join(
            f"C 1s,{i}" for i in range(30))))                    # words
        self.assertFalse(column_text.sniff(b"\x00\x01" * 500, ".txt"))

    def test_a_file_cut_mid_line_by_the_read_still_counts(self):
        be, y = peak(900)
        body = text(zip(be, y), ["BE_C1s", "CPS_C1s"]).encode()
        self.assertGreater(len(body), column_text.HEAD)
        self.assertTrue(column_text.sniff(body[:column_text.HEAD], ".asc"))


class Files(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, True)

    def write(self, name, body, enc="utf-8"):
        path = os.path.join(self.dir, name)
        with open(path, "wb") as fh:
            fh.write(body.encode(enc) if isinstance(body, str) else body)
        return path


class TestReader(Files):
    def test_load_file_picks_the_column_reader_last(self):
        be, y = peak()
        p = self.write("a.asc", text(zip(be, y), ["BE_Cl2p", "CPS_Cl2p"]))
        self.assertIs(readers.reader_for(p), column_text.ColumnTextFile)
        self.assertEqual(readers.supported_names()[-1],
                         "Column text (CSV / ASC / TXT)")
        doc = readers.load_file(p)
        (r,) = doc.regions
        self.assertEqual((r.name, r.energy_label, r.count_units),
                         ("Cl 2p", "Binding Energy", "counts/s"))
        self.assertEqual(r.energy, be)
        self.assertEqual(r.counts, y)
        self.assertAlmostEqual(r.step, 12 / 59, 3)
        self.assertEqual(r.source, "a.asc")
        self.assertIsNone(r.photon_energy)
        self.assertEqual(doc.warnings, [])
        self.assertEqual(doc.import_options["names"], ["Cl2p"])

    def test_intensity_first_without_a_header(self):
        be, y = peak()
        doc = readers.load_file(self.write("b.txt", text(zip(y, be))))
        r = doc.regions[0]
        self.assertEqual(r.energy, be)
        self.assertEqual(r.counts, y)
        self.assertEqual(r.count_units, "a.u.")
        self.assertEqual(r.name, "b")            # the file name, not a guess

    def test_ascending_binding_energy_is_stored_high_to_low(self):
        be, y = peak()
        doc = readers.load_file(self.write(
            "c.csv", text(zip(be[::-1], y[::-1]), sep=",")))
        r = doc.regions[0]
        self.assertEqual(r.energy, be)
        self.assertEqual(r.counts, y)
        self.assertIn("opposite order", r.note)

    def test_kinetic_energy_is_converted_when_the_photon_energy_is_given(self):
        hv = 1486.6
        be, y = peak()
        ke = [hv - b for b in be][::-1]
        body = text(zip(ke, y[::-1]), ["KE_C 1s", "CPS_C 1s"])
        p = self.write("d.asc", body)
        doc = readers.load_file(p, options={"photon_energy": hv})
        r = doc.regions[0]
        self.assertEqual(r.energy_label, "Binding Energy")
        for a, b in zip(r.energy, be):
            self.assertAlmostEqual(a, b, 6)
        self.assertEqual(r.counts, y)
        self.assertEqual(r.photon_energy, hv)
        self.assertIn("converted", r.note)
        self.assertEqual(r.name, "C 1s")

    def test_kinetic_energy_without_it_stays_kinetic_with_a_warning(self):
        ke = [1180 + i * 0.1 for i in range(40)]
        doc = readers.load_file(self.write(
            "e.asc", text(zip(ke, range(40)), ["KE_C 1s", "CPS_C 1s"])))
        r = doc.regions[0]
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertEqual(r.energy, ke)
        self.assertTrue(any("photon energy" in w for w in doc.warnings))

    def test_what_the_user_enters_reaches_the_regions(self):
        be, y = peak()
        p = self.write("f.csv", text(zip(be, y), ["Binding Energy (eV)",
                                                  "c/s"], sep=","))
        doc = readers.load_file(p, options={
            "sample": "Probe 1", "photon_energy": 1486.6, "anode": "Al Kα",
            "pass_energy": 20, "names": ["C1s"], "units": "counts"})
        r = doc.regions[0]
        self.assertEqual((r.sample, r.photon_energy, r.anode, r.pass_energy,
                          r.name, r.count_units),
                         ("Probe 1", 1486.6, "Al Kα", 20.0, "C 1s", "counts"))

    def test_options_override_the_guess_column_by_column(self):
        be, y = peak(40)
        p = self.write("g.txt", text(zip(be, y)))
        doc = readers.load_file(p, options={"energy_col": 1,
                                            "intensity_cols": [0]})
        # forcing the wrong columns reads the counts as the energy
        self.assertEqual(sorted(doc.regions[0].counts), sorted(be))

    def test_several_intensity_columns_make_several_regions(self):
        be, y = peak(30)
        body = text([(a, b, b * 2) for a, b in zip(be, y)],
                    ["BE", "Sample A", "Sample B"])
        doc = readers.load_file(self.write("h.csv", body.replace("\t", ",")))
        self.assertEqual([r.name for r in doc.regions],
                         ["Sample A", "Sample B"])
        self.assertEqual(doc.regions[1].counts, [b * 2 for b in y])

    def test_stray_lines_are_ignored_and_reported(self):
        be, y = peak(30)
        p = self.write("i.asc", text(zip(be, y), ["BE", "Counts/s"], sep=",",
                                     post=["65406.7,", "65298.3,"]))
        doc = readers.load_file(p)
        self.assertEqual(len(doc.regions[0].energy), 30)
        self.assertTrue(any("2 line(s) after the table" in w
                            for w in doc.warnings))

    def test_a_preamble_goes_to_the_comments(self):
        be, y = peak(30)
        p = self.write("j.txt", text(zip(be, y), ["BE", "CPS"],
                                     pre=["Sample A, Al Ka"]))
        self.assertEqual(readers.load_file(p).instrument["Comments"],
                         "Sample A, Al Ka")

    def test_a_survey_width_axis_is_a_survey(self):
        be = [1300 - i * 1.0 for i in range(1301)]
        doc = readers.load_file(self.write(
            "k.asc", text(zip(be, [100 + i for i in range(1301)]),
                          ["BE_Gen", "CPS_Gen"])))
        self.assertEqual(doc.regions[0].name, "Survey")
        self.assertTrue(doc.regions[0].is_survey)

    def test_it_is_a_normal_spectrum_file(self):
        be, y = peak()
        doc = readers.load_file(self.write("l.asc", text(zip(be, y),
                                                         ["BE_O1s", "CPS"])))
        self.assertIsNotNone(doc.tree)
        self.assertTrue(doc.regions[0].decodable)
        self.assertEqual(doc.format_name, "Column text (CSV / ASC / TXT)")

    def test_a_non_table_is_not_taken_and_says_what_it_is(self):
        p = self.write("m.csv", "name,value\n" + "\n".join(
            f"k{i},v{i}" for i in range(30)))
        with self.assertRaises(readers.UnsupportedFormat) as cm:
            readers.load_file(p)
        self.assertIn("text file starting", str(cm.exception))

    def test_loading_directly_without_a_table_says_so(self):
        p = self.write("n.csv", "a,b\n1,2\n")
        with self.assertRaises(ValueError):
            column_text.ColumnTextFile().load(p)


class TestStoredOptions(Files):
    def test_the_options_round_trip_through_the_annotations(self):
        be, y = peak()
        p = self.write("o.asc", text(zip(be, y), ["BE_C1s", "CPS_C1s"]))
        doc = readers.load_file(p, options={"sample": "S", "photon_energy":
                                            1486.6, "pass_energy": 20})
        ann = annotations.Annotations()
        self.assertTrue(ann.is_empty())
        ann.set_import("f1", doc.import_options)
        self.assertFalse(ann.is_empty())
        back = annotations.Annotations.from_json(ann.to_json())
        self.assertEqual(back.import_for("f1"), doc.import_options)
        again = readers.load_file(p, options=back.import_for("f1"))
        self.assertEqual(again.regions[0].energy, doc.regions[0].energy)
        self.assertEqual(again.regions[0].sample, "S")
        self.assertEqual(again.regions[0].photon_energy, 1486.6)
        ann.set_import("f1", None)
        self.assertTrue(ann.is_empty())

    def test_bad_stored_options_are_repaired_or_dropped(self):
        a = annotations.Annotations.from_json(
            {"imports": {"f1": {"axis": "no", "photon_energy": "x"},
                         "f2": "junk", "f3": {}}})
        self.assertEqual(sorted(a.imports), ["f1"])
        self.assertEqual(a.import_for("f1")["axis"], "BE")
        self.assertEqual(a.import_for("zzz"), {})


@unittest.skipUnless(os.environ.get("XPS_COLUMN_CORPUS"), "no corpus folder")
class TestCorpus(unittest.TestCase):
    def test_every_column_file_loads(self):
        folder = os.environ["XPS_COLUMN_CORPUS"]
        seen = 0
        for d, _s, files in os.walk(folder):
            for f in files:
                p = os.path.join(d, f)
                if not f.lower().endswith(ct.EXTS):
                    continue
                try:
                    cls = readers.reader_for(p)
                except readers.UnsupportedFormat:
                    continue
                if cls is not column_text.ColumnTextFile:
                    continue
                doc = readers.load_file(p)
                self.assertTrue(doc.regions, p)
                for r in doc.regions:
                    self.assertTrue(r.decodable and len(r.energy) >= 8, p)
                    self.assertEqual(len(r.energy), len(r.counts), p)
                seen += 1
        self.assertGreater(seen, 0)


if __name__ == "__main__":
    unittest.main()
