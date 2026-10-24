"""Old-format (converted) Avantage files and the scan count of a SnapMap.

The shared logic is tested through synthetic ``.avg`` dumps (helpers of
``test_thermo_kinds``); the ``.vgd`` property mapping through a hand-made
property set. ``XPS_OLD_VGD=<folder of old .VGD files>`` (the MGCO3.DATA of
Avantage's sample data) checks the real files.

Run:  python -m unittest discover tests
"""

import os
import struct
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import txfn  # noqa: E402
from readers import load_file  # noqa: E402
from readers import thermo_vgd  # noqa: E402
from test_thermo_kinds import Tmp, dump, snapmap_text  # noqa: E402

AXES = [(1188.0, 0.5, 4, "ENERGY", "LINEAR", "E", "eV", "Energy")]


def scan_text(extra=""):
    text = dump("C1s Scan", AXES, [(0, 3, 1)], [((), [], [1, 2, 3, 4])])
    return text.replace("DS_ANPROPID_PASS", extra + "DS_ANPROPID_PASS", 1)


def prop(line):
    return line + "\n"


class TestScanCount(Tmp):
    def test_a_snapmap_takes_the_dumps_periods(self):
        # the dump's PERIODS equals the .vgd's id 0x132 on all 613 pairs
        text = snapmap_text().replace(
            "DS_ANPROPID_PASS",
            prop("DS_ACPROPID_PERIODS : VT_I4 = 1") + "DS_ANPROPID_PASS", 1)
        f = load_file(self.write("map.avg", text))
        self.assertEqual(f.kind, "map")
        self.assertEqual(f.regions[0].extra["n_scans"], 1)

    def test_a_map_with_only_the_unverified_count_stays_unknown(self):
        text = snapmap_text().replace(
            "DS_ANPROPID_PASS",
            prop("DS_ACPROPID_PERIODS : VT_I4 = 2")
            + prop("DS_ACPROPID_PERIODS_UNVERIFIED : VT_BOOL = True")
            + "DS_ANPROPID_PASS", 1)
        f = load_file(self.write("map.avg", text))
        self.assertNotIn("n_scans", f.regions[0].extra)

    def test_an_ordinary_scan_is_unchanged(self):
        f = load_file(self.write("c.avg", scan_text(
            prop("DS_ACPROPID_PERIODS : VT_I4 = 4"))))
        self.assertEqual(f.regions[0].extra["n_scans"], 4)


class TestOldLayout(Tmp):
    START = prop("DS_ACPROPID_START_TIME : VT_DATE = 18/8/1998   14:50:12")

    def test_the_zone_is_utc_only_for_the_modern_layout(self):
        modern = load_file(self.write("m.avg", scan_text(self.START)))
        self.assertEqual(modern.regions[0].extra["tz"], "UTC")
        old = load_file(self.write("o.avg", scan_text(
            self.START + prop("DS_EXT_SOFTWARE_VERSION : VT_BSTR = "
                              "' v4.63 : Build 02770'"))))
        x = old.regions[0].extra
        self.assertEqual(x["t_start"], "1998-08-18 14:50:12")
        self.assertEqual(x["tz"], "")

    def test_a_non_monochromated_source_is_named(self):
        text = scan_text(prop("DS_SOPROPID_STYPE : VT_I4 = 10")).replace(
            "DS_SOPROPID_MONO                            : VT_BOOL = True",
            "DS_SOPROPID_MONO                            : VT_BOOL = False")
        f = load_file(self.write("n.avg", text))
        self.assertEqual(f.regions[0].anode, "Al K\u03b1 (non-mono)")
        self.assertEqual(f.instrument["X-ray source"],
                         "Al K\u03b1, non-monochromated (1486.68 eV)")

    def test_the_monochromated_source_is_as_before(self):
        f = load_file(self.write("m.avg", scan_text()))
        self.assertEqual(f.regions[0].anode, "Al K\u03b1 (mono)")
        self.assertEqual(f.instrument["X-ray source"],
                         "Al K\u03b1, monochromated (1486.68 eV)")

    def test_an_unknown_unmonochromated_source_is_not_named(self):
        text = scan_text().replace(
            "DS_SOPROPID_MONO                            : VT_BOOL = True",
            "DS_SOPROPID_MONO                            : VT_BOOL = False")
        f = load_file(self.write("n.avg", text))
        self.assertEqual(f.regions[0].anode, "")

    def test_provenance_rows_name_files_not_paths(self):
        text = scan_text(
            prop("DS_EXT_SOFTWARE_VERSION : VT_BSTR = ' v4.63 : Build 02770'")
            + prop("DS_EXT_SOURCE_FILE : VT_BSTR = 'F:\\\\MGCO3.DATA\\\\C1s.VGD'")
            + prop("DS_EXT_AUDIT_CREATED : VT_BSTR = 'Created from Eclipse "
                   "File : C:\\\\Program Files\\\\Thermo\\\\Mgco3\\\\C1S.DTS'"))
        f = load_file(self.write("p.avg", text))
        md = f.region_metadata(f.regions[0])
        self.assertEqual(md["Avantage version"], "v4.63 : Build 02770")
        self.assertEqual(md["Original file"], "C1s.VGD")
        self.assertEqual(md["Provenance"],
                         "Converted from Eclipse file C1S.DTS")
        for v in md.values():
            self.assertNotIn("Program Files", str(v))

    def test_a_modern_file_has_no_provenance_rows(self):
        f = load_file(self.write("m.avg", scan_text()))
        md = f.region_metadata(f.regions[0])
        for k in ("Avantage version", "Original file", "Provenance"):
            self.assertNotIn(k, md)

    def test_coefficients_without_the_removed_flag_are_never_applied(self):
        text = scan_text(prop("DS_ANPROPID_TXFN_COEFFS : VT_I4 = 4")
                         + prop("DS_ANPROPID_TXFN_COEFF[0] : VT_R4 = 3.867640")
                         + prop("DS_ANPROPID_TXFN_COEFF[1] : VT_R4 = -0.075012")
                         + prop("DS_ANPROPID_TXFN_COEFF[2] : VT_R4 = 0.003691")
                         + prop("DS_ANPROPID_TXFN_COEFF[3] : VT_R4 = -0.045752"))
        f = load_file(self.write("t.avg", text))
        r = f.regions[0]
        self.assertIsNone(r.extra["txfn"]["removed"])
        self.assertIsNone(txfn.casa_form(r))
        self.assertNotIn("Transmission function", f.region_metadata(r))


def property_set(props):
    """A minimal OLE property set: ``[(id, vt, value)]`` with int32 (3),
    float32 (4) and unicode string (31) values."""
    body = []
    for pid, vt, v in props:
        if vt == 3:
            data = struct.pack("<i", v)
        elif vt == 4:
            data = struct.pack("<f", v)
        else:
            data = struct.pack("<I", len(v) + 1) + (v + "\0").encode("utf-16le")
        data = struct.pack("<I", vt) + data
        body.append((pid, data + b"\0" * (-len(data) % 4)))
    head = 8 + 8 * len(body)
    offsets, pos = [], head
    for _pid, data in body:
        offsets.append(pos)
        pos += len(data)
    sec = struct.pack("<II", pos, len(body))
    sec += b"".join(struct.pack("<II", pid, off)
                    for (pid, _d), off in zip(body, offsets))
    sec += b"".join(d for _p, d in body)
    return b"\0" * 24 + struct.pack("<I", 1) + b"\0" * 16 \
        + struct.pack("<I", 48) + sec


class FakeOle:
    def __init__(self, streams):
        self.streams = streams

    def listdir(self):
        return list(self.streams)

    def read(self, name):
        return self.streams[name]


class DS:
    def __init__(self):
        self.props = {}


class TestVgdProperties(unittest.TestCase):
    def read(self, props, extra=None):
        streams = {"\x05Props": property_set(props)}
        streams.update(extra or {})
        ds = DS()
        thermo_vgd.ThermoVgdFile()._read_props(FakeOle(streams), ds)
        return ds.props

    def test_0x132_is_the_scan_count_and_beats_0x148(self):
        p = self.read([(0x132, 3, 1), (0x148, 3, 2)])
        self.assertEqual(p["DS_ACPROPID_PERIODS"], 1)
        self.assertNotIn("DS_ACPROPID_PERIODS_UNVERIFIED", p)

    def test_without_0x132_the_count_is_flagged_unverified(self):
        p = self.read([(0x148, 3, 2)])
        self.assertEqual(p["DS_ACPROPID_PERIODS"], 2)
        self.assertTrue(p["DS_ACPROPID_PERIODS_UNVERIFIED"])

    def test_old_layout_strings_and_the_audit_trail(self):
        text = ("Created from Eclipse File : C:\\Thermo\\C1S.DTS\0"
                .encode("utf-16le"))
        p = self.read([(0x12, 31, " v4.63 : Build 02770"),
                       (0x10, 31, "F:\\MGCO3.DATA\\C1s.VGD"),
                       (0x74, 31, "Al K-alpha")],
                      {"AuditTrail/00000000": b"\x01\x02" + text})
        self.assertEqual(p["DS_EXT_SOFTWARE_VERSION"].strip(),
                         "v4.63 : Build 02770")
        self.assertEqual(p["DS_EXT_SOURCE_NAME"], "Al K-alpha")
        self.assertEqual(p["DS_EXT_AUDIT_CREATED"],
                         "Created from Eclipse File : C:\\Thermo\\C1S.DTS")


@unittest.skipUnless(os.environ.get("XPS_OLD_VGD"),
                     "set XPS_OLD_VGD=<folder of old .VGD files>")
class TestRealOldFiles(unittest.TestCase):
    def test_every_file_reads_as_a_converted_old_one(self):
        folder = os.environ["XPS_OLD_VGD"]
        names = sorted(n for n in os.listdir(folder)
                       if n.lower().endswith(".vgd"))
        self.assertTrue(names)
        for n in names:
            f = load_file(os.path.join(folder, n))
            r = f.regions[0]
            md = f.region_metadata(r)
            self.assertEqual(r.anode, "Al K\u03b1 (non-mono)", n)
            self.assertIn(r.extra.get("n_scans"), (1, 2, 4), n)
            self.assertEqual(r.extra["tz"], "", n)
            self.assertIn("Eclipse", md["Provenance"], n)
            self.assertNotIn("Transmission function", md, n)
            self.assertEqual(f.warnings, [], n)


if __name__ == "__main__":
    unittest.main()
