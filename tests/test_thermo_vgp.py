"""Avantage project (``.vgp``) decoder: the layout of the objects, their links
and, on the real projects, the fit rebuilt from the decoded numbers against the
counts of the same project's ``.VGD``.

The synthetic projects are written with ``tests/olewriter.py`` from the layout
documented in ``readers/thermo_vgp.py``. The real ones are the 129 projects of
``C:\\AvantageSystem\\Knowledge`` (``XPS_VGP_CORPUS=<folder>`` to use another).

Run:  python -m unittest tests.test_thermo_vgp
"""

import os
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import olewriter  # noqa: E402
import readers  # noqa: E402
from readers import thermo_vgp as V  # noqa: E402
from readers.base import Region, SpectrumFile  # noqa: E402
from readers.thermo_experiment import looks_like_experiment  # noqa: E402

CORPUS = os.environ.get("XPS_VGP_CORPUS",
                        os.path.join("C:\\", "AvantageSystem", "Knowledge"))

GUID_F = bytes.fromhex("141a0bb94fa5d211852a00104b053b4c")
cs = olewriter.cstring


def class_tag(name, schema=3):
    return b"\xff\xff" + struct.pack("<HH", schema, len(name)) + name.encode()


def fitted_stream(peaks, count=None):
    """``peaks``: dicts with name, ke, fwhm, mix, area, height, tail, ver."""
    out = struct.pack("<ii", 1, len(peaks) if count is None else count)
    out += class_tag("CFittedPeakItem") + GUID_F
    for k, p in enumerate(peaks):
        if k:
            out += b"\x01\x80" + GUID_F
        ver = p.get("ver", 1)
        off = 40 if ver == 1 else 56
        rec = struct.pack("<iii", 2, ver, 3).ljust(off, b"\x00")
        name = p["name"]
        raw = (name + "\x00").encode("utf-16le") if name else b""
        rec += struct.pack("<I", len(raw)) + raw
        tail = p.get("tail", (0.0, 1.0, 0.0))
        rec += struct.pack("<8d", p["area"], p["ke"], p["height"], p["fwhm"],
                           p["mix"], *tail)
        rec += b"\x00" * 48                   # the parameter blocks, unread
        out += rec
    return out


def background_stream(bgs):
    """``bgs``: dicts with start, end, avg, width, code."""
    out = struct.pack("<ii", 1, len(bgs)) + class_tag("CBackgroundItem")
    for k, b in enumerate(bgs):
        if k:
            out += b"\x01\x80"
        out += struct.pack("<ddiidi", b["start"], b["end"], 1, b["avg"],
                           b["width"], b["code"])
        if b["code"] == 1:
            out += struct.pack("<idi", 50, 0.01, 1)
    return out + b"\x00\x00\x00\x00"


def data_space_stream(paths):
    out = struct.pack("<ii", 1, len(paths)) + class_tag("CDataSpaceItem", 1)
    for k, p in enumerate(paths):
        out += cs(p) + cs("{%08X-0000-0000-0000-000000000000}" % k) + b"\x01\x80"
    return out


def xrefs(*tables):
    out = b""
    for t in tables:
        out += struct.pack("<ii", 0, len(t))
        for a, b in t:
            out += struct.pack("<ii", a, b)
    return out


GUID_P = bytes.fromhex("5ae90d01e568d211851a006097c81436")


def peaks_stream(defs):
    """``defs``: dicts with ``lib``, ``name``, ``rsf`` and optional ``groups``
    (1 or 2): the layout of a CPeakItem record read from real projects (the
    sensitivity factor 37 bytes after the name, 53 for two groups, then -1.0)."""
    out = struct.pack("<ii", 1, len(defs)) + class_tag("CPeakItem") + GUID_P
    for k, d in enumerate(defs):
        if k:
            out += b"\x01\x80" + GUID_P
        lib = (d["lib"] + "\x00").encode("utf-16le") if d["lib"] else b""
        name = (d["name"] + "\x00").encode("utf-16le")
        groups = d.get("groups", 1)
        at = 37 if groups == 1 else 53
        tail = struct.pack("<i", groups).ljust(at, b"\x00")
        tail += struct.pack("<dd", d["rsf"], -1.0).ljust(102 if groups == 1
                                                          else 150, b"\x00")[:]
        out += (struct.pack("<ii", 6, len(lib)) + lib + struct.pack("<i", 0)
                + struct.pack("<i", len(name)) + name + tail)
    return out


PEAKS = [
    dict(name="C1s C-C", ke=1202.2, fwhm=1.3, mix=0.23, area=5566.0,
         height=94.0),
    dict(name="", ke=1200.9, fwhm=1.3, mix=0.23, area=680.0, height=11.5,
         ver=2),                                        # empty name, long prefix
    dict(name="Pd3d5 metal", ke=1151.3, fwhm=0.87, mix=0.26, area=3398.0,
         height=7719.0, tail=(1e-17, 0.367, 0.0428)),
]


def make_project(peaks=PEAKS, bgs=None, paths=None, xr_peak_ds=None,
                 count=None, defs=None):
    paths = paths or [".\\Demo.DATA\\C1s Scan.VGD", ".\\Demo.DATA\\Pd.VGD",
                      ".\\Demo.DATA\\Unfitted.VGD"]
    bgs = bgs or [dict(start=1195.0, end=1206.0, avg=10, width=0.05, code=3),
                  dict(start=1140.0, end=1157.0, avg=5, width=0.1, code=1)]
    xr_peak_ds = xr_peak_ds or [(0, 0), (1, 0), (2, 1)]
    return olewriter.build({
        "Contents": cs("Processing Document") + b"\x01\x00\x00\x00"
                    + cs("ALTHERMO1"),
        "Data Spaces": {"Contents": data_space_stream(paths),
                        "XRefs": xrefs([(0, 0)], [(0, 0)], [(0, 0), (1, 1)])},
        "Fitted Peaks": {"Contents": fitted_stream(peaks, count),
                         "XRefs": xrefs(xr_peak_ds)},
        "Backgrounds": {"Contents": background_stream(bgs),
                        "XRefs": xrefs([(0, 0), (0, 1), (1, 2)],
                                       [(0, 0), (1, 1)])},
        "Peaks": ({"Contents": peaks_stream(defs),
                   "XRefs": xrefs([(i, i) for i in range(len(defs))])}
                  if defs else
                  {"Contents": b"\x01\x00\x00\x00\x00\x00\x00\x00",
                   "XRefs": b""}),
    })


class TestLayout(unittest.TestCase):
    def setUp(self):
        self.prj = V.parse(make_project())

    def test_library_and_data_spaces(self):
        self.assertEqual(self.prj.library, "ALTHERMO1")
        self.assertEqual([d.path for d in self.prj.data_spaces],
                         [".\\Demo.DATA\\C1s Scan.VGD", ".\\Demo.DATA\\Pd.VGD",
                          ".\\Demo.DATA\\Unfitted.VGD"])
        self.assertTrue(self.prj.data_spaces[0].guid.startswith("{"))

    def test_fitted_peak_fields(self):
        p = self.prj.peaks[0]
        self.assertEqual((p.name, p.ke, p.fwhm, p.mix, p.area, p.height),
                         ("C1s C-C", 1202.2, 1.3, 0.23, 5566.0, 94.0))
        self.assertEqual(p.tail, (0.0, 1.0, 0.0))
        self.assertTrue(p.symmetric)

    def test_an_empty_name_and_the_long_prefix_are_read(self):
        p = self.prj.peaks[1]
        self.assertEqual(p.name, "")
        self.assertEqual((p.ke, p.area, p.height), (1200.9, 680.0, 11.5))

    def test_tail_terms_are_kept_and_make_a_peak_asymmetric(self):
        p = self.prj.peaks[2]
        self.assertEqual(p.tail, (1e-17, 0.367, 0.0428))
        self.assertFalse(p.symmetric)

    def test_links(self):
        self.assertEqual([p.data_space for p in self.prj.peaks], [0, 0, 1])
        self.assertEqual([p.background for p in self.prj.peaks], [0, 0, 1])
        self.assertEqual(self.prj.fitted_data_spaces(), [0, 1])
        self.assertEqual(len(self.prj.peaks_of(0)), 2)
        self.assertEqual(self.prj.background_of(1).index, 1)
        self.assertIsNone(self.prj.background_of(2))
        self.assertEqual(self.prj.notes, [])

    def test_background_forms(self):
        a, b = self.prj.backgrounds
        self.assertEqual((a.kind, a.start_ke, a.end_ke, a.avg_points,
                          a.width_ev), ("Smart", 1195.0, 1206.0, 10, 0.05))
        self.assertEqual(len(a.groups), 1)
        self.assertEqual(b.kind, "Shirley")
        self.assertEqual(len(b.groups), 2)         # the Shirley record's extra group
        lin = V.parse(make_project(bgs=[
            dict(start=1.0, end=2.0, avg=5, width=0.1, code=0),
            dict(start=1.0, end=2.0, avg=5, width=0.1, code=0)]))
        self.assertEqual(lin.backgrounds[0].kind, "Linear")

    def test_an_unknown_background_code_is_said(self):
        prj = V.parse(make_project(bgs=[
            dict(start=1.0, end=2.0, avg=5, width=0.1, code=7),
            dict(start=1.0, end=2.0, avg=5, width=0.1, code=7)]))
        self.assertEqual(prj.backgrounds[0].kind, "type 7")
        self.assertTrue(any("not known" in n for n in prj.notes))

    def test_a_peak_with_no_spectrum_is_said(self):
        prj = V.parse(make_project(xr_peak_ds=[(0, 0), (1, 0)]))
        self.assertIsNone(prj.peaks[2].data_space)
        self.assertTrue(any("no spectrum" in n for n in prj.notes))

    def test_a_project_without_fits_is_read(self):
        prj = V.parse(make_project(peaks=[], xr_peak_ds=[(0, 0)]))
        self.assertEqual(prj.peaks, [])
        self.assertEqual(prj.fitted_data_spaces(), [])

    def test_path_resolves_below_the_project_folder(self):
        ref = self.prj.data_spaces[0]
        self.assertEqual(ref.resolve(os.path.join("x", "proj")),
                         os.path.join("x", "proj", "Demo.DATA", "C1s Scan.VGD"))


class TestRefused(unittest.TestCase):
    def test_not_an_ole_file(self):
        with self.assertRaises(V.VgpError):
            V.parse(b"hello world, this is not a project")

    def test_an_ole_file_that_is_not_a_project(self):
        with self.assertRaises(V.VgpError):
            V.parse(olewriter.build({"Contents": b"\x00" * 16}))

    def test_a_peak_count_that_disagrees_with_the_records(self):
        with self.assertRaises(V.VgpError) as cm:
            V.parse(make_project(count=5))
        self.assertIn("5 peaks", str(cm.exception))

    def test_unread_bytes_after_the_backgrounds(self):
        raw = olewriter.build({
            "Data Spaces": {"Contents": data_space_stream(["a.VGD"]),
                            "XRefs": b""},
            "Fitted Peaks": {"Contents": struct.pack("<ii", 1, 0),
                             "XRefs": b""},
            "Backgrounds": {"Contents": background_stream(
                [dict(start=1.0, end=2.0, avg=5, width=0.1, code=0)])
                + b"\x01" * 16, "XRefs": b""}})
        with self.assertRaises(V.VgpError):
            V.parse(raw)

    def test_a_data_space_count_that_disagrees_with_the_paths(self):
        bad = struct.pack("<ii", 1, 3) + class_tag("CDataSpaceItem", 1) \
            + cs("a.VGD")
        raw = olewriter.build({
            "Data Spaces": {"Contents": bad, "XRefs": b""},
            "Fitted Peaks": {"Contents": struct.pack("<ii", 1, 0),
                             "XRefs": b""}})
        with self.assertRaises(V.VgpError):
            V.parse(raw)


# -- the shape ---------------------------------------------------------------------
class TestAvantageShape(unittest.TestCase):
    def test_parse_and_flags(self):
        import lineshapes as ls
        p = ls.parse_shape("AV(23,0,1,0,0.05)")
        self.assertEqual((p["kind"], p["mix"], p["tail_terms"], p["step"]),
                         ("AV", 23.0, (0.0, 1.0, 0.0), 0.05))
        self.assertEqual(ls.parse_shape("AVP(23,0,1,0,0.05)")["kind"], "AVP")
        for s in ("AV(23,0,1,0,0.05)", "AVP(23,0,1,0,0.05)"):
            self.assertFalse(ls.is_exact(s))                  # a reconstruction
            self.assertEqual(ls.check_shape(s), [])
        self.assertTrue(ls.check_shape("AV(120,0,1,0,0.05)"))
        self.assertTrue(ls.check_shape("AV()"))
        self.assertTrue(ls.check_shape("AV(23,0,1.5,0,0.05)"))   # TM above 1
        self.assertTrue(ls.check_shape("AV(23,0,0.5,0,-1)"))     # negative step

    def test_a_tail_that_cannot_be_drawn_is_said(self):
        import lineshapes as ls
        self.assertEqual(ls.shape_notes("AV(23,0,1,0,0.05)"), [])
        self.assertEqual(ls.shape_notes("AV(26,1e-17,0.367,0.0428,0.05)"), [])
        self.assertTrue(ls.shape_notes("AV(26,1e-17,0.367,0.0428,0)"))  # no step
        self.assertEqual(ls.shape_notes("GL(30)"), [])

    def test_the_sum_and_the_product_are_the_manuals_functions(self):
        import numpy as np
        import lineshapes as ls
        ke = np.linspace(1190.0, 1214.0, 2401)
        s = ls.component_curve(ke, "AV(23,0,1,0,0.01)", 1202.0, 1.3, 100.0)
        b = ls.component_curve(ke, "SGL(23)", 1202.0, 1.3, 100.0)
        # the sum is SGL (which is the program's table, so a part in a
        # thousand apart, with the tail area measured over another range)
        self.assertTrue(np.allclose(s, b, rtol=3e-3, atol=3e-3 * float(b.max())))
        # the product, 1 / ([1 + M t^2] exp{(1 - M) ln2 t^2}), t = 2 dx / FWHM
        raw = ls._raw_values(np.array([1202.0, 1202.65, 1203.9]),
                             ls.parse_shape("AVP(40,0,1,0,0.01)"), 1202.0, 1.3)
        for x, v in zip((0.0, 0.65, 1.9), raw):
            t = 2 * x / 1.3
            want = 1.0 / ((1 + 0.4 * t * t) * np.exp(0.6 * np.log(2) * t * t))
            self.assertAlmostEqual(float(v), want, places=9)

    def test_the_tail_is_the_manuals_on_the_high_binding_energy_side_only(self):
        import math
        import numpy as np
        import lineshapes as ls
        # Gaussian (mix 0), FWHM 1, centre 1000: 2 eV away the core is 1.5e-5
        # of the peak, so the value is the tail T = TM CT + (1 - TM) exp(-ET dx)
        # with dx in channels of 0.1 eV
        x = np.array([998.0, 1002.0])                  # low KE (high BE), high KE
        for terms, want in (((0.0, 0.0, 0.05), math.exp(-0.05 * 20)),
                            ((0.5, 0.4, 0.05), 0.4 * 0.5 + 0.6 * math.exp(-1.0)),
                            ((0.0, 1.0, 0.05), 0.0)):
            shape = "AV(0,%g,%g,%g,0.1)" % terms
            v = ls._raw_values(x, ls.parse_shape(shape), 1000.0, 1.0)
            self.assertAlmostEqual(float(v[0]), want, places=3)
            self.assertLess(float(v[1]), 1e-4)         # no tail on the other side

    def test_the_height_is_the_one_asked_for_with_or_without_a_tail(self):
        import numpy as np
        import lineshapes as ls
        ke = np.linspace(1100.0, 1300.0, 40001)
        for shape in ("AV(0,0,1,0,0.05)", "AV(23,0,1,0,0.05)",
                      "AVP(60,0,1,0,0.05)", "AV(26,0,0.367,0.0428,0.05)",
                      "AVP(40,0.1,0.5,0.05,0.1)"):
            area = ls.area_for_height(shape, 1.3, 1000.0)
            curve = ls.component_curve(ke, shape, 1200.0, 1.3, area)
            self.assertAlmostEqual(float(curve.max()) / 1000.0, 1.0, delta=0.01,
                                   msg=shape)

    def test_window_integral_is_the_area_inside_the_window(self):
        import lineshapes as ls
        whole = ls.window_integral("AVP(30,0,1,0,0.05)", 0.0, 1.0, -60.0, 60.0)
        half = ls.window_integral("AVP(30,0,1,0,0.05)", 0.0, 1.0, -60.0, 0.0)
        self.assertAlmostEqual(half / whole, 0.5, places=2)   # symmetric


# -- the project as a session --------------------------------------------------------
HV, DWELL = 1486.68, 0.025
SPECS = {"C1s Scan": ("C 1s", 1190.0, 1210.0), "Pd": ("Pd 3d", 1130.0, 1160.0),
         "Unfitted": ("O 1s", 940.0, 960.0)}
BGS = [dict(start=1195.0, end=1206.0, avg=10, width=0.05, code=0),
       dict(start=1140.0, end=1157.0, avg=5, width=0.1, code=3)]


def _shape(x, fwhm, mix):
    import numpy as np
    t = x / fwhm
    return (1 - mix) * np.exp(-4 * np.log(2) * t * t) + mix / (1 + 4 * t * t)


def _peak_shape(p, mode, step=0.1):
    """The shape string of a demo peak dict (``p``) in ``mode`` (AV / AVP)."""
    return V.shape_of(V.FittedPeak(
        index=0, name=p["name"], ke=p["ke"], fwhm=p["fwhm"], mix=p["mix"],
        area=p["area"], height=p["height"],
        tail=p.get("tail", (0.0, 1.0, 0.0))), mode, step)


def consistent(mode, scale=1.0):
    """The demo peaks with the areas Avantage would store for them: the drawn
    shape's own area over each peak's background window (``scale`` fakes a
    project that stores areas in another unit)."""
    import lineshapes as ls
    out = []
    for p in PEAKS:
        win = (1195.0, 1206.0) if p["ke"] > 1190 else (1140.0, 1157.0)
        q = dict(p)
        q["area"] = scale * p["height"] / DWELL * ls.window_integral(
            _peak_shape(p, mode), p["ke"], p["fwhm"], *win)
        out.append(q)
    return out


class StubVgd(SpectrumFile):
    """Stands in for ``ThermoVgdFile``: one synthetic spectrum per file name,
    a flat 1000 counts/s plus the fitted peaks of the demo project drawn in
    ``MODE`` (the Avantage sum ``AV`` or product ``AVP``, tail included).
    ``TWO`` gives the Pd file a second trace (a "Compare" overlay): ``True``
    an identical one, ``"second"`` a first trace that is only the background,
    so that only the second one is described by the fit. ``DROP`` leaves the
    named peak out of the data (a project whose fit is incomplete)."""
    format_name = "stub"
    MODE = "AV"
    TWO = False
    DROP = None

    def load(self, path):
        import numpy as np
        import lineshapes as ls
        self.path = path
        base = os.path.splitext(os.path.basename(path))[0]
        name, lo, hi = SPECS[base]
        ke = np.arange(lo, hi, 0.1)
        y = np.full(len(ke), 1000.0)
        for p in PEAKS:
            if lo < p["ke"] < hi and p["name"] != self.DROP:
                sp = ls.parse_shape(_peak_shape(p, self.MODE))
                y += p["height"] / DWELL * ls._raw_values(ke, sp, p["ke"],
                                                          p["fwhm"])
        self.kind = "spectrum"
        self.regions = []
        for i in range(2 if (self.TWO and base == "Pd") else 1):
            flat = self.TWO == "second" and i == 0
            self.regions.append(Region(
                name=name, index=i, offset=i, decodable=True,
                energy=list(HV - ke),
                counts=list((np.full(len(ke), 1000.0) if flat else y) * DWELL),
                photon_energy=HV, dwell=DWELL, extra={"n_scans": 1}))
        return self._finish()


class TestProjectFile(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="vgp_")
        self.path = os.path.join(self.dir, "Demo.vgp")
        self._write(consistent("AV"))
        os.makedirs(os.path.join(self.dir, "Demo.DATA"))
        for n in ("C1s Scan", "Pd", "Unfitted"):
            with open(os.path.join(self.dir, "Demo.DATA", n + ".VGD"), "wb"):
                pass
        self._real = V.ThermoVgdFile
        V.ThermoVgdFile = StubVgd
        StubVgd.MODE, StubVgd.TWO = "AV", False

    def _write(self, peaks):
        with open(self.path, "wb") as fh:
            fh.write(make_project(peaks=peaks, bgs=BGS))

    def tearDown(self):
        V.ThermoVgdFile = self._real
        StubVgd.MODE, StubVgd.TWO, StubVgd.DROP = "AV", False, None
        import shutil
        shutil.rmtree(self.dir, ignore_errors=True)

    def _doc(self):
        return V.ThermoProjectFile().load(self.path)

    def test_spectra_come_in_the_projects_order_under_its_name(self):
        doc = self._doc()
        self.assertEqual([r.name for r in doc.regions],
                         ["C 1s", "Pd 3d", "O 1s"])
        self.assertEqual({r.sample for r in doc.regions}, {"Demo"})
        self.assertEqual(doc.experiment["name"], "Demo")
        self.assertEqual(doc.errors, [])

    def test_a_spectrum_the_project_fits_has_its_fit_and_the_others_none(self):
        c1s, pd, o1s = self._doc().regions
        self.assertEqual([c.name for c in c1s.fit.components],
                         ["C1s C-C", "Peak 2"])        # the empty name is numbered
        self.assertEqual(len(pd.fit.components), 1)
        self.assertIsNone(o1s.fit)
        self.assertEqual(c1s.extra["fit_source"], "Avantage")
        reg = c1s.fit.regions[0]
        self.assertEqual((reg.background, reg.start_ke, reg.end_ke, reg.avg),
                         ("Linear", 1195.0, 1206.0, 10))
        self.assertEqual(reg.rsf, 0.0)                  # nothing is guessed
        self.assertTrue(all(c.rsf == 0.0 for c in c1s.fit.components))

    def test_positions_and_areas_are_in_the_fit_codes_units(self):
        c1s = self._doc().regions[0]
        c = c1s.fit.components[0]
        import lineshapes as ls
        self.assertEqual((c.pos_ke, c.fwhm), (1202.2, 1.3))
        self.assertEqual(c.shape, "AV(23,0,1,0,0.1)")        # step = the data's
        self.assertAlmostEqual(
            c.area, ls.area_for_height(c.shape, 1.3, 94.0 / DWELL), places=6)

    def test_the_rebuilt_fit_lands_on_the_spectrum(self):
        import casafit
        c1s = self._doc().regions[0]
        cv = casafit.curves(c1s.fit, c1s.energy, c1s.counts, HV, DWELL, 1)[0]
        self.assertTrue(cv.approximate)                 # never claimed exact
        self.assertLess(cv.residual_rms, 0.01)
        self.assertEqual(cv.notes, [])                  # symmetric, linear

    def test_a_tail_is_drawn_and_the_smart_background_is_reproduced(self):
        import casafit
        import numpy as np
        pd = self._doc().regions[1]
        self.assertTrue(pd.fit.components[0].shape.startswith("AV(26,"))
        cv = casafit.curves(pd.fit, pd.energy, pd.counts, HV, DWELL, 1)[0]
        self.assertTrue(cv.approximate)
        self.assertEqual(pd.fit.regions[0].background, "Smart")
        self.assertTrue(cv.background_known)
        self.assertFalse(any("Smart" in n for n in cv.notes))  # reproduced
        self.assertFalse(any("tail" in n for n in cv.notes))   # it is drawn
        # the tail is on the low-KE (high binding energy) side only
        v = np.array(cv.components[0][1])
        ke = HV - np.array(pd.energy)
        lo = v[np.argmin(np.abs(ke - (1151.3 - 2.0)))]
        hi = v[np.argmin(np.abs(ke - (1151.3 + 2.0)))]
        self.assertGreater(lo, 10 * hi)

    def test_the_product_is_read_from_the_stored_areas(self):
        self._write(consistent("AVP"))
        StubVgd.MODE = "AVP"
        c1s = self._doc().regions[0]
        self.assertTrue(c1s.fit.components[0].shape.startswith("AVP("))
        self.assertIn("product", c1s.fit.components[0].line)
        self._write(consistent("AV"))
        StubVgd.MODE = "AV"
        self.assertTrue(self._doc().regions[0].fit.components[0]
                        .shape.startswith("AV("))

    def test_areas_in_an_unknown_unit_are_drawn_as_the_sum_and_say_so(self):
        import casafit
        self._write(consistent("AVP", scale=0.03))
        c1s = self._doc().regions[0]
        self.assertTrue(c1s.fit.components[0].shape.startswith("AV("))
        cv = casafit.curves(c1s.fit, c1s.energy, c1s.counts, HV, DWELL, 1)[0]
        self.assertTrue(any("neither" in n for n in cv.notes))

    def test_of_several_identical_traces_the_first_carries_the_fit(self):
        StubVgd.TWO = True
        regs = self._doc().regions
        self.assertEqual([r.name for r in regs],
                         ["C 1s", "Pd 3d", "Pd 3d", "O 1s"])
        self.assertIsNotNone(regs[1].fit)
        self.assertIsNone(regs[2].fit)
        self.assertEqual(regs[1].extra["avantage_level"], 0)

    def test_the_fit_goes_to_the_trace_it_describes(self):
        import casafit
        StubVgd.TWO = "second"            # the first trace is only background
        regs = self._doc().regions
        self.assertIsNone(regs[1].fit)
        self.assertIsNotNone(regs[2].fit)
        self.assertEqual(regs[2].extra["avantage_level"], 1)
        r = regs[2]
        cv = casafit.curves(r.fit, r.energy, r.counts, HV, DWELL, 1)[0]
        self.assertLess(cv.residual_rms, 0.05)
        self.assertTrue(any("trace 2" in n for n in cv.notes))

    def test_a_fit_that_leaves_much_unexplained_says_so(self):
        import casafit
        c1s = self._doc().regions[0]
        cv = casafit.curves(c1s.fit, c1s.energy, c1s.counts, HV, DWELL, 1)[0]
        self.assertFalse(any("unexplained" in n for n in cv.notes))
        StubVgd.DROP = "C1s C-C"          # the data lack the project's main peak
        try:
            c1s = self._doc().regions[0]
        finally:
            StubVgd.DROP = None
        cv = casafit.curves(c1s.fit, c1s.energy, c1s.counts, HV, DWELL, 1)[0]
        self.assertGreater(cv.residual_rms, V.POOR_FIT)
        note = next(n for n in cv.notes if "unexplained" in n)
        self.assertIn("2 peaks", note)

    def test_members_name_the_project_and_its_spectra(self):
        doc = self._doc()
        rels = [rel for _p, rel in doc.members]
        self.assertEqual(rels[0], "Demo.vgp")
        self.assertEqual(sorted(rels[1:]),
                         sorted(os.path.join("Demo.DATA", n + ".VGD")
                                for n in ("C1s Scan", "Pd", "Unfitted")))

    def test_the_metadata_names_the_source_of_the_fit(self):
        doc = self._doc()
        md = doc.region_metadata(doc.regions[0])
        self.assertIn("Avantage fit", md)
        self.assertNotIn("CasaXPS fit", md)

    def test_a_missing_spectrum_is_reported_and_the_rest_loads(self):
        os.remove(os.path.join(self.dir, "Demo.DATA", "Pd.VGD"))
        doc = self._doc()
        self.assertEqual([r.name for r in doc.regions], ["C 1s", "O 1s"])
        self.assertEqual(len(doc.errors), 1)
        self.assertIn("Pd.VGD", doc.errors[0][0])

    def test_the_registry_knows_a_project_file_and_its_folder(self):
        self.assertIs(readers.reader_for(self.path), V.ThermoProjectFile)
        self.assertIs(readers.reader_for(self.dir), V.ThermoProjectFile)
        self.assertFalse(looks_like_experiment(self.dir))
        doc = readers.load_file(self.dir)               # a workbook's folder
        self.assertIsInstance(doc, V.ThermoProjectFile)
        self.assertEqual(len(doc.regions), 3)

    def test_sniff_by_extension_or_by_the_stream_names(self):
        with open(self.path, "rb") as fh:
            data = fh.read(4096)
        self.assertTrue(V.sniff(data, ".vgp"))
        self.assertTrue(V.sniff(data, ".dat"))          # renamed: by content
        self.assertFalse(V.sniff(b"not ole", ".vgp"))
        self.assertFalse(V.sniff(olewriter.build({"Contents": b"x"})[:4096],
                                 ".dat"))

    def test_a_file_that_is_not_a_project_says_so(self):
        bad = os.path.join(self.dir, "Bad.vgp")
        with open(bad, "wb") as fh:
            fh.write(olewriter.build({"Contents": b"\x00" * 8}))
        with self.assertRaises(ValueError) as cm:
            V.ThermoProjectFile().load(bad)
        self.assertIn("Avantage project", str(cm.exception))

    def test_a_folder_of_two_projects_is_not_one_project(self):
        with open(self.path, "rb") as fh:
            data = fh.read()
        with open(os.path.join(self.dir, "Other.vgp"), "wb") as fh:
            fh.write(data)
        self.assertIsNone(V.project_in_folder(self.dir))
        self.assertIs(readers.reader_for(self.dir), readers.ThermoExperiment)


def _corpus():
    if not os.path.isdir(CORPUS):
        return []
    return sorted(os.path.join(CORPUS, f) for f in os.listdir(CORPUS)
                  if f.lower().endswith(".vgp"))


@unittest.skipUnless(_corpus(), "no .vgp projects (XPS_VGP_CORPUS)")
class TestRealProjects(unittest.TestCase):
    """Every project decodes, its files are there, every peak sits inside the
    energy range of the spectrum it is linked to, and the fit rebuilt from the
    decoded numbers lands on that spectrum's counts at the right scale."""

    @classmethod
    def setUpClass(cls):
        import numpy as np
        import lineshapes
        cls.np, cls.ls, cls.readers = np, lineshapes, readers
        cls.projects = [(p, V.load(p)) for p in _corpus()]

    def test_every_project_decodes_and_its_spectra_exist(self):
        n_peaks = 0
        for path, prj in self.projects:
            for ref in prj.data_spaces:
                self.assertTrue(os.path.isfile(ref.resolve(CORPUS)),
                                f"{os.path.basename(path)}: {ref.path}")
            for p in prj.peaks:
                self.assertIsNotNone(p.data_space, (path, p.name))
                self.assertIsNotNone(p.background, (path, p.name))
            self.assertEqual(prj.notes, [], path)
            n_peaks += len(prj.peaks)
        self.assertGreater(n_peaks, 0)

    def _spectrum(self, prj, di, cache):
        if di not in cache:
            f = self.readers.load_file(prj.data_spaces[di].resolve(CORPUS))
            cache[di] = f.regions[0]
        return cache[di]

    def test_peaks_lie_inside_their_spectrum(self):
        for path, prj in self.projects:
            cache = {}
            for p in prj.peaks:
                r = self._spectrum(prj, p.data_space, cache)
                lo = r.photon_energy - max(r.energy)
                hi = r.photon_energy - min(r.energy)
                self.assertTrue(lo - 1.5 < p.ke < hi + 1.5,
                                f"{os.path.basename(path)}: {p.name} {p.ke}")

    def test_the_rebuilt_fit_lands_on_the_data_at_the_right_scale(self):
        np, ls = self.np, self.ls

        def shape(x, fwhm, mix):            # Gaussian/Lorentzian sum, mix = L share
            t = x / fwhm
            return ((1 - mix) * np.exp(-4 * np.log(2) * t * t)
                    + mix / (1 + 4 * t * t))

        rms, scale = [], []
        for path, prj in self.projects:
            cache = {}
            for di in prj.fitted_data_spaces():
                peaks = prj.peaks_of(di)
                bg = prj.background_of(di)
                if bg is None or not all(p.symmetric for p in peaks):
                    continue                    # the tail terms are not modelled
                r = self._spectrum(prj, di, cache)
                dwell, scans = r.dwell_and_scans()
                ke = r.photon_energy - np.array(r.energy)
                order = np.argsort(ke)
                ke = ke[order]
                y = np.array(r.counts)[order] / (dwell * scans)   # counts/s
                m = (ke >= bg.start_ke) & (ke <= bg.end_ke)
                if m.sum() < 10:
                    continue
                kk, yy = ke[m], y[m]
                base = (ls.linear_bg(yy, bg.avg_points) if bg.code == 0
                        else ls.shirley(yy, bg.avg_points))
                model = sum(p.height / (dwell * scans)
                            * shape(kk - p.ke, p.fwhm, p.mix) for p in peaks)
                rms.append(float(np.sqrt(np.mean((yy - base - model) ** 2))
                                 / (yy.max() - yy.min())))
                scale.append(float(np.dot(model, yy - base)
                                   / np.dot(model, model)))
        self.assertGreaterEqual(len(rms), 20)
        rms, scale = np.array(rms), np.array(scale)
        # measured on the 129 Knowledge projects: median 1.2 % of the range,
        # 97 % of the 67 spectra under 5 %, every scale within 15 % of 1
        self.assertLess(float(np.median(rms)), 0.02)
        self.assertGreater(float(np.mean(rms < 0.05)), 0.9)
        self.assertLess(float(np.max(np.abs(scale - 1))), 0.15)
        self.assertLess(abs(float(np.median(scale)) - 1), 0.02)

    def test_every_project_opens_as_a_session_with_its_fits(self):
        import casafit
        import quant
        np = self.np
        sym, tails, n_fits, modes, poor = [], [], 0, set(), []
        for path, prj in self.projects:
            doc = readers.load_file(path)
            self.assertIsInstance(doc, V.ThermoProjectFile, path)
            self.assertEqual(doc.errors, [], path)
            fitted = [r for r in doc.regions if r.fit is not None]
            n_fits += len(fitted)
            self.assertEqual(len(fitted) > 0, len(prj.peaks) > 0, path)
            for r in fitted:
                dwell, scans = r.dwell_and_scans()
                cv = casafit.curves(r.fit, r.energy, r.counts,
                                    r.photon_energy, dwell, scans)
                self.assertTrue(cv, path)
                self.assertTrue(cv[0].approximate)       # never claimed exact
                self.assertTrue(quant.fit_rows(r), path)
                if cv[0].residual_rms > V.POOR_FIT:
                    poor.append((r.source, any("unexplained" in n
                                               for n in cv[0].notes)))
                has_tail = any(
                    not (abs(t[0]) < 1e-9 and abs(t[1] - 1.0) < 1e-9)
                    for t in (self.ls.parse_shape(c.shape)["tail_terms"]
                              for c in r.fit.components))
                (tails if has_tail else sym).append(cv[0].residual_rms)
                modes.add(r.fit.components[0].shape.split("(")[0])
        self.assertGreater(n_fits, 100)
        self.assertEqual(modes, {"AV", "AVP"})            # both are found
        # measured on the 129 projects (114 fitted spectra): no tail median
        # 1.1 % of the range (p90 2.8 %, max 5.7 %); with the manual's tail
        # median 1.3 % (p90 2.8 %). One spectrum does not match its fit:
        # native Al oxide on Al's al2p_1, whose project lists no oxide peak
        # for it; it says so (the Cd native oxide file's fit sits on its
        # second trace and is read as the product shape, 7 %)
        self.assertLess(float(np.median(sym)), 0.02)
        self.assertLess(float(np.max(sym)), 0.08)
        self.assertLess(float(np.median(tails)), 0.025)
        self.assertLess(float(np.percentile(tails, 90)), 0.05)
        self.assertEqual([src for src, said in poor], ["al2p_1.VGD"])
        self.assertTrue(all(said for _src, said in poor))

    def test_a_two_level_file_gets_its_fit_on_the_level_it_describes(self):
        path = os.path.join(CORPUS, "Cd native oxide on metal.vgp")
        if not os.path.isfile(path):
            self.skipTest("Cd native oxide on metal.vgp not in the corpus")
        doc = readers.load_file(path)
        fitted = [r for r in doc.regions if r.fit is not None]
        self.assertEqual(len(fitted), 1)
        self.assertEqual(fitted[0].extra["avantage_level"], 1)
        self.assertEqual(fitted[0].fit.components[0].shape.split("(")[0], "AVP")
        import casafit
        r = fitted[0]
        dwell, scans = r.dwell_and_scans()
        cv = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                            dwell, scans)[0]
        self.assertLess(cv.residual_rms, 0.09)


class TestPeakDefinitions(unittest.TestCase):
    DEFS = [dict(lib="ALTHERMO1", name="C1s C-C", rsf=1.0),
            dict(lib="", name="Satellite", rsf=1.0, groups=2),
            dict(lib="ALTHERMO1", name="Pd3d5 metal", rsf=5.4, groups=2)]

    def setUp(self):
        self.prj = V.parse(make_project(defs=self.DEFS))

    def test_definitions_are_read(self):
        d = self.prj.peak_defs
        self.assertEqual([(x.library, x.name, x.groups) for x in d],
                         [("ALTHERMO1", "C1s C-C", 1), ("", "Satellite", 2),
                          ("ALTHERMO1", "Pd3d5 metal", 2)])

    def test_the_library_rsf_is_read_whatever_the_group_count(self):
        self.assertEqual([x.rsf for x in self.prj.peak_defs],
                         [1.0, None, 5.4])         # a hand-made peak: no RSF

    def test_each_fitted_peak_takes_its_definitions_rsf(self):
        self.assertEqual([(p.rsf, p.library) for p in self.prj.peaks],
                         [(1.0, "ALTHERMO1"), (None, ""),
                          (5.4, "ALTHERMO1")])

    def test_an_unreadable_layout_leaves_the_rsf_unknown(self):
        bad = make_project(defs=[dict(lib="ALTHERMO1", name="C1s", rsf=-3.0),
                                 dict(lib="ALTHERMO1", name="O1s", rsf=2.9),
                                 dict(lib="ALTHERMO1", name="N1s", rsf=1.8)])
        self.assertEqual([p.rsf for p in V.parse(bad).peaks],
                         [None, 2.9, 1.8])

    def test_a_project_without_definitions_has_none(self):
        prj = V.parse(make_project())
        self.assertEqual(prj.peak_defs, [])
        self.assertEqual([p.rsf for p in prj.peaks], [None] * 3)

    def test_the_components_carry_the_rsf_into_quantification(self):
        import quant
        import casafit
        from readers import thermo_vgp as vgp
        region = Region(name="C 1s", index=0, offset=0,
                        energy=[292.0 - 0.1 * i for i in range(121)],
                        counts=[100.0 + 5.0 * (i % 3) for i in range(121)],
                        decodable=True, photon_energy=1486.68, dwell=0.1,
                        step=0.1, pass_energy=20.0)
        region.extra["n_scans"] = 1
        fit = vgp.fit_for(self.prj.peaks_of(0)[:1], None, region, 0.1, 1)
        self.assertEqual([c.rsf for c in fit.components], [1.0])
        self.assertIn("KE^0.6", fit.regions[0].note)           # said, not applied
        fit2 = vgp.fit_for([self.prj.peaks[1]], None, region, 0.1, 1)
        self.assertEqual([c.rsf for c in fit2.components], [0.0])
        self.assertNotIn("KE^0.6", fit2.regions[0].note)
        self.assertIsNotNone(quant)
        self.assertIsNotNone(casafit)


@unittest.skipUnless(_corpus() and os.environ.get("XPS_ALTHERMO1_JSON"),
                     "set XPS_VGP_CORPUS and XPS_ALTHERMO1_JSON (ALTHERMO1 "
                     "dumped from QuantLib01 .mdb: lookup, rsf per row)")
class TestRsfAgainstTheLibrary(unittest.TestCase):
    def test_every_library_peak_has_the_librarys_sensitivity(self):
        import collections
        import json
        with open(os.environ["XPS_ALTHERMO1_JSON"],
                  encoding="utf-8-sig") as fh:
            rows = json.load(fh)
        lib = collections.defaultdict(set)
        for r in rows:
            lib[str(r["lookup"]).lower()].add(r["rsf"])
        n = 0
        for f in _corpus():
            for d in V.load(f).peak_defs:
                self.assertIsNotNone(d.rsf if d.library else 1.0, f)
                if d.library != "ALTHERMO1":
                    self.assertIsNone(d.rsf)
                    continue
                vals = lib.get((d.name.split() or [""])[0].lower())
                if not vals:
                    continue                 # a doublet entered as one peak
                n += 1
                self.assertTrue(any(abs(d.rsf - v) <= 2e-3 * max(1, v)
                                    for v in vals), (f, d.name, d.rsf))
        self.assertGreater(n, 400)


class TestSmartBackground(unittest.TestCase):
    """``lineshapes.smart_bg``: Shirley that never exceeds the data."""

    def peak(self, n=121, base=100.0, slope=0.0):
        import numpy as np
        x = np.linspace(0.0, 12.0, n)
        return x, base + slope * x + 800.0 * np.exp(-((x - 6.0) / 0.8) ** 2)

    def test_it_never_exceeds_the_data(self):
        import numpy as np
        import lineshapes as ls
        rng = np.random.default_rng(3)
        for slope in (0.0, 25.0, -25.0):
            _x, y = self.peak(slope=slope)
            y = y + rng.normal(0.0, 4.0, len(y))
            b = ls.smart_bg(y, avg=5)
            self.assertTrue(np.all(b <= y + 1e-9))

    def test_it_is_shirley_where_shirley_stays_below_the_data(self):
        import numpy as np
        import lineshapes as ls
        _x, y = self.peak()
        np.testing.assert_allclose(ls.smart_bg(y), ls.shirley(y), atol=1e-9)

    def test_it_follows_the_data_where_shirley_would_rise_above_it(self):
        import numpy as np
        import lineshapes as ls
        _x, y = self.peak(slope=-40.0)          # the low-KE end is the high one
        y[:8] = y[:8] - 150.0                   # a dip: the end level is too high
        plain = ls.shirley(y, avg=8)
        smart = ls.smart_bg(y, avg=8)
        self.assertTrue(np.any(plain > y + 1.0))
        self.assertTrue(np.all(smart <= y + 1e-9))

    def test_background_dispatch_knows_smart(self):
        import lineshapes as ls
        _x, y = self.peak()
        import numpy as np
        np.testing.assert_allclose(ls.background("Smart", y, 1, params=(0, 0)),
                                   ls.smart_bg(y))


if __name__ == "__main__":
    unittest.main()