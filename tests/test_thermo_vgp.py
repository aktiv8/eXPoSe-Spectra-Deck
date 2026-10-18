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
from readers import thermo_vgp as V  # noqa: E402

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


PEAKS = [
    dict(name="C1s C-C", ke=1202.2, fwhm=1.3, mix=0.23, area=5566.0,
         height=94.0),
    dict(name="", ke=1200.9, fwhm=1.3, mix=0.23, area=680.0, height=11.5,
         ver=2),                                        # empty name, long prefix
    dict(name="Pd3d5 metal", ke=1151.3, fwhm=0.87, mix=0.26, area=3398.0,
         height=7719.0, tail=(1e-17, 0.367, 0.0428)),
]


def make_project(peaks=PEAKS, bgs=None, paths=None, xr_peak_ds=None,
                 count=None):
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
        "Peaks": {"Contents": b"\x01\x00\x00\x00\x00\x00\x00\x00",
                  "XRefs": b""},
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
        import readers
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


if __name__ == "__main__":
    unittest.main()
