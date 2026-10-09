"""The NeXus ``.nxs`` reader for beamline NEXAFS scans (Diamond B07 / GDA):
the layout read from a synthetic file, the numbers against the text exports of
four real scans, the axis (photon energy, never a binding energy), the methods
text, the VAMAS round trip and what the NeXus export does with them.

Real files: ``XPS_NEXAFS_CORPUS=<folder with b07-*.nxs and XY/b07-*_NEXAFS.dat>``.

Run:  python -m unittest discover tests
"""

import glob
import os
import re
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import annotations  # noqa: E402
import exporters  # noqa: E402
import methods  # noqa: E402
import nexafs  # noqa: E402
import readers  # noqa: E402
import viewdata  # noqa: E402
from readers import nexus_nexafs  # noqa: E402

try:
    import h5py
    import numpy as np
    HAVE_H5 = True
except ImportError:                                 # pragma: no cover
    HAVE_H5 = False

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False


def text_of(path, encoding="latin-1"):
    with open(path, encoding=encoding) as fh:
        return fh.read()


def write_scan(path, start=270.0, stop=280.0, step=0.5, channels=("ca18b",
               "ca35b"), t0="2022-03-01T13:43:41.080Z",
               t1="2022-03-01T13:58:04.586Z", units="A", jitter=0.002):
    """A small B07-style file: one group per channel, a scattered readback."""
    n = int(round((stop - start) / step)) + 1
    grid = start + step * np.arange(n)
    rng = np.random.RandomState(3)
    readback = grid + rng.uniform(-jitter, jitter, n)
    with h5py.File(path, "w") as f:
        e = f.create_group("entry1")
        e.attrs["NX_class"] = b"NXentry"
        e["scan_command"] = (f"scan pgm_energy {start:g} {stop:g} {step:g} "
                             + " ".join(f"{c} 0.2" for c in channels)
                             + " ring_current").encode()
        e["start_time"], e["end_time"] = t0.encode(), t1.encode()
        e["program_name"] = b"GDA 9.25.0pre"
        e["experiment_identifier"] = b"si30483-1"
        e["user01/username"] = b"joy55749"
        e["instrument/name"] = b"b07"
        e["instrument/source/name"] = b"DLS"
        e["instrument/source/type"] = b"Synchrotron X-Ray Source"
        ring = 300.0 - 0.01 * np.arange(n)
        e["instrument/ring_current/ring_current"] = ring
        for i, c in enumerate(channels):
            g = e.create_group(c)
            ds = g.create_dataset(c, data=-1e-10 * (i + 1) * (1 + grid / 300))
            if units:
                ds.attrs["units"] = units.encode()
            gain = g.create_dataset("gain", data=np.full(n, 1e9))
            gain.attrs["units"] = b"V/A"
            pe = g.create_dataset("pgm_energy", data=readback)
            pe.attrs["primary"] = b"1"
            e[f"instrument/{c}/count_time"] = np.array([0.2])
            e[f"instrument/{c}/mode"] = np.array([b"Low Noise"], dtype=object)
        e["before_scan/pgm_grating/pgm_grating"] = b"600 l/mm Au"
        e["before_scan/pgm_mirror/pgm_mirror"] = b"Rh Stripe"
        e["before_scan/pgm_cff/pgm_cff"] = 2.25
    return readback


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestSyntheticFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        cls.readback = write_scan(cls.path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def load(self):
        return readers.load_file(self.path)

    def test_one_region_per_channel_named_by_detector(self):
        f = self.load()
        self.assertIsInstance(f, nexus_nexafs.NexusNexafsFile)
        self.assertEqual([r.name for r in f.regions], ["ca18b", "ca35b"])
        r = f.regions[0]
        self.assertEqual((r.technique, r.energy_label, r.energy_units),
                         ("NEXAFS", "Photon Energy", "eV"))
        self.assertEqual((r.count_label, r.count_units), ("Current", "A"))
        self.assertIsNone(r.photon_energy)

    def test_the_readback_is_kept_as_read_not_regenerated(self):
        r = self.load().regions[0]
        self.assertEqual(r.energy, [float(x) for x in self.readback])
        self.assertAlmostEqual(r.step, 0.5, places=2)
        self.assertEqual(r.n_points, 21)

    def test_a_unit_the_file_does_not_state_is_not_invented(self):
        p = os.path.join(self.dir, "nounits.nxs")
        write_scan(p, units="")
        r = readers.load_file(p).regions[0]
        self.assertEqual((r.count_label, r.count_units), ("Signal", "arb."))

    def test_times_dwell_and_instrument(self):
        f = self.load()
        r = f.regions[0]
        self.assertEqual(r.extra["t_start"], "2022-03-01 13:43:41")
        self.assertEqual(r.extra["tz"], "UTC")
        self.assertEqual(r.dwell, 0.2)
        self.assertEqual(f.instrument["Instrument"],
                         "B07 (Diamond Light Source)")
        self.assertEqual(f.instrument["Operator"], "joy55749")

    def test_an_offset_of_whole_hours_is_converted_to_utc(self):
        wall, utc = nexus_nexafs.parse_time("2022-07-26T02:38:40.707+01")
        self.assertEqual((wall, utc), ("2022-07-26 02:38:40",
                                       "2022-07-26 01:38:40"))
        self.assertEqual(nexus_nexafs.parse_time("2022-03-01T13:43:41.080Z"),
                         ("2022-03-01 13:43:41", "2022-03-01 13:43:41"))
        self.assertEqual(nexus_nexafs.parse_time("2022-03-01 13:43:41"),
                         ("2022-03-01 13:43:41", ""))
        self.assertEqual(nexus_nexafs.parse_time("soon"), ("", ""))

    def test_a_wide_scan_is_not_a_survey(self):
        p = os.path.join(self.dir, "wide.nxs")
        write_scan(p, 1500, 1800, 1.0, ("ca18b",))
        r = readers.load_file(p).regions[0]
        self.assertGreater(max(r.energy) - min(r.energy), 250)
        self.assertFalse(r.is_survey)

    def test_metadata_keys_only_when_the_file_has_them(self):
        md = self.load().region_metadata(self.load().regions[0])
        self.assertEqual(md["Technique"], "NEXAFS")
        self.assertEqual(md["Detector mode"], "Low Noise")
        self.assertEqual(md["Detector gain"], "1e+09 V/A")
        self.assertEqual(md["Monochromator grating"], "600 l/mm Au")
        self.assertIn("Photon energy start (eV)", md)
        self.assertNotIn("KE start (eV)", md)
        self.assertEqual(md["Photon energy (eV)"], "")
        self.assertEqual(md["Anode"], "")
        self.assertNotIn("Detector coupling", md)       # not written, not shown

    def test_it_is_not_taken_for_a_kfit_nor_a_kfit_for_it(self):
        head = b"\x89HDF\r\n\x1a\n" + b"\0" * 100
        self.assertTrue(nexus_nexafs.sniff(head, ".nxs"))
        self.assertFalse(nexus_nexafs.sniff(head, ".kfit"))
        self.assertIs(readers.reader_for(self.path),
                      nexus_nexafs.NexusNexafsFile)

    def test_hdf5_that_is_not_a_scan_says_so(self):
        p = os.path.join(self.dir, "other.nxs")
        with h5py.File(p, "w") as f:
            g = f.create_group("entry")
            g.attrs["NX_class"] = b"NXentry"
            g["title"] = b"nothing to see"
        with self.assertRaises(ValueError) as cm:
            readers.load_file(p)
        self.assertIn("holds no photon-energy scan", str(cm.exception))
        q = os.path.join(self.dir, "plain.nxs")
        with h5py.File(q, "w") as f:
            f["x"] = [1, 2, 3]
        with self.assertRaises(ValueError) as cm:
            readers.load_file(q)
        self.assertIn("no NXentry", str(cm.exception))

    def test_a_channel_of_the_wrong_length_is_left_out_with_a_note(self):
        p = os.path.join(self.dir, "short.nxs")
        write_scan(p)
        with h5py.File(p, "r+") as f:
            del f["entry1/ca35b/ca35b"]
            f["entry1/ca35b/ca35b"] = [1.0, 2.0]
        f = readers.load_file(p)
        self.assertEqual([r.name for r in f.regions], ["ca18b"])
        self.assertTrue(any("ca35b" in w for w in f.warnings))


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestAxisAndText(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)
        cls.f = readers.load_file(cls.path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_neither_energy_scale_converts_the_photon_energy(self):
        r = self.f.regions[0]
        for scale in viewdata.ENERGY_SCALES:
            a = viewdata.energy_axis(r, scale)
            self.assertEqual(a.x, r.energy)
            self.assertEqual(a.label, "Photon Energy")
            self.assertFalse(a.invert)

    def test_the_methods_text_describes_nexafs_not_xps(self):
        text = methods.generate(self.f.metadata_rows())
        self.assertIn("Near-edge X-ray absorption fine structure (NEXAFS)",
                      text)
        self.assertIn("B07 (Diamond Light Source)", text)
        self.assertIn("scanned from 270 to 280 eV", text)
        self.assertIn("step size of 0.5 eV", text)
        self.assertIn("dwell time of 0.2 s per point", text)
        self.assertIn("ca18b, ca35b", text)
        self.assertNotIn("photoelectron", text)
        self.assertNotIn("charge-corrected", text)
        self.assertNotIn("Anode", text)

    def test_next_to_xps_each_gets_its_own_sentence(self):
        xps = {"Technique": "XPS", "Region": "C 1s", "Instrument": "AXIS",
               "Photon energy (eV)": "1486.69", "Anode": "Al",
               "Pass energy (eV)": "20", "Step (eV)": "0.1"}
        text = methods.generate([xps] + self.f.metadata_rows())
        self.assertIn("X-ray photoelectron spectroscopy (XPS) measurements "
                      "were made on a AXIS spectrometer", text)
        self.assertIn("NEXAFS", text)
        self.assertEqual(text.count("B07"), 1)

    def test_csv_names_the_axis(self):
        p = os.path.join(self.dir, "x.csv")
        exporters.export_csv(self.f.regions, p)
        head = text_of(p, "utf-8-sig").splitlines()[0]
        self.assertIn("ca18b Photon Energy (eV)", head)
        self.assertIn("ca18b Current (A)", head)

    def test_nexus_export_leaves_them_out_with_a_reason(self):
        import nexus_export
        if not nexus_export.HAVE_H5PY:                  # pragma: no cover
            self.skipTest("h5py missing")
        with self.assertRaises(ValueError) as cm:
            nexus_export.export_nexus(self.f.regions,
                                      os.path.join(self.dir, "y.nxs"))
        self.assertIn("NXxps", str(cm.exception))


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestVamasRoundTrip(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)
        cls.f = readers.load_file(cls.path)
        cls.vms = os.path.join(cls.dir, "b07-1.vms")
        exporters.export_vamas(cls.f.regions, cls.vms,
                               metadata=cls.f.metadata_rows())
        cls.g = readers.load_file(cls.vms)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_armoury_convention_on_the_way_out(self):
        text = text_of(self.vms).splitlines()
        self.assertIn("AES dir", text)
        self.assertNotIn("XPS", text)
        self.assertTrue(any(ln.startswith("Technique : NEXAFS") for ln in text))

    def test_it_reads_back_as_nexafs_with_the_same_numbers(self):
        for a, b in zip(self.f.regions, self.g.regions):
            self.assertEqual((b.technique, b.energy_label),
                             ("NEXAFS", "Photon Energy"))
            self.assertIsNone(b.photon_energy)
            self.assertFalse(b.is_survey)
            self.assertEqual((b.name, b.sample, b.count_units),
                             (a.name, a.sample, "A"))
            for x, y in zip(a.counts, b.counts):         # currents of 1e-10 A
                self.assertAlmostEqual(x / y, 1.0, places=8)
            # a regular grid through a readback that scatters by 0.002 eV
            self.assertLess(max(abs(x - y) for x, y in
                                zip(a.energy, b.energy)), 0.01)

    def test_the_regularising_is_said(self):
        lines = self.g.regions[0].extra["comment_lines"]
        self.assertTrue(any("regularised" in ln for ln in lines))

    def test_a_casa_style_file_with_only_the_comment_line_is_recognised(self):
        # written without the metadata block: "Technique : NEXAFS" alone
        p = os.path.join(self.dir, "casa.vms")
        exporters.export_vamas(self.f.regions, p)
        text = text_of(p)
        self.assertIn("Technique : NEXAFS", text)
        self.assertNotIn("metadata ===", text)
        r = readers.load_file(p).regions[0]
        self.assertEqual((r.technique, r.energy_label),
                         ("NEXAFS", "Photon Energy"))
        self.assertIsNone(r.photon_energy)

    def test_a_plain_aes_dir_file_is_left_alone(self):
        from readers import Region
        r = Region("Cu LMM", 0, 0, technique="XPS", energy=[1.0, 2.0, 3.0],
                   counts=[5.0, 6.0, 7.0], energy_label="Kinetic Energy",
                   decodable=True)
        q = os.path.join(self.dir, "kin.vms")
        exporters.export_vamas([r], q)
        txt = text_of(q).replace("\r\n", "\n")
        self.assertEqual(txt.count("\nXPS\n"), 1)
        p = os.path.join(self.dir, "aes.vms")
        with open(p, "w", encoding="latin-1", newline="\r\n") as fh:
            fh.write(txt.replace("\nXPS\n", "\nAES dir\n"))
        g = readers.load_file(p).regions[0]
        self.assertNotEqual(g.technique, "NEXAFS")
        self.assertNotEqual(g.energy_label, "Photon Energy")


        self.assertNotEqual(g.technique, "NEXAFS")
        self.assertNotEqual(g.energy_label, "Photon Energy")


@unittest.skipUnless(HAVE_H5 and HAVE_MPL, "h5py / matplotlib / Tk missing")
class TestInTheWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror, ee.messagebox.askyesno)
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None
        ee.messagebox.askyesno = lambda *a, **k: True
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror, ee.messagebox.askyesno) = cls._boxes
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_it_opens_ticks_and_draws_on_a_photon_energy_axis(self):
        ws = ee.Workspace(self.root)
        ws._add_files([self.path])
        regions = ws.docs[0].regions
        ws.checked = {id(r) for r in regions}
        for scale in ("Binding", "Kinetic"):
            ws.scale_var.set(scale)
            ws._render()
            self.root.update_idletasks()
            labels = [ax.get_xlabel() for ax in ws.fig.axes]
            self.assertTrue(any("photon energy" in t.lower() for t in labels), labels)
            ax = ws.fig.axes[0]
            self.assertLess(ax.get_xlim()[0], ax.get_xlim()[1])   # not inverted

    def test_the_html_page_gets_the_axis_and_no_binding_flag(self):
        import htmlbrowser
        ws = ee.Workspace(self.root)
        ws._add_files([self.path])
        data = htmlbrowser.build_payload(ws.docs)
        specs = data["files"][0]["samples"][0]["spectra"] \
            if "samples" in data["files"][0] else None
        flat = []

        def walk(o):
            if isinstance(o, dict):
                if "elabel" in o:
                    flat.append(o)
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for v in o:
                    walk(v)
        walk(data)
        del specs
        self.assertEqual(len(flat), 2)
        for s in flat:
            self.assertEqual((s["elabel"], s["eunits"], s["binding"]),
                             ("Photon Energy", "eV", False))


class TestRingScalingMaths(unittest.TestCase):
    def test_each_point_is_scaled_to_the_mean_ring_current(self):
        counts, mean = nexafs.scale_to_mean([1.0, 1.0, 1.0], [300.0, 200.0,
                                                              100.0])
        self.assertEqual(mean, 200.0)
        self.assertEqual(counts, [200.0 / 300.0, 1.0, 2.0])

    def test_a_constant_ring_current_changes_nothing(self):
        counts, mean = nexafs.scale_to_mean([1.0, 2.0], [250.0, 250.0])
        self.assertEqual((counts, mean), ([1.0, 2.0], 250.0))

    def test_a_point_without_a_usable_ring_current_is_nan_not_guessed(self):
        counts, mean = nexafs.scale_to_mean(
            [1.0, 1.0, 1.0, 1.0], [300.0, 0.0, float("nan"), 100.0])
        self.assertEqual(mean, 200.0)                  # of the usable ones
        self.assertEqual(counts[0], 200.0 / 300.0)
        self.assertTrue(counts[1] != counts[1] and counts[2] != counts[2])
        self.assertEqual(counts[3], 2.0)

    def test_nothing_usable_gives_nothing(self):
        self.assertEqual(nexafs.scale_to_mean([1.0], [0.0]), (None, None))
        self.assertEqual(nexafs.scale_to_mean([1.0, 2.0], [5.0]), (None, None))

    def test_the_annotation_is_saved_and_read_back(self):
        a = annotations.Annotations()
        self.assertTrue(a.is_empty())
        a.nexafs_ring = True
        self.assertFalse(a.is_empty())
        b = annotations.Annotations.from_json(a.to_json())
        self.assertTrue(b.nexafs_ring)
        self.assertFalse(annotations.Annotations.from_json(
            {"nexafs_ring": "yes"}).nexafs_ring)       # only a real true
        self.assertNotIn("nexafs_ring", annotations.Annotations().to_json())


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestRingScalingOfAFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)
        cls.f = readers.load_file(cls.path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_the_reader_keeps_the_ring_current_of_every_point(self):
        r = self.f.regions[0]
        ring = nexafs.ring_points(r)
        self.assertEqual(len(ring), r.n_points)
        self.assertAlmostEqual(ring[0], 300.0)
        self.assertAlmostEqual(ring[-1], 300.0 - 0.01 * 20)

    def test_a_file_without_a_ring_current_has_nothing_to_scale(self):
        p = os.path.join(self.dir, "noring.nxs")
        write_scan(p)
        with h5py.File(p, "r+") as h:
            del h["entry1/instrument/ring_current"]
        r = readers.load_file(p).regions[0]
        self.assertIsNone(nexafs.ring_points(r))

    def test_the_metadata_and_methods_text_say_so_only_when_chosen(self):
        a = annotations.Annotations()
        self.f.annotations, self.f.file_id = a, "f1"
        plain = self.f.region_metadata(self.f.regions[0])
        self.assertNotIn("Normalisation", plain)
        self.assertIn("the spectra are not normalised",
                      methods.generate(self.f.metadata_rows()))
        a.nexafs_ring = True
        md = self.f.region_metadata(self.f.regions[0])
        self.assertIn("scaled to the mean ring current", md["Normalisation"])
        self.assertIn("unit is not recorded", md["Normalisation"])
        text = methods.generate(self.f.metadata_rows())
        self.assertIn("scaled to the mean ring current of its scan", text)
        self.assertNotIn("not normalised", text)
        self.f.annotations = None


@unittest.skipUnless(HAVE_H5 and HAVE_MPL, "h5py / matplotlib / Tk missing")
class TestRingScalingInTheWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls._boxes = (ee.messagebox.showinfo, ee.messagebox.showwarning,
                      ee.messagebox.showerror, ee.messagebox.askyesno)
        ee.messagebox.showinfo = ee.messagebox.showwarning = \
            ee.messagebox.showerror = lambda *a, **k: None
        ee.messagebox.askyesno = lambda *a, **k: True
        cls.dir = tempfile.mkdtemp()
        cls.path = os.path.join(cls.dir, "b07-1.nxs")
        write_scan(cls.path)

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror, ee.messagebox.askyesno) = cls._boxes
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def setUp(self):
        self.ws = ee.Workspace(self.root)
        self.ws._add_files([self.path])
        self.region = self.ws.docs[0].regions[0]
        self.ws.checked = {id(r) for r in self.ws.docs[0].regions}

    def choose(self, on):
        self.ws.nexafs_ring_var.set(on)
        self.ws.set_nexafs_ring()

    def test_off_by_default_and_the_data_are_what_the_file_holds(self):
        self.assertFalse(self.ws.ann.nexafs_ring)
        self.assertIs(self.ws._display(self.region), self.region)

    def test_the_choice_scales_the_drawn_copy_never_the_reader_output(self):
        original = list(self.region.counts)
        self.choose(True)
        shown = self.ws._display(self.region)
        self.assertEqual(self.region.counts, original)
        ring = nexafs.ring_points(self.region)
        mean = sum(ring) / len(ring)
        self.assertAlmostEqual(shown.counts[-1] / original[-1],
                               mean / ring[-1])
        self.assertEqual(shown.count_label, "Normalised current")
        self.assertEqual(shown.count_units, "A")
        self.assertAlmostEqual(shown.extra["ring_norm"], mean)
        self.choose(False)
        self.assertIs(self.ws._display(self.region), self.region)

    def test_it_draws_and_exports_the_scaled_values(self):
        self.choose(True)
        self.ws._render()
        self.root.update_idletasks()
        self.assertTrue(any("Normalised current" in ax.get_ylabel()
                            for ax in self.ws.fig.axes))
        p = os.path.join(self.dir, "scaled.csv")
        exporters.export_csv([self.ws._display(r)
                              for r in self.ws.docs[0].regions], p)
        head = text_of(p, "utf-8-sig").splitlines()[0]
        self.assertIn("Normalised current (A)", head)

    def test_the_menu_box_follows_the_annotation_and_a_saved_workbook(self):
        self.choose(True)
        self.assertTrue(self.ws.nexafs_ring_var.get())
        saved = self.ws.ann.to_json()
        self.ws.ann = annotations.Annotations()
        self.ws._ann_changed(relabel=False)
        self.assertFalse(self.ws.nexafs_ring_var.get())
        self.ws.ann = annotations.Annotations.from_json(saved)
        self.ws._ann_changed(relabel=False)
        self.assertTrue(self.ws.nexafs_ring_var.get())
        self.assertIsNot(self.ws._display(self.region), self.region)

    def test_a_spectrum_that_is_not_nexafs_is_untouched(self):
        from readers import Region
        r = Region("C 1s", 0, 0, energy=[285.0, 284.0], counts=[1.0, 2.0],
                   decodable=True)
        r.extra["ring_current_points"] = [300.0, 100.0]     # never used
        self.ws.ann.nexafs_ring = True
        self.assertIsNone(self.ws._ring_scaled(r))


@unittest.skipUnless(HAVE_H5 and os.environ.get("XPS_NEXAFS_CORPUS"),
                     "set XPS_NEXAFS_CORPUS to the folder of b07-*.nxs files")
class TestRealFiles(unittest.TestCase):
    """Each channel against the text export of the same scan (the .dat
    carries the values to 8 digits: 5e-8 relative)."""

    @classmethod
    def setUpClass(cls):
        cls.root = os.environ["XPS_NEXAFS_CORPUS"]
        cls.files = sorted(glob.glob(os.path.join(cls.root, "*.nxs")))

    def test_there_are_files(self):
        self.assertGreaterEqual(len(self.files), 1)

    def test_every_channel_equals_the_text_export(self):
        for path in self.files:
            stem = os.path.splitext(os.path.basename(path))[0]
            dat = os.path.join(self.root, "XY", f"{stem}_NEXAFS.dat")
            if not os.path.exists(dat):
                continue
            table = np.genfromtxt(dat, names=True)
            f = readers.load_file(path)
            self.assertEqual([r.name for r in f.regions],
                             ["ca18b", "ca35b", "ca36b"], stem)
            for r in f.regions:
                self.assertEqual(r.n_points, len(table), stem)
                np.testing.assert_allclose(r.energy, table["pgm_energy"],
                                           atol=1e-4, rtol=0)
                np.testing.assert_allclose(r.counts, table[r.name],
                                           rtol=1e-7, atol=1e-30)

    def test_the_range_matches_the_scan_command(self):
        for path in self.files:
            f = readers.load_file(path)
            cmd = f.regions[0].extra["nexafs"]["scan_command"]
            lo, hi, step = (float(x) for x in re.match(
                r"scan pgm_energy (\S+) (\S+) (\S+)", cmd).groups())
            r = f.regions[0]
            self.assertEqual(r.n_points, round((hi - lo) / step) + 1)
            self.assertAlmostEqual(r.energy[0], lo, delta=0.05)
            self.assertAlmostEqual(r.energy[-1], hi, delta=0.05)
            self.assertAlmostEqual(r.step, step, delta=0.002)

    def test_times_instrument_and_metadata(self):
        for path in self.files:
            f = readers.load_file(path)
            md = f.region_metadata(f.regions[0])
            self.assertTrue(md["Run started"].endswith("UTC"), path)
            self.assertTrue(md["Run finished"].endswith("UTC"), path)
            self.assertEqual(md["Instrument"], "B07 (Diamond Light Source)")
            self.assertEqual(md["Monochromator grating"], "600 l/mm Au")

    def test_a_kfit_is_not_claimed(self):
        self.assertFalse(nexus_nexafs.sniff(b"\x89HDF\r\n\x1a\n", ".kfit"))

    def test_the_whole_folder_loads_and_exports_to_vamas_and_back(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        for path in self.files:
            f = readers.load_file(path)
            p = os.path.join(d, "x.vms")
            exporters.export_vamas(f.regions, p, metadata=f.metadata_rows())
            g = readers.load_file(p)
            for a, b in zip(f.regions, g.regions):
                self.assertEqual(b.technique, "NEXAFS")
                np.testing.assert_allclose(b.counts, a.counts, rtol=1e-8)
                self.assertLess(np.max(np.abs(np.array(b.energy)
                                              - np.array(a.energy))), 0.06)


if __name__ == "__main__":
    unittest.main()
