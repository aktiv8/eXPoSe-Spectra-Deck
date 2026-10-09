"""The NeXus ``.nxs`` reader for analyser (NXxps) spectra.

No instrument-written XPS ``.nxs`` file exists on the author's machine (the 200
there are all Diamond NEXAFS scans), so what is checked is (1) a round trip of
every field ``nexus_export`` writes, on synthetic regions and, with
``XPS_NXXPS_CORPUS=<folder of spectra files>``, on every real file of a folder
(Kratos, Avantage, VAMAS from CasaXPS, ...), and (2) a hand-written file in the
layout the NeXus definitions name (other group names, no metadata note), which
is what a converter such as pynxtools-xps produces. Run:

    python -m unittest discover tests
"""

import glob
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import nexus_export  # noqa: E402
import readers  # noqa: E402
from readers import Region, nexus_nexafs, nexus_xps  # noqa: E402

try:
    import h5py
    import numpy as np
    HAVE_H5 = True
except ImportError:                                 # pragma: no cover
    HAVE_H5 = False

from test_nexus_nexafs import write_scan  # noqa: E402

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False


def region(name="C 1s", sample="PET", **kw):
    n = 21
    energy = [292.0 - 0.1 * i for i in range(n)]       # binding, high to low
    r = Region(name, 0, 0, energy=energy,
               counts=[100.0 + i * i for i in range(n)],
               energy_label="Binding Energy", decodable=True, sample=sample,
               photon_energy=1486.69, pass_energy=20.0, dwell=0.1, step=0.1,
               lens_mode="Hybrid", anode="Al Kα", date="2026-08-24 10:00:00")
    r.extra.update(n_scans=3, t_start="2026-08-24 10:00:00",
                   t_end="2026-08-24 10:05:00", tz="UTC",
                   analyser_mode="Constant analyser energy (CAE)")
    r.conditions.update({"X-ray Power": "225 W", "Anode voltage (kV)": "15",
                         "Emission current (mA)": "15", "X-ray spot (µm)": "300"})
    r.pos_x, r.pos_y = 12.3456789, -4.5
    r.tf_ke = [1000.0, 1100.0, 1200.0, 1300.0]
    r.tf_values = [1.0, 0.9, 0.8, 0.7]
    for k, v in kw.items():
        setattr(r, k, v)
    return r


def write(path, regions, instrument=None, **kw):
    mds = [{"Region": r.name, "Sample": r.sample, "Technique": r.technique,
            "Dwell (s)": "" if not r.dwell else f"{r.dwell:g}"}
           for r in regions]
    nexus_export.export_nexus(regions, path, metadata=mds,
                              instrument=instrument or {}, **kw)


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestRoundTrip(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def trip(self, regions, instrument=None, name="x.nxs"):
        p = os.path.join(self.dir, name)
        write(p, regions, instrument)
        return readers.load_file(p)

    def test_the_registry_gives_one_reader_for_both_kinds(self):
        p = os.path.join(self.dir, "x.nxs")
        write(p, [region()])
        self.assertIs(readers.reader_for(p), nexus_xps.NexusFile)
        self.assertTrue(issubclass(nexus_xps.NexusFile,
                                   nexus_nexafs.NexusNexafsFile))
        f = readers.load_file(p)
        self.assertEqual(f.format_name, "NeXus NXxps (.nxs)")

    def test_the_numbers_come_back_exactly(self):
        a = region()
        b = self.trip([a]).regions[0]
        self.assertEqual((b.name, b.sample), ("C 1s", "PET"))
        self.assertEqual(b.energy, a.energy)
        self.assertEqual(b.counts, a.counts)
        self.assertEqual((b.energy_label, b.energy_units, b.technique),
                         ("Binding Energy", "eV", "XPS"))
        self.assertEqual((b.photon_energy, b.pass_energy, b.step),
                         (1486.69, 20.0, 0.1))
        self.assertEqual(b.dwell_and_scans(), (0.1, 3))
        self.assertEqual((b.pos_x, b.pos_y), (12.3456789, -4.5))
        self.assertEqual(b.lens_mode, "Hybrid")
        self.assertEqual(b.extra["analyser_mode"],
                         "Constant analyser energy (CAE)")

    def test_source_conditions_and_times(self):
        b = self.trip([region()]).regions[0]
        self.assertEqual(b.conditions["X-ray Power"], "225 W")
        self.assertEqual(b.conditions["Anode voltage (kV)"], "15")
        self.assertEqual(b.conditions["Emission current (mA)"], "15")
        self.assertEqual(b.conditions["X-ray spot (µm)"], "300")
        self.assertEqual((b.extra["t_start"], b.extra["t_end"], b.extra["tz"]),
                         ("2026-08-24 10:00:00", "2026-08-24 10:05:00", "UTC"))
        self.assertEqual(b.anode, "Al Kα")

    def test_the_transmission_function_comes_back(self):
        a = region()
        b = self.trip([a]).regions[0]
        ta, tb = a.transmission(), b.transmission()
        self.assertEqual(len(ta), len(tb))
        for x, y in zip(ta, tb):
            self.assertAlmostEqual(x, y, places=9)

    def test_the_order_of_the_spectra_is_the_order_written(self):
        names = ["Survey", "C 1s", "O 1s", "Cl2p 2", "Cl 2p", "Au 4f"]
        regions = [region(n) for n in names]
        f = self.trip(regions)
        self.assertEqual([r.name for r in f.regions],
                         ["Survey", "C 1s", "O 1s", "Cl2p 2", "Cl 2p",
                          "Au 4f"])

    def test_a_correction_already_in_the_axis_is_not_added_again(self):
        r = region(shift_applied=0.5, calibration_shift=0.5)
        b = self.trip([r]).regions[0]
        self.assertEqual(b.calibration_shift, 0.0)
        self.assertEqual(b.extra["nexus_referencing"],
                         {"offset": 0.5, "applied": True})

    def test_a_correction_the_file_still_has_to_add_is_kept_as_such(self):
        r = region(calibration_shift=-1.25)
        r.extra["casa_calib"] = {"measured": 285.4, "assigned": 284.8,
                                 "shift": -1.25, "inherited": False}
        b = self.trip([r]).regions[0]
        self.assertEqual(b.calibration_shift, -1.25)
        self.assertEqual(b.extra["casa_calib"]["measured"], 285.4)
        self.assertEqual(b.extra["casa_calib"]["assigned"], 284.8)

    def test_depth_profile_fields(self):
        r = region(etch_level=4, etch_time=120.0)
        b = self.trip([r]).regions[0]
        self.assertEqual((b.etch_level, b.etch_time), (4, 120.0))

    def test_a_kinetic_axis_without_photon_energy_stays_kinetic(self):
        r = Region("ISS", 0, 0, technique="ISS", energy=[100.0, 500.0, 933.0],
                   counts=[1.0, 2.0, 3.0], energy_label="Kinetic Energy",
                   decodable=True, sample="Au")
        b = self.trip([r]).regions[0]
        self.assertEqual((b.energy_label, b.technique, b.energy),
                         ("Kinetic Energy", "ISS", [100.0, 500.0, 933.0]))
        self.assertIsNone(b.photon_energy)
        self.assertTrue(b.is_iss)

    def test_instrument_and_user(self):
        f = self.trip([region()], {"Instrument": "AXIS Supra",
                                   "Operator": "dm",
                                   "Institution": "Cardiff"})
        self.assertEqual(f.instrument["Instrument"], "AXIS Supra")
        self.assertEqual(f.instrument["Operator"], "dm")
        self.assertEqual(f.instrument["Institution"], "Cardiff")

    def test_the_whole_metadata_is_the_same_afterwards(self):
        a = region()
        f = self.trip([a])
        b = f.regions[0]
        md_a = {"Region": a.name, "Sample": a.sample}
        md_b = f.region_metadata(b)
        for k, v in md_a.items():
            self.assertEqual(md_b[k], v)
        self.assertEqual(md_b["Photon energy (eV)"], "1486.69")
        self.assertEqual(md_b["Pass energy (eV)"], "20")
        self.assertEqual(md_b["Scans"], "3")
        self.assertEqual(md_b["Run started"], "2026-08-24 10:00:00 UTC")

    def test_the_kratos_dwell_summed_over_sweeps_is_restored(self):
        a = region(dwell=0.6)           # as .experiment holds it: all sweeps
        a.extra["dwell_total"] = True
        p = os.path.join(self.dir, "kratos.nxs")
        mds = [{"Region": "C 1s", "Sample": "PET", "Dwell (s)": "0.6"}]
        nexus_export.export_nexus([a], p, metadata=mds)
        b = readers.load_file(p).regions[0]
        self.assertTrue(b.extra["dwell_total"])
        self.assertAlmostEqual(b.dwell, 0.6)
        self.assertEqual(b.dwell_and_scans()[1], 3)
        self.assertAlmostEqual(
            b.dwell_and_scans()[0] * b.dwell_and_scans()[1], 0.6)

    def test_a_fit_is_not_rebuilt_and_the_user_is_told(self):
        p = os.path.join(self.dir, "fit.nxs")
        write(p, [region()])
        with h5py.File(p, "r+") as h:
            g = h["C_1s"] if "C_1s" in h else next(iter(h.values()))
            fit = g.create_group("fit_C_1s")
            fit.attrs["NX_class"] = "NXfit"
        f = readers.load_file(p)
        self.assertEqual(len(f.regions), 1)
        self.assertTrue(any("NXfit" in w for w in f.warnings))


@unittest.skipUnless(HAVE_H5, "h5py not installed")
class TestForeignLayout(unittest.TestCase):
    """What a converter other than ours writes: different group names, no
    metadata note, kinetic energy with the photon energy in an NXbeam."""

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def build(self, path, axis_type="kinetic", n=11, signal_dims=1):
        ke = 1190.0 + np.arange(n) * 0.2
        with h5py.File(path, "w") as f:
            e = f.create_group("Cu_2p_scan")
            e.attrs["NX_class"] = "NXentry"
            e["definition"] = "NXmpes"
            e["title"] = "Cu 2p3/2 on Si"
            e["start_time"] = "2025-03-01T10:00:00+01:00"
            e["end_time"] = "2025-03-01T10:20:00+01:00"
            ins = e.create_group("instrument")
            ins.attrs["NX_class"] = "NXinstrument"
            ana = ins.create_group("analyser")
            ana.attrs["NX_class"] = "NXelectronanalyzer"
            col = ana.create_group("col")
            col.attrs["NX_class"] = "NXcollectioncolumn"
            col["lens_mode"] = "Transmission"
            ed = ana.create_group("disp")
            ed.attrs["NX_class"] = "NXenergydispersion"
            ed["pass_energy"] = 40.0
            ed["energy_scan_mode"] = "fixed_analyzer_transmission"
            beam = ins.create_group("xray")
            beam.attrs["NX_class"] = "NXbeam"
            beam["incident_energy"] = 1486.6
            src = ins.create_group("tube")
            src.attrs["NX_class"] = "NXsource"
            src["name"] = "Al Kα"
            src["power"] = 150.0
            smp = e.create_group("specimen")
            smp.attrs["NX_class"] = "NXsample"
            smp["name"] = "Cu foil"
            d = e.create_group("spectrum")
            d.attrs["NX_class"] = "NXdata"
            shape = (n,) if signal_dims == 1 else (n, 3)
            d["intensity"] = np.ones(shape) * 5.0
            d["energy"] = ke
            d["energy"].attrs["units"] = "eV"
            if axis_type:
                d["energy"].attrs["type"] = axis_type
            d.attrs["signal"] = "intensity"
            d.attrs["axes"] = "energy"
        return ke

    def test_names_are_found_by_class_not_by_group_name(self):
        p = os.path.join(self.dir, "a.nxs")
        ke = self.build(p)
        f = readers.load_file(p)
        r = f.regions[0]
        self.assertEqual((r.name, r.sample), ("Cu 2p3/2 on Si", "Cu foil"))
        self.assertEqual((r.photon_energy, r.pass_energy, r.lens_mode,
                          r.anode), (1486.6, 40.0, "Transmission", "Al Kα"))
        self.assertEqual(r.extra["analyser_mode"],
                         "Constant analyser energy (CAE)")
        self.assertEqual(r.conditions["X-ray Power"], "150 W")
        self.assertEqual(r.counts, [5.0] * 11)

    def test_a_kinetic_axis_with_a_photon_energy_becomes_binding(self):
        p = os.path.join(self.dir, "b.nxs")
        ke = self.build(p)
        r = readers.load_file(p).regions[0]
        self.assertEqual(r.energy_label, "Binding Energy")
        for k, b in zip(ke, r.energy):
            self.assertAlmostEqual(1486.6 - k, b)

    def test_a_binding_axis_is_kept(self):
        p = os.path.join(self.dir, "c.nxs")
        ke = self.build(p, axis_type="binding")
        r = readers.load_file(p).regions[0]
        self.assertEqual(r.energy_label, "Binding Energy")
        self.assertEqual(r.energy, [float(x) for x in ke])

    def test_an_axis_without_a_type_keeps_its_own_name_not_a_guess(self):
        p = os.path.join(self.dir, "d.nxs")
        self.build(p, axis_type=None)
        r = readers.load_file(p).regions[0]
        self.assertEqual(r.energy_label, "Energy")

    def test_times_with_an_offset_become_a_utc_run_window(self):
        p = os.path.join(self.dir, "e.nxs")
        self.build(p)
        r = readers.load_file(p).regions[0]
        self.assertEqual((r.extra["t_start"], r.extra["t_end"], r.extra["tz"]),
                         ("2025-03-01 09:00:00", "2025-03-01 09:20:00", "UTC"))
        self.assertEqual(r.date, "2025-03-01 10:00:00")

    def test_a_two_dimensional_signal_is_left_out_with_a_reason(self):
        p = os.path.join(self.dir, "f.nxs")
        self.build(p, signal_dims=2)
        with self.assertRaises(ValueError):
            readers.load_file(p)               # nothing else in the file
        h = h5py.File(p, "r")
        try:
            f = nexus_xps.NexusFile()
            f.path = p
            entry = h["Cu_2p_scan"]
            f._add_other("Cu_2p_scan", entry, "stem")
        finally:
            h.close()
        self.assertTrue(any("2 dimensions" in w for w in f.warnings))
        self.assertEqual(f.regions, [])

    def test_mixed_with_a_nexafs_scan_both_are_read(self):
        p = os.path.join(self.dir, "mixed.nxs")
        self.build(p)
        q = os.path.join(self.dir, "scan.nxs")
        write_scan(q)
        with h5py.File(p, "r+") as dst, h5py.File(q, "r") as src:
            src.copy("entry1", dst)
        f = readers.load_file(p)
        self.assertEqual(f.format_name, "NeXus (.nxs)")
        kinds = sorted({r.technique for r in f.regions})
        self.assertEqual(kinds, ["NEXAFS", "XPS"])

    def test_an_entry_that_is_neither_is_noted(self):
        p = os.path.join(self.dir, "extra.nxs")
        self.build(p)
        with h5py.File(p, "r+") as h:
            g = h.create_group("misc")
            g.attrs["NX_class"] = "NXentry"
            g["title"] = "calibration notes"
        f = readers.load_file(p)
        self.assertEqual(len(f.regions), 1)
        self.assertTrue(any("left out" in w for w in f.warnings))

    def test_a_file_with_nothing_to_read_says_what_it_looks_for(self):
        p = os.path.join(self.dir, "none.nxs")
        with h5py.File(p, "w") as h:
            g = h.create_group("entry")
            g.attrs["NX_class"] = "NXentry"
            g["title"] = "nothing"
        with self.assertRaises(ValueError) as cm:
            readers.load_file(p)
        self.assertIn("NXxps spectrum or photon-energy scan", str(cm.exception))


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
        cls.path = os.path.join(cls.dir, "spectra.nxs")
        write(cls.path, [region("C 1s"), region("O 1s")])

    @classmethod
    def tearDownClass(cls):
        (ee.messagebox.showinfo, ee.messagebox.showwarning,
         ee.messagebox.showerror, ee.messagebox.askyesno) = cls._boxes
        cls.root.destroy()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_it_opens_and_draws_as_binding_energy_spectra(self):
        ws = ee.Workspace(self.root)
        ws._add_files([self.path])
        regions = ws.docs[0].regions
        self.assertEqual([r.name for r in regions], ["C 1s", "O 1s"])
        ws.checked = {id(r) for r in regions}
        ws._render()
        self.root.update_idletasks()
        self.assertTrue(any("inding" in ax.get_xlabel()
                            for ax in ws.fig.axes))
        self.assertTrue(any("NeXus" in n and "*.nxs" in p
                            for n, p in readers.supported_patterns()))


@unittest.skipUnless(HAVE_H5 and os.environ.get("XPS_NXXPS_CORPUS"),
                     "set XPS_NXXPS_CORPUS to a folder of spectra files")
class TestRealFilesRoundTrip(unittest.TestCase):
    """Every spectrum file of a folder, exported to NeXus and read back."""

    def test_every_file_comes_back_the_same(self):
        root = os.environ["XPS_NXXPS_CORPUS"]
        files = [p for p in glob.glob(os.path.join(root, "**", "*"),
                                      recursive=True)
                 if os.path.isfile(p) and os.path.splitext(p)[1].lower() in
                 (".vms", ".vamas", ".avg", ".vgd", ".kal", ".dset",
                  ".experiment", ".spe", ".kfit")]
        self.assertTrue(files)
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        done = 0
        for src in files[:60]:
            try:
                f = readers.load_file(src)
            except Exception:                         # noqa: BLE001
                continue
            regs = [r for r in f.regions if r.decodable and r.counts
                    and not r.is_nexafs]
            if not regs:
                continue
            mds = [f.region_metadata(r) for r in regs]
            p = os.path.join(d, "x.nxs")
            nexus_export.export_nexus(regs, p, metadata=mds,
                                      instrument=f.instrument)
            g = readers.load_file(p)
            self.assertEqual(len(g.regions), len(regs), src)
            for a, b in zip(regs, g.regions):
                self.assertEqual(a.name, b.name, src)
                self.assertEqual(a.sample, b.sample, src)
                np.testing.assert_allclose(b.energy, a.energy, rtol=1e-12)
                np.testing.assert_allclose(b.counts, a.counts, rtol=1e-12)
                self.assertEqual(a.photon_energy, b.photon_energy, src)
                self.assertEqual(a.pass_energy, b.pass_energy, src)
                self.assertEqual(a.technique, b.technique, src)
                ta, tb = a.transmission(), b.transmission()
                self.assertEqual(ta is None, tb is None, src)
                if ta:
                    np.testing.assert_allclose(tb, ta, rtol=1e-9)
            done += 1
        self.assertGreater(done, 0)


if __name__ == "__main__":
    unittest.main()
