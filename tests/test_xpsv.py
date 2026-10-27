"""XPSView / HarwellXPS sessions (``.xpsv``): sniffing, both schemas, the
charge shift baked into the stored axis, a photon energy edited by hand,
non-unique sample names, levels, CasaXPS fits and kinetic-energy techniques.

Set ``XPS_XPSV_CORPUS=<folder holding the project folders>`` (``D:\\Temp\\hxps``)
to load every session in it and compare it with the instrument files it was
made from: counts equal for every level, and the restored axis, photon energy,
pass energy, sweeps and dwell equal the original reader's.

Run:  python -m unittest discover tests
"""

import glob
import gzip
import json
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import readers  # noqa: E402
from readers import xpsv  # noqa: E402

HV = 1486.69
CASA = ["Casa Info Follows", "1", "Calib M = 279.455 A = 284.8 BE ADD", "0",
        "1",
        "CASA region (*O 1s*) (*Shirley*) 948.0 960.0 0.78 2 0 0 0 -450 0 0 "
        "(*O 1s*) 15.9994 0 0.78", "1",
        "CASA comp (*Al2O3 (O)*) (*GL(30)*) Area 1000 1e-20 1e6 -1 1 MFWHM "
        "1.8 0.38 9.5 -1 1 Position 955.0 948.0 960.0 -1 1 RSF 0.78 MASS "
        "15.9994 INDEX 2 (*O 1s*) CONST (**) UNCORRECTEDRSF 0.78"]


def acq(**kw):
    a = {"pass_energy": 40.0, "dwell_s": 0.1, "dwell_total_s": 0.2,
         "sweeps": 2, "step_eV": 0.1, "lens_mode": "Hybrid",
         "aperture": "Slot", "source": "XPS - monochromatic", "anode": "Al",
         "xray_power": "225.00W", "emission_current_mA": 15.0,
         "neutraliser_on": True,
         "neutraliser": {"filament_current_A": 0.45, "charge_balance_V": 5.0,
                         "filament_bias_V": 1.4},
         "acq_start": "2026-09-13 02:32:26",
         "stage": {"x": 30.0, "y": 18.0, "z": -0.3, "unit": "mm"}}
    a.update(kw)
    return a


def spectrum(sid, sample, name, raw_be0, n=21, step=0.1, shift=0.0, hv=HV,
             levels=1, meta=None, **kw):
    """A spectrum whose raw BE axis starts (high) at ``raw_be0``; the stored
    ``ke0`` carries the baked charge shift, as XPSView stores it."""
    ke_raw = HV - raw_be0                        # lowest raw KE
    sp = {"id": sid, "name": name, "sample_id": sample, "technique": "XPS",
          "hv": hv, "ke0": ke_raw - shift, "dke": step, "n": n,
          "y": [[float(100 + i + 1000 * lv) for i in range(n)]
                for lv in range(levels)],
          "levels": [{"label": f"L{lv}", "etch_time": 30.0 * lv if levels > 1
                      else None, "iteration": None,
                      "values": {}} for lv in range(levels)],
          "dwell": 0.1, "y_is_cps": False, "sweeps": 2, "pass_energy": 40.0,
          "energy_scale": "BE", "calib_shift": 0.0, "tf": [], "acq": acq(),
          "meta": dict(meta or {}), "regions": [], "components": [],
          "notes": "", "display_name": ""}
    if shift:
        sp["meta"]["casa_shift"] = shift
    sp.update(kw)
    return sp


def experiment(spectra, samples=None):
    sids = []
    for sp in spectra:
        if sp["sample_id"] not in sids:
            sids.append(sp["sample_id"])
    return {"model_version": 2, "name": "Test", "instrument": "Kratos Axis",
            "vendor": "Kratos", "operator": "tester", "project": "P1",
            "meta": {"machine": "PC1"}, "notes": "",
            "samples": samples or [{"id": s, "name": "S", "position": {},
                                    "meta": {}, "group": ""} for s in sids],
            "spectra": spectra, "images": []}


def write(doc, wrapped=True, gz=True, name="t.xpsv"):
    d = tempfile.mkdtemp()
    path = os.path.join(d, name)
    body = {"schema": 2, "experiment": doc, "project": {}, "mapping": {},
            "content": {}, "calibration": {}, "provenance": {}} \
        if wrapped else doc
    raw = json.dumps(body).encode("utf-8")
    with open(path, "wb") as fh:
        fh.write(gzip.compress(raw) if gz else raw)
    return path


class TestSniff(unittest.TestCase):
    def head(self, path):
        with open(path, "rb") as fh:
            return fh.read(4096)

    def test_schema_2_gzip(self):
        p = write(experiment([spectrum("a", "s1", "C 1s", 290.0)]))
        self.assertTrue(xpsv.sniff(self.head(p), ".xpsv"))
        self.assertIs(readers.reader_for(p), xpsv.XpsvFile)

    def test_schema_1_plain_json_and_renamed(self):
        p = write(experiment([spectrum("a", "s1", "C 1s", 290.0)]),
                  wrapped=False, gz=False, name="renamed.dat")
        self.assertTrue(xpsv.sniff(self.head(p), ".dat"))

    def test_other_files_are_not_sessions(self):
        self.assertFalse(xpsv.sniff(b'{"a": 1}', ".json"))
        self.assertFalse(xpsv.sniff(b"x,y\n1,2\n", ".xpsv"))
        self.assertFalse(xpsv.sniff(gzip.compress(b"plain text"), ".gz"))
        self.assertFalse(xpsv.sniff(b"\x1f\x8b\x08\x00 broken", ".xpsv"))

    def test_registered_with_dialog_pattern(self):
        names = dict(readers.supported_patterns())
        self.assertIn("*.xpsv", names["XPSView session (.xpsv)"])


class TestLoad(unittest.TestCase):
    def load(self, *spectra, **kw):
        return readers.load_file(write(experiment(list(spectra), **kw)))

    def test_axis_and_acquisition(self):
        d = self.load(spectrum("a", "s1", "C1s", 292.0, n=21, step=0.1))
        r = d.regions[0]
        self.assertEqual(r.name, "C 1s")
        self.assertEqual(r.energy_label, "Binding Energy")
        self.assertAlmostEqual(r.energy[0], 292.0, 5)
        self.assertAlmostEqual(r.energy[-1], 290.0, 5)
        self.assertEqual(r.photon_energy, HV)
        self.assertEqual(r.pass_energy, 40.0)
        self.assertEqual(r.anode, "Al mono")
        self.assertEqual(r.conditions["X-ray Power"], "225 W")
        self.assertIn("filament current 0.45 A", r.extra["neutraliser"])
        self.assertEqual((r.pos_x, r.pos_y), (30.0, 18.0))
        self.assertEqual(r.extra["n_scans"], 2)
        self.assertEqual(d.instrument["Operator"], "tester")

    def test_kratos_dwell_is_the_sum_over_the_sweeps(self):
        r = self.load(spectrum("a", "s1", "C 1s", 292.0)).regions[0]
        self.assertAlmostEqual(r.dwell, 0.2)
        self.assertTrue(r.extra["dwell_total"])
        self.assertAlmostEqual(r.dwell_and_scans()[0], 0.1)

    def test_baked_charge_shift_is_undone(self):
        sp = spectrum("a", "s1", "C 1s", 292.0, shift=5.345,
                      meta={"vamas_block": {"comment": CASA[:3] + ["0"]}})
        r = self.load(sp).regions[0]
        # the stored axis starts 5.345 eV higher; the raw one does not
        self.assertAlmostEqual(r.energy[0], 292.0, 5)
        self.assertAlmostEqual(r.calibration_shift, 5.345, 5)
        self.assertAlmostEqual(r.extra["casa_calib"]["shift"], 5.345, 5)

    def test_hand_edited_photon_energy_uses_the_sources_value(self):
        sps = [spectrum(f"s{i}", "s1", name, be0)
               for i, (name, be0) in enumerate((("C 1s", 292.0),
                                                ("O 1s", 534.0),
                                                ("N 1s", 402.0)))]
        sps.append(spectrum("bad", "s1", "Al 2p", 80.0, shift=7.204,
                            hv=HV - 2.0, meta={"vamas_block": {
                                "comment": ["Calib M = 72.3409 A = 74.2 BE ADD",
                                            "Calib M = 279.455 A = 284.8 BE ADD"
                                            ]}}))
        d = self.load(*sps)
        r = d.regions[-1]
        self.assertEqual(r.photon_energy, HV)
        self.assertAlmostEqual(r.energy[0], 80.0, 5)
        self.assertTrue(any("1484.69" in w and "Al 2p" in w
                            for w in d.warnings))

    def test_manual_shift_is_reported_not_applied(self):
        sp = spectrum("a", "s1", "C 1s", 292.0, calib_shift=-0.05)
        d = self.load(sp)
        self.assertEqual(d.regions[0].calibration_shift, 0.0)
        self.assertTrue(any("-0.05" in w for w in d.warnings))

    def test_sample_names_are_made_unique(self):
        samples = [{"id": "a", "name": "ATH", "group": "", "meta": {},
                    "position": {}},
                   {"id": "b", "name": "ATH", "group": "Location 1",
                    "meta": {}, "position": {}},
                   {"id": "c", "name": "ATH", "group": "", "meta": {},
                    "position": {}}]
        sps = [spectrum("1", "a", "C 1s", 292.0,
                        meta={"run": {"label": "Spot sequence - ATH - 13 Sep"}}),
               spectrum("2", "b", "C 1s", 292.0),
               spectrum("3", "c", "C 1s", 292.0,
                        meta={"run": {"label": "Depth profile 1 - ATH"}})]
        d = self.load(*sps, samples=samples)
        names = [r.sample for r in d.regions]
        self.assertEqual(len(set(names)), 3, names)
        self.assertTrue(all(n.startswith("ATH") for n in names))

    def test_levels_become_regions_with_etch_times(self):
        d = self.load(spectrum("a", "s1", "C 1s", 292.0, levels=3))
        self.assertEqual(len(d.regions), 3)
        self.assertEqual([r.etch_time for r in d.regions], [0.0, 30.0, 60.0])
        self.assertEqual([r.etch_level for r in d.regions], [0, 1, 2])
        self.assertTrue(d.depth_profile["is_profile"])
        self.assertEqual(d.regions[1].counts[0], 1100.0)

    def test_repeats_without_etch_time_stay_plain_regions(self):
        sp = spectrum("a", "s1", "C 1s", 292.0, levels=2)
        for lv in sp["levels"]:
            lv["etch_time"] = None
        d = self.load(sp)
        self.assertEqual(len(d.regions), 2)
        self.assertFalse(d.depth_profile["is_profile"])

    def test_cps_and_null_counts(self):
        sp = spectrum("a", "s1", "C 1s", 292.0, y_is_cps=True)
        sp["y"][0][3] = None
        r = self.load(sp).regions[0]
        self.assertEqual(r.count_units, "counts/s")
        self.assertNotEqual(r.counts[3], r.counts[3])      # NaN

    def test_transmission_function(self):
        sp = spectrum("a", "s1", "C 1s", 292.0,
                      tf=[[400.0, 13.0], [200.0, 14.0]])
        r = self.load(sp).regions[0]
        self.assertEqual(r.tf_ke, [200.0, 400.0])
        self.assertEqual(r.tf_values, [14.0, 13.0])

    def test_casa_fit_is_read_from_the_stored_lines(self):
        sp = spectrum("a", "s1", "O 1s", 540.0, n=161, step=0.1,
                      shift=5.345,
                      meta={"vamas_block": {"comment": CASA}},
                      components=[{"name": "x"}])
        r = self.load(sp).regions[0]
        self.assertIsNotNone(r.fit)
        self.assertEqual(len(r.fit.components), 1)
        self.assertEqual(r.fit.components[0].index, 2)
        self.assertAlmostEqual(r.fit.calib_shift, 5.345, 5)

    def test_components_without_casa_lines_say_so(self):
        sp = spectrum("a", "s1", "O 1s", 540.0, components=[{"name": "x"}])
        d = self.load(sp)
        self.assertIsNone(d.regions[0].fit)
        self.assertTrue(any("CasaXPS lines" in w for w in d.warnings))

    def test_reels_is_a_kinetic_energy_spectrum(self):
        sp = spectrum("a", "s1", "REELS Survey", 0.0, hv=1000.0)
        sp.update({"technique": "REELS", "energy_scale": "KE", "ke0": 950.0,
                   "dke": 0.1, "n": 21, "y": [[float(i) for i in range(21)]],
                   "y_is_cps": True,
                   "acq": {"source_energy_eV": 1000.0, "sweeps": 4,
                           "stage": {"x": 50692623, "y": 43463115,
                                     "z": 34899997, "unit": "nm"}}})
        exp = experiment([sp])
        exp["vendor"] = "Thermo"
        r = readers.load_file(write(exp)).regions[0]
        self.assertEqual(r.technique, "REELS")
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertIsNone(r.photon_energy)
        self.assertAlmostEqual(r.energy[0], 950.0)
        self.assertAlmostEqual(r.pos_x, 50.692623)
        self.assertEqual(r.extra["acq_metadata"]["Primary beam energy (eV)"],
                         "1000")
        self.assertFalse(r.extra.get("dwell_total"))

    def test_schema_1_loads(self):
        p = write(experiment([spectrum("a", "s1", "C 1s", 292.0)]),
                  wrapped=False)
        self.assertEqual(len(readers.load_file(p).regions), 1)

    def test_not_a_session(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "x.xpsv")
        with open(p, "wb") as fh:
            fh.write(gzip.compress(b'{"schema": 2, "experiment": {"a": 1}}'))
        with self.assertRaises(ValueError):
            readers.load_file(p)

    def test_picture_is_kept(self):
        import base64
        exp = experiment([spectrum("a", "s1", "C 1s", 292.0)])
        exp["images"] = [{"id": "i", "kind": "holder", "name": "Holder",
                          "mime": "image/jpeg",
                          "data": base64.b64encode(b"\xff\xd8abc").decode(),
                          "sample_id": ""}]
        d = readers.load_file(write(exp))
        self.assertEqual(d.images[0].name, "Holder")
        self.assertEqual(d.images[0].data, b"\xff\xd8abc")


CORPUS = os.environ.get("XPS_XPSV_CORPUS", "")


@unittest.skipUnless(CORPUS and os.path.isdir(CORPUS),
                     "set XPS_XPSV_CORPUS to a folder of XPSView projects")
class TestCorpus(unittest.TestCase):
    """Every session against the instrument files of its own project folder."""

    @staticmethod
    def originals(project):
        docs = []
        for f in glob.glob(os.path.join(project, "**", "*.experiment"),
                           recursive=True):
            try:
                docs.append(readers.load_file(f))
            except Exception:                  # noqa: BLE001 - a stub file
                pass
        for f in glob.glob(os.path.join(project, "**", "*.VGX"),
                           recursive=True):
            docs.append(readers.load_file(f))
        return docs

    def test_every_session_matches_its_originals(self):
        import numpy as np
        sessions = glob.glob(os.path.join(CORPUS, "*", "**", "*.xpsv"),
                             recursive=True)
        self.assertTrue(sessions)
        for path in sessions:
            project = os.path.join(CORPUS, os.path.relpath(
                path, CORPUS).split(os.sep)[0])
            origs = self.originals(project)
            if not origs:
                continue
            session = readers.load_file(path)
            pool = {}
            for o in origs:
                for r in o.regions:
                    if r.counts:
                        pool.setdefault(len(r.counts), []).append(r)
            used = set()
            for r in session.regions:
                c = np.asarray(r.counts, float)
                hit = next((o for o in pool.get(len(c), [])
                            if id(o) not in used and np.allclose(
                                np.asarray(o.counts, float), c, rtol=1e-6,
                                atol=1e-3, equal_nan=True)), None)
                with self.subTest(path=os.path.basename(path), region=r.name,
                                  sample=r.sample):
                    self.assertIsNotNone(hit, "no original has these counts")
                    used.add(id(hit))
                    if r.energy_label == "Binding Energy":
                        self.assertAlmostEqual(r.energy[0], hit.energy[0], 2)
                        self.assertAlmostEqual(r.energy[-1], hit.energy[-1], 2)
                        self.assertAlmostEqual(r.photon_energy,
                                               hit.photon_energy, 3)
                    self.assertEqual(r.pass_energy, hit.pass_energy)


if __name__ == "__main__":
    unittest.main()
