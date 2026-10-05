"""Kratos stigmatic imaging maps: the reader (``.kal`` and ``.dset``), the
cube geometry, and the pure helpers the series viewer uses (``kratosmap``).
Synthetic files only; no instrument data needed.

Run:  python -m unittest tests.test_kratosmap
"""

import math
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import kratosmap  # noqa: E402
from readers import load_file  # noqa: E402
from test_kratos_dset import (block, dset, r_dbl, r_f32, r_int,  # noqa: E402
                              r_str)

NX = NY = 4
FSD_MM = 0.236
STEP = 2 / (NX - 1)                          # -1 .. +1 in NX points


def pixels_for(k):
    """A 4 x 4 image whose values identify the map ``k``."""
    return [float(k * 100 + i) for i in range(NX * NY)]


def kal_map(name, k, ke=1402.69, z_m=0.000775, x_m=0.039638, y_m=0.000088,
            when="17/09/08 08:48:12", ref="F_REFER_TO_XRAY_MONO_AL"):
    vals = ", ".join(f"{v:g}" for v in pixels_for(k))
    el, tr = name.split()
    return "\n".join([
        f"Object name               = {name}/{k}",
        "   1 Technique                 = F_XPS",
        "   2 Scan type                 = F_MAPPING",
        "   7 Dwell time                = 120 seconds",
        f"  12 Ordinate values           = {{{vals}}}",
        f"  25 # points per line in map  = {NX}",
        f"  26 # lines in map            = {NY}",
        f"  37 Acquisition name          = {name}",
        "  42 Pass energy               = 160 eV",
        f"  51 Map energy/mass           = {ke} eV",
        "  99 # Sweeps completed        = 1",
        f" 151 Date Acquired             = {when}",
        f"3003 step size x coord         = {STEP:.6g}",
        f"3004 step size y coord         = {STEP:.6g}",
        f"3005 Full Scale Deflection X   = {FSD_MM} mm",
        f"3006 Full Scale Deflection Y   = {FSD_MM} mm",
        f"3080 Xray Reference Energy     = {ref}",
        f"3088 Stage X Position          = {x_m} m",
        "3089 Stage Y Position          = 8.8e-05 m",
        f"3090 Stage Z Position          = {z_m} m",
        f"3113 Chemical symbol or formula = {el}",
        f"3114 Transition or charge state = {tr}",
        ""])


def kal_position(name, k):
    return "\n".join([
        f"Object name               = Sample Position/{k}",
        "  37 Acquisition name        = Sample Position",
        "  38 State change type       = F_POSITION",
        f"3087 Stage Position Name     = {name}",
        "3088 Stage X Position        = 0.039638 m", ""])


def kal_text(*objects):
    return "Dataset filename          = t.dset\n" + "\n".join(objects)


class Tmp(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, name, data):
        p = os.path.join(self.dir.name, name)
        with open(p, "wb" if isinstance(data, bytes) else "w") as fh:
            fh.write(data)
        return p


class TestKalMaps(Tmp):
    def load(self, *objects):
        return load_file(self.write("m.kal", kal_text(*objects)))

    def test_a_map_is_a_region_with_a_one_channel_cube_and_no_spectrum(self):
        f = self.load(kal_position("Grid on Tape", 1), kal_map("Au 4f", 2))
        (r,) = f.regions
        cube = r.extra["cube"]
        self.assertFalse(r.decodable)
        self.assertIsNone(r.counts)
        self.assertEqual((cube.nx, cube.ny, cube.n_energy), (NX, NY, 1))
        self.assertEqual(list(cube.array3d()[:, :, 0].ravel()), pixels_for(2))
        self.assertEqual(r.name, "Au 4f")
        self.assertEqual(r.technique, "XPS imaging")
        self.assertEqual(r.sample, "Grid on Tape")        # the preceding position
        self.assertEqual(r.extra["acq_mode"], "Stigmatic map")
        self.assertEqual(r.extra["position_name"], "Grid on Tape")
        self.assertEqual(r.extra["stage_z_um"], 775.0)
        self.assertEqual((r.pos_x, r.pos_y), (39.638, 0.088))
        self.assertEqual((r.pass_energy, r.dwell), (160.0, 120.0))
        self.assertEqual(r.date, "2017-09-08 08:48:12")

    def test_energy_is_kinetic_and_the_cube_channel_is_binding_energy(self):
        r = self.load(kal_map("Au 4f", 1)).regions[0]
        self.assertEqual(r.extra["map_ke"], 1402.69)
        self.assertAlmostEqual(r.extra["cube"].energy[0], 1486.6 - 1402.69, 6)
        self.assertEqual(r.energy_label, "Binding Energy")
        self.assertAlmostEqual(r.photon_energy, 1486.6)

    def test_geometry_is_centred_and_uses_the_measured_scale(self):
        cube = self.load(kal_map("Au 4f", 1)).regions[0].extra["cube"]
        want = float(f"{float(f'{STEP:.6g}') * FSD_MM * 1000 * kratosmap.MEASURED_SCALE:.6g}")
        self.assertAlmostEqual(cube.dx, want, 6)
        self.assertEqual(cube.dy, cube.dx)
        self.assertAlmostEqual(cube.x_of(0), -cube.x_of(NX - 1), 9)   # centred
        self.assertAlmostEqual(cube.y_of(0), -cube.y_of(NY - 1), 9)
        self.assertEqual((cube.stage_x_mm, cube.stage_y_mm), (39.638, 0.088))

    def test_the_region_says_the_scale_is_approximate(self):
        r = self.load(kal_map("Au 4f", 1)).regions[0]
        self.assertIn("approximate", r.note)
        self.assertIn("not applied", r.note)

    def test_no_position_object_gives_one_shared_sample(self):
        f = self.load(kal_map("Cu 2p", 1), kal_map("Au 4f", 2))
        self.assertEqual({r.sample for r in f.regions}, {"Maps"})

    def test_reference_none_or_missing_gives_kinetic_energy(self):
        f = self.load(kal_map("Au 4f", 1, ref="F_REFER_TO_NONE"))
        r = f.regions[0]
        self.assertEqual(r.energy_label, "Kinetic Energy")
        self.assertEqual(r.extra["cube"].energy, [1402.69])
        self.assertEqual(f.warnings, [])
        text = kal_map("Au 4f", 1).replace(
            "3080 Xray Reference Energy     = F_REFER_TO_XRAY_MONO_AL\n", "")
        g = self.load(text, kal_map("Cu 2p", 2).replace(
            "3080 Xray Reference Energy     = F_REFER_TO_XRAY_MONO_AL\n", ""))
        self.assertEqual(len(g.warnings), 1)                  # one for the file
        self.assertIn("X-ray energy unknown for 2 maps", g.warnings[0])

    def test_a_map_with_the_wrong_pixel_count_is_skipped_with_one_warning(self):
        bad = kal_map("Au 4f", 1).replace(f"{NX}\n", "5\n", 1)   # 5 x 4 != 16
        f = self.load(bad, kal_map("Cu 2p", 2))
        self.assertEqual([r.name for r in f.regions], ["Cu 2p"])
        self.assertEqual(len(f.warnings), 1)
        self.assertIn("not read", f.warnings[0])

    def test_tree_row_says_map_and_shows_the_energy_and_size(self):
        f = self.load(kal_position("P", 1), kal_map("Au 4f", 2))
        row = f.tree.children[0]
        self.assertIn("(Map)", row.label)
        self.assertEqual(row.cols[0], "83.9 eV")
        self.assertEqual(row.cols[1], f"{NX}x{NY}")

    def test_metadata_counts_the_image_time_and_names_the_mode(self):
        f = self.load(kal_map("Au 4f", 1))
        md = f.region_metadata(f.regions[0])
        self.assertEqual(md["Technique"], "XPS imaging")
        self.assertEqual(md["Acquisition mode"], "Stigmatic map")
        self.assertEqual(md["Counting time"], "2 min")
        self.assertEqual(md["Pass energy (eV)"], "160")

    def test_methods_text_keeps_images_apart_from_spectra(self):
        import methods
        f = self.load(kal_map("Au 4f", 1), kal_map("Cu 2p", 2, ke=554.69))
        text = methods.generate([f.region_metadata(r) for r in f.regions])
        self.assertIn("Stigmatic images (2: Au 4f, Cu 2p) were recorded", text)
        self.assertIn("an acquisition time of 120 s per image", text)
        self.assertIn("comprises 2 images", text)
        self.assertNotIn("High-resolution spectra", text)


class TestDsetMaps(Tmp):
    @staticmethod
    def dset_map(k, name="Au 4f", ordinal=1):
        el, tr = name.split()
        return block(
            ordinal, r_int(1, 3), r_int(2, 2), r_dbl(7, 120.0),
            r_f32(12, pixels_for(k)), r_int(25, NX), r_int(26, NY),
            r_str(37, name), r_dbl(42, 160.0), r_dbl(51, 1402.69), r_int(99, 1),
            r_str(151, "17/09/08 08:48:12"), r_dbl(3003, STEP),
            r_dbl(3004, STEP), r_dbl(3005, FSD_MM), r_dbl(3006, FSD_MM),
            r_int(3080, 8), r_dbl(3088, 0.039638), r_dbl(3089, 8.8e-05),
            r_dbl(3090, 0.000775), r_str(3113, el), r_str(3114, tr))

    def test_dset_map_equals_the_kal_map(self):
        # the .dset holds full doubles, the .kal six digits: the cubes agree
        a = load_file(self.write("a.dset", dset(self.dset_map(2))))
        b = load_file(self.write("b.kal", kal_text(kal_map("Au 4f", 2))))
        ra, rb = a.regions[0], b.regions[0]
        self.assertEqual(ra.extra["cube"], rb.extra["cube"])
        for key in ("map_ke", "stage_z_um", "acq_mode"):
            self.assertEqual(ra.extra[key], rb.extra[key])
        self.assertEqual((ra.name, ra.date, ra.pass_energy, ra.dwell),
                         (rb.name, rb.date, rb.pass_energy, rb.dwell))
        self.assertEqual(a.warnings, [])           # F_MAPPING is a known enum


class TestHelpers(unittest.TestCase):
    @staticmethod
    def frames_of(specs):
        """Frames from ``(name, ke, x_mm, y_mm, z_um, time)`` tuples, built
        through the reader so the cubes are real."""
        objs = []
        for k, (name, ke, x, z, when) in enumerate(specs, 1):
            objs.append(kal_map(name, k, ke=ke, x_m=x / 1000, z_m=z / 1e6,
                                when=when))
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "s.kal")
            with open(p, "w") as fh:
                fh.write(kal_text(*objs))
            f = load_file(p)
        return kratosmap.frames(f.regions)

    def test_focus_series_opens_on_the_position_and_plots_against_z(self):
        fr = self.frames_of([
            ("Au 4f", 1402.69, 10.0, 900, "17/09/08 09:00:00"),
            ("Au 4f", 1402.69, 10.0, 930, "17/09/08 09:02:00"),
            ("Au 4f", 1402.69, 10.0, 960, "17/09/08 09:04:00"),
            ("Au 4f", 1402.69, 12.0, 960, "17/09/08 09:10:00")])
        self.assertEqual(kratosmap.default_filter(fr, fr[0]), "position")
        same = kratosmap.select(fr, fr[0], "position")
        self.assertEqual(len(same), 3)
        label, xs = kratosmap.series_axis(same)
        self.assertEqual(label, "Stage Z (um)")
        self.assertEqual(xs, [900.0, 930.0, 960.0])

    def test_time_series_plots_minutes_and_multi_element_opens_on_all(self):
        t = self.frames_of([
            ("Cu 2p", 537.74, 10.0, 1000, "17/09/11 12:26:00"),
            ("Cu 2p", 537.74, 10.0, 1000, "17/09/11 12:46:00")])
        self.assertEqual(kratosmap.default_filter(t, t[0]), "energy")
        label, xs = kratosmap.series_axis(t)
        self.assertEqual((label, xs), ("Minutes from first", [0.0, 20.0]))
        m = self.frames_of([("Cu 2p", 554.69, 10, 700, "17/09/11 09:24:36"),
                            ("O 1s", 956.69, 10, 700, "17/09/11 09:24:36")])
        self.assertEqual(kratosmap.default_filter(m, m[0]), "all")
        self.assertEqual(kratosmap.series_axis(m)[0], "Frame")
        self.assertEqual(m[1].be, round(1486.6 - 956.69, 4))
        self.assertIn("BE 529.91 eV", m[1].label)

    def test_frames_are_in_acquisition_order_and_ignore_spectra(self):
        fr = self.frames_of([("Au 4f", 1402.69, 1, 1, "17/09/08 09:00:00"),
                             ("Cu 2p", 554.69, 1, 1, "17/09/08 09:01:00")])
        self.assertEqual([f.index for f in fr], [0, 1])
        self.assertEqual([f.region.name for f in fr], ["Au 4f", "Cu 2p"])

    def test_blur_keeps_the_total_and_does_nothing_at_zero(self):
        import numpy as np
        img = np.zeros((20, 20))
        img[10, 10] = 100.0
        self.assertTrue(np.array_equal(kratosmap.blur(img, 0), img))
        b = kratosmap.blur(img, 1.5)
        self.assertAlmostEqual(b.sum(), 100.0, 6)
        self.assertLess(b.max(), 100.0)
        self.assertEqual(np.unravel_index(b.argmax(), b.shape), (10, 10))

    def test_focus_metric_is_higher_for_the_sharper_image(self):
        import numpy as np
        yy, xx = np.mgrid[0:64, 0:64]
        def ring(w):
            return 5 + 50 * np.exp(-((np.hypot(yy - 32, xx - 32) - 15) / w) ** 2)
        self.assertGreater(kratosmap.focus_metric(ring(1.5)),
                           kratosmap.focus_metric(ring(6.0)))
        self.assertEqual(kratosmap.focus_metric(np.zeros((8, 8))), 0.0)
        # brightness alone does not change it
        self.assertAlmostEqual(kratosmap.focus_metric(ring(3.0)),
                               kratosmap.focus_metric(10 * ring(3.0)), 9)

    def test_roi_means_use_the_mask_and_give_nan_for_a_wrong_one(self):
        import numpy as np
        fr = self.frames_of([("Au 4f", 1402.69, 1, 1, "17/09/08 09:00:00")])
        whole = kratosmap.roi_means(fr)
        self.assertAlmostEqual(whole[0], sum(pixels_for(1)) / 16)
        m = np.zeros((NY, NX), bool)
        m[0, :] = True
        self.assertAlmostEqual(kratosmap.roi_means(fr, m)[0],
                               sum(pixels_for(1)[:NX]) / NX)
        self.assertTrue(math.isnan(kratosmap.roi_means(fr, np.zeros((2, 2), bool))[0]))

    def test_scale_constant_matches_the_two_registered_pairs(self):
        # 1.728 and 1.730 um/pixel were measured; nominal is step * FSD
        nominal = (2 / 255) * 236.0
        for measured in (1.728, 1.730):
            self.assertAlmostEqual(nominal * kratosmap.MEASURED_SCALE,
                                   measured, delta=0.005)


if __name__ == "__main__":
    unittest.main()
