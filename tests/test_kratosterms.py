import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import kratosterms  # noqa: E402
from kratosterms import friendly  # noqa: E402


class TestFriendly(unittest.TestCase):
    def test_every_lens_constant_with_or_without_prefix(self):
        for key, text in kratosterms.LENS_MODES.items():
            self.assertEqual(friendly(key), text)
            self.assertEqual(friendly("F_" + key), text)
            self.assertEqual(friendly(key.lower()), text)

    def test_fields_of_view_are_named(self):
        self.assertEqual(friendly("F_MHSA_LOW_MAGN"),
                         "Low magnification (FoV1)")
        self.assertEqual(friendly("MHSA_MEDIUM_MAGN"),
                         "Medium magnification (FoV2)")
        self.assertEqual(friendly("F_MHSA_HIGH_MAGN"),
                         "High magnification (FoV3)")

    def test_iss_is_recognised_whatever_the_wording_is(self):
        from unittest import mock
        for s in ("F_HSA_LENS_ISS", "HSA_LENS_ISS", "ISS",
                  kratosterms.LENS_MODES["HSA_LENS_ISS"]):
            self.assertTrue(kratosterms.is_iss_lens(s), s)
        for s in ("", None, "Hybrid", "HSA_LENS_HYBRID",
                  kratosterms.LENS_MODES["HSA_LENS_HYBRID"]):
            self.assertFalse(kratosterms.is_iss_lens(s), s)
        words = dict(kratosterms.LENS_MODES, HSA_LENS_ISS="Ion scattering")
        with mock.patch.object(kratosterms, "LENS_MODES", words):
            self.assertEqual(friendly("F_HSA_LENS_ISS"), "Ion scattering")
            self.assertTrue(kratosterms.is_iss_lens("Ion scattering"))
            self.assertTrue(kratosterms.is_iss_lens("HSA_LENS_ISS"))

    def test_unknown_constant_is_tidied_not_blank(self):
        self.assertEqual(friendly("F_HSA_LENS_NEW_THING"), "New thing")

    def test_other_vendors_text_is_untouched(self):
        for s in ("Hybrid", "Large Area", "Transmission", ""):
            self.assertEqual(friendly(s), s)
        self.assertEqual(friendly(None), "")


class TestVamasRoute(unittest.TestCase):
    def test_harwell_comment_reads_friendly(self):
        from readers.vamas import VamasFile
        kv = {"lens mode": "HSA_LENS_HYBRID"}
        self.assertEqual(
            friendly(VamasFile._lookup(VamasFile.__new__(VamasFile),
                                       "Lens mode", kv)),
            kratosterms.LENS_MODES["HSA_LENS_HYBRID"])


if __name__ == "__main__":
    unittest.main()
