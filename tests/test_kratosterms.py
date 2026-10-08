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

    def test_magnification(self):
        self.assertEqual(friendly("F_MHSA_MEDIUM_MAGN"), "Medium magnification")
        self.assertEqual(friendly("MHSA_HIGH_MAGN"), "High magnification")

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
                                       "Lens mode", kv)), "Hybrid")


if __name__ == "__main__":
    unittest.main()
