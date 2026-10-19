"""Portable mode: settings beside the program only when asked."""

import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import appinfo  # noqa: E402


class TestDataDir(unittest.TestCase):
    def test_default_is_home(self):
        with tempfile.TemporaryDirectory() as base:
            self.assertFalse(appinfo.portable(base))
            self.assertEqual(appinfo.data_dir(base, home="H"), "H")
            self.assertFalse(os.path.exists(os.path.join(base, "data")))

    def test_marker_gives_data_folder_beside_program(self):
        with tempfile.TemporaryDirectory() as base:
            open(os.path.join(base, appinfo.PORTABLE_MARKER), "w").close()
            self.assertTrue(appinfo.portable(base))
            got = appinfo.data_dir(base, home="H")
            self.assertEqual(got, os.path.join(base, "data"))
            self.assertTrue(os.path.isdir(got))

    def test_source_run_is_not_frozen(self):
        self.assertFalse(appinfo.frozen())
        self.assertEqual(appinfo.app_dir(), ROOT)


class TestBuildFiles(unittest.TestCase):
    def test_build_files_exist(self):
        for rel in ("build_exe.py", "requirements-build.txt",
                    os.path.join("packaging", "spectradeck.spec"),
                    "run_portable.bat"):
            self.assertTrue(os.path.isfile(os.path.join(ROOT, rel)), rel)

    def test_numpy_is_a_stated_requirement(self):
        with open(os.path.join(ROOT, "requirements.txt")) as fh:
            self.assertIn("numpy", fh.read())


if __name__ == "__main__":
    unittest.main()
