"""Every Python file of the project compiles without a SyntaxWarning.

An invalid escape sequence (a Windows path such as ``D:\\Temp`` written in a
plain docstring) only warns when a module is compiled, i.e. on the first
launch after an edit, then disappears because the cached ``.pyc`` hides it. A
valid-but-unintended one (``\\f`` in ``...\\for claude files``) is worse: it
silently becomes a control character. Write such paths in a raw string.

Run:  python -m unittest discover tests
"""

import glob
import os
import unittest
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestCompileClean(unittest.TestCase):
    def test_no_file_warns_when_compiled(self):
        files = sorted(
            glob.glob(os.path.join(ROOT, "*.py"))
            + glob.glob(os.path.join(ROOT, "readers", "*.py"))
            + glob.glob(os.path.join(ROOT, "tests", "*.py"))
            + glob.glob(os.path.join(ROOT, "manual", "*.py")))
        self.assertGreater(len(files), 100)
        bad = []
        for path in files:
            with open(path, encoding="utf-8") as fh:
                src = fh.read()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                compile(src, path, "exec")
            bad += [f"{os.path.relpath(path, ROOT)}: {w.message}"
                    for w in caught
                    if issubclass(w.category, (SyntaxWarning,
                                               DeprecationWarning))]
        self.assertEqual(bad, [])


if __name__ == "__main__":
    unittest.main()
