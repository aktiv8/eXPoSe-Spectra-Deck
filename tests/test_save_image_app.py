"""Save plot image / Save this panel as image in the real window (needs a
display and matplotlib; skipped otherwise): PNG, JPEG and SVG are written,
and one panel's image differs from the whole page's.

Run:  python -m unittest discover tests
"""

import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL and ee.HAVE_PIL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

from test_panel_views_app import doc, reg  # noqa: E402


@unittest.skipUnless(HAVE_MPL, "matplotlib / Pillow / Tk not available")
class TestSaveImage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        cls._rc = matplotlib.rcParams.copy()
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()
        cls.docs = {
            "a.vms": doc("a.vms", [reg("C 1s", "A", 286.0)]),
            "b.vms": doc("b.vms", [reg("O 1s", "B", 288.0)])}
        cls._load = ee.load_file
        ee.load_file = lambda p: cls.docs[os.path.basename(p)]
        cls._load_config = ee.load_config
        ee.load_config = lambda: {}
        cls.ws = ee.Workspace(cls.root)
        cls.ws.group_var.set("Element, per sample")
        cls.ws.panels_var.set("Auto")
        for p in cls.docs:
            cls.ws._add_file(p, file_id=p[:-4], refresh=False)
        cls.ws._finish_adding()
        cls.ws.checked = {id(r) for d in cls.ws.docs for r in d.regions}
        cls._dialog = ee.workbook_ui.ImageOptionsDialog
        ee.workbook_ui.ImageOptionsDialog = (
            lambda root, app, go: go("Paper (white)", 6, 4, 60))
        cls._ask = ee.filedialog.asksaveasfilename
        cls._tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        ee.workbook_ui.ImageOptionsDialog = cls._dialog
        ee.filedialog.asksaveasfilename = cls._ask
        ee.load_file = cls._load
        ee.load_config = cls._load_config
        cls.root.destroy()
        cls._tmp.cleanup()
        import matplotlib
        matplotlib.rcParams.update(cls._rc)

    def save(self, name, panel=None):
        path = os.path.join(self._tmp.name, name)
        ee.filedialog.asksaveasfilename = lambda **kw: path
        if panel is None:
            self.ws.save_plot_image()
        else:
            self.ws.save_panel_image(panel)
        self.assertTrue(os.path.isfile(path), name)
        return path

    def test_formats_are_written(self):
        from PIL import Image
        for name, fmt in (("p.png", "PNG"), ("p.jpg", "JPEG")):
            with Image.open(self.save(name)) as im:
                self.assertEqual(im.format, fmt)
        with open(self.save("p.svg"), encoding="utf-8") as f:
            self.assertIn("<svg", f.read(2000))

    def test_one_panel_is_not_the_page(self):
        from PIL import Image
        keys = [k for k, _ in self.ws._groups()]
        self.assertEqual(len(keys), 2)
        page = Image.open(self.save("page.png")).tobytes()
        one = Image.open(self.save("one.png", keys[0])).tobytes()
        self.assertNotEqual(page, one)


if __name__ == "__main__":
    unittest.main()
