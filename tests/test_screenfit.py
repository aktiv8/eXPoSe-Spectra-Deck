"""Dialogs stay on their monitor: sized to what the monitor can show, centred
over their window, pulled back inside (screenfit.py).

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest
from types import SimpleNamespace
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import screenfit as sf  # noqa: E402

try:
    import tkinter as tk
except Exception:                                   # pragma: no cover
    tk = None

LAPTOP = (0, 0, 1366, 728)              # 1366 x 768 less the taskbar
FULLHD = (0, 0, 1920, 1032)


def inside(bounds, w, h, x, y):
    left, top, right, bottom = bounds
    return (x >= left and y >= top and x + w + sf.DECOR_X <= right
            and y + h + sf.DECOR_Y <= bottom)


class TestFitGeometry(unittest.TestCase):
    def test_a_dialog_that_fits_is_centred_over_its_window(self):
        w, h, x, y = sf.fit_geometry(FULLHD, (100, 100, 1000, 600), 500, 300)
        self.assertEqual((w, h), (500, 300))
        self.assertEqual(x, 100 + (1000 - 500) // 2)
        self.assertEqual(y, 100 + (600 - 300) // 3)

    def test_a_window_in_the_corner_does_not_take_the_dialog_off_screen(self):
        # the main window at the bottom right, as in a real session
        for size in ((960, 660), (640, 760), (520, 680), (380, 260)):
            w, h, x, y = sf.fit_geometry(FULLHD, (920, 432, 1000, 600), *size)
            self.assertTrue(inside(FULLHD, w, h, x, y), (size, x, y))

    def test_a_dialog_taller_than_the_screen_is_cut_to_it(self):
        w, h, x, y = sf.fit_geometry(LAPTOP, (0, 0, 1366, 700), 920, 760)
        self.assertLess(h, 760)
        self.assertEqual(h, 728 - sf.DECOR_Y)
        self.assertEqual(y, 0)
        self.assertTrue(inside(LAPTOP, w, h, x, y))

    def test_a_dialog_wider_than_the_screen_is_cut_to_it(self):
        w, h, x, y = sf.fit_geometry(LAPTOP, (0, 0, 1366, 700), 1500, 300)
        self.assertEqual(w, 1366 - sf.DECOR_X)
        self.assertTrue(inside(LAPTOP, w, h, x, y))

    def test_a_small_dialog_is_never_enlarged(self):
        w, h, *_ = sf.fit_geometry(LAPTOP, (0, 0, 1366, 700), 360, 200)
        self.assertEqual((w, h), (360, 200))

    def test_a_monitor_with_negative_coordinates(self):
        left = (-1920, 0, 0, 1032)                     # left of the primary
        w, h, x, y = sf.fit_geometry(left, (-500, 600, 480, 400), 960, 660)
        self.assertTrue(inside(left, w, h, x, y), (x, y))
        above = (0, -1080, 1920, -48)                  # above the primary
        w, h, x, y = sf.fit_geometry(above, (1500, -300, 400, 250), 520, 680)
        self.assertTrue(inside(above, w, h, x, y), (x, y))


@unittest.skipUnless(tk, "Tk not available")
class TestPlaceDialog(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def test_the_geometry_it_sets_is_inside_the_monitor(self):
        win = tk.Toplevel(self.root)
        self.addCleanup(win.destroy)
        sw = self.root.winfo_screenwidth()
        win.geometry(f"900x500+{sw - 900}+300")
        win.update()
        dlg = tk.Toplevel(self.root)
        self.addCleanup(dlg.destroy)
        area = sf.monitor_bounds(win, win.winfo_rootx() + 10,
                                 win.winfo_rooty() + 10, work=True)
        geom = sf.place_dialog(dlg, win, 960, 5000)
        w, h, x, y = (int(v) for v in geom.replace("x", "+").split("+"))
        self.assertTrue(inside(area, w, h, x, y), geom)
        self.assertLess(h, 5000)                       # cut to the monitor

    def test_the_work_area_is_inside_the_monitor(self):
        mon = sf.monitor_bounds(self.root, 10, 10)
        work = sf.monitor_bounds(self.root, 10, 10, work=True)
        self.assertTrue(mon[0] <= work[0] and mon[1] <= work[1]
                        and work[2] <= mon[2] and work[3] <= mon[3])

    def test_the_shared_dialog_setup_uses_it(self):
        import workbook_ui
        app = SimpleNamespace(
            root=self.root, palette={"bg": "#ffffff"},
            themes=SimpleNamespace(recolor_tk=lambda w: None))
        dlg = tk.Toplevel(self.root)
        self.addCleanup(dlg.destroy)
        with mock.patch.object(sf, "monitor_bounds",
                               lambda *_a, **_k: LAPTOP),                 mock.patch.object(dlg, "geometry") as geometry:
            workbook_ui._finish(dlg, app, 920, 760)
        geom = geometry.call_args[0][0]               # what it asked for
        w, h, x, y = (int(v) for v in geom.replace("x", "+").split("+"))
        self.assertEqual(h, 728 - sf.DECOR_Y)
        self.assertTrue(inside(LAPTOP, w, h, x, y), geom)


if __name__ == "__main__":
    unittest.main()
