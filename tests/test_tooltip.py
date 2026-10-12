"""Where a hover hint appears: under its widget, but never past the right edge
of the window / screen or off the bottom (the ribbon's About, CasaXPS and
KherveFitting buttons sit at the far right).

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_APP = True
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_APP = None, None, False

PALETTE = {"border": "#888888", "panel": "#ffffff", "fg": "#000000"}
LONG = ("Open CasaXPS\nC:\\Program Files\\…\\CasaXPS 64bit\\CasaXPS.exe\n"
        "Right-click: another copy, forget or hide")


@unittest.skipUnless(HAVE_APP, "Tk not available")
class TestTooltipPlacement(unittest.TestCase):
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

    def window(self, x, y, w=500, h=200):
        win = tk.Toplevel(self.root)
        self.addCleanup(win.destroy)
        win.geometry(f"{w}x{h}+{x}+{y}")
        return win

    def button(self, win, **pack):
        b = tk.Button(win, text="B")
        b.pack(**pack)
        win.update()
        return b

    def shown(self, button):
        tip = ee.Tooltip(button, LONG, lambda: PALETTE)
        tip._show()
        self.addCleanup(tip._hide)
        tip._tip.update_idletasks()
        return (tip._tip.winfo_x(), tip._tip.winfo_y(),
                tip._tip.winfo_width(), tip._tip.winfo_height())

    def test_a_button_at_the_right_edge_keeps_its_tip_inside_the_window(self):
        win = self.window(100, 100)
        b = self.button(win, side="right", anchor="n")
        x, y, w, h = self.shown(b)
        self.assertGreater(w, 150)                       # a real, wide tip
        self.assertLessEqual(x + w, win.winfo_rootx() + win.winfo_width())
        self.assertGreaterEqual(x, win.winfo_rootx())
        self.assertGreaterEqual(y, b.winfo_rooty() + b.winfo_height())

    def test_a_window_at_the_screen_edge_keeps_the_tip_on_the_screen(self):
        sw = self.root.winfo_screenwidth()
        win = self.window(sw - 260, 100, w=260)
        b = self.button(win, side="right", anchor="n")
        x, _y, w, _h = self.shown(b)
        self.assertLessEqual(x + w, sw)

    def test_a_button_in_the_middle_is_under_its_own_left_edge(self):
        win = self.window(100, 100, w=900)
        b = self.button(win, side="left", anchor="n")
        x, y, _w, _h = self.shown(b)
        self.assertEqual(x, b.winfo_rootx())
        self.assertEqual(y, b.winfo_rooty() + b.winfo_height() + 4)

    def test_no_room_below_puts_the_tip_above_the_button(self):
        sh = self.root.winfo_screenheight()
        win = self.window(100, sh - 80, h=70)
        b = self.button(win, side="left", anchor="s")
        _x, y, _w, h = self.shown(b)
        self.assertLess(y + h, b.winfo_rooty())
        self.assertLessEqual(y + h, sh)

    def fake(self, wx, wy, win_x, win_w, wh=24):
        """A widget at (wx, wy) in a window at win_x, as _place sees it."""
        top = type("Top", (), {
            "winfo_rootx": lambda s: win_x, "winfo_width": lambda s: win_w})()
        return type("W", (), {
            "winfo_toplevel": lambda s: top, "winfo_rootx": lambda s: wx,
            "winfo_rooty": lambda s: wy, "winfo_height": lambda s: wh})()

    def place(self, widget, monitor, w=214, h=40):
        tip = ee.Tooltip.__new__(ee.Tooltip)
        tip.widget = widget
        with mock.patch.object(ee, "_monitor_bounds",
                               lambda *_a: monitor):
            return tip._place(w, h)

    def test_a_second_monitor_to_the_right_is_not_the_primary_one(self):
        # window on a monitor spanning 1920..3840, button at its right edge
        x, y = self.place(self.fake(3800, 10, 2900, 940), (1920, 0, 3840, 1080))
        self.assertGreaterEqual(x, 1920)
        self.assertLessEqual(x + 214, 3840)
        self.assertEqual(x, 2900 + 940 - 214 - 4)           # window's edge

    def test_a_monitor_left_of_the_primary_has_negative_coordinates(self):
        x, y = self.place(self.fake(-60, 10, -1000, 960), (-1920, 0, 0, 1080))
        self.assertGreaterEqual(x, -1920)
        self.assertLessEqual(x + 214, 0)
        x, y = self.place(self.fake(-1900, 10, -1920, 960),
                          (-1920, 0, 0, 1080))
        self.assertEqual(x, -1900)                         # nothing to fix

    def test_a_monitor_above_the_primary_flips_up_at_its_own_bottom(self):
        x, y = self.place(self.fake(100, -50, 0, 900), (0, -1080, 1920, 0),
                          h=40)
        self.assertEqual(y, -50 - 40 - 4)                  # above the button
        self.assertGreaterEqual(y, -1080)

    def test_the_monitor_lookup_contains_the_point_it_was_given(self):
        left, top, right, bottom = ee._monitor_bounds(self.root, 10, 10)
        self.assertLess(left, right)
        self.assertLess(top, bottom)
        self.assertTrue(left <= 10 < right and top <= 10 < bottom)

    def test_the_text_can_change_without_a_second_hint(self):
        win = self.window(100, 100)
        b = self.button(win)
        tip = ee.Tooltip(b, "one", lambda: PALETTE)
        tip.set_text("two")
        self.assertEqual(tip.text, "two")


if __name__ == "__main__":
    unittest.main()
