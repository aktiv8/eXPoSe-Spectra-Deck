"""Where a hover hint appears: under its widget, but never past the right edge
of the window / screen or off the bottom (the ribbon's About, CasaXPS and
KherveFitting buttons sit at the far right).

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest

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

    def test_the_text_can_change_without_a_second_hint(self):
        win = self.window(100, 100)
        b = self.button(win)
        tip = ee.Tooltip(b, "one", lambda: PALETTE)
        tip.set_text("two")
        self.assertEqual(tip.text, "two")


if __name__ == "__main__":
    unittest.main()
