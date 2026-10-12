"""Keeping windows on the monitor they are on: where a hover hint or a dialog
may go, whatever the screen size and however many monitors there are.

``monitor_bounds`` asks Windows for the monitor holding a screen point (its
usable area, taskbar left out, with ``work=True``) and falls back to the
virtual screen Tk reports. ``fit_geometry`` is the pure rule for a dialog:
make it no bigger than the monitor can show, centre it over the window it
belongs to, and pull it back inside the monitor. ``place_dialog`` applies it.
Only ``monitor_bounds`` and ``place_dialog`` touch Tk.
"""

from __future__ import annotations

import sys

# a window's frame around the client area that ``geometry`` sets: the title bar
# above, the borders at the sides and below
DECOR_Y = 40
DECOR_X = 16


def monitor_bounds(widget, x, y, work=False):
    """``(left, top, right, bottom)`` of the monitor holding the screen point
    ``(x, y)`` (the nearest one if it is on none); with ``work`` the part the
    taskbar leaves free. Windows asks the system, which knows every monitor;
    elsewhere, and if that fails, the whole virtual screen Tk reports."""
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes

            class MonitorInfo(ctypes.Structure):
                _fields_ = [("cbSize", wintypes.DWORD),
                            ("rcMonitor", wintypes.RECT),
                            ("rcWork", wintypes.RECT),
                            ("dwFlags", wintypes.DWORD)]
            user32 = ctypes.windll.user32
            user32.MonitorFromPoint.argtypes = [wintypes.POINT, wintypes.DWORD]
            user32.MonitorFromPoint.restype = wintypes.HANDLE
            user32.GetMonitorInfoW.argtypes = [wintypes.HANDLE,
                                               ctypes.POINTER(MonitorInfo)]
            mon = user32.MonitorFromPoint(wintypes.POINT(int(x), int(y)), 2)
            info = MonitorInfo()
            info.cbSize = ctypes.sizeof(MonitorInfo)
            if mon and user32.GetMonitorInfoW(mon, ctypes.byref(info)):
                r = info.rcWork if work else info.rcMonitor
                return r.left, r.top, r.right, r.bottom
        except Exception:
            pass
    vx, vy = widget.winfo_vrootx(), widget.winfo_vrooty()
    return (vx, vy, vx + max(widget.winfo_vrootwidth(), 1),
            vy + max(widget.winfo_vrootheight(), 1))


def fit_geometry(bounds, anchor, width, height):
    """``(w, h, x, y)`` for a ``width`` x ``height`` dialog belonging to the
    window ``anchor`` = ``(x, y, w, h)``, on the monitor area ``bounds``.

    The size is cut to what the area can show (frame included); the position
    centres the dialog across the anchor and a third of the way down, then is
    pulled back so the whole frame is inside the area."""
    left, top, right, bottom = bounds
    ax, ay, aw, ah = anchor
    w = max(1, min(width, right - left - DECOR_X))
    h = max(1, min(height, bottom - top - DECOR_Y))
    x = ax + max(0, (aw - w) // 2)
    y = ay + max(0, (ah - h) // 3)
    x = max(left, min(x, right - w - DECOR_X))
    y = max(top, min(y, bottom - h - DECOR_Y))
    return w, h, x, y


def place_dialog(dialog, anchor, width, height):
    """Size and place ``dialog`` over the window ``anchor`` (see
    ``fit_geometry``); returns the geometry string it set."""
    anchor.update_idletasks()
    ax, ay = anchor.winfo_rootx(), anchor.winfo_rooty()
    aw, ah = anchor.winfo_width(), anchor.winfo_height()
    bounds = monitor_bounds(anchor, ax + aw // 2, ay + ah // 2, work=True)
    w, h, x, y = fit_geometry(bounds, (ax, ay, aw, ah), width, height)
    geom = f"{w}x{h}+{x}+{y}"
    dialog.geometry(geom)
    return geom
