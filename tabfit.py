"""Keep the labels of a ``ttk.Notebook`` readable.

A notebook never wraps or scrolls its tab strip: when the tabs are wider than
the notebook, Tk squeezes them, eating their padding and then cutting the
text. ``fit_tabs`` measures what the labels need and, in order, (1) switches
the notebook to the compact tab style, (2) widens its dialog (never past the
monitor). A notebook in a pane that cannot grow keeps the compact style and
short labels; the full name is available through ``tab_names``.
"""

import tkinter as tk
from tkinter import ttk
import tkinter.font as tkfont

import screenfit

COMPACT = "Compact.TNotebook"
MARGIN = 10                    # slack for the strip's own edges


def _padding(nb, style):
    """Left + right padding of a tab in ``style`` (pixels)."""
    try:
        raw = ttk.Style(nb).lookup(style + ".Tab", "padding") or ""
        nums = [int(float(v)) for v in str(raw).replace("{", "").replace(
            "}", "").split()]
    except (tk.TclError, ValueError):
        nums = []
    if len(nums) == 1:
        return 2 * nums[0]
    if len(nums) >= 2:
        return 2 * nums[0]
    return 24


def needed_width(nb):
    """Pixels the visible tab labels need in the notebook's present style."""
    try:
        style = str(nb.cget("style")) or "TNotebook"
        name = str(ttk.Style(nb).lookup(style + ".Tab", "font")
                   or "TkDefaultFont")
        font = tkfont.Font(root=nb, name=name, exists=True)
    except tk.TclError:
        return 0
    pad = _padding(nb, style)
    total = MARGIN
    for i in range(len(nb.tabs())):
        if nb.tab(i, "state") == "hidden":
            continue
        total += font.measure(nb.tab(i, "text")) + pad
    return total


def is_clipped(nb):
    """True when the notebook is mapped and its tabs need more room than it
    has (Tk is then squeezing them)."""
    width = nb.winfo_width()
    return width > 20 and needed_width(nb) > width


def tab_names(nb):
    """The full name of each tab, for hints (``nb._tab_names`` if the owner
    set them, else the labels)."""
    return getattr(nb, "_tab_names", None) or [
        nb.tab(i, "text") for i in range(len(nb.tabs()))]


def fit_tabs(nb, window=None):
    """Make the tabs of ``nb`` fit. ``window`` is the dialog that may grow."""
    def check(_event=None):
        try:
            if not nb.winfo_exists() or nb.winfo_width() <= 20:
                return
            if is_clipped(nb) and str(nb.cget("style")) != COMPACT:
                nb.configure(style=COMPACT)
            if is_clipped(nb) and window is not None:
                grow(window, nb)
        except tk.TclError:
            pass

    nb.bind("<Configure>", check, add="+")
    nb.bind("<Map>", check, add="+")
    return check


def grow(window, nb):
    """Widen ``window`` by what the tabs of ``nb`` lack, within its monitor."""
    deficit = needed_width(nb) - nb.winfo_width()
    if deficit <= 0:
        return False
    window.update_idletasks()
    w, h = window.winfo_width(), window.winfo_height()
    x, y = window.winfo_x(), window.winfo_y()
    bounds = monitor_bounds(window)
    new_w = min(w + deficit, bounds[2] - bounds[0])
    if new_w <= w:
        return False
    x = max(bounds[0], min(x, bounds[2] - new_w))
    window.geometry(f"{new_w}x{h}+{x}+{y}")
    return True


def monitor_bounds(window):
    try:
        return screenfit.monitor_bounds(
            window, window.winfo_rootx(), window.winfo_rooty(), work=True)
    except Exception:                                   # pragma: no cover
        return (0, 0, window.winfo_screenwidth(), window.winfo_screenheight())
