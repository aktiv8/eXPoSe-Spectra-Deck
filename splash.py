"""The start-up splash screen and the window icon.

``begin`` is the very first thing the application script does (before the slow
imports such as matplotlib): it registers the bundled fonts, makes the one Tk
root the program will use, keeps it hidden and, unless switched off, shows the
splash on top. ``status`` updates the line under the picture while the rest
loads; ``finish`` shows the main window and closes the splash once it has been
up for ``MIN_SECONDS`` (a click or Esc closes it sooner).

The picture is whatever the user put in ``assets/`` (see
``appinfo.image_path``); without one, or if it cannot be read, the splash is a
plain card with the name. Nothing here may stop the program from starting, so
every step is best effort. Turn it off with ``--no-splash`` or the Help menu
(``show_splash`` in the settings file).
"""

from __future__ import annotations

import json
import os
import sys
import time
import tkinter as tk

import appinfo

MIN_SECONDS = 1.2
MAX_SCREEN_FRACTION = 0.6           # the picture never fills more than this
CARD_BG, CARD_FG, CARD_MUTED = "#12161A", "#D7DEE4", "#8E9AA5"
CARD_ACCENT, CARD_EDGE = "#4CC2E0", "#2A333B"
CONFIG_NAME = ".spectradeck_config.json"

_S = {"root": None, "win": None, "status": None, "photo": None, "t0": 0.0,
      "icon": None}


def wanted(argv=None, config_path=None):
    """Whether to show the splash: not with ``--no-splash`` and not when the
    settings file says ``show_splash: false``."""
    argv = sys.argv if argv is None else argv
    if "--no-splash" in argv:
        return False
    path = config_path or os.path.join(os.path.expanduser("~"), CONFIG_NAME)
    try:
        with open(path, encoding="utf-8") as fh:
            return bool(json.load(fh).get("show_splash", True))
    except Exception:
        return True


def dpi_aware_wanted(argv=None, config_path=None):
    """Whether to tell Windows this program draws at the screen's real pixel
    density (sharp text on a 125-200 % display instead of a stretched bitmap).
    Off unless asked: ``--dpi-aware`` or ``dpi_aware: true`` in the settings
    file (Help menu); ``--no-dpi-aware`` wins. Pixel-sized parts of the layout
    were designed at 100 %, so this stays a choice until it has been looked at
    on every screen the program is used on."""
    argv = sys.argv if argv is None else argv
    if "--no-dpi-aware" in argv:
        return False
    if "--dpi-aware" in argv:
        return True
    path = config_path or os.path.join(os.path.expanduser("~"), CONFIG_NAME)
    try:
        with open(path, encoding="utf-8") as fh:
            return bool(json.load(fh).get("dpi_aware", False))
    except Exception:
        return False


def dpi_plan(platform=None):
    """The Windows calls ``make_dpi_aware`` tries, in order (``[]`` elsewhere)."""
    platform = sys.platform if platform is None else platform
    if not platform.startswith("win"):
        return []
    return ["shcore.SetProcessDpiAwareness(1)", "user32.SetProcessDPIAware()"]


def make_dpi_aware(platform=None):
    """Declare the process DPI-aware (system level). Call before the first Tk
    window exists. Returns the call that worked, or ``""``; never raises."""
    for call in dpi_plan(platform):
        try:
            import ctypes
            if call.startswith("shcore"):
                if ctypes.windll.shcore.SetProcessDpiAwareness(1) in (0, ):
                    return call
            elif ctypes.windll.user32.SetProcessDPIAware():
                return call
        except Exception:                   # noqa: BLE001 - best effort
            continue
    return ""


def create_root():
    """The program's Tk root (drag-and-drop capable when tkinterdnd2 is
    installed)."""
    try:
        from tkinterdnd2 import TkinterDnD
        return TkinterDnD.Tk()
    except Exception:
        return tk.Tk()


def begin(argv=None):
    """Prelude of the application script: fonts, the hidden root and the
    splash. Returns the root."""
    set_app_id()
    if dpi_aware_wanted(argv):
        make_dpi_aware()                    # before any window exists
    try:
        import fonts
        fonts.register_process_fonts()      # before Tk enumerates fonts
    except Exception:
        pass
    root = create_root()
    root.withdraw()
    _S["root"] = root
    if wanted(argv):
        try:
            show(root)
        except Exception:
            _S["win"] = None
    return root


def take_root():
    """The root made by ``begin`` (once), or None."""
    root, _S["root"] = _S["root"], None
    return root


def set_app_id():
    """Give Windows a stable identity so the taskbar shows the program's own
    icon instead of Python's."""
    if not sys.platform.startswith("win"):
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
            "aktiv8.SpectraDeck")
    except Exception:
        pass


def _load_picture(root, path):
    """A Tk image of ``path`` no larger than the screen fraction, or None."""
    if not path:
        return None
    limit_w = int(root.winfo_screenwidth() * MAX_SCREEN_FRACTION)
    limit_h = int(root.winfo_screenheight() * MAX_SCREEN_FRACTION)
    try:
        from PIL import Image, ImageTk
        im = Image.open(path).convert("RGBA")
        scale = min(1.0, limit_w / im.width, limit_h / im.height)
        if scale < 1.0:
            im = im.resize((max(1, int(im.width * scale)),
                            max(1, int(im.height * scale))), Image.LANCZOS)
        return ImageTk.PhotoImage(im, master=root)
    except Exception:
        pass
    try:                                     # no Pillow: PNG / GIF only
        img = tk.PhotoImage(file=path, master=root)
        k = max(1, -(-img.width() // limit_w), -(-img.height() // limit_h))
        return img.subsample(k) if k > 1 else img
    except Exception:
        return None


def show(root):
    """Show the splash window (a borderless, centred, topmost Toplevel)."""
    win = tk.Toplevel(root)
    win.withdraw()
    win.overrideredirect(True)
    try:
        win.attributes("-topmost", True)
    except tk.TclError:
        pass
    card = tk.Frame(win, bg=CARD_BG, highlightthickness=1,
                    highlightbackground=CARD_EDGE)
    card.pack(fill="both", expand=True)
    photo = _load_picture(root, appinfo.image_path("splash"))
    _S["photo"] = photo
    if photo is not None:
        tk.Label(card, image=photo, bg=CARD_BG, bd=0).pack()
    family = "IBM Plex Sans"
    tk.Label(card, text=appinfo.NAME, bg=CARD_BG, fg=CARD_FG,
             font=(family, 20, "bold")).pack(padx=40, pady=(18, 0))
    tk.Label(card, text=f"Version {appinfo.VERSION}", bg=CARD_BG,
             fg=CARD_MUTED, font=(family, 10)).pack()
    status = tk.Label(card, text="Starting…", bg=CARD_BG, fg=CARD_ACCENT,
                      font=(family, 9))
    status.pack(padx=40, pady=(10, 18))
    _S.update(win=win, status=status, t0=time.monotonic())
    win.update_idletasks()
    w, h = win.winfo_reqwidth(), win.winfo_reqheight()
    x = max(0, (root.winfo_screenwidth() - w) // 2)
    y = max(0, (root.winfo_screenheight() - h) // 2)
    win.geometry(f"{w}x{h}+{x}+{y}")
    win.deiconify()
    for widget in (win, card):
        widget.bind("<Button-1>", lambda e: _destroy())
    win.bind("<Escape>", lambda e: _destroy())
    win.lift()
    win.update()
    return win


def is_open():
    return _S["win"] is not None


def status(text):
    """Update the line under the picture (no-op without a splash)."""
    label, win = _S["status"], _S["win"]
    if label is None or win is None:
        return
    try:
        label.configure(text=text)
        win.update()
    except tk.TclError:
        _S["win"] = _S["status"] = None


def _destroy():
    win = _S["win"]
    _S["win"] = _S["status"] = _S["photo"] = None
    if win is not None:
        try:
            win.destroy()
        except tk.TclError:
            pass


def finish(root):
    """Show the main window and close the splash once it has been up for
    ``MIN_SECONDS``. The main window is mapped straight away (behind the
    splash) so its layout can settle."""
    root.deiconify()
    if _S["win"] is None:
        return
    wait = MIN_SECONDS - (time.monotonic() - _S["t0"])
    if wait > 0:
        root.after(int(wait * 1000), _destroy)
    else:
        _destroy()


def set_window_icon(root):
    """Use the supplied picture, shrunk, as the window icon."""
    path = appinfo.image_path("splash")
    if not path:
        return False
    try:
        from PIL import Image, ImageTk
        im = Image.open(path).convert("RGBA")
        side = max(im.size)
        square = Image.new("RGBA", (side, side), (0, 0, 0, 0))
        square.paste(im, ((side - im.width) // 2, (side - im.height) // 2))
        _S["icon"] = ImageTk.PhotoImage(square.resize((64, 64),
                                                      Image.LANCZOS),
                                        master=root)
        root.iconphoto(True, _S["icon"])
        return True
    except Exception:
        return False
