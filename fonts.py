"""Bundled typefaces (all SIL Open Font License, the licence beside each in
``assets/fonts``) for Tk, matplotlib and the PDF report.

IBM Plex Sans is the default and the only face the Tk window uses. Tk cannot
load web fonts, so its files are registered with the operating system *for
this process only* (call :func:`register_process_fonts` before ``tk.Tk()``).
The other faces in :data:`CATALOG` are for plots (matplotlib) and, where
``Face.report`` is set, for the PDF text. Every step is best effort: if
anything fails, the app keeps the system font and the plots keep matplotlib's
default.
"""

from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass

FAMILY = "IBM Plex Sans"
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "assets", "fonts")
FILES = ("IBMPlexSans-Regular.ttf", "IBMPlexSans-Bold.ttf")


@dataclass(frozen=True)
class Face:
    """One bundled family: static TrueType files below ``FONT_DIR`` (variable
    and CFF fonts are not used: reportlab embeds neither, and matplotlib would
    get one weight), the licence beside them, and whether it has the Greek
    letters (α of "Al Kα") a PDF's text needs, since reportlab does not fall
    back to another font the way matplotlib does."""
    regular: str
    bold: str
    licence: str
    report: bool = True


# name = the family name inside the files (what matplotlib calls it)
CATALOG = {
    FAMILY: Face(FILES[0], FILES[1], "OFL.txt"),
    "Source Sans 3": Face("SourceSans3/SourceSans3-Regular.ttf",
                          "SourceSans3/SourceSans3-Bold.ttf",
                          "SourceSans3/OFL.txt"),
    "Inter": Face("Inter/Inter-Regular.ttf", "Inter/Inter-Bold.ttf",
                  "Inter/OFL.txt"),
    "IBM Plex Serif": Face("IBMPlexSerif/IBMPlexSerif-Regular.ttf",
                           "IBMPlexSerif/IBMPlexSerif-Bold.ttf",
                           "IBMPlexSerif/OFL.txt"),
    "STIX Two Text": Face("STIXTwoText/STIXTwoText-Regular.ttf",
                          "STIXTwoText/STIXTwoText-Bold.ttf",
                          "STIXTwoText/OFL.txt"),
    # no Greek: plots fall back to DejaVu Sans for α, a PDF cannot
    "IBM Plex Mono": Face("IBMPlexMono/IBMPlexMono-Regular.ttf",
                          "IBMPlexMono/IBMPlexMono-Bold.ttf",
                          "IBMPlexMono/OFL.txt", report=False),
}
FAMILIES = tuple(CATALOG)
REPORT_FAMILIES = tuple(n for n, f in CATALOG.items() if f.report)


def family_paths(family):
    """``(regular, bold)`` absolute paths of a catalog family, or ``()`` when
    it is unknown or either file is missing."""
    face = CATALOG.get(family)
    if face is None:
        return ()
    paths = tuple(os.path.join(FONT_DIR, f) for f in (face.regular, face.bold))
    return paths if all(os.path.isfile(p) for p in paths) else ()

# Type scale (points): caption / body / strong / section / panel title
SIZE = {"caption": 8, "body": 9, "strong": 9, "section": 10, "title": 10}

_TK_NAMED = ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont",
             "TkCaptionFont", "TkSmallCaptionFont", "TkIconFont",
             "TkTooltipFont")


def font_paths():
    return [p for p in (os.path.join(FONT_DIR, f) for f in FILES)
            if os.path.isfile(p)]


def register_process_fonts() -> bool:
    """Make the bundled font files visible to Tk. True if all were added."""
    paths = font_paths()
    if len(paths) != len(FILES):
        return False
    try:
        if sys.platform.startswith("win"):
            add = ctypes.windll.gdi32.AddFontResourceExW
            add.argtypes = [ctypes.c_wchar_p, ctypes.c_uint, ctypes.c_void_p]
            return all(add(p, 0x10, None) for p in paths)      # FR_PRIVATE
        if sys.platform == "darwin":
            ct = ctypes.cdll.LoadLibrary(
                "/System/Library/Frameworks/CoreText.framework/CoreText")
            cf = ctypes.cdll.LoadLibrary(
                "/System/Library/Frameworks/CoreFoundation.framework/"
                "CoreFoundation")
            make = cf.CFURLCreateFromFileSystemRepresentation
            make.restype = ctypes.c_void_p
            make.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_long,
                             ctypes.c_bool]
            reg = ct.CTFontManagerRegisterFontsForURL
            reg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
            reg.restype = ctypes.c_bool
            ok = True
            for p in paths:
                raw = os.fsencode(p)
                ok &= bool(reg(make(None, raw, len(raw), False), 1, None))
            return ok
        fc = ctypes.cdll.LoadLibrary("libfontconfig.so.1")
        fc.FcConfigAppFontAddFile.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        return all(fc.FcConfigAppFontAddFile(None, os.fsencode(p))
                   for p in paths)
    except Exception:                       # noqa: BLE001 - never block start-up
        return False


def apply_tk_fonts(root) -> str:
    """Point Tk's named fonts at the bundled family. Returns the family now in
    use (the system default's family if the bundled one isn't available)."""
    import tkinter.font as tkfont
    default = tkfont.nametofont("TkDefaultFont").actual("family")
    try:
        if FAMILY not in set(tkfont.families(root)):
            return default
        for name in _TK_NAMED:
            f = tkfont.nametofont(name)
            f.configure(family=FAMILY, size=SIZE["body"])
        tkfont.nametofont("TkCaptionFont").configure(weight="bold")
        return FAMILY
    except Exception:                       # noqa: BLE001
        return default


def register_matplotlib() -> str | None:
    """Add every bundled font file to matplotlib; returns the default family
    (IBM Plex Sans) or None."""
    try:
        from matplotlib import font_manager
        paths = list(font_paths())
        for fam in FAMILIES:
            paths += [p for p in family_paths(fam) if p not in paths]
        for p in paths:
            font_manager.fontManager.addfont(p)
        names = {f.name for f in font_manager.fontManager.ttflist}
        return FAMILY if FAMILY in names else None
    except Exception:                       # noqa: BLE001
        return None
