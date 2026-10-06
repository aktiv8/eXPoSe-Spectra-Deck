"""Type and colour of the PDF report (Tk-free; reportlab only inside functions).

One ``Look`` says which fonts and which accent colour the report is set in.
The typeface is the bundled IBM Plex Sans (registered with reportlab on first
use; it has the "α" of "Al Kα" and the "µ" / "°" that the built-in Helvetica
lacks) and the accent is the one chosen for the cover (``reportspec``). Without
the font files everything falls back to Helvetica. The text and the fill of a
table header use ``ink`` (the accent darkened), so a pale accent such as amber
still reads on white; ``accent`` itself is for rules and bars.
"""

from __future__ import annotations

from dataclasses import dataclass

import covers
import fonts

REGULAR, BOLD = "IBMPlexSans", "IBMPlexSans-Bold"
FALLBACK = ("Helvetica", "Helvetica-Bold")
TEXT = "#222222"
MUTED = "#6B7785"
RULE = "#C8D0D8"

_registered: dict = {}          # family -> (regular, bold) reportlab names


def _names(family):
    base = "".join(family.split())
    return base, base + "-Bold"


def resolve_family(family=fonts.FAMILY):
    """The family the report is really set in: ``family`` when it may be used
    in a PDF and its files load, else IBM Plex Sans, else '' (Helvetica)."""
    for fam in (family, fonts.FAMILY):
        if fam in fonts.REPORT_FAMILIES and register_fonts(fam)[0] != FALLBACK[0]:
            return fam
    return ""


def register_fonts(family=fonts.FAMILY):
    """``(regular, bold)`` reportlab font names of a bundled family when its
    files are there and load, else Helvetica. Registered once per process and
    family. Asking for a family that is not in the catalog (or not meant for
    a PDF) gives the same answer as IBM Plex Sans."""
    if family not in fonts.REPORT_FAMILIES:
        family = fonts.FAMILY
    if family in _registered:
        return _registered[family]
    _registered[family] = FALLBACK
    paths = fonts.family_paths(family)
    if paths:
        reg_name, bold_name = _names(family)
        try:
            from reportlab.pdfbase import pdfmetrics
            from reportlab.pdfbase.ttfonts import TTFont
            pdfmetrics.registerFont(TTFont(reg_name, paths[0]))
            pdfmetrics.registerFont(TTFont(bold_name, paths[1]))
            # no italic files: <i> keeps the upright face rather than failing
            pdfmetrics.registerFontFamily(reg_name, normal=reg_name,
                                          bold=bold_name, italic=reg_name,
                                          boldItalic=bold_name)
            _registered[family] = (reg_name, bold_name)
        except Exception:                    # noqa: BLE001 - never stop a report
            _registered[family] = FALLBACK
    return _registered[family]


def font_file(bold=False, family=fonts.FAMILY):
    """The regular (or, ``bold``, bold) typeface's file, for PyMuPDF (page
    footers, divider pages); '' = none (Helvetica)."""
    fam = resolve_family(family)
    if not fam:
        return ""
    return fonts.family_paths(fam)[1 if bold else 0]


def page_size(option="a4"):
    """``(portrait, figure_landscape_in)`` for the PDF's ``page`` option
    ("a4" or "letter", see ``reportspec.OPTIONS``; anything else is "a4"):
    the reportlab page size (points, portrait) for the text-flow document,
    and the matching landscape size in inches a figure or camera/SnapMap
    page should use so every page in the report is the same size."""
    from reportlab.lib.pagesizes import A4, LETTER
    from reportlab.lib.units import inch
    portrait = LETTER if option == "letter" else A4
    return portrait, (portrait[1] / inch, portrait[0] / inch)


@dataclass(frozen=True)
class Look:
    accent: str = covers.DEFAULT_ACCENT
    font: str = FALLBACK[0]
    bold: str = FALLBACK[1]
    family: str = ""             # bundled family behind font/bold ('' = none)

    @property
    def ink(self):
        """The accent darkened: text and table headers."""
        return covers.mix(self.accent, "#000000", 0.2)

    def tint(self, t):
        """The accent faded towards white (0 = the accent, 1 = white)."""
        return covers.tint(self.accent, t)


def look(accent="", family=fonts.FAMILY):
    """The ``Look`` for an accent colour ('' = the default one) and a bundled
    family (see ``resolve_family`` for what an unusable one becomes)."""
    fam = resolve_family(family)
    font, bold = register_fonts(fam) if fam else FALLBACK
    return Look(covers.valid_accent(accent) or covers.DEFAULT_ACCENT, font,
                bold, fam)


def rgb(hex_colour):
    """A '#RRGGBB' colour as the 0-1 tuple PyMuPDF takes."""
    h = hex_colour.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def styles(lk):
    """The paragraph styles of the report, as a dict of reportlab styles."""
    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_LEFT

    text = colors.HexColor(TEXT)
    ink = colors.HexColor(lk.ink)
    base = ParagraphStyle("base", fontName=lk.font, fontSize=10, leading=14.5,
                          textColor=text, alignment=TA_LEFT)

    def make(name, **kw):
        return ParagraphStyle(name, parent=base, **kw)

    return {
        "body": make("body"),
        "small": make("small", fontSize=8, leading=10),
        "title": make("title", fontName=lk.bold, fontSize=26, leading=30,
                      textColor=ink, spaceAfter=6),
        "h1": make("h1", fontName=lk.bold, fontSize=16, leading=20,
                   textColor=ink, spaceBefore=16, spaceAfter=6,
                   keepWithNext=1),
        "h2": make("h2", fontName=lk.bold, fontSize=13, leading=16,
                   textColor=ink, spaceBefore=10, spaceAfter=4,
                   keepWithNext=1),
        "label": make("label", fontName=lk.bold, fontSize=8.5, leading=12,
                      textColor=colors.HexColor(MUTED)),
        "toc1": make("toc1", fontName=lk.bold, fontSize=11, leading=15,
                     textColor=ink),
        "toc2": make("toc2", fontSize=9.5, leading=13, leftIndent=14,
                     textColor=colors.HexColor("#444444")),
    }
