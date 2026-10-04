"""Cover pictures for the report and the slides.

Tk-free. A *cover* is the small dict kept in the report spec
(``reportspec``): ``{"design": id, "image": path, "accent": "#RRGGBB"}``.
Designs are

* built in and drawn here from one accent colour (so they match whatever
  colour the report uses): ``ribbon`` (stacked, fading peaks), ``band`` (a solid
  band with fine rules and a silhouette), ``minimal`` (a rule and whitespace),
  ``data`` (the spectrum you are reporting on, as a hero line; the ribbon when
  there is none) and ``grid`` (a peak map);
* a picture in ``assets/covers/`` (``file:<name>``), or any picture the user
  browsed to (``image`` with the path in ``cover["image"]``): cropped to the
  shape of the space from the middle, never stretched, and stored as JPEG so a
  photograph does not make the report heavy;
* ``none``: text only.

``art(cover, "pdf" | "pptx" | "docx", data)`` gives the banner and the size to
place it at (slides and Word). The PDF's cover is a full page: ``page_art``
draws the design (or crops the picture) to the whole page and the text sits on
a clear panel (``panel_span``) at ``cover["zone"]``: "top", "middle" or
"bottom". Without matplotlib (built-ins) or Pillow (pictures) there is no art
and the cover is text only; nothing raises.
"""

from __future__ import annotations

import hashlib
import io
import os
from dataclasses import dataclass

import appinfo

DEFAULT_ACCENT = "#2C3E50"
ACCENTS = (("Navy", "#2C3E50"), ("Teal", "#1F7A8C"), ("Green", "#2E7D5B"),
           ("Crimson", "#A23B4A"), ("Amber", "#B7791F"), ("Slate", "#5B6470"))
DEFAULT_DESIGN = "ribbon"

BUILTIN = (("ribbon", "Spectrum ribbon"), ("band", "Band"),
           ("minimal", "Minimal"), ("data", "Your data"),
           ("grid", "Peak map"))
NAMES = dict(BUILTIN)

# where the art goes and how big: (width, height) in mm for the PDF (inside the
# page margins), in inches for a slide (a strip along the bottom) and for the
# Word document (a banner at the top, the same 4:1 shape as the PDF's, sized
# to its default content width); the resolution it is drawn at is chosen so a
# photograph stays light
SIZES = {"pdf": (180.0, 45.0, "mm", 200), "pptx": (13.333, 1.55, "in", 150),
         "docx": (6.5, 1.625, "in", 200)}
THUMB = (240, 60, 100)                       # pixels and dpi of a preview

# the PDF's full-page cover: the page in mm, the resolution the picture is drawn
# at (about 1240 x 1754 px on A4: a photograph stays near 400 KB as JPEG), the
# zones the text panel can sit in and the part of the page (fractions of its
# height, from the bottom) each one covers
PAGE_MM = {"a4": (210.0, 297.0), "letter": (215.9, 279.4)}
PAGE_DPI = 150
ZONES = ("top", "middle", "bottom")
DEFAULT_ZONE = "bottom"
_PANEL = {"top": (0.62, 1.0), "middle": (0.31, 0.69), "bottom": (0.0, 0.38)}
PAGE_THUMB = (90, 127, 45)                  # a page preview for the picker

try:
    from PIL import Image
    HAVE_PIL = True
except ImportError:                                      # pragma: no cover
    HAVE_PIL = False
try:
    import matplotlib                                     # noqa: F401
    HAVE_MPL = True
except ImportError:                                      # pragma: no cover
    HAVE_MPL = False


@dataclass
class Cover:
    id: str                 # what goes in cover["design"]
    name: str               # what the picker shows
    kind: str               # none | builtin | file


@dataclass
class Art:
    data: bytes | None      # the picture; None: no art (text-only cover)
    width: float = 0.0      # size to place it at, in ``unit``
    height: float = 0.0
    unit: str = "mm"
    note: str = ""          # why there is none, or what was substituted
    fmt: str = "png"        # "png" or "jpeg"


# -- colour --------------------------------------------------------------------
def valid_accent(text):
    """``text`` as '#RRGGBB' (upper case), or '' when it is not one."""
    t = str(text or "").strip()
    if len(t) == 7 and t[0] == "#":
        try:
            int(t[1:], 16)
            return t.upper()
        except ValueError:
            pass
    return ""


def _rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def mix(a, b, t):
    """The colour ``t`` of the way from ``a`` to ``b`` (both '#RRGGBB')."""
    ra, rb = _rgb(a), _rgb(b)
    return "#%02X%02X%02X" % tuple(round(x + (y - x) * t)
                                   for x, y in zip(ra, rb))


def tint(accent, t):
    """``accent`` faded towards white: 0 = the accent, 1 = white."""
    return mix(accent, "#FFFFFF", t)


# -- the list ------------------------------------------------------------------
def list_covers(folder=None):
    """Everything the picker can offer, in order: none, the built-in designs
    and the pictures in ``assets/covers``."""
    out = [Cover("none", "No picture", "none")]
    out += [Cover(i, n, "builtin") for i, n in BUILTIN]
    for name, _path in appinfo.cover_images(folder):
        out.append(Cover("file:" + name, os.path.splitext(name)[0], "file"))
    return out


def resolve(cover, folder=None):
    """``(design, picture path or None, note)`` for a cover dict: the design
    actually used (``none`` / a built-in / ``picture``) and, for a picture, its
    file. A picture that is gone or a design that is unknown says so in the
    note and falls back to text only."""
    cover = cover or {}
    design = str(cover.get("design") or "none")
    if design == "none" or design in NAMES:
        return design, None, ""
    if design == "image":
        path = str(cover.get("image") or "")
        return ("picture", path, "") if os.path.isfile(path) else (
            "none", None, f"the cover picture {path or '(none chosen)'} "
                          "was not found")
    if design.startswith("file:"):
        name = design[5:]
        for n, path in appinfo.cover_images(folder):
            if n == name:
                return "picture", path, ""
        return "none", None, f"the cover picture {name} was not found"
    return "none", None, f"unknown cover design {design}"


# -- drawing -------------------------------------------------------------------
def _gauss(x, c, w):
    import numpy as np
    return np.exp(-0.5 * ((x - c) / w) ** 2)


def _tall(w, h):
    """True for a page-shaped picture (a cover page), False for a banner."""
    return h > 1.3 * w


def panel_span(zone):
    """``(low, high)``: the part of the page height (0 = bottom) the text
    panel covers for ``zone``."""
    return _PANEL.get(zone, _PANEL[DEFAULT_ZONE])


def _free_span(zone):
    """The part of a page the panel leaves clear (the larger one, for the
    middle zone, the upper)."""
    lo, hi = panel_span(zone)
    if zone == "bottom":
        return hi + 0.04, 0.97
    if zone == "top":
        return 0.03, lo - 0.04
    return hi + 0.03, 0.97


def _draw_ribbon(ax, w, h, accent, _data, zone=None):
    import numpy as np
    x = np.linspace(0.0, 1.0, 600)
    # peak positions, widths and heights of each layer: fixed, so the design is
    # the same every time
    layers = [((0.16, 0.030, 0.9), (0.34, 0.020, 0.5), (0.61, 0.026, 0.7),
               (0.83, 0.018, 0.4)),
              ((0.12, 0.028, 0.6), (0.29, 0.022, 0.9), (0.55, 0.020, 0.5),
               (0.78, 0.030, 0.8)),
              ((0.21, 0.025, 0.8), (0.42, 0.030, 0.6), (0.66, 0.022, 0.9),
               (0.90, 0.020, 0.5)),
              ((0.10, 0.020, 0.5), (0.37, 0.026, 0.8), (0.58, 0.030, 0.6),
               (0.74, 0.020, 0.7)),
              ((0.18, 0.030, 0.7), (0.47, 0.020, 0.9), (0.70, 0.028, 0.5),
               (0.87, 0.022, 0.6))]
    n = len(layers)
    if _tall(w, h):
        # a page: the layers spread over the whole height and the peaks keep
        # the height they have on the banner (a share of the width, not of
        # the page height)
        lo, hi = 0.04, 0.96
        step, amp = (hi - lo) / n, 0.12 * w / h
    else:
        lo, step, amp = 0.10, 0.145, 0.26
    for k in range(n - 1, -1, -1):                # back (top, pale) to front
        base = lo + step * k
        y = sum(a * _gauss(x, c, s) for c, s, a in layers[k]) * amp
        ax.fill_between(x, base, base + y, color=tint(accent, 0.12 + 0.13 * k),
                        linewidth=0, zorder=2 + (n - k) * 2)
        ax.plot(x, base + y, color=tint(accent, 0.05 + 0.10 * k),
                linewidth=0.9, zorder=3 + (n - k) * 2)


def _draw_band(ax, w, h, accent, _data, zone=None):
    import numpy as np
    ax.add_patch(matplotlib_rect(0, 0, 1, 1, accent))
    tall = _tall(w, h)
    rules = (0.22, 0.42, 0.62, 0.82) if not tall else [
        i / 10 for i in range(1, 10)]
    for y in rules:
        ax.plot([0, 1], [y, y], color=tint(accent, 0.82), linewidth=0.6,
                alpha=0.35, zorder=3)
    x = np.linspace(0.0, 1.0, 500)
    if tall:
        # the silhouette stands on the clear part of the page
        lo, hi = _free_span(zone or DEFAULT_ZONE)
        base, amp = (hi if zone == "top" else lo + 0.02), 0.34 * 0.5 * w / h
        base = lo + 0.02 if zone != "top" else lo + 0.02
    else:
        base, amp = 0.05, 0.34
    y = base + amp * (0.9 * _gauss(x, 0.30, 0.035)
                      + 0.6 * _gauss(x, 0.52, 0.03)
                      + 1.0 * _gauss(x, 0.74, 0.04)
                      + 0.45 * _gauss(x, 0.90, 0.025))
    ax.fill_between(x, base - (0.0 if tall else 0.05), y,
                    color=tint(accent, 0.55), alpha=0.55, linewidth=0, zorder=4)
    ax.plot(x, y, color=tint(accent, 0.75), linewidth=1.0, zorder=5)


def _draw_minimal(ax, w, h, accent, _data, zone=None):
    import numpy as np
    if _tall(w, h):
        lo, hi = _free_span(zone or DEFAULT_ZONE)
        y0, amp = (lo + hi) / 2.0, 0.22 * 0.25 * w / h
    else:
        y0, amp = 0.5, 0.22
    ax.plot([0.0, 1.0], [y0, y0], color=tint(accent, 0.75), linewidth=0.8)
    ax.plot([0.0, 0.13], [y0, y0], color=accent, linewidth=4.0,
            solid_capstyle="butt")
    x = np.linspace(0.84, 1.0, 120)
    y = y0 + amp * (0.8 * _gauss(x, 0.90, 0.012) + _gauss(x, 0.95, 0.010))
    ax.plot(x, y, color=accent, linewidth=1.2)


def _draw_grid(ax, w, h, accent, _data, zone=None):
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    nx = 96
    ny = max(6, round(nx * h / w))
    xs = np.linspace(0, 1, nx)[None, :]
    ys = np.linspace(0, 1, ny)[:, None]
    field = np.zeros((ny, nx))
    blobs = ((0.18, 0.60, 0.07, 0.30, 0.8), (0.42, 0.35, 0.10, 0.35, 1.0),
             (0.66, 0.65, 0.06, 0.28, 0.7), (0.85, 0.40, 0.08, 0.32, 0.9))
    if _tall(w, h):
        # a page: the same peaks repeated in rows, as roundish blobs, each
        # row shifted a little so they do not line up
        rows = 6
        blobs = tuple(((cx + 0.07 * r) % 1.0, (r + 0.5 + 0.3 * (cy - 0.5)) / rows,
                       sx, 0.45 * sy / rows * 2.2, a * (0.55 + 0.45 * ((r * 7) % 5) / 4))
                      for r in range(rows) for cx, cy, sx, sy, a in blobs)
    for cx, cy, sx, sy, a in blobs:
        field += a * np.exp(-0.5 * (((xs - cx) / sx) ** 2
                                    + ((ys - cy) / sy) ** 2))
    rng = np.random.default_rng(7)                # fixed: the same every time
    field += 0.05 * rng.random((ny, nx))
    cmap = LinearSegmentedColormap.from_list("cover", [tint(accent, 0.94),
                                                       tint(accent, 0.45),
                                                       accent])
    ax.imshow(field, cmap=cmap, extent=(0, 1, 0, 1), origin="lower",
              interpolation="nearest", aspect="auto", vmin=0,
              vmax=float(field.max()))


def _draw_data(ax, w, h, accent, data, zone=None):
    import numpy as np
    try:
        e = np.asarray(data[0], dtype=float)
        c = np.asarray(data[1], dtype=float)
        ok = e.size == c.size and e.size >= 3 and np.isfinite(c).all()
    except (TypeError, ValueError, IndexError):
        ok = False
    if not ok:
        return _draw_ribbon(ax, w, h, accent, None, zone)
    step = max(1, e.size // 900)
    e, c = e[::step], c[::step]
    span = float(e.max() - e.min()) or 1.0
    x = 0.03 + 0.94 * (float(e.max()) - e) / span      # binding energy: high at left
    order = np.argsort(x)
    x, c = x[order], c[order]
    lo, hi = float(c.min()), float(c.max())
    norm = (c - lo) / ((hi - lo) or 1.0)
    if _tall(w, h):
        # the spectrum stands in the clear part of the page, as tall as the
        # banner's is wide-proportioned (a share of the width)
        flo, fhi = _free_span(zone or DEFAULT_ZONE)
        base = flo + 0.10 * (fhi - flo)
        y = base + min(0.72 * (fhi - flo), 0.42 * w / h) * norm
    else:
        base = 0.14
        y = base + 0.72 * norm
    ax.fill_between(x, base, y, color=tint(accent, 0.72), linewidth=0, zorder=2)
    ax.plot(x, y, color=accent, linewidth=1.3, zorder=3)
    ax.plot([0.03, 0.97], [base, base], color=tint(accent, 0.5), linewidth=0.8,
            zorder=3)


def matplotlib_rect(x, y, w, h, colour):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, facecolor=colour, edgecolor="none",
                     zorder=1)


_DRAW = {"ribbon": _draw_ribbon, "band": _draw_band, "minimal": _draw_minimal,
         "grid": _draw_grid, "data": _draw_data}


def _render_builtin(design, px, dpi, accent, data, zone=None):
    if not HAVE_MPL:
        return None
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        w, h = px
        fig = Figure(figsize=(w / dpi, h / dpi), dpi=dpi, facecolor="white")
        FigureCanvasAgg(fig)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        _DRAW[design](ax, w, h, accent, data, zone)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
        return buf.getvalue()
    except Exception:                                    # noqa: BLE001
        return None                # art is a nicety: never stop a report


def _render_picture(path, px, fmt="jpeg"):
    """The picture at ``path`` scaled to cover ``px`` and cropped to it from
    the middle, as JPEG (or PNG); None when it cannot be read."""
    if not HAVE_PIL:
        return None
    try:
        with Image.open(path) as im:
            im = im.convert("RGB")
            w, h = px
            scale = max(w / im.width, h / im.height)
            size = (max(w, round(im.width * scale)),
                    max(h, round(im.height * scale)))
            im = im.resize(size, Image.LANCZOS)
            left, top = (size[0] - w) // 2, (size[1] - h) // 2
            im = im.crop((left, top, left + w, top + h))
            buf = io.BytesIO()
            if fmt == "jpeg":
                im.save(buf, format="JPEG", quality=90, optimize=True)
            else:
                im.save(buf, format="PNG")
            return buf.getvalue()
    except Exception:                                    # noqa: BLE001
        return None


_cache = {}


def _data_key(data):
    if not data:
        return ""
    try:
        e, c = list(data[0]), list(data[1])
        stride = max(1, len(e) // 64)
        return hashlib.sha1(repr((len(e), e[::stride], c[::stride]))
                            .encode()).hexdigest()
    except Exception:                                    # noqa: BLE001
        return ""


def render(cover, px, dpi, data=None, folder=None, jpeg=False, zone=None):
    """``(picture bytes or None, note)`` for ``cover`` at ``px`` = (w, h)
    pixels. A picture of the user's is JPEG when ``jpeg`` (report and slides)
    and PNG otherwise (the picker's previews, which Tk can show); a built-in
    design is always PNG."""
    design, path, note = resolve(cover, folder)
    if design == "none":
        return None, note
    accent = valid_accent((cover or {}).get("accent")) or DEFAULT_ACCENT
    mtime = ""
    if path:
        try:
            mtime = f"{os.path.getmtime(path):.0f}:{os.path.getsize(path)}"
        except OSError:
            mtime = ""
    key = (design, path, mtime, px, dpi, accent, bool(jpeg),
           _data_key(data) if design == "data" else "", zone)
    if key in _cache:
        return _cache[key], note
    if design == "picture":
        png = _render_picture(path, px, "jpeg" if jpeg else "png")
        if png is None:
            return None, note or f"the picture {os.path.basename(path)} " \
                                 "could not be read"
    else:
        png = _render_builtin(design, px, dpi, accent, data, zone)
        if png is None:
            return None, note or "the cover picture needs matplotlib"
    if len(_cache) > 48:
        _cache.pop(next(iter(_cache)))
    _cache[key] = png
    return png, note


def art(cover, kind, data=None, folder=None):
    """The cover picture for the ``"pdf"`` or ``"pptx"`` cover, with the size
    (in ``unit``) to place it at. ``data`` = ``(energy, counts)`` for the
    ``data`` design. ``Art.png`` is None when there is nothing to show."""
    width, height, unit, dpi = SIZES[kind]
    inches = (width / 25.4, height / 25.4) if unit == "mm" else (width, height)
    px = (round(inches[0] * dpi), round(inches[1] * dpi))
    data_, note = render(cover, px, dpi, data, folder, jpeg=True)
    fmt = "jpeg" if data_ and data_[:2] == b"\xff\xd8" else "png"
    return Art(data_, width, height, unit, note, fmt)


def page_art(cover, page="a4", data=None, folder=None):
    """The cover picture for the PDF's full-page cover: the design drawn (or
    the picture cropped from its middle) to the whole ``page`` ("a4" or
    "letter"), with the part opposite the text panel kept for the artwork
    (``cover["zone"]``). ``Art.width`` / ``height`` are the page in mm;
    ``Art.data`` is None when there is nothing to show (``Art.note`` says why)."""
    width, height = PAGE_MM.get(page, PAGE_MM["a4"])
    px = (round(width / 25.4 * PAGE_DPI), round(height / 25.4 * PAGE_DPI))
    zone = (cover or {}).get("zone")
    zone = zone if zone in ZONES else DEFAULT_ZONE
    data_, note = render(cover, px, PAGE_DPI, data, folder, jpeg=True,
                         zone=zone)
    fmt = "jpeg" if data_ and data_[:2] == b"\xff\xd8" else "png"
    return Art(data_, width, height, "mm", note, fmt)


def page_thumbnail(cover, data=None, folder=None):
    """A small portrait PNG (about 100 x 141) of the full-page cover for the
    picker, the text panel shown as a pale band with grey lines at the chosen
    zone; None for no picture."""
    w, h, dpi = PAGE_THUMB
    zone = (cover or {}).get("zone")
    zone = zone if zone in ZONES else DEFAULT_ZONE
    png, _note = render(cover, (w, h), dpi, data, folder, zone=zone)
    if png is None or not HAVE_PIL:
        return png
    try:
        im = Image.open(io.BytesIO(png)).convert("RGBA")
        lo, hi = panel_span(zone)
        y0, y1 = round(h * (1 - hi)), round(h * (1 - lo))
        over = Image.new("RGBA", im.size, (0, 0, 0, 0))
        from PIL import ImageDraw
        d = ImageDraw.Draw(over)
        d.rectangle((0, y0, w, y1), fill=(255, 255, 255, 215))
        top = y0 + (y1 - y0) // 4
        d.rectangle((8, top, 8 + w // 3, top + 2),
                    fill=_rgb(valid_accent((cover or {}).get("accent"))
                              or DEFAULT_ACCENT) + (255,))
        for k, frac in enumerate((0.8, 0.55, 0.65)):
            yy = top + 9 + k * 6
            d.rectangle((8, yy, 8 + round((w - 16) * frac), yy + 2),
                        fill=(150, 156, 164, 255))
        out = io.BytesIO()
        Image.alpha_composite(im, over).convert("RGB").save(out, format="PNG")
        return out.getvalue()
    except Exception:                                    # noqa: BLE001
        return png


def thumbnail(cover, data=None, folder=None):
    """A small PNG (240 x 60) of ``cover`` for the picker, or None."""
    w, h, dpi = THUMB
    return render(cover, (w, h), dpi, data, folder)[0]
