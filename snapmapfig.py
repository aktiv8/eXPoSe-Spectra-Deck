"""A SnapMap view as data, and the one function that draws it.

Tk-free (matplotlib and numpy are imported inside the drawing functions). A
*view* is a plain JSON-able dict, so it can be stored as a workbook figure and
drawn again for the dialog, a picture file, the PDF, the slides, Word and the
HTML page with the same result::

    {"kind": "snapmap", "ref": <workbook.region_ref of the map's region>,
     "window": [lo, hi], "background": bool, "scale": "Viridis",
     "overlay": {"on": bool, "alpha": 0.65},
     "rois": [{"shape": "rect" | "ellipse", "rect": [x0, y0, x1, y1],
               "label": str}, ...],           # µm, the map's own frame
     "smooth": {"method": "None", "strength": 0.5},
     "show_spectrum": bool, "scalebar": bool,
     "mode": "single" | "grid" | "composite" | "ratio",
     "layers": [{"ref": <region_ref>, "window": [lo, hi] | None}, ...]}

``mode`` other than single lays several elements of one map file on the one
grid: *grid* = every layer as its own small map, *composite* = the first three
as red, green and blue, *ratio* = the first divided by the second.

Nothing here edits the cube: the window, background and ROI are applied on the
way out, as everywhere else in the app.
"""

from __future__ import annotations

import math

import snapmap
import snapshot

KIND = "snapmap"

DEFAULTS = {
    "kind": KIND, "ref": None, "window": None, "background": False,
    "scale": "Viridis", "overlay": {"on": False, "alpha": 0.65}, "rois": [],
    "smooth": {"method": "None", "strength": 0.5}, "show_spectrum": True,
    "scalebar": True, "tiles": [], "mode": "single", "layers": [],
}
MODES = ("single", "grid", "composite", "ratio")
MODE_NAMES = {"single": "Single map", "grid": "Element grid",
              "composite": "RGB composite", "ratio": "Ratio map"}
MAX_LAYERS = 12
RATIO_FLOOR = 0.05      # a denominator under this share of its 99th percentile is blank


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) \
        and math.isfinite(v)


def sanitise(spec) -> dict:
    """Any input as a valid view: unknown keys dropped, bad values reset."""
    spec = spec if isinstance(spec, dict) else {}
    out = {"kind": KIND}
    ref = spec.get("ref")
    out["ref"] = dict(ref) if isinstance(ref, dict) else None
    w = spec.get("window")
    out["window"] = ([float(w[0]), float(w[1])]
                     if isinstance(w, (list, tuple)) and len(w) == 2
                     and all(_num(v) for v in w) else None)
    out["background"] = bool(spec.get("background", False))
    out["scale"] = (spec["scale"] if isinstance(spec.get("scale"), str)
                    and spec["scale"] else DEFAULTS["scale"])
    ov = spec.get("overlay") if isinstance(spec.get("overlay"), dict) else {}
    a = ov.get("alpha", 0.65)
    out["overlay"] = {"on": bool(ov.get("on", False)),
                      "alpha": float(min(1.0, max(0.1, a))) if _num(a) else 0.65}
    out["rois"] = sanitise_rois(spec.get("rois"), spec.get("roi"))
    sm = spec.get("smooth") if isinstance(spec.get("smooth"), dict) else {}
    s = sm.get("strength", 0.5)
    out["smooth"] = {"method": sm.get("method") if isinstance(
        sm.get("method"), str) else "None",
        "strength": float(min(1.0, max(0.0, s))) if _num(s) else 0.5}
    out["show_spectrum"] = bool(spec.get("show_spectrum", True))
    out["scalebar"] = bool(spec.get("scalebar", True))
    out["mode"] = spec.get("mode") if spec.get("mode") in MODES else "single"
    out["layers"] = sanitise_layers(spec.get("layers"))
    tl = spec.get("tiles")
    out["tiles"] = ([dict(t) for t in tl if isinstance(t, dict)]
                    if isinstance(tl, (list, tuple)) else [])
    return out


def sanitise_layers(layers) -> list:
    """``[{"ref", "window"}]`` with a reference dict each, a window of two
    numbers or None (the map's own default), at most ``MAX_LAYERS``."""
    out = []
    for it in layers if isinstance(layers, (list, tuple)) else []:
        ref = it.get("ref") if isinstance(it, dict) else None
        if not isinstance(ref, dict):
            continue
        w = it.get("window")
        out.append({"ref": dict(ref),
                    "window": ([float(w[0]), float(w[1])]
                               if isinstance(w, (list, tuple)) and len(w) == 2
                               and all(_num(v) for v in w) else None)})
        if len(out) == MAX_LAYERS:
            break
    return out


MAX_ROIS = 6
ROI_COLOURS = ("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7",
               "#0072B2")
SHAPES = ("rect", "ellipse")


def roi_colour(i: int) -> str:
    return ROI_COLOURS[i % len(ROI_COLOURS)]


def new_roi(rect, shape="rect", label="") -> dict:
    return {"shape": shape if shape in SHAPES else "rect",
            "rect": [float(v) for v in rect], "label": label}


def sanitise_rois(rois, legacy=None) -> list:
    """A valid list of areas: each ``{"shape", "rect", "label"}``; a
    ``legacy`` single rectangle (earlier views stored ``roi``) becomes the
    first one. At most ``MAX_ROIS``; labels default to "Area N"."""
    items = list(rois) if isinstance(rois, (list, tuple)) else []
    if not items and isinstance(legacy, (list, tuple)):
        items = [{"rect": legacy}]
    out = []
    for it in items:
        rect = it.get("rect") if isinstance(it, dict) else None
        if not (isinstance(rect, (list, tuple)) and len(rect) == 4
                and all(_num(v) for v in rect)):
            continue
        label = it.get("label")
        out.append(new_roi(rect, it.get("shape"),
                           label.strip() if isinstance(label, str)
                           and label.strip() else f"Area {len(out) + 1}"))
        if len(out) == MAX_ROIS:
            break
    return out


def is_view(state) -> bool:
    """True for the state of a saved figure that is a SnapMap view."""
    return isinstance(state, dict) and state.get("kind") == KIND


def roi_rect_of_pixel(cube, ix, iy):
    """The ROI rectangle (µm) that selects exactly one pixel."""
    hx, hy = abs(cube.dx) / 2, abs(cube.dy) / 2
    x, y = cube.x_of(ix), cube.y_of(iy)
    return [x - hx * 0.9, y - hy * 0.9, x + hx * 0.9, y + hy * 0.9]


def roi_mask(cube, roi):
    """Boolean ``(ny, nx)`` mask of one area (a ``new_roi`` dict, or a bare
    rectangle), or None when it holds no pixel. An ellipse is the one inscribed
    in its rectangle, taken by pixel centres."""
    if not roi:
        return None
    if not isinstance(roi, dict):
        roi = new_roi(roi)
    x0, y0, x1, y1 = roi["rect"]
    m = cube.rect_mask(x0, y0, x1, y1)
    if roi.get("shape") == "ellipse" and m.any():
        import numpy as np
        xs = cube.x0 + np.arange(cube.nx) * cube.dx
        ys = cube.y0 + np.arange(cube.ny) * cube.dy
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        rx, ry = abs(x1 - x0) / 2, abs(y1 - y0) / 2
        if rx > 0 and ry > 0:
            inside = (((xs[None, :] - cx) / rx) ** 2
                      + ((ys[:, None] - cy) / ry) ** 2) <= 1.0
            m = m & inside
        # a very small ellipse must still hold the pixel it was clicked on
        if not m.any():
            m = cube.rect_mask(cx - abs(cube.dx) / 2, cy - abs(cube.dy) / 2,
                               cx + abs(cube.dx) / 2, cy + abs(cube.dy) / 2)
    return m if m.any() else None


def roi_masks(cube, rois):
    """``[mask or None]`` for each area."""
    return [roi_mask(cube, r) for r in rois]


def roi_stats(cube, window, background, rois, image=None) -> list:
    """One row per area, from the map as shown (the counts summed over the
    window, or the given ``image``): ``label, shape, px, area_um2, mean, sum,
    sd`` and ``ratio``, the area's mean over the whole map's. An area with no
    pixel is left out; blank (NaN) pixels are not counted."""
    import numpy as np
    if image is None:
        image = cube.image(window[0], window[1], background=background)
    image = np.asarray(image, dtype=float)
    whole = float(np.nanmean(image)) if np.isfinite(image).any() else 0.0
    rows = []
    for r in rois:
        m = roi_mask(cube, r)
        if m is None:
            continue
        v = image[m]
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        mean = float(v.mean())
        rows.append({
            "label": r.get("label", ""), "shape": r.get("shape", "rect"),
            "px": int(v.size),
            "area_um2": float(v.size * abs(cube.dx * cube.dy)),
            "mean": mean, "sum": float(v.sum()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "ratio": mean / whole if whole else float("nan")})
    return rows


def nice_length(width_um: float) -> float:
    """A round scale-bar length (1, 2 or 5 x 10^n µm) of about a fifth of the
    map's width."""
    if width_um <= 0:
        return 1.0
    target = width_um / 5.0
    mag = 10 ** math.floor(math.log10(target))
    for m in (5, 2, 1):
        if m * mag <= target:
            return float(m * mag)
    return float(mag)


def mean_spectra(cube, masks, smooth):
    """``(whole, [spectrum or None per mask])``: the mean spectra per pixel of
    the whole map and of each area, smoothed for display when a method is
    chosen."""
    n_all = cube.nx * cube.ny
    whole = [v / n_all for v in cube.total()]
    rois = []
    for mask in masks:
        if mask is not None and mask.any():
            n = int(mask.sum())
            rois.append([v / n for v in cube.roi_spectrum(mask)])
        else:
            rois.append(None)
    method = (smooth or {}).get("method", "None")
    if method != "None":
        import smoothing
        strength = (smooth or {}).get("strength", 0.5)
        whole = smoothing.smooth(whole, method, strength).tolist()
        rois = [None if r is None else
                smoothing.smooth(r, method, strength).tolist() for r in rois]
    return whole, rois


def describe(spec, region=None) -> str:
    """One line saying what a view shows (speaker notes, figure list)."""
    s = sanitise(spec)
    name = region.name if region is not None else (s["ref"] or {}).get("name", "")
    parts = [f"SnapMap of {name}".strip()]
    if s["window"]:
        lo, hi = s["window"]
        parts.append(f"energy window {min(lo, hi):.1f}-{max(lo, hi):.1f}")
    if s["background"]:
        parts.append("background removed")
    parts.append(f"{s['scale']} colours")
    if s["mode"] != "single" and s["layers"]:
        parts.append(f"{MODE_NAMES[s['mode']].lower()} of "
                     f"{len(s['layers'])} elements")
    if s["rois"]:
        n = len(s["rois"])
        parts.append("1 area of interest marked" if n == 1
                     else f"{n} areas of interest marked")
    if s["overlay"]["on"]:
        parts.append("over the camera image")
    return "; ".join(parts) + "."


def draw_view(fig, rect, spec, region, pal, cmap, photo=None,
              unit_mark=None, tiles=None, layers=None):
    """Draw a view onto ``fig`` inside ``rect`` = ``(left, bottom, right, top)``
    figure fractions (None = the whole figure) and return its axes
    ``(ax_map, ax_spec or None)``.

    ``region`` carries the cube (``extra["cube"]``) and the energy labels;
    ``pal`` needs ``plot_bg``, ``plot_fg``, ``muted`` and ``accent``;
    ``photo`` is ``(PIL image, calibration)`` for the camera overlay.
    ``unit_mark`` is appended to the scale bar (e.g. "(approx.)").
    ``layers`` = ``[(region, window)]`` for the modes other than single.
    ``tiles`` are further regions whose maps are laid beside this one by their
    stage positions (one canvas, one window, one colour range); the camera
    overlay and the area are not used then."""
    from matplotlib.patches import Rectangle
    s = sanitise(spec)
    cube = region.extra["cube"]
    window = tuple(s["window"]) if s["window"] else snapmap.default_window(cube)
    if s["mode"] != "single" and layers and not tiles:
        out = draw_layers(fig, rect, s, region, layers, pal, cmap, unit_mark)
        if out is not None:
            return out
    l, b, r, t = rect if rect else (0.0, 0.0, 1.0, 1.0)
    w, h = r - l, t - b
    show_spec = s["show_spectrum"]
    if show_spec:
        ax_map = fig.add_axes([l + 0.07 * w, b + 0.12 * h, 0.36 * w, 0.76 * h])
        ax_spec = fig.add_axes([l + 0.58 * w, b + 0.12 * h, 0.39 * w, 0.76 * h])
    else:
        ax_map = fig.add_axes([l + 0.10 * w, b + 0.10 * h, 0.72 * w, 0.78 * h])
        ax_spec = None
    for ax in (ax_map, ax_spec):
        if ax is None:
            continue
        ax.set_facecolor(pal["plot_bg"])
        ax.tick_params(colors=pal["muted"], labelsize=8)
        for sp in ax.spines.values():
            sp.set_color(pal["muted"])

    left, right, bottom, top = cube.extent()
    mos = None
    if tiles:
        mos = mosaic_image([region.extra["cube"]]
                           + [t.extra["cube"] for t in tiles], window,
                           s["background"])
        if mos is not None:
            left, right, bottom, top = mos.extent
    use_photo = s["overlay"]["on"] and photo is not None and mos is None
    if use_photo:
        img, cal = photo
        cx, cy = snapshot.stage_to_pixel(cal, cube.stage_x_mm, cube.stage_y_mm)
        ux, uy = cal["um_per_px_x"], cal["um_per_px_y"]
        pw, ph = img.size
        ax_map.imshow(img, extent=[-cx * ux, (pw - cx) * ux,
                                   (ph - cy) * uy, -cy * uy], zorder=0)
    if mos is not None:
        image, extent = mos.image, mos.extent
        try:
            cmap = cmap.with_extremes(bad=pal["plot_bg"])
        except Exception:
            pass
    else:
        image = cube.image(window[0], window[1], background=s["background"])
        extent = cube.extent()
    vmin, vmax = snapmap.colour_range(image)
    im = ax_map.imshow(image, extent=extent, origin="upper", cmap=cmap,
                       interpolation="nearest", vmin=vmin, vmax=vmax, zorder=1,
                       alpha=s["overlay"]["alpha"] if use_photo else 1.0)
    if use_photo:
        padx, pady = (right - left) * 0.3, (bottom - top) * 0.3
        ax_map.set_xlim(left - padx, right + padx)
        ax_map.set_ylim(bottom + pady, top - pady)
    else:
        ax_map.set_xlim(left, right)
        ax_map.set_ylim(bottom, top)
    ax_map.set_xlabel("X (µm)", color=pal["plot_fg"], fontsize=9)
    ax_map.set_ylabel("Y (µm)", color=pal["plot_fg"], fontsize=9)
    lo, hi = window
    ax_map.set_title(f"{region.name}   {min(lo, hi):.1f}–{max(lo, hi):.1f} "
                     f"{region.energy_units}", fontsize=9, color=pal["plot_fg"])
    cb = fig.colorbar(im, ax=ax_map, fraction=0.046, pad=0.02)
    cb.ax.tick_params(colors=pal["muted"], labelsize=8)
    cb.set_label(cube.label or "Counts", color=pal["plot_fg"], fontsize=8)

    masks = []
    if mos is not None:
        for _i, fl, ft, fw, fh in mosaic_footprints(
                [cube] + [t.extra["cube"] for t in tiles]):
            ax_map.add_patch(Rectangle((fl, ft), fw, fh, fill=False, lw=0.8,
                                       ls=":", ec=pal["muted"], zorder=4))
    else:
        rois = s["rois"]
        masks = roi_masks(cube, rois)
        for k, (roi, m) in enumerate(zip(rois, masks)):
            if m is None:
                continue
            _outline_roi(ax_map, cube, roi, m, roi_colour(k))
    if s["scalebar"]:
        _scalebar(ax_map, abs(right - left), unit_mark)

    if ax_spec is not None:
        whole, spectra = mean_spectra(cube, masks, s["smooth"])
        if mos is not None:
            whole, _n = mean_spectrum_of(mos.used)
            spectra = []
        e = cube.energy
        shown = [(k, sp) for k, sp in enumerate(spectra) if sp is not None]
        ax_spec.plot(e, whole, color=pal["muted"] if shown else pal["accent"],
                     lw=1.4, label="whole map")
        for k, sp in shown:
            ax_spec.plot(e, sp, color=roi_colour(k), lw=1.6,
                         label=f"{s['rois'][k]['label']}, "
                               f"{int(masks[k].sum())} px")
        if shown:
            ax_spec.legend(fontsize=8, frameon=False, labelcolor=pal["plot_fg"])
        binding = "inding" in (region.energy_label or "")
        ax_spec.set_xlim((max(e), min(e)) if binding else (min(e), max(e)))
        ax_spec.axvspan(lo, hi, color=pal["accent"], alpha=0.14, lw=0)
        ax_spec.set_xlabel(f"{region.energy_label} ({region.energy_units})",
                           color=pal["plot_fg"], fontsize=9)
        ax_spec.set_ylabel("Counts per pixel", color=pal["plot_fg"], fontsize=9)
        ax_spec.set_title(region.name, fontsize=9, color=pal["plot_fg"])
    return ax_map, ax_spec


def _scalebar(ax, width_um, unit_mark=None):
    """A white bar of a round length near the bottom-left of the map."""
    length = nice_length(width_um)
    ymin, ymax = ax.get_ylim()
    by = ymin + (ymax - ymin) * 0.06
    xmin, xmax = ax.get_xlim()
    bx = xmin + (xmax - xmin) * 0.05
    ax.plot([bx, bx + length], [by, by], color="white", lw=3,
            solid_capstyle="butt", zorder=6, path_effects=_outline())
    txt = f"{length:g} µm" + (f" {unit_mark}" if unit_mark else "")
    ax.text(bx + length / 2, by + (ymax - ymin) * 0.025, txt, color="white",
            fontsize=8, ha="center", va="bottom", zorder=6,
            path_effects=_outline())


def _outline_roi(ax, cube, roi, mask, colour):
    """Draw one area on the map: its rectangle (or ellipse) snapped to the
    pixels it holds, and its label."""
    import numpy as np
    from matplotlib.patches import Ellipse, Rectangle
    rows, cols = np.flatnonzero(mask.any(axis=1)), np.flatnonzero(
        mask.any(axis=0))
    x0, x1 = cube.x_of(cols[0]) - cube.dx / 2, cube.x_of(cols[-1]) + cube.dx / 2
    y0, y1 = cube.y_of(rows[0]) - cube.dy / 2, cube.y_of(rows[-1]) + cube.dy / 2
    if roi.get("shape") == "ellipse":
        ax.add_patch(Ellipse(((x0 + x1) / 2, (y0 + y1) / 2), abs(x1 - x0),
                             abs(y1 - y0), fill=False, lw=1.8, ls="--",
                             ec=colour, zorder=5))
    else:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, lw=1.8,
                               ls="--", ec=colour, zorder=5))
    ax.text(min(x0, x1), min(y0, y1), " " + roi.get("label", ""), color=colour,
            fontsize=8, va="bottom", ha="left", zorder=6,
            path_effects=_outline())


def _outline():
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=1.6, foreground="black")]


# -- several maps laid side by side on one canvas ---------------------------------------
class Mosaic:
    """``image`` (NaN where no map reaches), its ``extent`` (µm, the common
    frame), the ``used`` cubes, ``notes`` and the ``gap`` flag."""

    def __init__(self, image, extent, used, notes, overlap_px):
        self.image, self.extent, self.used = image, extent, used
        self.notes, self.overlap_px = notes, overlap_px


def tile_offset(cube, ref_cube):
    """Where ``cube``'s centre lies from ``ref_cube``'s, in µm of the common
    (image) frame: image x runs against stage X, image y with stage Y, as for a
    camera picture (``snapshot.X_DIR`` / ``Y_DIR``)."""
    return (snapshot.X_DIR * (cube.stage_x_mm - ref_cube.stage_x_mm) * 1000.0,
            snapshot.Y_DIR * (cube.stage_y_mm - ref_cube.stage_y_mm) * 1000.0)


def mosaic_image(cubes, window, background=False):
    """Lay ``cubes`` (the first is the reference) on one canvas by their stage
    positions and sum each over the same energy ``window``.

    The canvas pixel is the finest of the maps'; each canvas pixel takes the
    pixel of a map that covers it (the mean where maps overlap) by nearest
    neighbour, never interpolating. A map with no stage position or a zero
    pixel size is left out, with a note."""
    import numpy as np
    notes, tiles = [], []
    ref = cubes[0]
    for c in cubes:
        if c.stage_x_mm is None or c.stage_y_mm is None:
            notes.append("a map with no stage position was left out")
            continue
        if not c.dx or not c.dy:
            notes.append("a map with no pixel size was left out")
            continue
        ox, oy = tile_offset(c, ref)
        l, r, b, t = c.extent()
        tiles.append((c, ox, oy, min(l, r) + ox, max(l, r) + ox,
                      min(t, b) + oy, max(t, b) + oy))
    if not tiles:
        return None
    gdx = min(abs(c.dx) for c, *_ in tiles)
    gdy = min(abs(c.dy) for c, *_ in tiles)
    x0 = min(t[3] for t in tiles)
    x1 = max(t[4] for t in tiles)
    y0 = min(t[5] for t in tiles)
    y1 = max(t[6] for t in tiles)
    nx = max(1, int(round((x1 - x0) / gdx)))
    ny = max(1, int(round((y1 - y0) / gdy)))
    gx = x0 + (np.arange(nx) + 0.5) * gdx
    gy = y0 + (np.arange(ny) + 0.5) * gdy
    total = np.zeros((ny, nx))
    count = np.zeros((ny, nx))
    for c, ox, oy, tl, tr, tt, tb in tiles:
        img = c.image(window[0], window[1], background=background)
        ix = np.floor((gx - tl) / abs(c.dx)).astype(int)
        iy = np.floor((gy - tt) / abs(c.dy)).astype(int)
        okx = (ix >= 0) & (ix < c.nx)
        oky = (iy >= 0) & (iy < c.ny)
        if not okx.any() or not oky.any():
            continue
        sub = img[np.ix_(np.clip(iy, 0, c.ny - 1), np.clip(ix, 0, c.nx - 1))]
        cover = oky[:, None] & okx[None, :]
        total += np.where(cover, sub, 0.0)
        count += cover
    overlap = int((count > 1).sum())
    out = np.where(count > 0, total / np.maximum(count, 1), np.nan)
    if overlap:
        notes.append(f"where maps overlap ({overlap} px) the mean is shown")
    if (count == 0).any():
        notes.append("blank where no map reaches")
    used = [t[0] for t in tiles]
    return Mosaic(out, (x0, x1, y1, y0), used, notes, overlap)


def mosaic_footprints(cubes):
    """``[(label_index, left, top, width, height)]`` of each map in the
    reference frame, for outlining the tiles."""
    out = []
    ref = cubes[0]
    for i, c in enumerate(cubes):
        if c.stage_x_mm is None or c.stage_y_mm is None:
            continue
        ox, oy = tile_offset(c, ref)
        l, r, b, t = c.extent()
        out.append((i, min(l, r) + ox, min(t, b) + oy, abs(r - l), abs(b - t)))
    return out


def mean_spectrum_of(cubes):
    """Mean counts per pixel over several maps when they share one energy
    axis, else the first map's."""
    import numpy as np
    first = cubes[0]
    same = [c for c in cubes if len(c.energy) == len(first.energy)
            and np.allclose(c.energy, first.energy)]
    n = sum(c.nx * c.ny for c in same)
    tot = np.zeros(len(first.energy))
    for c in same:
        tot += np.asarray(c.total())
    return (tot / n).tolist(), len(same)


# -- several elements of one map file on one grid -----------------------------------------
def same_grid(a, b) -> bool:
    """True when two maps have the very same pixels (what the elements of one
    SnapMap file have), so images can be combined pixel by pixel."""
    return (a.nx, a.ny) == (b.nx, b.ny) and all(
        abs(x - y) <= 1e-6 * max(1.0, abs(x))
        for x, y in ((a.x0, b.x0), (a.dx, b.dx), (a.y0, b.y0), (a.dy, b.dy)))


def ratio_image(num, den, floor=RATIO_FLOOR):
    """``num / den`` pixel by pixel, NaN where the denominator is not
    positive or is under ``floor`` of its own 99th percentile (a ratio of
    noise is not a map)."""
    import numpy as np
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    ok = den > 0
    if ok.any():
        ok &= den >= floor * float(np.percentile(den[ok], 99))
    out = np.full(num.shape, np.nan)
    out[ok] = num[ok] / den[ok]
    return out


def unit_range(image, low=1.0, high=99.0):
    """``image`` scaled 0..1 between its ``low`` and ``high`` percentiles."""
    import numpy as np
    vmin, vmax = snapmap.colour_range(image, low, high)
    return np.clip((np.asarray(image, dtype=float) - vmin) / (vmax - vmin),
                   0.0, 1.0)


def rgb_image(images):
    """``(ny, nx, 3)`` from up to three images (red, green, blue), each
    scaled on its own 1st-99th percentile range; a missing channel is 0."""
    import numpy as np
    first = np.asarray(images[0])
    out = np.zeros(first.shape + (3,))
    for k, im in enumerate(images[:3]):
        out[..., k] = unit_range(im)
    return out


def layer_image(layer, background):
    """``(cube, window, image)`` of one ``(region, window)`` layer."""
    rg, win = layer
    cube = rg.extra["cube"]
    win = tuple(win) if win else snapmap.default_window(cube)
    return cube, win, cube.image(win[0], win[1], background=background)


def mode_image(spec, region, layers):
    """The image a ratio view shows (numerator over denominator), else None.
    Layers that are not on the region's pixels are ignored."""
    s = sanitise(spec)
    cube = region.extra["cube"]
    use = [ly for ly in layers if same_grid(cube, ly[0].extra["cube"])]
    if s["mode"] != "ratio" or len(use) < 2:
        return None
    a = layer_image(use[0], s["background"])[2]
    b = layer_image(use[1], s["background"])[2]
    return ratio_image(a, b)


def _map_axes(fig, box, pal):
    ax = fig.add_axes(box)
    ax.set_facecolor(pal["plot_bg"])
    ax.tick_params(colors=pal["muted"], labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(pal["muted"])
    return ax


def _rois_on(ax, cube, rois):
    for k, (roi, m) in enumerate(zip(rois, roi_masks(cube, rois))):
        if m is not None:
            _outline_roi(ax, cube, roi, m, roi_colour(k))


def draw_layers(fig, rect, s, region, layers, pal, cmap, unit_mark=None):
    """The multi-element modes of :func:`draw_view`; returns ``(first axes,
    None)``, or None when there are too few usable layers (the caller then
    draws the single map)."""
    import numpy as np
    from matplotlib.lines import Line2D
    cube = region.extra["cube"]
    use = [ly for ly in layers if same_grid(cube, ly[0].extra["cube"])]
    mode = s["mode"]
    need = {"grid": 2, "composite": 2, "ratio": 2}[mode]
    if len(use) < need:
        return None
    use = use[:{"grid": MAX_LAYERS, "composite": 3, "ratio": 2}[mode]]
    l, b, r, t = rect if rect else (0.0, 0.0, 1.0, 1.0)
    w, h = r - l, t - b
    ext = cube.extent()
    left, right, bottom, top = ext
    bg = s["background"]
    rois = s["rois"]

    if mode == "grid":
        n = len(use)
        cols = int(np.ceil(np.sqrt(n)))
        rows = int(np.ceil(n / cols))
        cw, ch = 0.94 * w / cols, 0.90 * h / rows
        first = None
        for k, ly in enumerate(use):
            cx, cy = k % cols, k // cols
            box = [l + 0.05 * w + cx * cw, b + h * 0.95 - (cy + 1) * ch + 0.10 * ch,
                   cw * 0.72, ch * 0.66]
            ax = _map_axes(fig, box, pal)
            _c, win, img = layer_image(ly, bg)
            vmin, vmax = snapmap.colour_range(img)
            im = ax.imshow(img, extent=ext, origin="upper", cmap=cmap,
                           interpolation="nearest", vmin=vmin, vmax=vmax)
            ax.set_xlim(left, right)
            ax.set_ylim(bottom, top)
            ax.set_title(f"{ly[0].name}  {min(win):.1f}–{max(win):.1f}",
                         fontsize=8, color=pal["plot_fg"])
            ax.tick_params(labelsize=6)
            cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
            cb.ax.tick_params(colors=pal["muted"], labelsize=6)
            _rois_on(ax, cube, rois)
            if k == 0 and s["scalebar"]:
                _scalebar(ax, abs(right - left), unit_mark)
            first = first or ax
        return first, None

    ax = _map_axes(fig, [l + 0.10 * w, b + 0.10 * h, 0.72 * w, 0.78 * h], pal)
    ax.set_xlabel("X (µm)", color=pal["plot_fg"], fontsize=9)
    ax.set_ylabel("Y (µm)", color=pal["plot_fg"], fontsize=9)
    if mode == "composite":
        imgs, handles = [], []
        for k, ly in enumerate(use):
            _c, win, img = layer_image(ly, bg)
            imgs.append(img)
            handles.append(Line2D([0], [0], color=("#FF5050", "#40E060",
                                                   "#5080FF")[k], lw=4,
                                  label=f"{ly[0].name}  {min(win):.1f}–"
                                        f"{max(win):.1f}"))
        ax.imshow(rgb_image(imgs), extent=ext, origin="upper",
                  interpolation="nearest")
        ax.legend(handles=handles, loc="upper right", fontsize=8,
                  facecolor="black", edgecolor="none", framealpha=0.65,
                  labelcolor="white")
        ax.set_title("RGB composite  (each channel on its own 1–99 % range)",
                     fontsize=9, color=pal["plot_fg"])
    else:
        img = mode_image(s, region, use)
        try:
            cmap = cmap.with_extremes(bad=pal["plot_bg"])
        except Exception:
            pass
        vmin, vmax = snapmap.colour_range(img)
        im = ax.imshow(img, extent=ext, origin="upper", cmap=cmap,
                       interpolation="nearest", vmin=vmin, vmax=vmax)
        ax.set_title(f"{use[0][0].name} / {use[1][0].name}", fontsize=9,
                     color=pal["plot_fg"])
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(colors=pal["muted"], labelsize=8)
        cb.set_label("ratio of summed counts", color=pal["plot_fg"],
                     fontsize=8)
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    _rois_on(ax, cube, rois)
    if s["scalebar"]:
        _scalebar(ax, abs(right - left), unit_mark)
    return ax, None
