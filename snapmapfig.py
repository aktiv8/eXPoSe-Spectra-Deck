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
     "show_spectrum": bool, "scalebar": bool}

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
    "scalebar": True, "tiles": [],
}


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
    tl = spec.get("tiles")
    out["tiles"] = ([dict(t) for t in tl if isinstance(t, dict)]
                    if isinstance(tl, (list, tuple)) else [])
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


def roi_stats(cube, window, background, rois) -> list:
    """One row per area, from the map as shown (the counts summed over the
    window): ``label, shape, px, area_um2, mean, sum, sd`` and ``ratio``, the
    area's mean over the whole map's. An area with no pixel is left out."""
    import numpy as np
    image = cube.image(window[0], window[1], background=background)
    whole = float(np.mean(image)) if image.size else 0.0
    rows = []
    for r in rois:
        m = roi_mask(cube, r)
        if m is None:
            continue
        v = image[m]
        mean = float(v.mean())
        rows.append({
            "label": r.get("label", ""), "shape": r.get("shape", "rect"),
            "px": int(m.sum()),
            "area_um2": float(m.sum() * abs(cube.dx * cube.dy)),
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
    if s["rois"]:
        n = len(s["rois"])
        parts.append("1 area of interest marked" if n == 1
                     else f"{n} areas of interest marked")
    if s["overlay"]["on"]:
        parts.append("over the camera image")
    return "; ".join(parts) + "."


def draw_view(fig, rect, spec, region, pal, cmap, photo=None,
              unit_mark=None, tiles=None):
    """Draw a view onto ``fig`` inside ``rect`` = ``(left, bottom, right, top)``
    figure fractions (None = the whole figure) and return its axes
    ``(ax_map, ax_spec or None)``.

    ``region`` carries the cube (``extra["cube"]``) and the energy labels;
    ``pal`` needs ``plot_bg``, ``plot_fg``, ``muted`` and ``accent``;
    ``photo`` is ``(PIL image, calibration)`` for the camera overlay.
    ``unit_mark`` is appended to the scale bar (e.g. "(approx.)").
    ``tiles`` are further regions whose maps are laid beside this one by their
    stage positions (one canvas, one window, one colour range); the camera
    overlay and the area are not used then."""
    from matplotlib.patches import Rectangle
    s = sanitise(spec)
    cube = region.extra["cube"]
    window = tuple(s["window"]) if s["window"] else snapmap.default_window(cube)
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
        length = nice_length(abs(right - left))
        # the bar sits near the bottom edge of the map
        ymin, ymax = ax_map.get_ylim()
        by = ymin + (ymax - ymin) * 0.06
        xmin, xmax = ax_map.get_xlim()
        bx = xmin + (xmax - xmin) * 0.05
        ax_map.plot([bx, bx + length], [by, by], color="white", lw=3,
                    solid_capstyle="butt", zorder=6,
                    path_effects=_outline())
        txt = f"{length:g} µm" + (f" {unit_mark}" if unit_mark else "")
        ax_map.text(bx + length / 2, by + (ymax - ymin) * 0.025, txt,
                    color="white", fontsize=8, ha="center", va="bottom",
                    zorder=6, path_effects=_outline())

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
