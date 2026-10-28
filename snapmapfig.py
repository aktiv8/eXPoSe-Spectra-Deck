"""A SnapMap view as data, and the one function that draws it.

Tk-free (matplotlib and numpy are imported inside the drawing functions). A
*view* is a plain JSON-able dict, so it can be stored as a workbook figure and
drawn again for the dialog, a picture file, the PDF, the slides, Word and the
HTML page with the same result::

    {"kind": "snapmap", "ref": <workbook.region_ref of the map's region>,
     "window": [lo, hi], "background": bool, "scale": "Viridis",
     "overlay": {"on": bool, "alpha": 0.65},
     "roi": [x0, y0, x1, y1] | None,          # µm, the map's own frame
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
    "scale": "Viridis", "overlay": {"on": False, "alpha": 0.65}, "roi": None,
    "smooth": {"method": "None", "strength": 0.5}, "show_spectrum": True,
    "scalebar": True,
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
    roi = spec.get("roi")
    out["roi"] = ([float(v) for v in roi]
                  if isinstance(roi, (list, tuple)) and len(roi) == 4
                  and all(_num(v) for v in roi) else None)
    sm = spec.get("smooth") if isinstance(spec.get("smooth"), dict) else {}
    s = sm.get("strength", 0.5)
    out["smooth"] = {"method": sm.get("method") if isinstance(
        sm.get("method"), str) else "None",
        "strength": float(min(1.0, max(0.0, s))) if _num(s) else 0.5}
    out["show_spectrum"] = bool(spec.get("show_spectrum", True))
    out["scalebar"] = bool(spec.get("scalebar", True))
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
    """Boolean ``(ny, nx)`` mask of a view's ROI, or None."""
    if not roi:
        return None
    m = cube.rect_mask(*roi)
    return m if m.any() else None


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


def mean_spectra(cube, mask, smooth):
    """``(whole, roi)`` mean spectra per pixel, smoothed for display when a
    method is chosen."""
    n_all = cube.nx * cube.ny
    whole = [v / n_all for v in cube.total()]
    roi = None
    if mask is not None and mask.any():
        n = int(mask.sum())
        roi = [v / n for v in cube.roi_spectrum(mask)]
    method = (smooth or {}).get("method", "None")
    if method != "None":
        import smoothing
        strength = (smooth or {}).get("strength", 0.5)
        whole = smoothing.smooth(whole, method, strength).tolist()
        if roi is not None:
            roi = smoothing.smooth(roi, method, strength).tolist()
    return whole, roi


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
    if s["roi"]:
        parts.append("area of interest marked")
    if s["overlay"]["on"]:
        parts.append("over the camera image")
    return "; ".join(parts) + "."


def draw_view(fig, rect, spec, region, pal, cmap, photo=None,
              unit_mark=None):
    """Draw a view onto ``fig`` inside ``rect`` = ``(left, bottom, right, top)``
    figure fractions (None = the whole figure) and return its axes
    ``(ax_map, ax_spec or None)``.

    ``region`` carries the cube (``extra["cube"]``) and the energy labels;
    ``pal`` needs ``plot_bg``, ``plot_fg``, ``muted`` and ``accent``;
    ``photo`` is ``(PIL image, calibration)`` for the camera overlay.
    ``unit_mark`` is appended to the scale bar (e.g. "(approx.)")."""
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
    use_photo = s["overlay"]["on"] and photo is not None
    if use_photo:
        img, cal = photo
        cx, cy = snapshot.stage_to_pixel(cal, cube.stage_x_mm, cube.stage_y_mm)
        ux, uy = cal["um_per_px_x"], cal["um_per_px_y"]
        pw, ph = img.size
        ax_map.imshow(img, extent=[-cx * ux, (pw - cx) * ux,
                                   (ph - cy) * uy, -cy * uy], zorder=0)
    image = cube.image(window[0], window[1], background=s["background"])
    vmin, vmax = snapmap.colour_range(image)
    im = ax_map.imshow(image, extent=cube.extent(), origin="upper", cmap=cmap,
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

    mask = roi_mask(cube, s["roi"])
    if mask is not None:
        import numpy as np
        rows = np.flatnonzero(mask.any(axis=1))
        cols = np.flatnonzero(mask.any(axis=0))
        x0, x1 = cube.x_of(cols[0]) - cube.dx / 2, cube.x_of(cols[-1]) + cube.dx / 2
        y0, y1 = cube.y_of(rows[0]) - cube.dy / 2, cube.y_of(rows[-1]) + cube.dy / 2
        ax_map.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                                   lw=1.6, ls="--", ec=pal["accent"], zorder=5))
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
        whole, roi = mean_spectra(cube, mask, s["smooth"])
        e = cube.energy
        ax_spec.plot(e, whole, color=pal["muted"] if roi else pal["accent"],
                     lw=1.4, label="whole map")
        if roi:
            ax_spec.plot(e, roi, color=pal["accent"], lw=1.6,
                         label=f"area, {int(mask.sum())} px")
            ax_spec.legend(fontsize=8, frameon=False, labelcolor=pal["plot_fg"])
        binding = "inding" in (region.energy_label or "")
        ax_spec.set_xlim((max(e), min(e)) if binding else (min(e), max(e)))
        ax_spec.axvspan(lo, hi, color=pal["accent"], alpha=0.14, lw=0)
        ax_spec.set_xlabel(f"{region.energy_label} ({region.energy_units})",
                           color=pal["plot_fg"], fontsize=9)
        ax_spec.set_ylabel("Counts per pixel", color=pal["plot_fg"], fontsize=9)
        ax_spec.set_title(region.name, fontsize=9, color=pal["plot_fg"])
    return ax_map, ax_spec


def _outline():
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=1.6, foreground="black")]
