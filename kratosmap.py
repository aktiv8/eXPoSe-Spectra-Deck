"""Kratos stigmatic imaging maps: single-energy images from a ``.kal`` / ``.dset``.

Tk-free. An object with ``Scan type = F_MAPPING`` is one image (256 x 256 in
every file seen) at one kinetic energy (``Map energy/mass``); a dataset holds
many, in acquisition order: a focus series (the stage Z steps between images),
several elements at one position, or the same pair repeated over time. Each
becomes a :class:`snapmap.MapCube` with a single energy channel in
``Region.extra["cube"]`` (the region itself has no spectrum), so what the
SnapMap viewer needs from a cube (``image``, ``rect_mask``, ``extent``) works
unchanged; the series viewer (``imaging_ui``) works from :func:`frames`.

**The field of view is not in the file.** ``step size`` x ``Full Scale
Deflection`` is a nominal 1.851 um per pixel (the scan runs from -1 to +1 of
the deflection), but registering neighbouring images that the stage had moved
by a known amount gave 1.728 and 1.730 um/pixel (two independent pairs, r 0.89
and 0.88 against 0.98 for images at one position, imaging3 of the Kratos Axis
test set): 0.935 of the nominal. ``MEASURED_SCALE`` applies that and the maps
say so. It is a measurement on one instrument in one lens mode (MHSA medium
magnification), not a recorded value. Directions, from the same pairs: image
x runs *against* stage X and image y (rows, downward) *with* stage Y, like the
camera pictures, so the image is not flipped; a residual rotation of about
3.7 degrees between image and stage axes is not modelled.

numpy is imported inside the functions that need it, so the readers load
without it.
"""

from __future__ import annotations

import re
from array import array
from dataclasses import dataclass, field
from datetime import datetime

import snapmap

MEASURED_SCALE = 0.935
SCALE_NOTE = ("Pixel size is approximate: the file does not record the field "
              "of view, so it was calibrated from stage moves between images "
              "(about 1.73 um/pixel); rotation between image and stage axes "
              "is not applied.")
_NUM = re.compile(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?")


def _num(text, default=None):
    m = _NUM.search(text or "")
    return float(m.group(0)) if m else default


def _six(x):
    """6 significant digits: what a ``.kal`` prints, so a ``.dset`` (full
    doubles) gives exactly the same cube."""
    return float(f"{x:.6g}")


def is_map(o) -> bool:
    return "MAPPING" in o.get("Scan type", "").upper()


def build_cube(o, values):
    """``MapCube`` from one object's ``{field name: value text}`` and its pixel
    values (row-major ``[line][point]``). Raises ``ValueError`` when the pixel
    count is not ``# points per line`` x ``# lines``."""
    nx = int(_num(o.get("# points per line in map"), 0) or 0)
    ny = int(_num(o.get("# lines in map"), 0) or 0)
    if nx < 1 or ny < 1 or len(values) != nx * ny:
        raise ValueError(f"{len(values)} pixels for a {nx} x {ny} map")
    fsd_x = _num(o.get("Full Scale Deflection X"), 0.0) * 1000.0      # mm -> um
    fsd_y = _num(o.get("Full Scale Deflection Y"), 0.0) * 1000.0
    step_x = _six(_num(o.get("step size x coord"), 2.0 / max(nx - 1, 1)))
    step_y = _six(_num(o.get("step size y coord"), 2.0 / max(ny - 1, 1)))
    dx = _six(step_x * fsd_x * MEASURED_SCALE)
    dy = _six(step_y * fsd_y * MEASURED_SCALE)
    ke = _num(o.get("Map energy/mass"))
    cube = snapmap.MapCube(
        [ke if ke is not None else 0.0], nx, ny,
        -(nx - 1) / 2 * dx, dx, -(ny - 1) / 2 * dy, dy,
        array("f", values), "Counts")
    sx, sy = _num(o.get("Stage X Position")), _num(o.get("Stage Y Position"))
    if sx is not None and sy is not None:                              # m -> mm
        cube.stage_x_mm, cube.stage_y_mm = _six(sx * 1000), _six(sy * 1000)
    return cube


# -- the series viewer's view of a file's maps ----------------------------------
@dataclass
class Frame:
    region: object
    cube: object
    ke: float | None
    be: float | None
    z_um: float | None
    when: datetime | None
    position: str
    index: int = 0
    label: str = field(default="")


def _when(text):
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime((text or "").strip(), fmt)
        except ValueError:
            pass
    return None


def frames(regions):
    """The maps among ``regions`` (those with a one-channel cube) as
    :class:`Frame`s in acquisition order."""
    out = []
    for r in sorted((r for r in regions
                     if r.extra.get("cube") is not None
                     and r.extra["cube"].n_energy == 1), key=lambda r: r.index):
        cube = r.extra["cube"]
        ke = r.extra.get("map_ke")
        hv = r.photon_energy
        out.append(Frame(r, cube, ke,
                         round(hv - ke, 4) if hv and ke is not None else None,
                         r.extra.get("stage_z_um"), _when(r.date),
                         r.extra.get("position_name", ""), len(out)))
    for f in out:
        f.label = describe(f)
    return out


def describe(f: Frame) -> str:
    """One line naming a frame: element, energy, height, time."""
    bits = [f.region.name]
    if f.be is not None:
        bits.append(f"BE {f.be:.2f} eV")
    elif f.ke is not None:
        bits.append(f"KE {f.ke:.2f} eV")
    if f.z_um is not None:
        bits.append(f"Z {f.z_um:.0f} um")
    if f.when is not None:
        bits.append(f.when.strftime("%H:%M"))
    return "  ".join(bits)


def select(frames_, current: Frame, how: str):
    """The frames to step through: ``all``, those at ``current``'s energy
    (``energy``) or at its stage position (``position``: X and Y within 1 um)."""
    if how == "energy" and current.ke is not None:
        return [f for f in frames_ if f.ke == current.ke]
    if how == "position" and current.cube.stage_x_mm is not None:
        return [f for f in frames_ if f.cube.stage_x_mm is not None
                and abs(f.cube.stage_x_mm - current.cube.stage_x_mm) < 1e-3
                and abs(f.cube.stage_y_mm - current.cube.stage_y_mm) < 1e-3]
    return list(frames_)


def default_filter(frames_, current: Frame) -> str:
    """Which frames to open on: those at the same stage position when the
    stage height varies among them (a focus series), else those at the same
    energy (a repeated or time series), else all of them."""
    here = select(frames_, current, "position")
    zs = [f.z_um for f in here if f.z_um is not None]
    if len(here) > 1 and zs and max(zs) - min(zs) > 1.0:
        return "position"
    if len(select(frames_, current, "energy")) > 1:
        return "energy"
    return "all"


def series_axis(frames_):
    """``(label, values)`` for plotting a measure against the frames: the
    stage height when it varies by more than 1 um (a focus series), else the
    time of day in minutes from the first frame, else the frame number."""
    zs = [f.z_um for f in frames_]
    if all(z is not None for z in zs) and max(zs) - min(zs) > 1.0:
        return "Stage Z (um)", [float(z) for z in zs]
    ts = [f.when for f in frames_]
    if all(t is not None for t in ts) and ts[-1] != ts[0]:
        return "Minutes from first", [(t - ts[0]).total_seconds() / 60 for t in ts]
    return "Frame", [float(i + 1) for i in range(len(frames_))]


# -- numbers from the pixels -------------------------------------------------------
def pixels(frame: Frame):
    """``(ny, nx)`` float64 view of a frame's counts."""
    import numpy as np
    return frame.cube.array3d()[:, :, 0].astype(np.float64)


def blur(img, sigma: float):
    """Gaussian blur by ``sigma`` pixels (separable, edges mirrored); the
    image itself when ``sigma`` is not positive. scipy is not required."""
    import numpy as np
    a = np.asarray(img, dtype=np.float64)
    if not sigma or sigma <= 0:
        return a
    r = max(1, int(round(3 * sigma)))
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    for axis in (0, 1):
        pad = [(0, 0), (0, 0)]
        pad[axis] = (r, r)
        p = np.pad(a, pad, mode="reflect")
        n = a.shape[axis]
        acc = np.zeros_like(a)
        for j, w in enumerate(k):
            acc += w * (p[j:j + n, :] if axis == 0 else p[:, j:j + n])
        a = acc
    return a


def roi_means(frames_, mask=None):
    """Mean counts per pixel of each frame, over the ``(ny, nx)`` boolean
    ``mask`` (the whole image when None)."""
    out = []
    for f in frames_:
        a = pixels(f)
        if mask is not None:
            if mask.shape != a.shape or not mask.any():
                out.append(float("nan"))
                continue
            a = a[mask]
        out.append(float(a.mean()))
    return out


def focus_metric(img, sigma: float = 1.5) -> float:
    """How sharp an image is: the mean squared gradient of the blurred image
    over the mean counts squared, so brighter or longer images do not score
    higher for that reason alone. Peaks at best focus in a Z series."""
    import numpy as np
    b = blur(img, sigma)
    m = float(b.mean())
    if m <= 0:
        return 0.0
    gy, gx = np.gradient(b)
    return float((gx ** 2 + gy ** 2).mean() / (m * m))
