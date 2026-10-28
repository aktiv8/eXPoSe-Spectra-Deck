"""Peak fitting for a SnapMap: which chemical state is where.

Tk-free, numpy only (scipy is not installed). A *fit spec* is a plain dict::

    {"components": [{"name": "Si", "be": 99.4, "fwhm": 0.9, "shape": "GL(30)",
                     "doublet": {"split": 0.61, "ratio": 0.5} | None,
                     "free_be": True, "free_fwhm": True,
                     "be_range": 0.5, "fwhm_range": [0.4, 2.5]}, ...],
     "window": [lo, hi],              # binding energies the fit uses
     "background": "shirley" | "linear",
     "avg": 3,                        # channels averaged at each end of the window
     "bin": 1,                        # pixels averaged for the per-pixel step
     "shared_fwhm": False}

The method is variable projection: the positions and widths are refined once
on the *mean* spectrum of the whole map (or of an area) by a small
Levenberg-Marquardt, the amplitudes being solved by non-negative least squares
inside every step. The shapes are then **frozen** and each pixel needs only its
amplitudes, a linear problem solved for all pixels at once (:func:`fit_maps`),
so a 100 x 100 map takes seconds and a noisy pixel cannot wander off to a
different peak. The background is a pair of non-negative columns (its level at
each end of the window, joined by the Shirley step of the mean spectrum or by a
straight line), also solved per pixel.

An amplitude is the component's **counts summed over its channels** in that
pixel (what the channel counts add up to), not a rate and not corrected for
sensitivity: nothing here says atomic percent.

Nothing edits a cube; every result is recomputed from the spec.
"""

from __future__ import annotations

import math

DEFAULT_SHAPE = "GL(30)"
DEFAULT_AVG = 3              # end-point averaging of the Shirley step (CasaXPS "Avg Width")
MAX_AVG = 10
DOUBLET_RATIO = {"p": 0.5, "d": 2.0 / 3.0, "f": 3.0 / 4.0}   # 2nd line / 1st
DEFAULT_HV = 1486.6          # only matters for asymmetric shapes
FRACTION_FLOOR = 0.05        # a total under this share of its 99th percentile is blank
MAX_COMPONENTS = 6
MIN_FWHM, MAX_FWHM = 0.2, 6.0
DEFAULT_BE_RANGE = 0.5       # eV a position may move from its seed (charging, chemistry)
MAX_SEEDED_STATES = 3        # a level with more table states is seeded by its main line


def _np():
    import numpy as np
    return np


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) \
        and math.isfinite(v)


# -- the spec ------------------------------------------------------------------------
def new_component(name, be, fwhm=1.0, shape=DEFAULT_SHAPE, doublet=None,
                  free_be=True, free_fwhm=True, be_range=DEFAULT_BE_RANGE,
                  fwhm_range=(0.4, 2.5)) -> dict:
    return sanitise_component({
        "name": name, "be": be, "fwhm": fwhm, "shape": shape,
        "doublet": doublet, "free_be": free_be, "free_fwhm": free_fwhm,
        "be_range": be_range, "fwhm_range": list(fwhm_range)})


def doublet_for(orbital: str, split: float):
    """The ``doublet`` of a spin-orbit pair: ``orbital`` is the subshell
    letter ("p", "d", "f" or a line such as "2p"), ``split`` the separation in
    eV (the second line lies at higher binding energy), the area ratio the
    degeneracy ratio (1:2, 2:3, 3:4). None for an s level."""
    letter = next((c for c in orbital if c in DOUBLET_RATIO), None)
    if letter is None or not _num(split) or split <= 0:
        return None
    return {"split": float(split), "ratio": DOUBLET_RATIO[letter]}


def sanitise_component(c) -> dict | None:
    """A valid component or None when it has no usable position."""
    if not isinstance(c, dict) or not _num(c.get("be")):
        return None
    import lineshapes
    shape = c.get("shape") if isinstance(c.get("shape"), str) else DEFAULT_SHAPE
    try:
        if lineshapes.check_shape(shape):
            shape = DEFAULT_SHAPE
    except Exception:                                   # noqa: BLE001
        shape = DEFAULT_SHAPE
    fw = c.get("fwhm")
    fw = float(min(MAX_FWHM, max(MIN_FWHM, fw))) if _num(fw) else 1.0
    d = c.get("doublet")
    doublet = None
    if isinstance(d, dict) and _num(d.get("split")) and d["split"] > 0 \
            and _num(d.get("ratio")) and d["ratio"] > 0:
        doublet = {"split": float(d["split"]), "ratio": float(d["ratio"])}
    r = c.get("be_range")
    fr = c.get("fwhm_range")
    if not (isinstance(fr, (list, tuple)) and len(fr) == 2
            and all(_num(v) for v in fr) and fr[0] < fr[1]):
        fr = (0.4, 2.5)
    name = c.get("name")
    return {"name": name.strip() if isinstance(name, str) and name.strip()
            else "Component",
            "be": float(c["be"]), "fwhm": fw, "shape": shape,
            "doublet": doublet,
            "free_be": bool(c.get("free_be", True)),
            "free_fwhm": bool(c.get("free_fwhm", True)),
            "be_range": float(r) if _num(r) and r > 0 else DEFAULT_BE_RANGE,
            "fwhm_range": [float(max(MIN_FWHM, fr[0])),
                           float(min(MAX_FWHM, fr[1]))]}


def sanitise(spec) -> dict:
    spec = spec if isinstance(spec, dict) else {}
    comps = []
    for c in spec.get("components") if isinstance(
            spec.get("components"), (list, tuple)) else []:
        c = sanitise_component(c)
        if c is not None:
            comps.append(c)
        if len(comps) == MAX_COMPONENTS:
            break
    names = []
    for c in comps:                                   # names are unique keys
        base, n = c["name"], 2
        while c["name"] in names:
            c["name"] = f"{base} {n}"
            n += 1
        names.append(c["name"])
    w = spec.get("window")
    b = spec.get("bin")
    av = spec.get("avg")
    return {"components": comps,
            "window": ([float(min(w)), float(max(w))]
                       if isinstance(w, (list, tuple)) and len(w) == 2
                       and all(_num(v) for v in w) and w[0] != w[1] else None),
            "background": spec.get("background")
            if spec.get("background") in ("shirley", "linear") else "shirley",
            "avg": (int(av) if isinstance(av, (int, float))
                    and not isinstance(av, bool) and 1 <= av <= MAX_AVG
                    else DEFAULT_AVG),
            "bin": int(b) if b in (1, 2, 3) else 1,
            "shared_fwhm": bool(spec.get("shared_fwhm", False))}


def noise_gain(cube) -> float:
    """How much one count is worth in the map's values: ``variance = gain x
    value`` (1 for integer counts; about 100 for the Avantage maps, whose
    values are counts scaled by transmission and detector corrections; below 1
    for counts per second). Estimated from the channel-to-channel differences
    of every pixel, which for Poisson noise on a smooth spectrum satisfy
    ``var(x[c] - x[c+1]) = gain (m[c] + m[c+1])``; the median over channels is
    used. It only sets the weights of the fit and the scale of the reduced
    chi-square (so that it is about 1 for a good fit); the amplitudes do not
    depend on it."""
    np = _np()
    a = cube.array3d()
    if a.shape[2] < 4:
        return 1.0
    a = np.asarray(a, dtype=float)
    d = a[:, :, 1:] - a[:, :, :-1]
    s = a[:, :, 1:] + a[:, :, :-1]
    num = (d * d).sum(axis=(0, 1))
    den = s.sum(axis=(0, 1))
    ok = den > 0
    if not ok.any():
        return 1.0
    return float(max(np.median(num[ok] / den[ok]), 1e-6))


# -- seeds from the reference tables -----------------------------------------------------
_LOWER_J = {"p": ("3/2", "1/2"), "d": ("5/2", "3/2"), "f": ("7/2", "5/2")}


def doublet_of_level(core_level, lines=None):
    """The spin-orbit ``doublet`` of a core level ("Si 2p", "Pt 4f"): the
    separation from the app's own line table (``xpslines``: the two lines of
    the element) and the degeneracy area ratio; None for an s level or when
    the table lacks a line."""
    import xpslines
    from readers.base import canon_region_name
    name = canon_region_name(core_level)
    parts = name.split()
    if len(parts) != 2 or parts[1][-1] not in _LOWER_J:
        return None
    el, shell = parts[0], parts[1]
    lo_j, hi_j = _LOWER_J[shell[-1]]
    lines = lines if lines is not None else xpslines.load_lines()
    be = {e["line"]: xpslines.line_be(e) for e in lines
          if e.get("el") == el and e.get("line") in (shell + lo_j,
                                                      shell + hi_j)}
    if len(be) != 2:
        return None
    split = be[shell + hi_j] - be[shell + lo_j]
    return doublet_for(shell, split)


def library_states(core_level, states=None):
    """The chemical states the table lists for a core level, one entry per
    state name (several worked examples are averaged by their median), sorted
    by binding energy: ``[{"name", "be", "fwhm", "source"}]`` (for a picker;
    ``Si2p3/2 SiO2`` -> ``SiO2``)."""
    import chemstates
    from readers.base import canon_region_name
    want = canon_region_name(core_level)
    rows = states if states is not None else chemstates.load_states()
    groups = {}
    for r in rows:
        if canon_region_name(r["core_level"]) == want:
            groups.setdefault(r["state"], []).append(r)
    out = []
    for state, rs in groups.items():
        label = state.split(None, 1)
        bes = sorted(float(r["be"]) for r in rs)
        fw = sorted(float(r["fwhm"]) for r in rs if _num(r.get("fwhm")))
        out.append({"name": label[1] if len(label) > 1 else state,
                    "be": bes[len(bes) // 2],
                    "fwhm": fw[len(fw) // 2] if fw else 1.0,
                    "source": rs[0].get("source", "")})
    return sorted(out, key=lambda e: e["be"])


def main_line(core_level, lines=None):
    """``(label, be)`` of a core level's principal line in the app's line
    table (``Si 2p`` -> ``2p3/2`` 99.4; ``F 1s`` -> 686.0), or None."""
    import xpslines
    from readers.base import canon_region_name
    parts = canon_region_name(core_level).split()
    if len(parts) != 2:
        return None
    el, shell = parts[0], parts[1].split("/")[0]
    shell = shell[:2] if len(shell) > 2 and shell[2].isdigit() else shell
    lines = lines if lines is not None else xpslines.load_lines()
    cands = [e for e in lines if e.get("el") == el
             and (e.get("line") == shell or str(e.get("line", "")).startswith(
                 shell) and e.get("be") is not None)]
    if not cands:
        return None
    best = min(cands, key=lambda e: (e.get("rank", 9), e.get("be")))
    return best["line"], float(xpslines.line_be(best))


def peak_position(be, y, near, span=3.0):
    """The binding energy of the highest smoothed point of ``y`` within
    ``span`` eV of ``near`` (``near`` itself when there is no data there or
    the maximum lies at the edge of that window, i.e. the peak is elsewhere)."""
    np = _np()
    be = np.asarray(be, dtype=float)
    y = np.asarray(y, dtype=float)
    sel = np.abs(be - near) <= span
    if sel.sum() < 3:
        return float(near)
    ys = np.convolve(y[sel], np.ones(3) / 3.0, mode="same")
    k = int(np.argmax(ys[1:-1])) + 1
    if k <= 2 or k >= len(ys) - 3:         # the highest point is the window's
        return float(near)                 # edge: a slope, not a peak nearby
    return float(be[sel][k])


def seed_components(core_level, states=None, lines=None, be=None, y=None):
    """Starting components for any core level, to be refined by
    :func:`fit_mean` (positions may move ``DEFAULT_BE_RANGE`` eV).

    A level with a few states in the chemical-state table (Si 2p: Si and
    SiO2) is seeded with all of them. Any other level (C 1s has many states,
    Al 2p and F 1s none) gets **one** component at the principal line of the
    app's line table, moved to the measured maximum of the spectrum (``be``,
    ``y``) when one is given and lies within 3 eV; more are added by hand or
    from :func:`library_states`. A doublet follows from the line table (2p, 3d,
    4f). Elements the tables do not know give []."""
    doublet = doublet_of_level(core_level, lines)
    lib = library_states(core_level, states)
    if 0 < len(lib) <= MAX_SEEDED_STATES:
        return [new_component(e["name"], e["be"], e["fwhm"], doublet=doublet)
                for e in lib]
    line = main_line(core_level, lines)
    if line is None:
        return []
    pos = line[1]
    if be is not None and y is not None:
        pos = peak_position(be, y, pos)
    label = str(core_level).strip()
    return [new_component(label.split()[0] if label.split() else label, pos,
                          1.0, doublet=doublet)]


def suggest_window(be, comps, pad=3.5):
    """A fit window round the components (``pad`` eV beyond the outermost,
    widened by their allowed movement and doublets), cut to the data."""
    import numpy as np
    be = np.asarray(be, dtype=float)
    lo = min(c["be"] - c["be_range"] for c in comps) - pad
    hi = max(c["be"] + c["be_range"] + (c["doublet"]["split"] if c.get("doublet")
                                       else 0.0) for c in comps) + pad
    return [float(max(lo, be.min())), float(min(hi, be.max()))]


# -- solvers ------------------------------------------------------------------------------
def nnls(A, b, max_iter=None):
    """``x >= 0`` minimising ``|A x - b|`` (Lawson-Hanson active set)."""
    np = _np()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    n = A.shape[1]
    x = np.zeros(n)
    passive = np.zeros(n, dtype=bool)
    w = A.T @ (b - A @ x)
    max_iter = max_iter or 5 * n + 10
    it = 0
    while (~passive).any() and w[~passive].max() > 1e-10 * max(
            1.0, float(np.abs(w).max())) and it < max_iter:
        it += 1
        j = int(np.argmax(np.where(~passive, w, -np.inf)))
        passive[j] = True
        while True:
            z = np.zeros(n)
            z[passive] = np.linalg.lstsq(A[:, passive], b, rcond=None)[0]
            if (z[passive] > 0).all():
                break
            neg = passive & (z <= 0)
            alpha = np.min(x[neg] / (x[neg] - z[neg] + 1e-300))
            x = x + alpha * (z - x)
            passive &= x > 1e-12
            x[~passive] = 0.0
            if not passive.any():
                break
        x = z
        w = A.T @ (b - A @ x)
    return x


def nnls_many(A, Y):
    """Columns of ``Y`` (``n x m``) each solved as :func:`nnls`: one
    unconstrained solve for all, the active-set solve only for the columns
    that came out negative. Returns ``k x m``."""
    np = _np()
    X = np.linalg.lstsq(A, Y, rcond=None)[0]
    bad = np.flatnonzero((X < 0).any(axis=0))
    for j in bad:
        X[:, j] = nnls(A, Y[:, j])
    return X


# -- the model ---------------------------------------------------------------------------
def _eval_line(be, hv, shape, be0, fwhm):
    """A unit-area line (per eV) on the binding-energy grid ``be`` (any
    order): evaluated on ascending kinetic energy, as every shape is."""
    np = _np()
    import lineshapes
    ke = hv - np.asarray(be, dtype=float)
    order = np.argsort(ke)
    out = np.empty(len(ke))
    out[order] = lineshapes.component_curve(ke[order], shape, hv - be0, fwhm,
                                            1.0)
    return out


def component_profile(be, hv, comp, be0=None, fwhm=None):
    """One state as a unit-area profile (per eV) on ``be``; a doublet is two
    lines (the second ``split`` eV higher in binding energy, ``ratio`` of the
    first's area) normalised so the pair sums to one."""
    be0 = comp["be"] if be0 is None else be0
    fw = comp["fwhm"] if fwhm is None else fwhm
    d = comp.get("doublet")
    first = _eval_line(be, hv or DEFAULT_HV, comp["shape"], be0, fw)
    if not d:
        return first
    second = _eval_line(be, hv or DEFAULT_HV, comp["shape"],
                        be0 + d["split"], fw)
    return (first + d["ratio"] * second) / (1.0 + d["ratio"])


def _step_of(be):
    np = _np()
    be = np.asarray(be, dtype=float)
    return float(abs(np.diff(be)).mean()) if len(be) > 1 else 1.0


def shirley_step(be, y, avg=3, iters=100):
    """The Shirley step of a spectrum, ``0`` at the low-binding-energy end to
    ``1`` at the high end (the integrated peak signal over the integrated
    total, iterated), on the order of ``be``."""
    np = _np()
    be = np.asarray(be, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(-be)                    # high BE first, as the sum runs
    ys = y[order]
    n = len(ys)
    k = max(1, min(avg, n // 2))
    lo, hi = float(ys[:k].mean()), float(ys[-k:].mean())   # high BE / low BE
    b = np.linspace(lo, hi, n)
    s = np.zeros(n)
    for _ in range(iters):
        cum = np.cumsum(np.maximum(ys - b, 0.0))
        tot = cum[-1] if cum[-1] > 0 else 1.0
        s = cum / tot                             # 0 at high BE -> 1 at low BE
        nb = lo + (hi - lo) * s
        done = float(np.abs(nb - b).max()) < 1e-9 * max(1.0, abs(hi - lo))
        b = nb
        if done:
            break
    out = np.empty(n)
    out[order] = 1.0 - s                          # 0 at low BE -> 1 at high BE
    return out


def background_columns(kind, be, y_mean, avg=DEFAULT_AVG):
    """Two non-negative columns of the background, with their coefficients.

    **Shirley:** ``[1, step]``: a base level (the background at the low-binding
    energy end) and a **non-negative step** up to the high-BE side, the step
    being the Shirley step of ``y_mean``. The background can only rise towards
    high binding energy (what inelastic scattering does), so it cannot sink
    below a peak on the low-BE side and take its intensity (a low-BE end level
    free to exceed the data there did exactly that to a small elemental-Si
    peak). ``avg`` channels at each end of the window set the end levels the
    Shirley step is iterated between (as CasaXPS's average width); it shapes
    the step only. **Linear:** ``[1 - ramp, ramp]``: the levels at the two
    ends, either of which may be the higher (solved with the peaks, so ``avg``
    does not enter)."""
    np = _np()
    be = np.asarray(be, dtype=float)
    if kind == "shirley":
        s = shirley_step(be, y_mean, avg=avg)
        return np.stack([np.ones(len(be)), s], axis=1)
    span = float(be.max() - be.min()) or 1.0
    s = (be - be.min()) / span
    return np.stack([1.0 - s, s], axis=1)          # (n, 2): low-BE level, high-BE level


def design(be, hv, spec, params=None, y_mean=None):
    """``(A, names)``: one column per component (counts summed over the
    channels = 1) then the two background columns."""
    np = _np()
    be = np.asarray(be, dtype=float)
    step = _step_of(be)
    cols = []
    for i, c in enumerate(spec["components"]):
        be0, fw = (c["be"], c["fwhm"]) if params is None else params[i]
        cols.append(component_profile(be, hv, c, be0, fw) * step)
    cols = np.stack(cols, axis=1) if cols else np.zeros((len(be), 0))
    bg = background_columns(spec["background"], be,
                            y_mean if y_mean is not None else np.zeros(len(be)),
                            spec.get("avg", DEFAULT_AVG))
    return np.concatenate([cols, bg], axis=1)


# -- refinement on the mean spectrum -------------------------------------------------------
class Result:
    """The refined spec and what to draw: ``spec`` (positions and widths
    updated), ``amplitudes`` (counts per component), ``background`` levels,
    ``be`` / ``y`` (the window), ``components`` ``(k, n)``, ``envelope``,
    ``bg`` (the background curve), ``chi2red``, ``errors`` (1-sigma of each
    free position / width, approximate) and ``notes``."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def _free_vector(spec):
    """``[(component index, "be" | "fwhm")]`` for the free parameters, with a
    shared width counted once."""
    free = []
    shared_done = False
    for i, c in enumerate(spec["components"]):
        if c["free_be"]:
            free.append((i, "be"))
        if c["free_fwhm"]:
            if spec["shared_fwhm"]:
                if not shared_done:
                    free.append((i, "fwhm"))
                    shared_done = True
            else:
                free.append((i, "fwhm"))
    return free


def _params_from(theta, spec, free):
    np = _np()
    pars = [[c["be"], c["fwhm"]] for c in spec["components"]]
    shared = None
    for t, (i, what) in zip(theta, free):
        c = spec["components"][i]
        if what == "be":
            lo, hi = spec["_be0"][i] - c["be_range"], spec["_be0"][i] + c["be_range"]
            pars[i][0] = float(np.clip(t, lo, hi))
        else:
            lo, hi = c["fwhm_range"]
            v = float(np.clip(t, lo, hi))
            pars[i][1] = v
            if spec["shared_fwhm"]:
                shared = v
    if shared is not None:
        for p in pars:
            p[1] = shared
    return pars


def _model(be, hv, spec, params, y, w):
    """Weighted NNLS of ``y`` on the columns for ``params``."""
    np = _np()
    A = design(be, hv, spec, params, y)
    sw = np.sqrt(w)
    x = nnls(A * sw[:, None], y * sw)
    return A, x


def fit_mean(be, y, hv, spec, window=None, n_pix=1, gain=1.0,
             max_iter=60) -> Result:
    """Refine the free positions and widths of ``spec`` on one spectrum
    (``be``, ``y`` counts per channel): see :func:`fit_multi`. ``y`` is a mean
    over ``n_pix`` pixels, so its variance is ``gain * y / n_pix`` (the
    weights; ``gain`` from :func:`noise_gain`), which makes the reduced
    chi-square about 1 for a good fit."""
    return fit_multi(be, [(y, n_pix)], hv, spec, window=window, gain=gain,
                     max_iter=max_iter)


def fit_multi(be, spectra, hv, spec, window=None, gain=1.0,
              max_iter=60) -> Result:
    """Refine the free positions and widths of ``spec`` on one or several
    spectra at once: ``spectra`` is ``[(y, n_pix)]``, each a mean over
    ``n_pix`` pixels. Levenberg-Marquardt over only those parameters, which
    the spectra **share**; each spectrum has its own amplitudes and background
    levels (and its own Shirley step), solved by NNLS inside every step
    (variable projection). The result describes the **first** spectrum (the
    whole map's mean), its chi-square included; the others only pin the shapes
    (a state that is strong in a few pixels is invisible in the mean)."""
    np = _np()
    spec = sanitise(spec)
    if not spec["components"]:
        raise ValueError("no components to fit")
    be = np.asarray(be, dtype=float)
    win = window or spec["window"]
    sel = np.ones(len(be), dtype=bool)
    if win:
        sel = (be >= min(win)) & (be <= max(win))
    be = be[sel]
    if len(be) < 3 * len(spec["components"]) + 4:
        raise ValueError("too few channels in the fit window")
    ys, ws = [], []
    for y, n_pix in spectra:
        y = np.asarray(y, dtype=float)[sel]
        ys.append(y)
        ws.append(float(max(1, n_pix)) / (gain * np.maximum(y, gain)))
    y, w = ys[0], ws[0]
    spec = dict(spec, _be0=[c["be"] for c in spec["components"]])
    free = _free_vector(spec)
    theta = np.array([spec["components"][i]["be" if what == "be" else "fwhm"]
                      for i, what in free], dtype=float)

    def resid(th):
        pars = _params_from(th, spec, free)
        parts, models = [], []
        for yy, ww in zip(ys, ws):
            A, x = _model(be, hv, spec, pars, yy, ww)
            parts.append((yy - A @ x) * np.sqrt(ww))
            models.append((A, x))
        return np.concatenate(parts), pars, models

    # a coarse look along each free position first: a component whose seed is
    # off the peak gets no amplitude, so its gradient is zero and the
    # refinement below could never find it
    for j, (i, what) in enumerate(free):
        if what != "be":
            continue
        rng_ = spec["components"][i]["be_range"]
        r0 = resid(theta)[0]
        best = (float(r0 @ r0), theta[j])
        for cand in np.linspace(spec["_be0"][i] - rng_, spec["_be0"][i] + rng_,
                                max(5, int(round(2 * rng_ / 0.1)) + 1)):
            th2 = theta.copy()
            th2[j] = cand
            rr = resid(th2)[0]
            if float(rr @ rr) < best[0] - 1e-9:
                best = (float(rr @ rr), cand)
        theta[j] = best[1]
    r, pars, models = resid(theta)
    cost = float(r @ r)
    lam = 1e-2
    notes = []
    J = None
    for _it in range(max_iter if free else 0):
        J = np.empty((len(r), len(theta)))
        for j in range(len(theta)):
            h = 1e-3 * (1.0 if free[j][1] == "be" else 0.5)
            th2 = theta.copy()
            th2[j] += h
            J[:, j] = (resid(th2)[0] - r) / h
        g = J.T @ r
        H = J.T @ J
        improved = False
        for _try in range(12):
            try:
                step = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-12), -g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            th2 = theta + step
            r2, pars2, models2 = resid(th2)
            c2 = float(r2 @ r2)
            if c2 < cost:
                improved = abs(cost - c2) > 1e-9 * cost
                theta, r, pars, models, cost = th2, r2, pars2, models2, c2
                lam = max(lam / 3.0, 1e-9)
                break
            lam *= 4.0
        if not improved:
            break
    A, x = models[0]
    nc = A.shape[1]
    ndf_all = max(1, len(r) - len(free) - nc * len(ys))
    r1 = r[:len(be)]
    chi2red = float(r1 @ r1) / max(1, len(be) - len(free) - nc)
    errors = {}
    if free and J is not None:
        live = [k for k in range(len(free)) if np.abs(J[:, k]).max() > 0]
        try:
            if live:
                cov = np.linalg.inv(J[:, live].T @ J[:, live]) * (cost / ndf_all)
                for n, k in enumerate(live):
                    i, what = free[k]
                    errors[(spec["components"][i]["name"], what)] = float(
                        math.sqrt(max(cov[n, n], 0.0)))
        except np.linalg.LinAlgError:
            notes.append("the uncertainties could not be estimated "
                         "(the free parameters are not independent)")
    out = sanitise(dict(spec, components=[
        dict(c, be=pars[i][0], fwhm=pars[i][1])
        for i, c in enumerate(spec["components"])]))
    k = len(out["components"])
    comps = (A[:, :k] * x[:k]).T
    bgcurve = A[:, k:] @ x[k:]
    edge = [c for i, c in enumerate(spec["components"])
            if abs(pars[i][0] - spec["_be0"][i]) >= 0.98 * c["be_range"]
            and c["free_be"]]
    if edge:
        notes.append("position at its limit: " + ", ".join(
            c["name"] for c in edge))
    wide = [c["name"] for i, c in enumerate(spec["components"])
            if c["free_fwhm"] and (pars[i][1] <= c["fwhm_range"][0] * 1.001
                                   or pars[i][1] >= c["fwhm_range"][1] * 0.999)]
    if wide:
        notes.append("width at its limit: " + ", ".join(wide))
    unused = [c["name"] for n, c in enumerate(out["components"])
              if all(m[1][n] <= 0 for m in models)]
    if unused:
        notes.append("no intensity in: " + ", ".join(unused))
    return Result(spec=out, amplitudes=x[:k], levels=x[k:], be=be, y=y,
                  components=comps, bg=bgcurve, envelope=comps.sum(axis=0) + bgcurve,
                  chi2red=chi2red, errors=errors, notes=notes)


def enriched_spectra(cube, be, hv, spec, gain=None, share=0.03,
                     min_pixels=30, window=None):
    """``[(mean spectrum over all channels, n_pixels, component name)]``: for
    **minor** component (under half of the map's total in a first partition
    with the starting shapes), the mean spectrum of the pixels where it is
    most prominent (the top ``share`` of the map, at least ``min_pixels``, by
    its share of the total), left out when it is not prominent anywhere (under
    10 % in those pixels). The dominant state is left alone: the mean already
    shows it, and the pixels where it "wins" are arbitrary. Feeding
    these to :func:`fit_multi` lets a state that is strong in only a few
    pixels (a stripe of elemental Si on an oxidised surface) fix its own
    position and width, which the mean of the whole map cannot."""
    np = _np()
    maps = fit_maps(cube, be, hv, spec, gain=gain, window=window)
    total = maps.total
    ok = total > 0
    if ok.any():
        ok &= total >= FRACTION_FLOOR * float(np.percentile(total[ok], 99))
    data = cube.array3d()
    out = []
    grand = float(total.sum()) or 1.0
    for k, name in enumerate(maps.names):
        if ok.sum() < min_pixels:
            break
        if float(maps.areas[k].sum()) / grand >= 0.5:
            continue
        sh = np.where(ok, maps.areas[k] / np.where(total > 0, total, 1.0), -1.0)
        n = max(min_pixels, int(share * ok.sum()))
        idx = np.argsort(sh.ravel())[::-1][:n]
        if float(sh.ravel()[idx].mean()) < 0.10:
            continue
        mask = np.zeros(sh.size, dtype=bool)
        mask[idx] = True
        mask = mask.reshape(sh.shape)
        out.append((data[mask].mean(axis=0), int(mask.sum()), name))
    return out


# -- the whole map ----------------------------------------------------------------------------
class FitMaps:
    """Per-pixel results: ``names``, ``areas`` ``(k, ny, nx)`` (counts summed
    over the channels), ``background`` ``(ny, nx)`` (counts under the
    components), ``chi2`` (reduced, ``ny x nx``), ``total`` (sum of the
    areas), and ``notes``."""

    def __init__(self, names, areas, background, chi2, notes):
        self.names, self.areas = names, areas
        self.background, self.chi2 = background, chi2
        self.total = areas.sum(axis=0)
        self.notes = notes

    def area_of(self, name):
        return self.areas[self.names.index(name)]


def _boxcar(a, b):
    """Mean over ``b`` x ``b`` pixels (edges use what exists) of ``(ny, nx,
    channels)``."""
    np = _np()
    if b <= 1:
        return a
    out = np.zeros_like(a, dtype=float)
    cnt = np.zeros(a.shape[:2])
    ny, nx = a.shape[:2]
    h = b // 2
    for dy in range(-h, b - h):
        for dx in range(-h, b - h):
            ys = slice(max(0, dy), ny + min(0, dy))
            xs = slice(max(0, dx), nx + min(0, dx))
            yd = slice(max(0, -dy), ny + min(0, -dy))
            xd = slice(max(0, -dx), nx + min(0, -dx))
            out[yd, xd] += a[ys, xs]
            cnt[yd, xd] += 1
    return out / cnt[..., None]


def fit_maps(cube, be, hv, spec, y_mean=None, window=None,
             gain=None) -> FitMaps:
    """Frozen shapes, per-pixel amplitudes: every pixel of ``cube`` solved at
    once by weighted non-negative least squares. ``be`` is the cube's energy
    axis as shown (shift applied), ``y_mean`` the mean spectrum the Shirley
    step is taken from (default: the whole map's), ``window`` the fit window
    (default ``spec["window"]``, else everything), ``gain`` the noise gain
    (default: estimated from the cube, :func:`noise_gain`). Nothing is written
    to the cube."""
    np = _np()
    spec = sanitise(spec)
    if not spec["components"]:
        raise ValueError("no components to fit")
    be = np.asarray(be, dtype=float)
    data = cube.array3d()
    win = window or spec["window"]
    sel = np.ones(len(be), dtype=bool)
    if win:
        sel = (be >= min(win)) & (be <= max(win))
    be_w = be[sel]
    Y = np.asarray(data[:, :, sel], dtype=float)
    if spec["bin"] > 1:
        Y = _boxcar(Y, spec["bin"])
    ny, nx, nc = Y.shape
    if y_mean is None:
        y_mean = Y.reshape(-1, nc).mean(axis=0)
    else:
        y_mean = np.asarray(y_mean, dtype=float)
        if len(y_mean) != len(be):
            raise ValueError("the mean spectrum does not match the energy axis")
        y_mean = y_mean[sel]
    A = design(be_w, hv, spec, None, y_mean)
    g = noise_gain(cube) if gain is None else gain
    w = 1.0 / (g * np.maximum(y_mean, g))
    sw = np.sqrt(w)
    flat = Y.reshape(-1, nc)
    X = nnls_many(A * sw[:, None], (flat * sw).T)          # (k+2, npix)
    k = len(spec["components"])
    model = A @ X
    res = (flat.T - model) * sw[:, None]
    chi2 = (res * res).sum(axis=0) / max(1, nc - A.shape[1])
    bg = (A[:, k:] @ X[k:]).sum(axis=0)
    return FitMaps([c["name"] for c in spec["components"]],
                   X[:k].reshape(k, ny, nx), bg.reshape(ny, nx),
                   chi2.reshape(ny, nx),
                   ["spatially binned %d x %d" % (spec["bin"], spec["bin"])]
                   if spec["bin"] > 1 else [])


def fraction(maps, names, floor=FRACTION_FLOOR):
    """The share of the named components in the total at each pixel, NaN where
    the total is zero or under ``floor`` of its own 99th percentile (a
    fraction of noise is not a map)."""
    np = _np()
    num = sum(maps.area_of(n) for n in names) if names else np.zeros_like(
        maps.total)
    tot = maps.total
    out = np.full(tot.shape, np.nan)
    ok = tot > 0
    if ok.any():
        ok &= tot >= floor * float(np.percentile(tot[ok], 99))
    out[ok] = num[ok] / tot[ok]
    return out


def evaluate(be, y, hv, spec, window=None, n_pix=1, gain=1.0) -> Result:
    """The model of ``spec`` as it stands (no refinement) on one spectrum:
    amplitudes and background solved by NNLS for the frozen shapes. What the
    overlay draws for an area of a map; the same fields as :func:`fit_mean`."""
    np = _np()
    spec = sanitise(spec)
    if not spec["components"]:
        raise ValueError("no components to fit")
    be = np.asarray(be, dtype=float)
    y = np.asarray(y, dtype=float)
    win = window or spec["window"]
    if win:
        sel = (be >= min(win)) & (be <= max(win))
        be, y = be[sel], y[sel]
    if len(be) < 3 * len(spec["components"]) + 4:
        raise ValueError("too few channels in the fit window")
    w = float(max(1, n_pix)) / (gain * np.maximum(y, gain))
    A, x = _model(be, hv, spec, None, y, w)
    k = len(spec["components"])
    comps = (A[:, :k] * x[:k]).T
    bg = A[:, k:] @ x[k:]
    r = (y - A @ x) * np.sqrt(w)
    ndf = max(1, len(be) - A.shape[1])
    return Result(spec=spec, amplitudes=x[:k], levels=x[k:], be=be, y=y,
                  components=comps, bg=bg, envelope=comps.sum(axis=0) + bg,
                  chi2red=float(r @ r) / ndf, errors={}, notes=[])
