"""Per-panel view settings (Tk-free).

The page has one *default* look (the view row: view, normalise, offset, z axis,
reverse, fit layers). A panel may override any of those keys, so one page can
show a depth profile as a waterfall, a survey as a stack and a CasaXPS fit as a
single fitted spectrum. Overrides are keyed by the panel's group label and are
stored beside the look (``Workspace.panel_views``, ``state["panel_views"]``).
Page-level settings (grouping, energy scale, style, colours) are not here.
"""
import viewdata

VIEWS = ("Stack", "Waterfall 3D", "Heatmap", "Fit", "Fit stack")
FIT_STACK_MAX = 8          # fitted traces on one "Fit stack" panel (the trace
                           # window steps through a longer group)
FIT_VIEWS = ("Fit", "Fit stack")
SERIES = ("Waterfall 3D", "Heatmap")           # read in z order
NORMS = ("None", "Max = 1", "Area = 1", "At cursor", "Range")
FIT_LAYERS = ("components", "envelope", "background")
IDENT_LAYERS = ("secondary", "auger")
# What the user sees: the "secondary" tier is other elements' lines that could
# also explain a peak (possibilities), not secondary lines of the same element.
IDENT_LABELS = {"secondary": "Alternatives", "auger": "Auger"}
KEYS = ("view", "norm", "offset", "z_axis", "reverse", "fit_show",
        "fit_merge", "ident_show")
OFFSET_RANGE = (0.0, 3.0)


def is_series(view):
    return view in SERIES


def sanitise(d):
    """The valid part of one panel's override; unknown keys and bad values are
    dropped. ``fit_show`` may name only some layers."""
    out = {}
    if not isinstance(d, dict):
        return out
    if d.get("view") in VIEWS:
        out["view"] = d["view"]
    if d.get("norm") in NORMS:
        out["norm"] = d["norm"]
    if d.get("z_axis") in viewdata.Z_MODES:
        out["z_axis"] = d["z_axis"]
    if isinstance(d.get("reverse"), bool):
        out["reverse"] = d["reverse"]
    if isinstance(d.get("fit_merge"), bool):    # INDEX groups drawn as one
        out["fit_merge"] = d["fit_merge"]
    v = d.get("offset")
    if isinstance(v, (int, float)) and not isinstance(v, bool) and v == v:
        lo, hi = OFFSET_RANGE
        out["offset"] = round(max(lo, min(hi, float(v))), 4)
    fs = d.get("fit_show")
    if isinstance(fs, dict):
        fs = {k: fs[k] for k in FIT_LAYERS if isinstance(fs.get(k), bool)}
        if fs:
            out["fit_show"] = fs
    ids = d.get("ident_show")
    if isinstance(ids, dict):
        ids = {k: ids[k] for k in IDENT_LAYERS if isinstance(ids.get(k), bool)}
        if ids:
            out["ident_show"] = ids
    return out


def sanitise_all(views):
    """``{panel label: override}`` with every override cleaned; empty ones
    disappear."""
    out = {}
    if not isinstance(views, dict):
        return out
    for label, d in views.items():
        d = sanitise(d)
        if d:
            out[str(label)] = d
    return out


def sanitise_order(order):
    """A panel order as a list of distinct labels (anything else: empty)."""
    out = []
    if isinstance(order, (list, tuple)):
        for k in order:
            if isinstance(k, str) and k not in out:
                out.append(k)
    return out


def order_groups(groups, order):
    """``[(label, regions)]`` with the labels in ``order`` first, in that
    order; the others keep their own order after them (a label in ``order``
    that is not on the page is ignored)."""
    rank = {k: i for i, k in enumerate(order or [])}
    return sorted(groups, key=lambda g: rank.get(g[0], len(rank)))


def moved(keys, key, how):
    """``keys`` (the labels as shown) with ``key`` moved: ``how`` is
    "first", "earlier", "later" or "last"."""
    keys = list(keys)
    if key not in keys:
        return keys
    i = keys.index(key)
    keys.pop(i)
    j = {"first": 0, "earlier": max(0, i - 1), "later": i + 1,
         "last": len(keys)}[how]
    keys.insert(min(j, len(keys)), key)
    return keys


def moved_to(keys, key, target):
    """``keys`` with ``key`` taken out and put where ``target`` is."""
    keys = list(keys)
    if key not in keys or target not in keys or key == target:
        return keys
    j = keys.index(target)
    keys.remove(key)
    keys.insert(j, key)
    return keys


def resolve(defaults, override=None):
    """The full look of a panel: the page default with the override on top."""
    look = dict(defaults)
    look["fit_show"] = dict(defaults.get("fit_show") or {})
    look["ident_show"] = dict(defaults.get("ident_show") or {})
    ov = sanitise(override)
    fs = ov.pop("fit_show", {})
    ids = ov.pop("ident_show", {})
    look.update(ov)
    look["fit_show"].update(fs)
    look["ident_show"].update(ids)
    return look


def with_value(views, label, key, value):
    """A copy of ``views`` with one key of one panel's override set. For
    ``fit_show``/``ident_show`` pass ``(layer, bool)`` as the value."""
    out = {k: dict(v) for k, v in views.items()}
    cur = out.setdefault(label, {})
    if key in ("fit_show", "ident_show"):
        layer, on = value
        cur[key] = {**cur.get(key, {}), layer: bool(on)}
    else:
        cur[key] = value
    cleaned = sanitise_all(out)
    return cleaned


def reset(views, label=None):
    """A copy of ``views`` without one panel's override (or without any)."""
    if label is None:
        return {}
    return {k: dict(v) for k, v in views.items() if k != label}


def arrange(regs, look, z_order):
    """Order one group's spectra for its view. ``z_order(regs, mode)`` returns
    indices (``viewdata.z_sorted`` in the app). Series views read in z order,
    stacks and fits run in file order, optionally reversed."""
    if is_series(look["view"]):
        return [regs[i] for i in z_order(regs, look["z_axis"])]
    return list(regs)[::-1] if look["reverse"] else list(regs)


def describe(override, defaults=None):
    """A short phrase for tooltips and slide notes, e.g. "Waterfall 3D, by
    Etch level"; only what the override says (or what changes the default)."""
    ov = sanitise(override)
    if not ov:
        return ""
    look = resolve(defaults, ov) if defaults else ov
    view = look.get("view", "Stack")
    parts = [view]
    if is_series(view) and look.get("z_axis") not in (None, "Auto"):
        parts[0] += f", by {look['z_axis']}"
    if "norm" in ov and ov["norm"] != "None":
        parts.append(f"normalised {ov['norm']}")
    if "offset" in ov and view in ("Stack", "Fit", "Fit stack"):
        parts.append(f"offset {ov['offset']:g}×")
    if ov.get("reverse"):
        parts.append("reversed")
    if "fit_merge" in ov:
        parts.append("INDEX groups merged" if ov["fit_merge"]
                     else "INDEX groups drawn separately")
    if view in FIT_VIEWS:
        on = [k for k in FIT_LAYERS if look.get("fit_show", {}).get(k, True)]
        if len(on) < len(FIT_LAYERS):
            parts.append("showing " + (", ".join(on) if on else "no layers"))
    ident_on = [k for k in IDENT_LAYERS if look.get("ident_show", {}).get(k)]
    if ident_on:
        parts.append("other possible lines: " + ", ".join(IDENT_LABELS[k].lower() for k in ident_on))
    return ", ".join(parts)
