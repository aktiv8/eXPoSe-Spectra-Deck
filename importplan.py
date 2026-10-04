"""Deciding what to import when a folder holds the same data twice.

Thermo Avantage writes both ``.avg`` and ``.vgd`` for one acquisition, and
Kratos Vision2 a ``.dset`` that DumpDataset turns into a ``.kal``; in both
cases the data are identical, so loading both just doubles everything in the
tree. Pure functions here; the dialog lives in the app.

A *family* is the two extensions that can hold the same data. The pair's
choices are the two extensions without the dot, or ``"both"``.
"""

from __future__ import annotations

import os

AVG_VGD = (".avg", ".vgd")
KAL_DSET = (".kal", ".dset")
FAMILIES = (AVG_VGD, KAL_DSET)
CHOICES = ("avg", "vgd", "both")                 # the .avg / .vgd family


def choices(family=AVG_VGD):
    """The valid choices for ``family``: its two extensions, or ``both``."""
    return tuple(e.lstrip(".") for e in family) + ("both",)


def find_pairs(paths, family=AVG_VGD):
    """``[(first_path, second_path)]`` for files in the same folder that share
    a (case-insensitive) name stem and have both extensions of ``family``
    (``.avg`` then ``.vgd`` by default). Order follows the first appearance
    of each pair in ``paths``."""
    first, second = family
    groups = {}
    for p in paths:
        stem, ext = os.path.splitext(os.path.basename(p))
        ext = ext.lower()
        if ext in family:
            key = (os.path.normcase(os.path.dirname(os.path.abspath(p))),
                   stem.lower())
            groups.setdefault(key, {}).setdefault(ext, p)
    return [(g[first], g[second]) for g in groups.values()
            if first in g and second in g]


def apply_choice(paths, pairs, choice, family=AVG_VGD):
    """The paths to load after resolving ``pairs``.

    ``choice`` is one of ``choices(family)`` (``"avg"``, ``"vgd"`` or
    ``"both"`` for the default family), or a dict mapping a pair's first path
    to one of those (pairs missing from the dict load both). Files that are
    not part of a pair are always kept, in order."""
    valid = choices(family)
    keep_first, keep_second = valid[0], valid[1]
    drop = set()
    for first, second in pairs:
        c = choice.get(first, "both") if isinstance(choice, dict) else choice
        if c not in valid:
            raise ValueError(f"Unknown choice {c!r}")
        if c == keep_first:
            drop.add(second)
        elif c == keep_second:
            drop.add(first)
    return [p for p in paths if p not in drop]


def classify_paths(paths):
    """Split dropped or chosen paths into ``(workbooks, folders, files)``
    (workbooks are ``.xpscontainer`` files; the order of ``paths`` is kept)."""
    books, folders, files = [], [], []
    for p in paths:
        p = str(p)
        if os.path.isdir(p):
            folders.append(p)
        elif p.lower().endswith(".xpscontainer"):
            books.append(p)
        else:
            files.append(p)
    return books, folders, files
