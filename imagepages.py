"""Camera pictures, SnapMaps and Kratos imaging maps as pages, for the PDF
report and the slides.

Tk-free; matplotlib, numpy and Pillow are imported only when a page is drawn.
:func:`plan` says which pages there are (cheap, geometry only); each
:class:`Page` draws itself onto a matplotlib ``Figure`` given a rectangle, so the
report (an A4 PDF page) and the deck (a picture sized for a slide) share one
drawing:

* **camera sheets**: the calibrated camera pictures, several to a page, each
  with the analysis points that fall in it (its own point highlighted), the
  outline of any SnapMap taken there and a scale bar;
* **SnapMap pages**: one per map site: the camera picture taken there (with the
  map's footprint) beside a grid of element maps, each the counts in the window
  round that element's strongest peak, with its own colour bar;
* **SnapMap mosaic pages** (with ``mosaics``): where the SnapMap sites of one
  file lie side by side (their footprints overlap or nearly touch), one page
  lays each element's maps on one canvas by their stage positions
  (``snapmapfig.mosaic_image``): one energy window and one colour range for all
  the tiles of an element;
* **image-map pages**: the Kratos stigmatic images (``kratosmap``) of one
  stage position, in acquisition order, up to ``IMAGES_PER_PAGE`` to a page,
  each with its own colour scale, a scale bar marked approximate (the field of
  view is not in the file) and a colour bar.
"""

from __future__ import annotations

import io
import os
from dataclasses import dataclass, field
from types import SimpleNamespace

import kratosmap
import mosaic
import reportspec
import snapmap
import snapshot

PICTURE_PX = 640                 # pictures are drawn from a copy this wide at most
IMAGES_PER_PAGE = 12             # image maps on one page (4 across, 3 down)
IMAGE_COLUMNS = 4


@dataclass
class Picture:
    name: str
    sample: str
    parser: object
    blob: object
    points: dict                 # label -> (column, row) in the full-size picture
    outlines: dict               # label -> (left, top, width, height) in pixels
    _array: object = field(default=None, repr=False)
    key: str = ""                # what a report choice calls it (``item_key``)
    _sizes: dict = field(default_factory=dict, repr=False)

    @property
    def calib(self):
        return self.blob.calib

    def array(self, max_px=PICTURE_PX):
        """The picture as an RGB array, shrunk to at most ``max_px`` wide
        (default ``PICTURE_PX``: a PDF page or slide gains nothing from more;
        the mosaic asks for more). Each size is read once."""
        if self._array is not None and max_px == PICTURE_PX:
            return self._array
        if max_px not in self._sizes:
            import numpy as np
            from PIL import Image
            data = self.parser.extract_jpeg(self.blob)
            im = Image.open(io.BytesIO(data)).convert("RGB")
            if im.width > max_px:
                im = im.resize((max_px, round(im.height * max_px
                                              / im.width)), Image.LANCZOS)
            self._sizes[max_px] = np.asarray(im)
        return self._sizes[max_px]


@dataclass
class Site:
    title: str
    maps: list                   # [(shown region, cube)]
    camera: Picture | None = None
    key: str = ""
    doc: object = None


@dataclass
class ImageSet:
    """The Kratos imaging maps of one stage position (a ``Region.sample``)."""
    title: str
    frames: list                 # kratosmap.Frame, in acquisition order
    key: str = ""


@dataclass
class Page:
    kind: str                    # "camera", "mosaic", "mapmosaic", "maps", "images"
    title: str
    n_items: int
    _draw: object = field(default=None, repr=False)
    _notes: object = field(default=None, repr=False)

    def draw(self, fig, rect=(0.0, 0.0, 1.0, 1.0), cmap=None):
        """Draw onto ``fig`` inside ``rect`` = (left, bottom, right, top) in
        figure fractions. ``cmap`` is the colour map for SnapMaps."""
        self._draw(fig, rect, cmap)

    def notes(self):
        return self._notes() if self._notes else ""


# -- planning ------------------------------------------------------------------------
def item_key(kind, doc, name):
    """The stable id of one camera picture (``kind`` "cam", ``name`` the
    picture's) or SnapMap site ("map", the sample's own name) of a file: what
    a report choice lists as a child of the pictures section."""
    return f"{kind}:{reportspec.doc_key(doc)}/{name}"


def _is_snapmap(r) -> bool:
    """A region with a spectrum behind every pixel. A Kratos imaging map is a
    one-channel cube (an image at one energy) and has no report page yet."""
    cube = r.extra.get("cube")
    return cube is not None and cube.n_energy > 1


def _is_image_map(r) -> bool:
    """A single-energy image with no spectrum of its own (a Kratos stigmatic
    map): the test ``htmlbrowser.is_image_map`` makes for the data browser."""
    cube = r.extra.get("cube")
    return cube is not None and cube.n_energy == 1 and not r.decodable


def items(docs, label_of=None):
    """``[(key, label)]`` of every calibrated camera picture, SnapMap site and
    stage position with Kratos imaging maps the files hold, for the Report
    generator (cheap: no picture is read)."""
    label_of = label_of or (lambda p, s: s)
    out = []
    for p in docs:
        for blob in p.images:
            if snapshot.has_calibration(blob.calib):
                sample = label_of(p, blob.sample) if blob.sample else ""
                # the picture's own name usually says the sample already
                out.append((item_key("cam", p, blob.name), blob.name + (
                    f" – {sample}" if sample and sample not in blob.name
                    else "")))
        seen = []
        for r in p.regions:
            if _is_snapmap(r) and r.sample not in seen:
                seen.append(r.sample)
                out.append((item_key("map", p, r.sample),
                            f"SnapMap – {label_of(p, r.sample)}"))
        seen = []
        for r in p.regions:
            if _is_image_map(r) and r.sample not in seen:
                seen.append(r.sample)
                out.append((item_key("img", p, r.sample),
                            f"Image maps – {label_of(p, r.sample)}"))
    return out


def plan(docs, label_of=None, display=None, per_sheet=6, columns=3, skip=(),
         mosaics=False):
    """The pages for ``docs``: camera sheets first, then (with ``mosaics``) one
    page for each set of overlapping pictures stitched together (``mosaic``),
    then one page per SnapMap site. ``label_of(parser, sample)`` gives a
    sample's shown name and ``display(region)`` a region as it should appear
    (names, energy shift); ``skip`` holds the keys (``items``) of pictures and
    sites left out (a picture left out is not drawn beside its map either, nor
    stitched)."""
    label_of = label_of or (lambda p, s: s)
    display = display or (lambda r: r)
    skip = set(skip)
    pictures, sites, sets, image_sets, map_sets = [], [], [], [], []
    for p in docs:
        doc_sites = []
        by_sample = {}                              # sample -> its image maps
        for r in p.regions:
            if _is_image_map(r):
                by_sample.setdefault(r.sample, []).append(display(r))
        mine_sets = [ImageSet(label_of(p, sample), kratosmap.frames(rs),
                              item_key("img", p, sample))
                     for sample, rs in by_sample.items()]
        image_sets.append([s for s in mine_sets if s.key not in skip])
        positions = {label_of(p, k): xy
                     for k, xy in p.sample_positions().items()}
        cubes = {}                                  # sample -> its map regions
        for r in p.regions:
            if _is_snapmap(r):
                cubes.setdefault(r.sample, []).append(r)
        first = {label_of(p, s): rs[0].extra["cube"] for s, rs in cubes.items()}
        mine = []
        for blob in p.images:
            if not snapshot.has_calibration(blob.calib):
                continue
            points, outlines = snapshot.view_of(blob.calib, positions, first)
            mine.append(Picture(blob.name, label_of(p, blob.sample)
                                if blob.sample else "", p, blob, points,
                                outlines, key=item_key("cam", p, blob.name)))
        pictures += [m for m in mine if m.key not in skip]
        if mosaics:                       # per file: another holder, another frame
            sets += [(c, positions, first) for c in mosaic.clusters(
                [m for m in mine if m.key not in skip])]
        for sample, rs in cubes.items():
            if item_key("map", p, sample) in skip:
                continue
            cube = rs[0].extra["cube"]
            near = None
            if cube.stage_x_mm is not None:
                b = snapshot.nearest_image([m.blob for m in mine],
                                           cube.stage_x_mm, cube.stage_y_mm)
                near = next((m for m in mine if m.blob is b), None)
                if near is not None and near.key in skip:
                    near = None
            site = Site(label_of(p, sample),
                        [(display(r), r.extra["cube"]) for r in rs],
                        near, item_key("map", p, sample), p)
            sites.append(site)
            doc_sites.append(site)
        if mosaics:
            map_sets += snapmap_clusters(doc_sites)
    pages = []
    slots = -(-per_sheet // columns) * columns
    chunks = [pictures[i:i + per_sheet] for i in range(0, len(pictures),
                                                       per_sheet)]
    for n, chunk in enumerate(chunks, 1):
        title = "Camera pictures" + (f" ({n} of {len(chunks)})"
                                     if len(chunks) > 1 else "")
        pages.append(Page(
            "camera", title, len(chunk),
            lambda fig, rect, cmap, c=chunk: draw_sheet(fig, rect, c, columns,
                                                        slots),
            lambda c=chunk: _sheet_notes(c)))
    for n, (cluster, positions, first) in enumerate(sets, 1):
        span = sample_span([q.sample for q in cluster])
        title = "Camera mosaic" + (f" \u2013 {span}" if span else "") + (
            f" ({n} of {len(sets)})" if len(sets) > 1 else "")
        pages.append(Page(
            "mosaic", title, len(cluster),
            lambda fig, rect, cmap, c=cluster, ps=positions, fs=first:
                draw_mosaic(fig, rect, c, ps, fs),
            lambda c=cluster: _mosaic_notes(c)))
    for n, cluster in enumerate(map_sets, 1):
        span = sample_span([q.title for q in cluster])
        title = "SnapMap mosaic" + (f" \u2013 {span}" if span else "") + (
            f" ({n} of {len(map_sets)})" if len(map_sets) > 1 else "")
        pages.append(Page(
            "mapmosaic", title, len(cluster),
            lambda fig, rect, cmap, c=cluster: draw_snap_mosaic(fig, rect, c,
                                                                cmap),
            lambda c=cluster: _snap_mosaic_notes(c)))
    for site in sites:
        pages.append(Page(
            "maps", f"SnapMap – {site.title}", len(site.maps),
            lambda fig, rect, cmap, s=site: draw_site(fig, rect, s, cmap),
            lambda s=site: _site_notes(s)))
    for doc_sets in image_sets:
        chunks = _pack_images(doc_sets)
        many = len(doc_sets) > 1               # then each cell names its position
        for n, chunk in enumerate(chunks, 1):
            title = "Image maps" + (f" – {doc_sets[0].title}" if not many
                                    else "")
            title += f" ({n} of {len(chunks)})" if len(chunks) > 1 else ""
            entries = [(s.title if many else "", f) for s, f in chunk]
            pages.append(Page(
                "images", title, len(entries),
                lambda fig, rect, cmap, e=entries: draw_image_sheet(
                    fig, rect, e, cmap),
                lambda e=entries, c=chunk: _image_notes(e, c)))
    return pages


SNAP_GAP = 0.05          # sites this close (share of a map's width) count as touching
SNAP_MOSAIC_ELEMENTS = 6


def _footprint(site):
    """``(left, top, right, bottom)`` of a site's first map in the common
    image frame (µm), or None without a stage position."""
    cube = site.maps[0][1]
    if cube.stage_x_mm is None or cube.stage_y_mm is None:
        return None
    l, r, b, t = cube.extent()
    ox = snapshot.X_DIR * cube.stage_x_mm * 1000.0
    oy = snapshot.Y_DIR * cube.stage_y_mm * 1000.0
    return (min(l, r) + ox, min(t, b) + oy, max(l, r) + ox, max(t, b) + oy)


def snapmap_clusters(sites):
    """Sets of two or more SnapMap sites of one file whose footprints overlap
    or lie within ``SNAP_GAP`` of a map's width of each other, directly or
    through neighbours (stage positions only; nothing is read)."""
    boxes = [(s, _footprint(s)) for s in sites]
    boxes = [(s, f) for s, f in boxes if f is not None]
    parent = list(range(len(boxes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            gap = SNAP_GAP * max(a[2] - a[0], b[2] - b[0])
            if (a[0] - gap <= b[2] and b[0] - gap <= a[2]
                    and a[1] - gap <= b[3] and b[1] - gap <= a[3]):
                parent[find(i)] = find(j)
    groups = {}
    for i, (s, _f) in enumerate(boxes):
        groups.setdefault(find(i), []).append(s)
    return [g for g in groups.values() if len(g) > 1]


def _snap_elements(cluster):
    """``[(name, [(site, region, cube)])]``: the core levels found in at least
    two sites of the cluster, in the order first met, up to
    ``SNAP_MOSAIC_ELEMENTS``."""
    by_name = {}
    for site in cluster:
        seen = set()
        for region, cube in site.maps:
            if region.name in seen:
                continue
            seen.add(region.name)
            by_name.setdefault(region.name, []).append((site, region, cube))
    return [(n, v) for n, v in by_name.items() if len(v) > 1][
        :SNAP_MOSAIC_ELEMENTS]


def _snap_mosaics(cluster):
    """``[(name, window, Mosaic, [(site, region, cube)])]``."""
    import snapmapfig
    out = []
    for name, members in _snap_elements(cluster):
        cubes = [c for _s, _r, c in members]
        window = snapmap.default_window(cubes[0])
        mos = snapmapfig.mosaic_image(cubes, window)
        if mos is not None:
            out.append((name, window, mos, members))
    return out


def _snap_mosaic_notes(cluster):
    import snapmapfig
    rows = []
    for name, window, mos, members in _snap_mosaics(cluster):
        rows.append(f"{name}: counts summed over {window[0]:.1f}\u2013"
                    f"{window[1]:.1f} eV from {len(mos.used)} maps; "
                    + ("; ".join(mos.notes) + "." if mos.notes else ""))
    where = "; ".join(f"{s.title} (stage {s.maps[0][1].stage_x_mm:.3f}, "
                      f"{s.maps[0][1].stage_y_mm:.3f} mm)" for s in cluster)
    return (f"SnapMaps laid side by side by their stage positions only (no "
            f"rotation or image registration), one energy window and one "
            f"colour range for all the maps of an element: {where}.\n"
            + "\n".join(rows))


def draw_snap_mosaic(fig, rect, cluster, cmap=None):
    """One mosaic per element, in a grid, each with its tiles outlined and
    named, a colour bar and a scale bar."""
    import matplotlib
    import snapmapfig
    left, bottom, right, top = rect
    cmap = cmap or matplotlib.colormaps["viridis"]
    try:
        cmap = cmap.with_extremes(bad="#EDEDED")
    except Exception:                      # noqa: BLE001
        pass
    mosaics = _snap_mosaics(cluster)
    if not mosaics:
        return
    n = len(mosaics)
    ext = mosaics[0][2].extent
    aspect = abs(ext[1] - ext[0]) / max(abs(ext[2] - ext[3]), 1e-9)
    # a wide mosaic gets the page's width, a tall one is one of several columns
    cols = 1 if aspect > 1.5 else min(n, 4) if aspect < 0.7 else min(2, n)
    rows = -(-n // cols)
    gs = fig.add_gridspec(rows, cols, left=left + 0.05, right=right - 0.06,
                          top=top - 0.03, bottom=bottom + 0.03, wspace=0.35,
                          hspace=0.45)
    tags = _site_tags([m[0].title for m in mosaics[0][3]])
    for i, (name, window, mos, members) in enumerate(mosaics):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        vmin, vmax = snap_range(mos.image)
        im = ax.imshow(mos.image, extent=mos.extent, origin="upper", cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        cubes = [c for _s, _r, c in members]
        for (k, fl, ft, fw, fh), (site, _r, _c) in zip(
                snapmapfig.mosaic_footprints(cubes), members):
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((fl, ft), fw, fh, fill=False, lw=0.6,
                                   ls=":", ec="#444444"))
            ax.text(fl + fw * 0.02, ft + fh * 0.03,
                    tags.get(site.title, site.title), fontsize=5.5,
                    color="white", va="top", ha="left",
                    bbox=dict(boxstyle="square,pad=0.1", fc="#00000088",
                              ec="none"))
        ax.set_xlim(mos.extent[0], mos.extent[1])
        ax.set_ylim(mos.extent[2], mos.extent[3])
        ax.set_title(f"{name}   {window[0]:.1f}\u2013{window[1]:.1f} eV",
                     fontsize=7, pad=3)
        ax.tick_params(labelsize=5.5, length=2, pad=1)
        ax.set_xlabel("X (\u00b5m)", fontsize=6, labelpad=1)
        ax.set_ylabel("Y (\u00b5m)", fontsize=6, labelpad=1)
        # a colour bar as tall as the (equal-aspect) map, not as its cell
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        cax = make_axes_locatable(ax).append_axes("right", size="3%", pad=0.08)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=5.5, length=2, pad=1)


def _site_tags(titles):
    """Short labels for the tiles: the part of each title that differs from
    the others (``a`` / ``b`` for "Point 002a ..." and "Point 002b ..."), the
    whole title when they do not share a start and an end."""
    titles = list(dict.fromkeys(titles))
    if len(titles) < 2:
        return {t: t for t in titles}
    pre = os.path.commonprefix(titles)
    suf = os.path.commonprefix([t[::-1] for t in titles])[::-1]
    tags = {t: t[len(pre):len(t) - len(suf)].strip() for t in titles}
    if any(not v for v in tags.values()) or len(set(tags.values())) < len(tags):
        return {t: t for t in titles}
    return tags


def snap_range(image):
    return snapmap.colour_range(image)


def _pack_images(sets):
    """``[[(ImageSet, Frame)]]``: a file's image maps, position after
    position, ``IMAGES_PER_PAGE`` to a page. A position that fits on one page
    is not split across two; a longer one fills pages and goes on."""
    pages, cur = [], []
    for s in sets:
        items = [(s, f) for f in s.frames]
        if cur and len(items) <= IMAGES_PER_PAGE and len(cur) + len(
                items) > IMAGES_PER_PAGE:
            pages.append(cur)
            cur = []
        for it in items:
            if len(cur) == IMAGES_PER_PAGE:
                pages.append(cur)
                cur = []
            cur.append(it)
    if cur:
        pages.append(cur)
    return pages


@dataclass
class MosaicSet:
    """One set of overlapping pictures of a file, ready to be stitched."""
    parser: object
    cluster: list                # Pictures
    positions: dict              # shown sample -> stage position
    sites: dict                  # shown sample -> its first SnapMap cube
    key: str = ""

    @property
    def title(self):
        return mosaic_title(self.cluster)


def sample_span(names, limit=3):
    """The samples a set of pictures was taken for as one short phrase:
    ``Stained Pt #003a/b/c`` when the names differ only in a short ending,
    else the names joined, cut after ``limit`` ("+N more")."""
    names = list(dict.fromkeys(n for n in names if n))
    if len(names) <= 1:
        return "".join(names)
    pre = os.path.commonprefix(names)
    ends = [n[len(pre):] for n in names]
    if pre and all(0 < len(e) <= 3 for e in ends):
        return pre + "/".join(ends)
    text = ", ".join(names[:limit])
    return text + (f" +{len(names) - limit} more" if len(names) > limit else "")


def mosaic_title(cluster):
    """"Mosaic of 3 pictures: Stained Pt #003a/b/c" (the samples part only
    when the pictures say which sample they were taken for)."""
    span = sample_span([p.sample for p in cluster])
    return f"Mosaic of {len(cluster)} pictures" + (f": {span}" if span else "")


def mosaic_sets(docs, label_of=None, skip=(), keep="latest"):
    """The sets of overlapping calibrated pictures in ``docs`` (per file:
    another holder is another frame), minus the pictures whose keys are in
    ``skip``. ``keep`` says which picture stands for a stage position that was
    photographed more than once: "latest" (after the analysis) or "earliest"
    (before it). Nothing is stitched until ``mosaic_picture`` is called."""
    label_of = label_of or (lambda p, s: s)
    skip = set(skip)
    out = []
    for p in docs:
        positions = {label_of(p, k): xy
                     for k, xy in p.sample_positions().items()}
        first = {}
        for r in p.regions:
            if _is_snapmap(r):
                first.setdefault(label_of(p, r.sample), r.extra["cube"])
        mine = []
        for blob in p.images:
            if not snapshot.has_calibration(blob.calib):
                continue
            points, outlines = snapshot.view_of(blob.calib, positions, first)
            key = item_key("cam", p, blob.name)
            if key in skip:
                continue
            mine.append(Picture(blob.name, label_of(p, blob.sample)
                                if blob.sample else "", p, blob, points,
                                outlines, key=key))
        for n, c in enumerate(mosaic.clusters(mine, keep), 1):
            out.append(MosaicSet(p, c, positions, first,
                                 item_key("mosaic", p, str(n))))
    return out


def mosaic_picture(cluster, positions, sites, refine=True, gain=True):
    """``(Mosaic, Picture)``: the stitched pictures as one picture with its
    own calibration, so the analysis points and map outlines on it come from
    ``snapshot.view_of`` like on any camera picture. ``refine`` False places
    the pictures by stage position only, ``gain`` False leaves their
    brightness as it is."""
    m = mosaic.build(cluster, refine=refine, gain=gain)
    points, outlines = snapshot.view_of(m.calib, positions, sites)
    pic = Picture(mosaic_title(cluster), "", None,
                  SimpleNamespace(calib=m.calib), points, outlines,
                  _array=m.array)
    return m, pic


def draw_mosaic(fig, rect, cluster, positions, sites):
    """One stitched mosaic filling ``rect``."""
    _m, pic = mosaic_picture(cluster, positions, sites)
    draw_sheet(fig, rect, [pic], 1, 1)


def _mosaic_notes(cluster):
    m = mosaic.build(cluster)
    return "\n".join([f"Mosaic: {m.description()}",
                      "Pictures: " + ", ".join(m.names) + "."] + m.notes)


def available(docs) -> bool:
    """True when ``docs`` hold a calibrated camera picture, a SnapMap or a
    Kratos imaging map."""
    for p in docs:
        if any(snapshot.has_calibration(b.calib) for b in p.images):
            return True
        if any(_is_snapmap(r) or _is_image_map(r) for r in p.regions):
            return True
    return False


def window_of(region, cube):
    """The energy window a map is drawn for, in the region's shown energies
    (the first-choice window of ``snapmap.default_window`` plus any shift the
    shown region carries)."""
    lo, hi = snapmap.default_window(cube)
    shift = region.energy[0] - cube.energy[0] if region.energy else 0.0
    return lo + shift, hi + shift


def _sheet_notes(pictures):
    rows = []
    for pic in pictures:
        c = pic.calib
        rows.append(f"{pic.name}: stage {c['x_mm']:.3f}, {c['y_mm']:.3f} mm; "
                    f"field of view "
                    f"{c['width'] * c['um_per_px_x'] / 1000:.1f} x "
                    f"{c['height'] * c['um_per_px_y'] / 1000:.1f} mm.")
    return "\n".join(rows)


def _site_notes(site):
    rows = []
    for region, cube in site.maps:
        try:
            lo, hi = window_of(region, cube)
            rows.append(f"{region.name}: counts summed over {lo:.1f}–"
                        f"{hi:.1f} {region.energy_units or 'eV'}")
        except ImportError:
            rows.append(region.name)
    cube = site.maps[0][1]
    text = (f"SnapMap of {site.title}: {cube.nx} x {cube.ny} pixels of "
            f"{abs(cube.dx):g} µm. Each map has its own colour scale "
            "(1st to 99th percentile of its pixels).\n") + "\n".join(rows)
    if site.camera is not None:
        text += f"\nCamera picture: {site.camera.name}."
    return text


def _image_notes(entries, chunk):
    """``entries`` as drawn, ``chunk`` the same as ``(ImageSet, Frame)``."""
    cube = chunk[0][1].cube
    where = list(dict.fromkeys(s.title for s, _f in chunk))
    rows = [f"{head}: {f.label}" if head else f.label for head, f in entries]
    return (f"Image maps ({len(chunk)} on this page; "
            f"{', '.join(where)}): {cube.nx} x {cube.ny} pixels of "
            f"{abs(cube.dx):.2f} µm (approx.), in acquisition order. Each "
            "image has its own colour scale (1st to 99th percentile of its "
            f"pixels; counts per pixel). {kratosmap.SCALE_NOTE}\n"
            + "\n".join(rows))


# -- drawing -------------------------------------------------------------------------
def _draw_picture(ax, pic, title_size=7):
    """One camera picture with its markers, map outlines and scale bar."""
    import matplotlib.patheffects as pe
    from matplotlib.patches import Rectangle
    import plots
    c = pic.calib
    w, h = c["width"], c["height"]
    ax.imshow(pic.array(), extent=(0, w, h, 0))
    halo = [pe.withStroke(linewidth=2.5, foreground="#0B1116")]
    for label, (left, top, mw, mh) in pic.outlines.items():
        ax.add_patch(Rectangle((left, top), mw, mh, fill=False, lw=1.0,
                               ls="--", ec="#FFD23F", zorder=2))
        ax.text(max(left, 0) + 4, max(top, 0) + 4, label, color="#FFD23F",
                fontsize=5.5,
                va="top", ha="left", path_effects=halo, zorder=6)
    hot = {pic.sample} if pic.sample else set()
    # a map site is already shown by its outline: no second marker on top
    points = {k: v for k, v in pic.points.items() if k not in pic.outlines}
    marks = plots.draw_holder_markers(ax, points, hot, size=6)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_axis_off()
    ax.set_title(pic.name, fontsize=title_size, pad=2)
    fov = w * c["um_per_px_x"] / 1000.0
    bar_mm = 1.0 if fov >= 4 else 0.5
    bar = bar_mm * 1000.0 / c["um_per_px_x"]
    ax.plot([w * 0.04, w * 0.04 + bar], [h * 0.94] * 2, color="white", lw=2,
            solid_capstyle="butt", path_effects=halo, zorder=6)
    ax.text(w * 0.04 + bar / 2, h * 0.91, f"{bar_mm:g} mm", color="white",
            ha="center", va="bottom", fontsize=6, path_effects=halo, zorder=6)
    ax.apply_aspect()
    try:
        plots.place_marker_labels(ax, marks, hot)
    except Exception:                          # noqa: BLE001 - labels are a nicety
        pass


def draw_sheet(fig, rect, pictures, columns, slots):
    """``pictures`` in a grid ``columns`` wide with ``slots`` cells (so a
    part-filled last sheet keeps the picture size of the others)."""
    left, bottom, right, top = rect
    rows = -(-slots // columns)
    gs = fig.add_gridspec(rows, columns, left=left + 0.01, right=right - 0.01,
                          bottom=bottom, top=top, wspace=0.06, hspace=0.16)
    for i, pic in enumerate(pictures):
        _draw_picture(fig.add_subplot(gs[i // columns, i % columns]), pic)


def draw_site(fig, rect, site, cmap=None):
    """A SnapMap site: the camera picture (if any) and a grid of element maps."""
    import matplotlib
    left, bottom, right, top = rect
    cmap = cmap or matplotlib.colormaps["viridis"]
    n = len(site.maps)
    cols = min(4, n)
    rows = -(-n // cols)
    has_cam = site.camera is not None
    # the grid is as tall as its maps need (about 2 in a row), from the top
    height = min(top - bottom - 0.07, rows * 2.05 / fig.get_figheight())
    gs = fig.add_gridspec(
        rows, cols + (1 if has_cam else 0), left=left + 0.04, right=right - 0.02,
        top=top - 0.03, bottom=top - 0.03 - height, wspace=0.42, hspace=0.42,
        width_ratios=([1.9] if has_cam else []) + [1.0] * cols)
    if has_cam:
        cam_ax = fig.add_subplot(gs[:, 0])
        cam_ax.set_anchor("N")                # top-aligned with the first row
        _draw_picture(cam_ax, site.camera)
    for i, (region, cube) in enumerate(site.maps):
        ax = fig.add_subplot(gs[i // cols, i % cols + (1 if has_cam else 0)])
        lo, hi = snapmap.default_window(cube)
        img = cube.image(lo, hi)
        vmin, vmax = snapmap.colour_range(img)
        im = ax.imshow(img, extent=cube.extent(), cmap=cmap, vmin=vmin,
                       vmax=vmax, interpolation="nearest")
        shift = region.energy[0] - cube.energy[0] if region.energy else 0.0
        ax.set_title(f"{region.name}   {lo + shift:.1f}–{hi + shift:.1f} eV",
                     fontsize=7, pad=3)
        ax.tick_params(labelsize=5.5, length=2, pad=1)
        if i // cols == rows - 1 or i + cols >= n:
            ax.set_xlabel("X (µm)", fontsize=6, labelpad=1)
        if i % cols == 0:
            ax.set_ylabel("Y (µm)", fontsize=6, labelpad=1)
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
        cb.ax.tick_params(labelsize=5.5, length=2, pad=1)


def draw_image_sheet(fig, rect, entries, cmap=None):
    """Kratos imaging maps, ``entries`` = ``[(heading, Frame)]`` (the heading,
    the stage position, may be empty), in a grid ``IMAGE_COLUMNS`` wide with
    ``IMAGES_PER_PAGE`` cells (so a part-filled last page keeps the image size
    of the others); each has its own colour range, scale bar and colour bar."""
    import matplotlib
    import matplotlib.patheffects as pe
    left, bottom, right, top = rect
    cmap = cmap or matplotlib.colormaps["viridis"]
    rows = -(-IMAGES_PER_PAGE // IMAGE_COLUMNS)
    gs = fig.add_gridspec(rows, IMAGE_COLUMNS, left=left + 0.03,
                          right=right - 0.02, bottom=bottom + 0.05,
                          top=top - 0.03, wspace=0.28, hspace=0.30)
    halo = [pe.withStroke(linewidth=2.0, foreground="#0B1116")]
    for i, (heading, f) in enumerate(entries):
        ax = fig.add_subplot(gs[i // IMAGE_COLUMNS, i % IMAGE_COLUMNS])
        img = kratosmap.pixels(f)
        vmin, vmax = snapmap.colour_range(img)
        left_, right_, bottom_, top_ = f.cube.extent()
        im = ax.imshow(img, extent=(left_, right_, bottom_, top_), cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_axis_off()
        ax.set_title(f"{heading}\n{f.label}" if heading else f.label,
                     fontsize=6.5, pad=2)
        width = abs(right_ - left_)
        height = abs(bottom_ - top_)
        bar = kratosmap.scale_bar(width)
        x1 = min(left_, right_) + 0.05 * width
        y = max(bottom_, top_) - 0.06 * height
        ax.plot([x1, x1 + bar], [y, y], color="white", lw=2,
                solid_capstyle="butt", path_effects=halo)
        ax.text(x1, y - 0.02 * height, f"{bar} µm (approx.)", color="white",
                fontsize=5.5, ha="left", va="bottom", path_effects=halo)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=5, length=2, pad=1)
