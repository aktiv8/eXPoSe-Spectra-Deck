"""Thermo Avantage project (``.vgp``): the peak fits of an analysis session.

A ``.vgp`` is an OLE2 compound file (the container of a ``.vgd``; ``ole2.py``)
holding MFC ``CArchive`` serialisations of the objects of a processing
document, one storage per object type, each with a ``Contents`` stream and an
``XRefs`` stream (the links between objects):

* ``Data Spaces`` -- which ``.VGD`` spectra the project uses (a relative path
  into the ``<project>.DATA`` folder beside it, and a GUID), in the order the
  other objects number them;
* ``Fitted Peaks`` -- the fitted peaks (this is what a fit is made of);
* ``Backgrounds`` -- the background of each fitted spectrum;
* ``Peaks`` (the peak *definitions* with their constraints), ``Grid Cells``,
  ``Chem States``, ``Procedural Peaks``, ``Contents`` -- **not read**: the
  definitions are not needed to draw a fit, the others are layout or empty in
  every file seen.

Reverse-engineered from the 129 projects of ``C:\\AvantageSystem\\Knowledge``
(88 with fits, 480 fitted peaks, 270 ``.VGD`` references, all present) and
**validated against the data, not against a second format** (no Avantage
export of a project exists): the model rebuilt from the decoded numbers is
compared with the counts of the same project's ``.VGD`` (``tests/test_thermo_vgp.py``,
``XPS_VGP_CORPUS``). What the data established:

* **A fitted peak** is a record split by the class marker ``01 80`` + the
  class GUID (480 of 480 records; the count in the stream header agrees). After
  a short fixed prefix come the name (a 4-byte byte length, UTF-16, NUL
  included; empty in 23 of 480) and then, packed and **unaligned**, doubles:
  ``area`` (counts/s x eV), ``ke`` (the position, **kinetic** energy, as
  CasaXPS stores it), ``height`` (**counts**, not counts/s: divide by dwell x
  scans), ``fwhm`` (eV), ``mix`` (the Lorentzian share) and three tail terms.
  (The first two numbers read as height and area at first; swapping them put
  the model at the data's absolute scale -- median 1.003 over 55 spectra
  fitted with symmetric peaks -- and made area / (height x FWHM) constant
  within a spectrum to 0.2 %.)
* **Tail terms** ``tail`` = (CT, TM, ET), the constant tail height, tail
  mixing ratio and exponential tail of Avantage's tail function (Data
  Processing Reference Manual v5, p.73): ``T = TM CT + (1 - TM) exp(-ET dx)``,
  ``dx`` the distance from the peak centre in **channels** (data points), the
  peak being ``H [GL + (1 - GL) T]`` on the high-binding-energy side and
  ``H GL`` on the other. ``(0, 1, 0)`` is a symmetric peak. The manual gives
  the formula, not where the file keeps the numbers: that CT, TM, ET are the
  three stored values **in this order** (and the channel is the data step) was
  found by testing every assignment against the stored areas (median error
  2.6 % against 23 % for the next best; 72 % of the 39 tail peaks within
  3 %) and the data (tail spectra: median residual 3.5 % of the range without
  the tail, 1.5 % with it).
* **The base shape** ``GL`` (manual pp.71-72) is Avantage's Gaussian/Lorentzian
  *product* ``1 / ([1 + M t^2] exp{(1 - M) ln2 t^2})`` (``M`` = L/G mix / 100)
  or a *sum*; the manual lists product, sum and convolution (Voigt) forms and
  the file does not say which was used. Where the project's stored areas equal
  height x FWHM x (the shape's area) to 2 % for the product but not for the sum
  (20 of 114 spectra; the stored area matches to 4 digits there) the product is
  drawn; otherwise the sum (which, for these mixes, is also what the Voigt
  looks like, residual 1.2 % against 1.2 %). Six projects store areas in an
  unrecognised unit (ratios 0.03 to 0.5 for either shape) and are drawn as sum
  with a note.
* **Links** come from the ``XRefs`` streams: the first table of ``Fitted
  Peaks/XRefs`` pairs a fitted peak with its data space (the data space's index
  in ``Data Spaces``; a project can list spectra it does not fit), the first
  table of ``Backgrounds/XRefs`` pairs a background with the peaks it serves and
  the second with its data space. Every fitted peak's energy lies inside the
  energy range of the spectrum it links to.
* **A background** record: the two limits (KE), one end-point group (flag,
  averaging points, width in eV) and a type code. Code ``0`` fits as a linear
  background, ``1`` as Shirley (and carries a second group: 50, 0.01, flag --
  read as iterations and tolerance, not confirmed), ``3`` as Shirley too
  (Avantage's "Smart" by name; that is an inference, not read from the file).

Anything this module cannot place is left out and said in ``Project.notes``;
a file whose records do not parse as described raises :class:`VgpError`.

``ThermoProjectFile`` (below) opens a project as a session: its ``.VGD``
spectra with a :class:`casafit.Fit` for every fitted one. A fitted peak
becomes ``AV(mix %, CT, TM, ET, step)`` (the sum) or ``AVP(...)`` (the product)
of ``lineshapes.py``, tail included, always flagged approximate (the form the
project used is inferred, not read); the area is that of the drawn shape at the
stored height; no sensitivity factor is invented.
"""

from __future__ import annotations

import math
import os
import re
import struct
from dataclasses import dataclass, field

import casafit
import lineshapes

from .ole2 import OleFile, OleError, is_ole
from .thermo_experiment import ThermoExperiment, LoadCancelled
from .thermo_vgd import ThermoVgdFile


class VgpError(Exception):
    """The file is not a project this module understands."""


_TAG = b"\x01\x80"                       # "same class as before" object tag
BACKGROUND_KINDS = {0: "Linear", 1: "Shirley", 3: "Smart"}


@dataclass
class DataSpaceRef:
    """A spectrum the project uses: ``path`` as stored (``.\\<project>.DATA\\
    <name>.VGD``) and its ``guid``."""
    path: str
    guid: str = ""

    def resolve(self, project_dir: str) -> str:
        rel = self.path.replace("\\", os.sep).replace("/", os.sep)
        if rel.startswith("." + os.sep):
            rel = rel[2:]
        return os.path.normpath(os.path.join(project_dir, rel))


@dataclass
class FittedPeak:
    """One fitted peak. Units: ``ke`` in eV (kinetic), ``fwhm`` in eV,
    ``area`` in counts/s x eV, ``height`` in counts (divide by dwell x scans
    for counts/s), ``mix`` the Lorentzian share (0-1). ``tail`` is ``(CT, TM,
    ET)`` (see the module docstring and the ``ct`` / ``tm`` / ``et``
    properties). ``data_space`` and ``background`` index
    ``Project.data_spaces`` / ``Project.backgrounds`` (``None`` when the
    project does not say)."""
    index: int
    name: str
    ke: float
    fwhm: float
    mix: float
    area: float
    height: float
    tail: tuple
    data_space: int | None = None
    background: int | None = None

    @property
    def ct(self) -> float:
        """Constant tail height (1 = as high as the peak)."""
        return self.tail[0]

    @property
    def tm(self) -> float:
        """Tail mixing ratio (0 = exponential tail only, 1 = constant only)."""
        return self.tail[1]

    @property
    def et(self) -> float:
        """Exponential tail, per channel (0 = infinite, 1 = no tail)."""
        return self.tail[2]

    @property
    def symmetric(self) -> bool:
        """No tail: ``CT = 0`` and ``TM = 1`` (``T`` is then 0)."""
        return abs(self.tail[0]) < 1e-9 and abs(self.tail[1] - 1.0) < 1e-9


@dataclass
class BackgroundDef:
    """A background: its ``start_ke`` / ``end_ke`` limits (eV, kinetic), the
    end-point averaging (``avg_points`` points, ``width_ev`` eV, from the first
    group), the raw type ``code`` and the groups as stored."""
    index: int
    start_ke: float
    end_ke: float
    avg_points: int
    width_ev: float
    code: int
    groups: tuple = ()
    data_space: int | None = None

    @property
    def kind(self) -> str:
        return BACKGROUND_KINDS.get(self.code, f"type {self.code}")


@dataclass
class Project:
    path: str = ""
    library: str = ""
    data_spaces: list = field(default_factory=list)
    peaks: list = field(default_factory=list)
    backgrounds: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def peaks_of(self, data_space: int) -> list:
        return [p for p in self.peaks if p.data_space == data_space]

    def background_of(self, data_space: int):
        for b in self.backgrounds:
            if b.data_space == data_space:
                return b
        return None

    def fitted_data_spaces(self) -> list:
        return sorted({p.data_space for p in self.peaks
                       if p.data_space is not None})


# -- the archive primitives ----------------------------------------------------
def _stream(ole, name: str) -> bytes:
    try:
        return ole.read(name)
    except Exception:
        return b""


def _header(b: bytes):
    """``(version, object count)`` of a ``Contents`` stream, ``(0, 0)`` when
    it is too short to have one."""
    if len(b) < 8:
        return 0, 0
    return struct.unpack_from("<ii", b, 0)


def _class_tag(b: bytes, pos: int):
    """The class definition ``FF FF <schema:2> <length:2> <name>`` at ``pos``:
    ``(name, offset after it)``."""
    if b[pos:pos + 2] != b"\xff\xff" or len(b) < pos + 6:
        raise VgpError("an object stream does not start with a class tag")
    (nlen,) = struct.unpack_from("<H", b, pos + 4)
    name = b[pos + 6:pos + 6 + nlen].decode("latin-1")
    return name, pos + 6 + nlen


def _cstrings(b: bytes) -> list:
    """The MFC Unicode ``CString``s of a stream: ``FF FE FF <chars:1>`` then
    that many UTF-16 characters."""
    out = []
    for m in re.finditer(rb"\xff\xfe\xff", b):
        i = m.end()
        if i >= len(b):
            continue
        n = b[i]
        raw = b[i + 1:i + 1 + 2 * n]
        if len(raw) == 2 * n:
            try:
                out.append(raw.decode("utf-16le"))
            except UnicodeDecodeError:
                pass
    return out


def _tables(b: bytes) -> list:
    """The ``XRefs`` tables: repeated ``[flag, count, count x (a, b)]`` of
    int32; reading stops at the first block that does not fit."""
    w = struct.unpack("<%di" % (len(b) // 4), b[:len(b) // 4 * 4])
    out, i = [], 0
    while i + 2 <= len(w):
        n = w[i + 1]
        if n < 0 or i + 2 + 2 * n > len(w):
            break
        out.append([(w[i + 2 + 2 * k], w[i + 3 + 2 * k]) for k in range(n)])
        i += 2 + 2 * n
    return out


# -- the object streams ---------------------------------------------------------
def _fitted_records(b: bytes) -> list:
    ver, n = _header(b)
    if n == 0:
        return []
    _name, pos = _class_tag(b, 8)
    guid = b[pos:pos + 16]
    if len(guid) < 16:
        raise VgpError("Fitted Peaks is cut short")
    marker = _TAG + guid
    recs, start = [], pos + 16
    while True:
        nxt = b.find(marker, start)
        if nxt < 0:
            recs.append(b[start:])
            break
        recs.append(b[start:nxt])
        start = nxt + len(marker)
    if len(recs) != n:
        raise VgpError(f"Fitted Peaks lists {n} peaks but {len(recs)} "
                       "records were found")
    return recs


def _fitted_peak(index: int, rec: bytes) -> FittedPeak:
    # prefix: int32 2, int32 1 or 2 (the name sits at +40, or +56 when it is
    # 2: 23 of the 480 peaks); an empty name is a zero length, so the second
    # int decides which offset to try first
    hint = struct.unpack_from("<i", rec, 4)[0] if len(rec) >= 8 else 1
    for off in ((56, 40) if hint == 2 else (40, 56)):
        if len(rec) < off + 4:
            continue
        (n,) = struct.unpack_from("<I", rec, off)
        if n % 2 or n > 512 or len(rec) < off + 4 + n + 64:
            continue
        raw = rec[off + 4:off + 4 + n]
        if n and raw[-2:] != b"\x00\x00":
            continue
        try:
            name = raw.decode("utf-16le").rstrip("\x00")
        except UnicodeDecodeError:
            continue
        d = off + 4 + n
        v = [struct.unpack_from("<d", rec, d + 8 * k)[0] for k in range(8)]
        return FittedPeak(index=index, name=name, area=v[0], ke=v[1],
                          height=v[2], fwhm=v[3], mix=v[4],
                          tail=(v[5], v[6], v[7]))
    raise VgpError(f"fitted peak {index}: no readable name field")


def _backgrounds(b: bytes) -> list:
    ver, n = _header(b)
    if n == 0:
        return []
    _name, pos = _class_tag(b, 8)
    out = []
    for k in range(n):
        if k:
            if b[pos:pos + 2] != _TAG:
                raise VgpError(f"background {k} does not follow {k - 1}")
            pos += 2
        try:
            start, end = struct.unpack_from("<dd", b, pos)
            flag, avg, width, code = struct.unpack_from("<iidi", b, pos + 16)
            pos += 36
            groups = [(flag, avg, width)]
            if code == 1:            # Shirley: a second group follows
                it, tol, last = struct.unpack_from("<idi", b, pos)
                groups.append((it, tol, last))
                pos += 16
        except struct.error:
            raise VgpError("Backgrounds is cut short")
        out.append(BackgroundDef(index=k, start_ke=start, end_ke=end,
                                 avg_points=avg, width_ev=width, code=code,
                                 groups=tuple(groups)))
    if b[pos:].strip(b"\x00"):                   # (four zero bytes in real files)
        raise VgpError(f"Backgrounds has unread data after its {n} records")
    return out


def _data_spaces(b: bytes) -> list:
    ver, n = _header(b)
    strings = _cstrings(b)
    paths = [s for s in strings if s.upper().endswith(".VGD")]
    guids = [s for s in strings if s.startswith("{") and s.endswith("}")]
    if len(paths) != n:
        raise VgpError(f"Data Spaces lists {n} spectra but {len(paths)} "
                       "paths were found")
    return [DataSpaceRef(p, guids[i] if i < len(guids) else "")
            for i, p in enumerate(paths)]


# -- the project ------------------------------------------------------------------
def parse(data: bytes, path: str = "") -> Project:
    """Decode the bytes of a ``.vgp``. Raises :class:`VgpError`."""
    if not is_ole(data[:8]):
        raise VgpError("not an OLE2 compound file")
    try:
        ole = OleFile(data)
    except (OleError, struct.error, IndexError) as exc:
        raise VgpError(f"unreadable container: {exc}")
    names = {e if isinstance(e, str) else e[0] for e in ole.listdir()}
    if "Fitted Peaks/Contents" not in names or "Data Spaces/Contents" not in names:
        raise VgpError("not an Avantage project (no Fitted Peaks / Data Spaces)")
    prj = Project(path=path)
    libs = [s for s in _cstrings(_stream(ole, "Contents")) if s != "Processing Document"]
    prj.library = libs[0] if libs else ""
    prj.data_spaces = _data_spaces(_stream(ole, "Data Spaces/Contents"))
    prj.backgrounds = _backgrounds(_stream(ole, "Backgrounds/Contents"))
    prj.peaks = [_fitted_peak(i, r) for i, r in
                 enumerate(_fitted_records(_stream(ole, "Fitted Peaks/Contents")))]

    # links: peak -> data space, background -> peaks / data space
    ft = _tables(_stream(ole, "Fitted Peaks/XRefs"))
    if prj.peaks:
        for pk, ds in (ft[0] if ft else []):
            if 0 <= pk < len(prj.peaks) and 0 <= ds < len(prj.data_spaces):
                prj.peaks[pk].data_space = ds
    bt = _tables(_stream(ole, "Backgrounds/XRefs"))
    if prj.backgrounds:
        for bg, pk in (bt[0] if bt else []):
            if 0 <= bg < len(prj.backgrounds) and 0 <= pk < len(prj.peaks):
                prj.peaks[pk].background = bg
        for bg, ds in (bt[1] if len(bt) > 1 else []):
            if 0 <= bg < len(prj.backgrounds) and 0 <= ds < len(prj.data_spaces):
                prj.backgrounds[bg].data_space = ds
    unlinked = [p.name or f"#{p.index}" for p in prj.peaks if p.data_space is None]
    if unlinked:
        prj.notes.append("no spectrum is named for " + ", ".join(unlinked[:5])
                         + (" ..." if len(unlinked) > 5 else ""))
    unknown = sorted({b.code for b in prj.backgrounds
                      if b.code not in BACKGROUND_KINDS})
    if unknown:
        prj.notes.append("background type code(s) not known: "
                         + ", ".join(map(str, unknown)))
    return prj


def load(path: str) -> Project:
    """Read the project file at ``path`` (``.vgp``)."""
    from .base import read_bytes
    return parse(read_bytes(path), path)


# -- the project opened as a session ----------------------------------------------
SOURCE = "Avantage"
PRODUCT_TOL = 0.02        # stored area within this of the product's: product
MATCH_TOL = 0.05          # within this of either shape: the unit is understood
POOR_FIT = 0.10           # residual (share of the data range) said on the plot


def sniff(head: bytes, ext: str) -> bool:
    """An OLE2 file named ``.vgp`` or whose directory (it sits in the first
    4 KB of every project seen) lists the two streams a project has."""
    if not is_ole(head):
        return False
    return ext == ".vgp" or (
        "Fitted Peaks".encode("utf-16le") in head
        and "Data Spaces".encode("utf-16le") in head)


def project_in_folder(folder: str):
    """The one ``.vgp`` directly in ``folder`` (a project restored from a
    workbook, or opened as a folder), else None."""
    try:
        hits = [n for n in os.listdir(folder) if n.lower().endswith(".vgp")
                and not n.startswith(".")]
    except OSError:
        return None
    return os.path.join(folder, hits[0]) if len(hits) == 1 else None


def channel_step(region) -> float:
    """The width of one channel (eV): the spacing of the spectrum's points."""
    e = region.energy
    return abs(e[-1] - e[0]) / (len(e) - 1) if e and len(e) > 1 else 0.0


def shape_of(peak: FittedPeak, mode: str = "AV", step: float = 0.0) -> str:
    """The ``lineshapes.py`` shape string of a fitted peak: ``AV`` (the sum) or
    ``AVP`` (the product) of ``(mix %, CT, TM, ET, channel step)``."""
    t = peak.tail
    return (f"{mode}({100.0 * peak.mix:.6g},{t[0]:.6g},{t[1]:.6g},{t[2]:.6g},"
            f"{step:.6g})")


def residual_of(fit, region):
    """The reconstructed fit's rms residual as a share of the data range on
    ``region``, or None when it cannot be worked out."""
    dwell, scans = region.dwell_and_scans()
    cv = casafit.curves(fit, region.energy, region.counts,
                        region.photon_energy, dwell, scans)
    return cv[0].residual_rms if cv else None


def detect_mode(peaks, k: float, step: float, lo: float, hi: float):
    """``(mode, note)``: ``"AVP"`` (Avantage's product function) when the
    stored areas of ``peaks`` equal the product shape's area, window ``lo`` to
    ``hi`` (KE), to ``PRODUCT_TOL`` and the sum's do not fit better, else
    ``"AV"`` (the sum). ``k`` is dwell x scans (the heights are in counts).
    ``note`` is empty unless the stored areas match neither shape to
    ``MATCH_TOL`` -- they are then in a unit this module does not know and the
    mode could not be read (the sum is drawn)."""
    def ratio(mode):
        r = []
        for p in peaks:
            if not (p.area > 0 and p.height > 0 and p.fwhm > 0):
                continue
            shape = shape_of(p, mode, step)
            model = (p.height / k) * lineshapes.window_integral(
                shape, p.ke, p.fwhm, lo, hi)
            if model > 0:
                r.append(p.area / model)
        return sum(r) / len(r) if r else None
    rp, rs = ratio("AVP"), ratio("AV")
    if rp is None or rs is None:
        return "AV", ""
    dp, ds = abs(rp - 1.0), abs(rs - 1.0)
    if dp <= PRODUCT_TOL and dp < ds:
        return "AVP", ""
    if min(dp, ds) > MATCH_TOL:
        return "AV", ("the areas stored in the project match neither "
                      "Avantage's product nor its sum shape, so which was used "
                      "is not known: the sum is drawn")
    return "AV", ""


def fit_for(peaks, bg, region, dwell: float, scans: int):
    """The :class:`casafit.Fit` of one spectrum from its fitted ``peaks`` and
    ``bg`` (a :class:`BackgroundDef` or None). Nothing is guessed: the RSFs are
    left unknown (0), the project stores none that was read."""
    k = dwell * (scans or 1)
    step = channel_step(region)
    fit = casafit.Fit()
    if bg is None:
        reg = casafit.FitRegion(name=region.name, background="None",
                                rsf=0.0, line="Avantage project, no background")
        reg.start_ke = region.photon_energy - max(region.energy)
        reg.end_ke = region.photon_energy - min(region.energy)
    else:
        name = "Linear" if bg.code == 0 else "Shirley"
        note = ("" if bg.code in (0, 1) else
                f"the {bg.kind} background (type {bg.code}) is drawn as a "
                "Shirley background")
        reg = casafit.FitRegion(
            name=region.name, background=name, start_ke=bg.start_ke,
            end_ke=bg.end_ke, rsf=0.0, avg=max(1, int(bg.avg_points)),
            line=f"Avantage {bg.kind} background", note=note)
    mode, mode_note = detect_mode(peaks, k, step, reg.start_ke, reg.end_ke)
    if mode_note:
        reg.note = "; ".join(x for x in (reg.note, mode_note) if x)
    fit.regions.append(reg)
    for i, p in enumerate(peaks):
        shape = shape_of(p, mode, step)
        fit.components.append(casafit.FitComponent(
            name=p.name or f"Peak {i + 1}", shape=shape,
            area=lineshapes.area_for_height(shape, p.fwhm, p.height / k),
            fwhm=p.fwhm, pos_ke=p.ke, rsf=0.0, region=region.name,
            line=f"Avantage fitted peak, area {p.area:.6g} counts/s eV, "
                 f"{'product' if mode == 'AVP' else 'sum'} shape"))
    return fit


class ThermoProjectFile(ThermoExperiment):
    """An Avantage project: its ``.VGD`` spectra with the project's fits."""
    format_name = "Thermo Avantage project (.vgp)"

    def __init__(self):
        super().__init__()
        self.project = None             # the decoded Project
        self._order = {}

    def _rank_of(self):                 # the project's own order of spectra
        return lambda config, sample, scan: (self._order.get(scan, 10 ** 6),)

    def load(self, path: str, progress=None):
        path = os.path.abspath(path)
        vgp = project_in_folder(path) if os.path.isdir(path) else path
        if not vgp:
            raise ValueError("this folder does not hold exactly one Avantage "
                             "project (.vgp)")
        try:
            prj = load(vgp)
        except VgpError as exc:
            raise ValueError(f"{os.path.basename(vgp)} could not be read as an "
                             f"Avantage project ({exc})")
        self.path, self.project = path, prj
        self.root_dir = os.path.dirname(vgp)
        name = os.path.splitext(os.path.basename(vgp))[0]
        loaded, subs, seen = [], {}, set()
        n = len(prj.data_spaces)
        for i, ref in enumerate(prj.data_spaces):
            if progress and progress(i, n, os.path.basename(ref.path)):
                raise LoadCancelled()
            f = ref.resolve(self.root_dir)
            rel = os.path.relpath(f, self.root_dir)
            if not os.path.isfile(f):
                self.errors.append((rel, "the file is missing"))
                continue
            try:
                sub = ThermoVgdFile().load(f)
            except ValueError as exc:
                if "is empty" in str(exc):
                    self.skipped.append((rel, "aborted acquisition: no data"))
                else:
                    self.errors.append((rel, str(exc)))
                continue
            except Exception as exc:                    # noqa: BLE001
                self.errors.append((rel, str(exc) or type(exc).__name__))
                continue
            if sub.kind == "table":
                self.skipped.append((rel, "auto-height table, not a spectrum"))
                continue
            subs[i] = sub
            scan = os.path.splitext(os.path.basename(f))[0]
            self._order.setdefault(scan, i)
            loaded.append((f, rel, sub, name, ""))
            if f not in seen:
                seen.add(f)
                self.members.append((f, rel))
        self.members.insert(0, (vgp, os.path.basename(vgp)))
        self._attach_fits(prj, subs)
        self.warnings += [f"{name}: {note}" for note in prj.notes]
        self._merge(loaded)
        self.experiment = {"name": name, "project": "", "platter": "",
                           "date": ""}
        return self._finish()

    def _attach_fits(self, prj, subs):
        for di in prj.fitted_data_spaces():
            sub = subs.get(di)
            if sub is None:
                continue
            regions = [r for r in sub.regions if r.decodable and r.energy
                       and r.photon_energy]
            peaks = prj.peaks_of(di)
            bg = prj.background_of(di)
            # A .VGD can hold several traces (the two levels of a "Compare"
            # file): the file does not say which one a fit belongs to, so the
            # fit goes to the one it describes best (the first on a tie).
            cands = []
            for r in regions:
                if len(regions) == 1:
                    mine = peaks
                else:
                    lo = r.photon_energy - max(r.energy) - 1.5
                    hi = r.photon_energy - min(r.energy) + 1.5
                    mine = [p for p in peaks if lo < p.ke < hi]
                if not mine:
                    continue
                dwell, scans = r.dwell_and_scans()
                if not dwell:
                    self.warnings.append(
                        f"{r.name}: the dwell time is unknown, so the fit "
                        "cannot be drawn (its heights are in counts)")
                    continue
                fit = fit_for(mine, bg, r, dwell, scans)
                cands.append((r, mine, fit, residual_of(fit, r)))
            if not cands:
                continue
            best = min(range(len(cands)), key=lambda i: (
                cands[i][3] if cands[i][3] is not None else float("inf"), i))
            r, mine, fit, rms = cands[best]
            extra = []
            if best:
                extra.append(f"the project's fit matches trace {best + 1} of "
                             "this file; the first fits far worse")
            if rms is not None and rms > POOR_FIT:
                n = len(mine)
                extra.append(f"the project's fit leaves {100 * rms:.0f} % of "
                             "the range unexplained (the project lists "
                             f"{n} peak{'s' if n != 1 else ''} for this "
                             "spectrum)")
            if extra:
                reg = fit.regions[0]
                reg.note = "; ".join(x for x in [reg.note] + extra if x)
            r.fit = fit
            r.extra["fit_source"] = SOURCE
            r.extra["avantage_level"] = best
            r.extra["avantage_peaks"] = [
                {"name": p.name, "area": p.area, "height_counts": p.height,
                 "tail": list(p.tail)} for p in mine]
