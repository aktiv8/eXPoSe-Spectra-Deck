"""Kratos Vision2 ``.dset`` reader (the binary dataset the instrument writes).

``DumpDataset`` turns a ``.dset`` into the ``.kal`` text that
:mod:`readers.kratos_kal` reads. This module decodes the ``.dset`` itself and
rebuilds the same ``{field name: value text}`` objects, so the field logic
(anode, aperture, neutraliser, sweeps, transmission, ...) stays in one place
and a ``.dset`` gives exactly the regions of its ``.kal``.

Layout (reverse-engineered from 36 ``.dset`` / ``.kal`` pairs, 325 objects;
every number below was checked against the ``.kal`` of the same data):

* everything is **big-endian 32-bit words**; the file starts
  ``00 00 DE 01 00 00 00 02`` and a 0x0B10-byte index (object names with
  duplicated offsets) that this reader does not need;
* the data blocks follow back to back, ``next = offset + 4 + size``, and
  normally end exactly at the end of the file; word 5 of the index is their
  number. A transmission-function run can leave a section of its own after
  them (two of 18 files seen; its header is not block-shaped): it is skipped,
  with a note, only when every indexed block was read, so a file whose chain
  stops early (``NPL_Trans``: 65 of 183) is still refused. A block is ``size``, the object's
  ordinal (the ``/N`` of the ``.kal``'s ``Object name``), ``0``, ``0``, then
  records ``id, value`` with ascending ids, closed by ``0, 0``;
* **the ids are the numbers the ``.kal`` prints** (``12`` Ordinate values,
  ``3`` Spectrum scan start, ``3087`` Stage Position Name ...);
* values: an enumeration or integer is one word, a float64 two words, a
  string a length (characters + NUL) then one UTF-32 character per word,
  an array a count then that many int32 / **float32** (the counts) /
  float64 (the transmission function) values;
* id ``5587`` (the transmission function object) is a container: three
  header words, its own records, then its own ``0, 0``.

The meaning of an *enumeration* number (which anode, which lens mode) is a
Kratos constant that is only known from pairs seen so far: an unseen number
is left out and named in a warning, never guessed (an unknown anode then
gives a kinetic-energy axis, as when the file does not record the source).
"""

from __future__ import annotations

import os
import struct

from .base import read_bytes
from .kratos_kal import KratosKalFile

MAGIC = b"\x00\x00\xde\x01\x00\x00\x00\x02"
FIRST_BLOCK = 0x0B10          # all 36 reference files; verified, see _chain
INDEX_COUNT = 20              # word 5 of the index: how many data blocks there are

# id -> (kal field name, type, unit text)
#   e enumeration  i int32  d float64  s string
#   ai int32 array  af float32 array  ad float64 array  c container
IDS = {
    1: ('Technique', 'e', ''),
    2: ('Scan type', 'e', ''),
    3: ('Spectrum scan start', 'd', 'eV'),
    4: ('Spectrum scan step size', 'd', 'eV'),
    5: ('Abscissa label', 's', ''),
    6: ('Abscissa units', 's', ''),
    7: ('Dwell time', 'd', 'seconds'),
    8: ('Ordinate label', 's', ''),
    9: ('Ordinate units', 's', ''),
    12: ('Ordinate values', 'af', ''),
    25: ('# points per line in map', 'i', ''),
    26: ('# lines in map', 'i', ''),
    37: ('Acquisition name', 's', ''),
    38: ('State change type', 'e', ''),
    42: ('Pass energy', 'd', 'eV'),
    51: ('Map energy/mass', 'd', 'eV'),
    59: ('Loop index', 'i', ''),
    99: ('# Sweeps completed', 'i', ''),
    151: ('Date Acquired', 's', ''),
    3001: ('start x coord', 'd', ''),
    3002: ('start y coord', 'd', ''),
    3003: ('step size x coord', 'd', ''),
    3004: ('step size y coord', 'd', ''),
    3005: ('Full Scale Deflection X', 'd', 'mm'),
    3006: ('Full Scale Deflection Y', 'd', 'mm'),
    3009: ('Real time display data index', 'i', ''),
    3011: ('Energy/line scan # steps', 'i', ''),
    3013: ('Enabled detectors', 'ai', ''),
    3016: ('Save data flag', 'e', ''),
    3017: ('Sum sweeps flag', 'e', ''),
    3045: ('Type Of Analyser', 'e', ''),
    3046: ('Analyser Scan Mode', 'e', ''),
    3047: ('Analyser Pass Energy', 'e', ''),
    3049: ('HSA Lens Mode', 'e', ''),
    3050: ('MHSA Lens Mode', 'e', ''),
    3065: ('Charge Neutraliser Charge Balance', 'd', ''),
    3066: ('Charge Neutraliser Filament Current', 'd', ''),
    3067: ('Charge Neutraliser Filament Bias', 'd', ''),
    3068: ('Magnet Lens Trim Coil', 'd', ''),
    3070: ('Neutraliser Switch State', 'e', ''),
    3073: ('Tuning Mode Flag', 'i', 'Flag'),
    3074: ('Reference Energy', 'd', 'eV'),
    3075: ('Max - Min Counts at Reference Energy', 'd', 'Counts'),
    3076: ('Half-Width at Reference Energy', 'd', 'eV'),
    3077: ('Peak Area at Reference Energy', 'd', ''),
    3078: ('Intensity at Reference Energy', 'd', 'Counts'),
    3080: ('Xray Reference Energy', 'e', ''),
    3087: ('Stage Position Name', 's', ''),
    3088: ('Stage X Position', 'd', 'm'),
    3089: ('Stage Y Position', 'd', 'm'),
    3083: ('# lines of map completed', 'i', ''),
    3090: ('Stage Z Position', 'd', 'm'),
    3091: ('Stage X Rotation', 'd', 'm'),
    3096: ('pixel square per point', 'i', ''),
    3097: ('Map type', 'e', ''),
    3102: ('Active Analyser Type', 'e', ''),
    3103: ('Manually Acquired Data Saved', 'i', ''),
    3105: ('Analyser Energy for required energy at key detector', 'd', 'eV'),
    3106: ('Stigmatic Image Shift X', 'd', 'm'),
    3107: ('Stigmatic Image Shift Y', 'd', 'm'),
    3112: ('Work function', 'd', ''),
    3113: ('Chemical symbol or formula', 's', ''),
    3114: ('Transition or charge state', 's', ''),
    3190: ('NICPU X-ray Gun Filament', 'e', ''),
    3191: ('NICPU X-ray Gun Active Anode Material', 'e', ''),
    3192: ('NICPU X-ray Gun Emission Current', 'd', 'A'),
    3193: ('NICPU X-ray Gun Anode HT Voltage', 'd', 'V'),
    3194: ('NICPU X-ray Gun Focus Voltage', 'd', 'V'),
    3201: ('Magnet Lens Trim Coil 2', 'd', ''),
    3202: ('Magnet Lens Trim Coil 3', 'd', ''),
    3244: ('Origin of platen', 'i', ''),
    3245: ('Type of platen', 'i', ''),
    3247: ('NICPU X-ray Gun Suppressor or Bias Voltage', 'd', 'V'),
    3250: ('Auto Z coordinate', 'i', ''),
    3254: ('Auto Z ordinate', 'i', ''),
    3255: ('Z Step size (mm) for auto z', 'd', ''),
    3256: ('Number of steps in auto z', 'i', ''),
    3261: ('This an index to a row in either position tables', 'i', ''),
    3278: ('NICPU Ion Gun Standby Delay Time', 'd', 'seconds'),
    3282: ('Descriptor for aperture size used in acquisition', 's', ''),
    3283: ('Descriptor for iris position used in acquisition', 's', ''),
    3290: ('Vision Software Version', 's', ''),
    # ion-gun "PAH" records of an etch block: ids far above the others
    # (0x2020_0CE0 ...), still ascending inside the block
    0x20200CE0: ('NICPU Ion Gun PSU PAH Gun Mode', 'i', ''),
    0x20200CE7: ('NICPU Ion Gun PSU PAH Action On Completion', 'e', ''),
    0x20500CDB: ('NICPU Ion Gun PSU PAH Wien Voltage', 'd', 'V'),
    0x20500CDC: ('NICPU Ion Gun PSU PAH Beam Bend Voltage', 'd', 'V'),
    0x20500CE3: ('NICPU Ion Gun PSU PAH Beam Monitor Current', 'd', 'A'),
    5587: ('Transmission Function Object (ke,t)', 'c', ''),
    5615: ('Transmission Function Kinetic Energy', 'ad', ''),
    5616: ('Transmission Function Value', 'ad', ''),
}

# id -> {number: Kratos constant}, every pair seen in the reference files
ENUMS = {
    1: {1: 'F_ISS', 3: 'F_XPS'},
    2: {0: 'F_SPECTRUM', 2: 'F_MAPPING'},
    38: {0: 'F_ETCH', 1: 'F_POSITION', 2: 'F_COUNTER', 9: 'F_ION_GUN_GAS',
         10: 'F_DELAY'},
    3016: {1: 'F_SAVE'},
    3017: {1: 'F_SUM'},
    3045: {0: 'F_HEMISPHERICAL', 1: 'F_MIRROR_HEMISPHERICAL'},
    3046: {0: 'F_FAT'},
    3047: {0: 'F_FAT_PASS_ENERGY_5_EV', 1: 'F_FAT_PASS_ENERGY_10_EV',
           2: 'F_FAT_PASS_ENERGY_20_EV', 3: 'F_FAT_PASS_ENERGY_40_EV',
           4: 'F_FAT_PASS_ENERGY_80_EV', 5: 'F_FAT_PASS_ENERGY_160_EV',
           6: 'F_FAT_PASS_ENERGY_320_EV'},
    3049: {0: 'F_HSA_LENS_HYBRID', 3: 'F_HSA_LENS_ELECTROSTATIC',
           6: 'F_HSA_LENS_ISS'},
    3050: {0: 'F_MHSA_LOW_MAGN', 1: 'F_MHSA_MEDIUM_MAGN', 2: 'F_MHSA_HIGH_MAGN'},
    3070: {3: 'F_NEUTRALISER_MANUAL_SETTINGS'},
    3080: {2: 'F_REFER_TO_NONE', 8: 'F_REFER_TO_XRAY_MONO_AL'},
    0x20200CE7: {0: 'F_ION_GUN_STANDBY'},
    3097: {1: 'F_MAP_STIGMATIC'},
    3102: {0: 'F_HEMISPHERICAL', 1: 'F_MIRROR_HEMISPHERICAL'},
    3190: {0: 'F_NICPU_XRAY_PSU_FILAMENT_MONO_1'},
    3191: {2: 'F_NICPU_XRAY_ANODE_STD_MONO'},
}

_MAX_ID = 0x2000              # ids seen are < 6000; used only to guess skips


def sniff(head: bytes, ext: str) -> bool:
    return head[:len(MAGIC)] == MAGIC


def is_index_only(path: str) -> bool:
    """A ``.dset`` that is just the index (the experiment-level file that
    sits beside the numbered datasets): nothing to load."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(len(MAGIC))
        return head == MAGIC and os.path.getsize(path) <= FIRST_BLOCK
    except OSError:
        return False


class _Unsynced(ValueError):
    """A block cannot be decoded (unknown id of unknown shape, truncation)."""


def _chain(d: bytes):
    """``([(offset, size)], end)``: the data blocks that chain from the first
    one, and the offset where the chain stops (``len(d)`` when it ends
    exactly at the end of the file)."""
    out, pos = [], FIRST_BLOCK
    while pos + 4 <= len(d):
        size = struct.unpack_from(">I", d, pos)[0]
        if size < 16 or size % 4 or pos + 4 + size > len(d):
            break
        out.append((pos, size))
        pos += 4 + size
    return out, pos


def _block_shaped(d: bytes, pos: int) -> bool:
    """A data block's header is ``size, ordinal, 0, 0``. What stops the chain
    in a truncated file is such a header with too little behind it; the
    trailing section of a transmission-function file is not (``size, 1024,
    25792, 0`` ...)."""
    if pos + 16 > len(d):
        return False
    _size, ordinal, z1, z2 = struct.unpack_from(">4I", d, pos)
    return ordinal != 0 and z1 == 0 and z2 == 0


def _is_text(w, pos, n):
    """w[pos:pos+n] is n-1 printable characters then a NUL."""
    if n < 1 or pos + n > len(w) or w[pos + n - 1] != 0:
        return False
    return all(9 <= c < 0x110000 and c not in (0xFFFE, 0xFFFF)
               for c in w[pos:pos + n - 1])


def _plausible(w, p, after):
    """w[p] can start the next record (or the ``0, 0`` that closes a block)."""
    if p + 1 < len(w) and w[p] == 0 and w[p + 1] == 0:
        return True
    # ids ascend; past _MAX_ID only the known high-numbered ones are accepted
    return p < len(w) and after < w[p] and (w[p] < _MAX_ID or w[p] in IDS)


def _guess(w, pos):
    """Words taken by the record at ``pos`` whose id this reader does not
    know: the first of string, int32, float64, int32 array, float64 array
    that lands on a plausible next record."""
    i = w[pos]
    n = w[pos + 1] if pos + 1 < len(w) else 0
    cands = []
    if 1 < n < 4096 and _is_text(w, pos + 2, n):
        cands.append(pos + 2 + n)
    cands.append(pos + 2)
    cands.append(pos + 3)
    if 0 < n < 1_000_000:
        cands.append(pos + 2 + n)
        cands.append(pos + 2 + 2 * n)
    for p in cands:
        if p <= len(w) and _plausible(w, p, i):
            return p
    return None


def _fmt(x):
    return repr(int(x)) if x == int(x) and abs(x) < 1e15 else f"{x:.15g}"


def _fmt32(x):
    return str(int(x)) if x == int(x) else f"{x:.9g}"


class KratosDsetFile(KratosKalFile):
    format_name = "Kratos Vision (.dset)"

    def load(self, path: str):
        self.path = path
        d = read_bytes(path)
        if d[:len(MAGIC)] != MAGIC:
            raise ValueError("not a Kratos Vision2 .dset file")
        if len(d) <= FIRST_BLOCK:
            raise ValueError(
                "this Vision2 dataset file holds only its index, no data "
                "blocks; open the .dset that holds the spectra")
        chain, end = _chain(d)
        listed = struct.unpack_from(">I", d, INDEX_COUNT)[0]
        if end != len(d) and not (
                chain and len(chain) == listed and not _block_shaped(d, end)):
            raise ValueError(
                "the .dset data blocks do not chain to the end of the file "
                f"({len(chain)} of the {listed} objects the index lists could "
                "be read; damaged, truncated or a layout not seen yet)")
        if end != len(d):
            # every indexed object is read; what follows is a section of its
            # own (the processing data of a transmission-function run)
            self.warnings.append(
                f"{len(d) - end:,} bytes after the last data block are not "
                f"part of the {listed} objects the index lists and were not "
                "read; all of those objects were loaded.")
        self._dataset = os.path.basename(path)
        unknown = {}
        objects = []
        for off, size in chain:
            w = struct.unpack_from(">%dI" % (size // 4), d, off + 4)
            try:
                fields = self._block(w, unknown)
            except (_Unsynced, struct.error, IndexError) as e:
                self.warnings.append(
                    f"Block at byte {off}: not decoded ({e}).")
                continue
            acq = fields.get("Acquisition name") or "object"
            fields["_object"] = f"{acq}/{w[0]}"
            objects.append(fields)
        for (name, n), count in sorted(unknown.items()):
            self.warnings.append(
                f"{name} = {n} (in {count} object{'s' if count != 1 else ''}) "
                "is a Kratos constant this reader has not met; the setting "
                "is left out.")
        return self._from_objects(objects)

    # ------------------------------------------------------------------
    def _block(self, w, unknown):
        """One block's words -> ``{field name: kal-style value text}``."""
        out = {}
        self._records(w, 3, out, unknown)
        return out

    def _records(self, w, pos, out, unknown):
        """Decode records from ``pos`` up to the closing ``0, 0``; returns
        the position after it."""
        n_w = len(w)
        while True:
            if pos + 1 >= n_w:
                raise _Unsynced("no closing marker")
            i = w[pos]
            if i == 0:
                return pos + 2
            spec = IDS.get(i)
            if spec is None:
                nxt = _guess(w, pos)
                if nxt is None:
                    raise _Unsynced(f"unknown field {i} of unknown shape")
                pos = nxt
                continue
            name, typ, unit = spec
            if typ == "c":
                pos = self._records(w, pos + 4, out, unknown)
                continue
            if typ in ("e", "i"):
                v = w[pos + 1]
                pos += 2
                if typ == "e":
                    text = ENUMS.get(i, {}).get(v)
                    if text is None:
                        unknown[(name, v)] = unknown.get((name, v), 0) + 1
                        continue
                    out[name] = text
                else:
                    v = struct.unpack(">i", struct.pack(">I", v))[0]
                    out[name] = str(v) + (f" {unit}" if unit else "")
                continue
            if typ == "d":
                x = struct.unpack(">d", struct.pack(">II", w[pos + 1],
                                                    w[pos + 2]))[0]
                out[name] = _fmt(x) + (f" {unit}" if unit else "")
                pos += 3
                continue
            n = w[pos + 1]
            if typ == "s":
                if not _is_text(w, pos + 2, n):
                    raise _Unsynced(f"{name}: malformed string")
                out[name] = "".join(chr(c) for c in w[pos + 2:pos + 1 + n])
                pos += 2 + n
                continue
            words = 2 * n if typ == "ad" else n
            if pos + 2 + words > n_w:
                raise _Unsynced(f"{name}: array past the end of the block")
            raw = struct.pack(">%dI" % words, *w[pos + 2:pos + 2 + words])
            if typ == "af":
                vals = [_fmt32(x) for x in struct.unpack(">%df" % n, raw)]
            elif typ == "ai":
                vals = [str(x) for x in struct.unpack(">%di" % n, raw)]
            else:
                vals = [_fmt(x) for x in struct.unpack(">%dd" % n, raw)]
            out[name] = "{" + ", ".join(vals) + "}"
            pos += 2 + words
