"""PHI / ULVAC-PHI MultiPak ``.spe`` reader.

Layout (verified against PHI's own VAMAS exports of the same data, for
float64 files from a PHI Model 5000 and float32 files from a Quantes): an
ASCII header ``SOFH ... EOFH`` of ``Key: value`` lines, then a binary section:
a 16-byte preamble (version, n regions, total block-header bytes, preamble
size), one 96-byte block header per region, then each region's intensities.
Values are counts per second (PHI's own VAMAS export leaves them as stored).

The 96-byte block header says how its data are stored (all 32-bit integers
except the text): ``+20`` number of points, ``+56`` the unit (``c/s``),
``+72`` the **data type** (``f8`` = float64, ``f4`` = float32), ``+76`` the
data length in bytes and ``+80`` the data offset from the start of the binary
section. The numbers are little-endian, except in files that Casa Software's
CasaXPS wrote for Thermo data, which are big-endian (the preamble tells).

``SpectralRegDef: idx area name Z npts step start end start end time PE ...``
gives the binding-energy axis (start, negative step). The header also holds
``SpectralRegDefFull:``, ``SpectralRegDef2:`` and similar lines that are not
region definitions and are ignored.
"""

from __future__ import annotations

import math
import re
import struct

from .base import (Region, SpectrumFile, analyser_mode_name, canon_region_name,
                   kv_from_lines, read_bytes)

PREAMBLE = 16
BLOCK_HEADER = 96            # bytes per region, all stored before the data
MAX_COUNTS = 1e12            # larger (or not finite) means a misread layout
_ITEM = {b"f8": ("d", 8), b"f4": ("f", 4)}      # block-header data-type tag
_DEF_RE = re.compile(r"^SpectralRegDef\s*:")
_SURVEY_NAMES = {"su1s", "xpssurvey"}            # PHI's names for a survey
_SUBSHELL_RE = re.compile(r"^([A-Z][a-z]?)(\d[spdf])(\d)$")
_NUMBER_RE = re.compile(r"^\d+(?:\.\d+)?$")
_SOURCE_RE = re.compile(r"^([A-Z][a-z]?)\s+\d")


def sniff(head: bytes, ext: str) -> bool:
    return head[:4] == b"SOFH"


def _num(s, default=None):
    m = re.search(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s or "")
    return float(m.group(0)) if m else default


def _photon_source(text):
    """(photon energy in eV or None, anode text) from ``XraySource``.

    The energy is the whole token that is a plain number, so ``KA1486 1486.68``
    gives 1486.68 and not the 1486 inside "KA1486". The anode is only named
    when the string starts with an element symbol and a number (``Al 1486.6
    mono`` -> "Al mono", as PHI's VAMAS export labels it); a source given as
    an energy only, or in another label, has no anode here."""
    text = (text or "").strip()
    hv = next((float(t) for t in text.split() if _NUMBER_RE.match(t)), None)
    m = _SOURCE_RE.match(text)
    anode = ""
    if m:
        anode = m.group(1) + (" mono" if re.search(r"\bmono", text, re.I)
                              else "")
    return hv, anode


def _region_name(name: str) -> str:
    if name.lower() in _SURVEY_NAMES:
        return "Survey"
    m = _SUBSHELL_RE.match(name)
    if m:                       # "Zr2p3" -> "Zr 2p3", as PHI's VAMAS export
        return f"{m.group(1)} {m.group(2)}{m.group(3)}"
    return canon_region_name(name)


def _date(text: str) -> str:
    m = re.match(r"^\s*(\d{4})\s+(\d{1,2})\s+(\d{1,2})\b", text or "")
    if m:
        return f"{m.group(1)}-{int(m.group(2)):02d}-{int(m.group(3)):02d}"
    return (text or "").strip()


def _spot_um(text):
    """The X-ray beam diameter in µm from ``100.0 um``; None for a unit that
    is not recognised (nothing is assumed)."""
    parts = (text or "").split()
    value = _num(parts[0]) if parts else None
    unit = parts[1].lower() if len(parts) > 1 else ""
    scale = {"um": 1.0, "µm": 1.0, "μm": 1.0, "mm": 1000.0}.get(unit)
    return value * scale if value and scale else None


class PhiSpeFile(SpectrumFile):
    format_name = "PHI MultiPak (.spe)"

    def load(self, path: str):
        self.path = path
        data = read_bytes(path)
        end = data.find(b"EOFH")
        if end < 0:
            raise ValueError("PHI .spe header terminator (EOFH) not found")
        header = data[:end].decode("latin-1").replace("\r", "").split("\n")
        pos = end + 4
        while data[pos:pos + 1] in (b"\r", b"\n"):
            pos += 1
        kv = kv_from_lines([l for l in header if ":" in l])
        defs = [line.split(":", 1)[1].split()
                for line in header if _DEF_RE.match(line)]
        if not defs:
            raise ValueError("no SpectralRegDef entries in the .spe header")
        n_areas = int(_num(kv.get("nospatialarea"), 1) or 1)
        n_regs = int(_num(kv.get("nospectralreg"), len(defs)) or len(defs))
        defs = defs[:n_regs]
        npts_all = [int(d[4]) for d in defs]
        tail = len(data) - pos

        bo = self._byte_order(data, pos, len(defs))
        layout = self._layout(data, pos, bo, npts_all)
        if layout is None:                  # block headers say nothing useful
            expect = (PREAMBLE + BLOCK_HEADER * len(defs)
                      + sum(n * 8 for n in npts_all))
            if tail < expect:
                raise ValueError(f"binary section is {tail} bytes; expected "
                                 f"{expect} for {len(defs)} region(s)")
            cursor = pos + PREAMBLE + BLOCK_HEADER * len(defs)
            layout = []
            for n in npts_all:
                layout.append((cursor, "d", 8))
                cursor += n * 8
        used = max(off + n * size for (off, _c, size), n
                   in zip(layout, npts_all)) - pos
        if tail > used:
            self.warnings.append(
                f"{tail - used} unexpected extra byte(s) after the spectra "
                f"(file has {n_areas} spatial area(s); only the first is read).")

        hv, anode = _photon_source(kv.get("xraysource"))
        power = _num(kv.get("xraypower"))
        wf = _num(kv.get("analyserworkfcn"))
        sample = ""
        if n_areas > 1:
            self.warnings.append("Several spatial areas: showing area 1 only.")
        date = _date(kv.get("acqfiledate") or kv.get("filedate") or "")
        desc = kv.get("filedesc", "")
        mode = analyser_mode_name(kv.get("analysermode"))
        spot = _spot_um(kv.get("xraybeamdiameter"))
        cal = kv.get("intensitycalcoeff", "").split()
        angle = _num(kv.get("sourceanalyserangle"))

        for idx, (d, (off, code, size)) in enumerate(zip(defs, layout)):
            name, npts = d[2], int(d[4])
            step, start = float(d[5]), float(d[6])
            dwell = _num(d[10] if len(d) > 10 else None)
            pe = _num(d[11] if len(d) > 11 else None)
            counts = list(struct.unpack_from(f"{bo}{npts}{code}", data, off))
            top = max((abs(c) for c in counts), default=0.0)
            if not math.isfinite(top) or top > MAX_COUNTS:
                raise ValueError(
                    f"region {name}: the values are not plausible counts "
                    f"(largest {top:g}); the byte order or data type of this "
                    "file is not understood")
            energy = [start + i * step for i in range(npts)]
            reg = Region(
                name=_region_name(name), index=idx, offset=idx,
                technique=kv.get("technique", "XPS") or "XPS",
                energy=energy, counts=counts, energy_label="Binding Energy",
                energy_units="eV", count_label="Intensity",
                count_units="counts/s",
                decodable=True, sample=sample, photon_energy=hv,
                pass_energy=pe, dwell=dwell, step=abs(step) or None,
                date=date, anode=anode)
            if power:
                reg.conditions["X-ray Power"] = f"{power:g} W"
            if spot:
                reg.conditions["X-ray spot (µm)"] = f"{spot:g}"
            if mode:
                reg.extra["analyser_mode"] = mode
            rows = {}
            if len(cal) >= 2:
                rows["Intensity calibration (a, b)"] = (
                    f"{cal[0]}, {cal[1]} (recorded, not applied)")
            if angle is not None:       # not the take-off angle: see CLAUDE.md
                rows["Source–analyser angle (°)"] = f"{angle:g}"
            if rows:
                reg.extra["acq_metadata"] = rows
            reg.extra["header"] = kv
            self.regions.append(reg)

        source = ""
        if anode:
            source = f"{anode} ({hv:g} eV)" if hv else anode
        version = kv.get("softwareversion", "")
        instr = {"Instrument": kv.get("instrumentmodel", ""),
                 "Operator": kv.get("operator", ""),
                 "Institution": kv.get("institution", ""),
                 "X-ray source": source}
        if wf:
            instr["Work function (eV)"] = f"{wf:g}"
        instr["Acquisition software"] = (
            "CasaXPS (MultiPak-format output)" if "casaxps" in desc.lower()
            else "PHI MultiPak" + (f" (software version {version})"
                                   if version else ""))
        ne = _num(kv.get("neutralizerenergy"))
        if ne:
            instr["Charge neutraliser"] = (
                f"on, {kv.get('neutralizerenergy', '')}, "
                f"{kv.get('neutralizercurrent', '')}").strip(", ")
        sc = _num(kv.get("sputtercurrent"), 0)
        if sc:
            instr["Ion gun / sputtering"] = (
                f"{kv.get('sputterion', '')} {kv.get('sputterenergy', '')} "
                f"{kv.get('sputtercurrent', '')}").strip()
        if desc:
            instr["Sample description"] = desc
        if kv.get("comments"):
            instr["Comments"] = kv["comments"]
        import sputter
        self.sputter_hint = sputter.from_text(" ".join(
            str(kv.get(k, "")) for k in ("sputterion", "sputterenergy",
                                         "sputtercurrent", "sputterraster",
                                         "sputterrastersize")))
        self.instrument = {k: v for k, v in instr.items() if v}
        return self._finish()

    @staticmethod
    def _byte_order(data, pos, n_regs):
        """'<' or '>': the order in which the preamble holds the number of
        regions. A preamble that matches neither is read as little-endian, as
        before (the length and plausibility checks then decide)."""
        if len(data) >= pos + 8:
            if struct.unpack_from("<i", data, pos + 4)[0] == n_regs:
                return "<"
            if struct.unpack_from(">i", data, pos + 4)[0] == n_regs:
                return ">"
        return "<"

    @staticmethod
    def _layout(data, pos, bo, npts_all):
        """[(absolute offset, struct code, item size)] per region from the
        block headers, or None when they do not describe the data (no type
        tag, a length that is not points x item size, data past the end of
        the file), in which case the caller uses the contiguous float64
        layout. A type tag this reader does not know is an error."""
        out = []
        for k, npts in enumerate(npts_all):
            h = pos + PREAMBLE + BLOCK_HEADER * k
            if h + BLOCK_HEADER > len(data):
                return None
            tag = data[h + 72:h + 76].split(b"\0")[0]
            if not tag:
                return None
            if tag not in _ITEM:
                raise ValueError(
                    f"unsupported data type {tag.decode('latin-1')!r} in the "
                    "PHI .spe block header (float32 'f4' and float64 'f8' "
                    "are read)")
            code, size = _ITEM[tag]
            n_hdr = struct.unpack_from(f"{bo}i", data, h + 20)[0]
            length = struct.unpack_from(f"{bo}i", data, h + 76)[0]
            off = struct.unpack_from(f"{bo}i", data, h + 80)[0]
            if (n_hdr != npts or length != npts * size or off < 0
                    or pos + off + length > len(data)):
                return None
            out.append((pos + off, code, size))
        return out
