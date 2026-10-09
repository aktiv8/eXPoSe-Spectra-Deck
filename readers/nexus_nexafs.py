"""NeXus ``.nxs`` reader for beamline NEXAFS scans (Diamond B07 / GDA).

Checked on four B07 files and the text exports (``b07-*_NEXAFS.dat``) of the
same scans: ``entry1`` holds one group per detector channel (``ca18b``,
``ca35b``, ``ca36b``) with the current in A (``<ch>/<ch>``, attribute
``units``), its amplifier ``gain`` and the monochromator readback
``pgm_energy`` (the x axis, attribute ``primary``); the scan command, times,
user name, program version, ring current per point and ~130 beamline readbacks
(``before_scan``) sit beside them. Every channel becomes one ``Region``:
photon energy (eV) against current (A).

Nothing is guessed. The readback is *not* a regular grid (it scatters by
0.001-0.02 eV around the set points), so it is kept as read. Which channel is
total electron yield or fluorescence yield, the polarisation and the sample are
not in the files, so channels are named by detector, the gain is recorded but
not applied (the currents are already in A) and the data are not normalised to
the ring current. Only HDF5 files with the extension ``.nxs`` are claimed (the
``.kfit`` reader accepts any HDF5 signature, so this one is registered before
it). An analyser (NXxps) ``.nxs`` is not read: no beamline file of that kind is
available to check against.
"""

from __future__ import annotations

import datetime as _dt
import io
import os
import re

import timing
from .base import Region, SpectrumFile, read_bytes

_HDF5_SIG = b"\x89HDF\r\n\x1a\n"
ABBREVIATIONS = {"DLS": "Diamond Light Source"}   # the source name as recorded
_SKIP = ("instrument", "before_scan", "user")      # groups that are not channels


def sniff(head: bytes, ext: str) -> bool:
    return ext == ".nxs" and head[:8] == _HDF5_SIG


def _h5py():
    try:
        import h5py
        return h5py
    except ImportError:
        raise ValueError(
            "NeXus .nxs files need the h5py package (pip install h5py, or "
            "run: python launch.py --reinstall).")


def _text(v):
    """A str from whatever h5py hands back for a string (bytes, 0-d or 1-element
    arrays); '' for anything else."""
    try:
        import numpy as np
        if isinstance(v, np.ndarray):
            v = v.reshape(-1)[0] if v.size else ""
    except ImportError:
        pass
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace").strip()
    return str(v).strip() if isinstance(v, str) else ""


def _read(group, path, default=None):
    """Dataset ``path`` below ``group`` as a Python/numpy value (None if absent)."""
    try:
        return group[path][()]
    except (KeyError, ValueError, TypeError):
        return default


def _string(group, path):
    return _text(_read(group, path, ""))


def _attr(ds, name):
    try:
        return _text(ds.attrs.get(name, ""))
    except Exception:                       # noqa: BLE001 - a bad attribute only
        return ""


def _number(v):
    try:
        import numpy as np
        return float(np.asarray(v).reshape(-1)[0])
    except (TypeError, ValueError, IndexError, ImportError):
        return None


_TIME_RE = re.compile(
    r"^(\d{4}-\d\d-\d\d)[T ](\d\d:\d\d:\d\d)(?:\.\d+)?\s*"
    r"(Z|[+-]\d\d(?::?\d\d)?)?$")


def parse_time(text):
    """``(wall 'YYYY-MM-DD HH:MM:SS', utc 'YYYY-MM-DD HH:MM:SS' or '')`` from a
    GDA time such as ``2022-03-01T13:43:41.080Z`` or ``...02:38:40.707+01``
    (an offset of whole hours is written without minutes). The wall time is
    what the file says; the UTC time is only given when the file states the
    offset. ('', '') for text that is not a time."""
    m = _TIME_RE.match((text or "").strip())
    if not m:
        return "", ""
    wall = f"{m.group(1)} {m.group(2)}"
    zone = m.group(3)
    if not zone:
        return wall, ""
    if zone == "Z":
        return wall, wall
    sign = -1 if zone[0] == "-" else 1
    digits = zone[1:].replace(":", "")
    minutes = int(digits[:2]) * 60 + (int(digits[2:4]) if len(digits) > 2 else 0)
    when = timing.parse_ts(wall) - sign * _dt.timedelta(minutes=minutes)
    return wall, when.strftime("%Y-%m-%d %H:%M:%S")


def _channels(entry):
    """``[(name, group)]`` of the detector channels of an entry: groups that
    hold a dataset of their own name (the current)."""
    h5py = _h5py()
    out = []
    for name, grp in entry.items():
        if (not isinstance(grp, h5py.Group) or name.startswith(_SKIP)
                or not isinstance(grp.get(name), h5py.Dataset)):
            continue
        out.append((name, grp))
    return out


def _is_scan(entry):
    command = _string(entry, "scan_command")
    return command.lower().startswith("scan pgm_energy") or any(
        "pgm_energy" in g for _n, g in _channels(entry))


class NexusNexafsFile(SpectrumFile):
    format_name = "NeXus NEXAFS (.nxs)"

    def load(self, path: str):
        self.path = path
        h5py = _h5py()
        try:
            fh = h5py.File(io.BytesIO(read_bytes(path)), "r")
        except OSError as exc:
            raise ValueError(f"{os.path.basename(path)} could not be read as "
                             f"HDF5: {exc}") from exc
        with fh as f:
            entries = [(n, g) for n, g in f.items()
                       if isinstance(g, h5py.Group)
                       and _attr(g, "NX_class") == "NXentry"]
            if not entries:
                raise ValueError(
                    f"{os.path.basename(path)} has no NXentry: it is HDF5 but "
                    "not NeXus.")
            scans = [(n, g) for n, g in entries if _is_scan(g)]
            if not scans:
                raise ValueError(
                    f"{os.path.basename(path)} is NeXus but holds no "
                    "photon-energy scan this reader recognises (NEXAFS scans "
                    "of Diamond beamlines only).")
            stem = os.path.splitext(os.path.basename(path))[0]
            for n, g in scans:
                self._add_entry(n, g, stem if len(scans) == 1
                                else f"{stem} {n}")
        if not self.regions:
            raise ValueError(f"{os.path.basename(path)}: no detector channel "
                             "matches the length of the photon-energy scan.")
        return self._finish()

    # -- one NXentry --------------------------------------------------------
    def _add_entry(self, name, entry, sample):
        channels = _channels(entry)
        energy = self._axis(entry, channels)
        if energy is None:
            self.warnings.append(f"{name}: no photon-energy readback found.")
            return
        n = len(energy)
        info = self._info(entry, name)
        ring = _read(entry, "instrument/ring_current/ring_current")
        if ring is None:
            for _c, g in channels:
                if "ring_current" in g:
                    ring = g["ring_current"][()]
                    break
        wall0, utc0 = parse_time(_string(entry, "start_time"))
        wall1, utc1 = parse_time(_string(entry, "end_time"))
        if not self.instrument:
            self.instrument = self._instrument(entry)
        # the readback scatters around the set points, so the step is the
        # mean one (the median of single steps is biased by that scatter)
        step = abs(energy[-1] - energy[0]) / (n - 1) if n > 1 else None
        group = f"{os.path.basename(self.path)}:{name}"
        for ch, grp in channels:
            ds = grp[ch]
            if ds.ndim != 1 or ds.shape[0] != n:
                self.warnings.append(
                    f"{ch} has {ds.shape} values for {n} photon energies: "
                    "left out.")
                continue
            units = _attr(ds, "units")
            reg = Region(
                name=ch, index=len(self.regions), offset=len(self.regions),
                technique="NEXAFS", energy=list(energy),
                counts=[float(v) for v in ds[()]],
                energy_label="Photon Energy", energy_units="eV",
                count_label="Current" if units == "A" else "Signal",
                count_units=units or "arb.", decodable=True,
                sample=sample, photon_energy=None,
                dwell=_number(_read(entry, f"instrument/{ch}/count_time")),
                step=step, date=wall0)
            reg.extra["n_scans"] = 1
            reg.extra["acq_mode"] = "Scan"
            if utc0 or wall0:
                reg.extra["t_start"] = utc0 or wall0
                if utc1 or wall1:
                    reg.extra["t_end"] = utc1 or wall1
                if utc0:
                    reg.extra["tz"] = "UTC"
            reg.extra["nexafs_group"] = group
            reg.extra["acq_metadata"] = self._rows(entry, ch, grp, ring, info)
            reg.extra["nexafs"] = dict(
                info, channel=ch,
                scan_command=_string(entry, "scan_command"))
            self.regions.append(reg)

    @staticmethod
    def _axis(entry, channels):
        """The monochromator readback (the dataset marked ``primary``, else the
        first ``pgm_energy`` found) as a list of floats, or None."""
        found = []
        for _c, g in channels:
            if "pgm_energy" in g:
                ds = g["pgm_energy"]
                found.append((0 if _attr(ds, "primary") == "1" else 1, ds))
        if not found:
            ds = _read(entry, "instrument/pgm_energy/pgm_energy")
            return [float(v) for v in ds] if ds is not None else None
        found.sort(key=lambda t: t[0])
        return [float(v) for v in found[0][1][()]]

    @staticmethod
    def _instrument(entry):
        beamline = _string(entry, "instrument/name")
        source = _string(entry, "instrument/source/name")
        source = ABBREVIATIONS.get(source, source)
        return {k: v for k, v in {
            "Instrument": (f"{beamline.upper()} ({source})" if beamline and source
                           else beamline.upper()),
            "Operator": _string(entry, "user01/username"),
            "Acquisition software": _string(entry, "program_name"),
            "X-ray source": _string(entry, "instrument/source/type"),
            "Experiment ID": _string(entry, "experiment_identifier"),
        }.items() if v}

    @staticmethod
    def _info(entry, name):
        """What is kept per scan for the reports (plain strings/numbers)."""
        before = entry.get("before_scan")

        def readback(key):
            if before is None or key not in before:
                return None
            g = before[key]
            for sub in (key, "value"):
                if sub in g:
                    return g[sub][()]
            return None

        out = {"entry": name, "scan_id": _string(entry, "scan_identifier"),
               "entry_id": _string(entry, "entry_identifier")}
        for key, label in (("pgm_grating", "grating"), ("pgm_mirror", "mirror")):
            out[label] = _text(readback(key))
        out["cff"] = _number(readback("pgm_cff"))
        out["slit_x"] = _number(readback("s1b_xsize"))
        out["slit_y"] = _number(readback("s1b_ysize"))
        return out

    @staticmethod
    def _rows(entry, ch, grp, ring, info):
        """Extra metadata rows (shown in Details, written to VAMAS comments)."""
        rows = {}
        mode = _string(entry, f"instrument/{ch}/mode")
        if mode:
            rows["Detector mode"] = mode
        coupling = _string(entry, f"instrument/{ch}/coupling")
        if coupling:
            rows["Detector coupling"] = coupling
        gain = grp.get("gain")
        if gain is not None and gain.shape and gain.shape[0]:
            vals = sorted({float(v) for v in gain[()]})
            unit = _attr(gain, "units")
            rows["Detector gain"] = (", ".join(f"{v:g}" for v in vals)
                                     + (f" {unit}" if unit else ""))
        if ring is not None and len(ring):
            mean = sum(float(v) for v in ring) / len(ring)
            rows["Ring current"] = (f"{mean:.1f} ({float(ring[0]):.1f} to "
                                    f"{float(ring[-1]):.1f}; unit not recorded)")
        if info.get("grating"):
            rows["Monochromator grating"] = info["grating"]
        if info.get("mirror"):
            rows["Monochromator mirror"] = info["mirror"]
        if info.get("cff") is not None:
            rows["Monochromator c_ff"] = f"{info['cff']:.3f}"
        if info.get("slit_x") is not None and info.get("slit_y") is not None:
            rows["Slit s1b (x × y)"] = (f"{info['slit_x']:g} × "
                                        f"{info['slit_y']:g} (unit not recorded)")
        cmd = _string(entry, "scan_command")
        if cmd:
            rows["Scan command"] = cmd
        return rows
