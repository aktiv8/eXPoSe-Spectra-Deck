#!/usr/bin/env python3
"""
build_exe.py - build the portable Windows copy of eXPoSe SpectraDeck.

    python build_exe.py            build dist/SpectraDeck/ and a zip
    python build_exe.py --clean    rebuild the build environment first
    python build_exe.py --no-zip   leave the folder, skip the zip

It makes an isolated build environment (./.buildvenv, never ./.venv), installs
requirements-build.txt into it, and runs PyInstaller in one-folder mode. The
result is a folder holding SpectraDeck.exe and an _internal folder with
the Python runtime and every library as real .dll / .pyd files: unzip it
anywhere and double-click the exe. No installer, no registry, no Python needed.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import venv
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD_VENV = os.path.join(HERE, ".buildvenv")
REQS = os.path.join(HERE, "requirements-build.txt")
SPEC = os.path.join(HERE, "packaging", "spectradeck.spec")
ICO = os.path.join(HERE, "packaging", "app.ico")
DIST = os.path.join(HERE, "dist")
WORK = os.path.join(HERE, "build")
NAME = "SpectraDeck"

PORTABLE_TXT = """eXPoSe SpectraDeck - portable copy
==================================

Double-click SpectraDeck.exe. Nothing is installed and nothing needs
Python: the interpreter and all libraries are in the _internal folder (keep
it next to the exe).

Settings: by default they are kept in your user folder (as the installed
program does). To keep them beside this copy instead - for a USB stick - add
an empty file called portable.txt next to the exe; a data folder is then
created for the settings, calibration and log.

Windows may say "Windows protected your PC" the first time (the exe is not
code-signed): choose More info, then Run anyway.
"""

ICON_CODE = """
import sys
from PIL import Image
im = Image.open(sys.argv[1]).convert('RGBA')
s = min(im.size)
l, t = (im.width - s) // 2, (im.height - s) // 2
im = im.crop((l, t, l + s, t + s))
im.save(sys.argv[2], sizes=[(16, 16), (32, 32), (48, 48), (64, 64),
                            (128, 128), (256, 256)])
"""


def log(msg):
    print(f"[build] {msg}", flush=True)


def venv_python():
    return os.path.join(BUILD_VENV, "Scripts" if os.name == "nt" else "bin",
                        "python.exe" if os.name == "nt" else "python")


def ensure_env(clean):
    if clean and os.path.exists(BUILD_VENV):
        shutil.rmtree(BUILD_VENV)
    if not os.path.exists(venv_python()):
        log("Creating the build environment in ./.buildvenv ...")
        venv.create(BUILD_VENV, with_pip=True, clear=True)
    log("Installing build requirements ...")
    subprocess.run([venv_python(), "-m", "pip", "install", "--upgrade", "pip"],
                   check=False)
    subprocess.run([venv_python(), "-m", "pip", "install", "-r", REQS],
                   check=True)


def make_icon():
    """packaging/app.ico from assets/splash.png (best effort)."""
    src = os.path.join(HERE, "assets", "splash.png")
    if os.path.isfile(src):
        subprocess.run([venv_python(), "-c", ICON_CODE, src, ICO], check=False)


def version():
    sys.path.insert(0, HERE)
    import appinfo
    return appinfo.VERSION


def build():
    for path in (WORK, os.path.join(DIST, NAME)):
        if os.path.exists(path):
            shutil.rmtree(path)
    log("Running PyInstaller ...")
    subprocess.run([venv_python(), "-m", "PyInstaller", SPEC, "--noconfirm",
                    "--distpath", DIST, "--workpath", WORK], check=True,
                   cwd=HERE)
    out = os.path.join(DIST, NAME)
    with open(os.path.join(out, "PORTABLE.txt"), "w", encoding="utf-8") as fh:
        fh.write(PORTABLE_TXT)
    readme = os.path.join(HERE, "README.md")
    if os.path.isfile(readme):
        shutil.copy2(readme, out)
    return out


def zip_folder(folder):
    target = os.path.join(DIST, f"{NAME}-{version()}-win64.zip")
    log(f"Zipping to {target} ...")
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(folder):
            for f in files:
                full = os.path.join(root, f)
                zf.write(full, os.path.relpath(full, DIST))
    return target


def main():
    args = set(sys.argv[1:])
    ensure_env("--clean" in args)
    make_icon()
    out = build()
    log(f"Built {out}")
    if "--no-zip" not in args:
        log(f"Zip: {zip_folder(out)}")


if __name__ == "__main__":
    main()
