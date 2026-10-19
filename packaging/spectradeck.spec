# PyInstaller spec for the portable one-folder build (see build_exe.py).
# -*- mode: python ; coding: utf-8 -*-
import os

from PyInstaller.utils.hooks import collect_all

ROOT = os.path.dirname(SPECPATH)
NAME = "SpectraDeck"

datas, binaries, hiddenimports = [], [], []
for pkg in ("tkinterdnd2", "h5py", "pymupdf", "fitz", "reportlab", "pptx",
            "docx", "lxml"):
    try:
        d, b, h = collect_all(pkg)
    except Exception:           # an optional library that is not installed
        continue
    datas += d
    binaries += b
    hiddenimports += h

# our own data, placed so every module's os.path.dirname(__file__) lookup works
for folder in ("assets", "viewer"):
    datas.append((os.path.join(ROOT, folder), folder))

hiddenimports += ["matplotlib.backends.backend_tkagg",
                  "matplotlib.backends.backend_agg", "PIL.ImageTk"]

a = Analysis(
    [os.path.join(ROOT, "spectradeck.py")],
    pathex=[ROOT],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["tests", "manual", "pytest", "IPython", "notebook", "scipy",
              "matplotlib.backends.backend_qt5agg",
              "matplotlib.backends.backend_qtagg",
              "matplotlib.backends.backend_wxagg",
              "matplotlib.backends.backend_gtk3agg",
              "matplotlib.backends.backend_webagg"],
    noarchive=False,
)
pyz = PYZ(a.pure)
icon = os.path.join(SPECPATH, "app.ico")
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name=NAME,
    console=False,
    icon=icon if os.path.isfile(icon) else None,
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name=NAME, upx=False)
