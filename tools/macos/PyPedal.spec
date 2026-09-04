# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the PyPedal macOS application bundle.

Build with ``python tools/macos/build_app.py``. Output stays outside the
repository. Bundle identifier: org.pypedal.PyPedal

Version strings are read from pyproject.toml so a patch rebuild does not
need a separate spec edit.

Analysis, PYZ, EXE, COLLECT, and BUNDLE are provided by PyInstaller when
it executes this spec.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

SPECDIR = Path(SPECPATH).resolve()  # noqa: F821
sys.path.insert(0, str(SPECDIR))
from bundle_config import (  # noqa: E402
    BUNDLE_IDENTIFIER,
    BUNDLE_NAME,
    drop_collected_path,
    icns_path,
    info_plist,
    keep_submodule,
    read_project_version,
    repository_root,
)

REPO = repository_root()
VERSION = read_project_version(REPO)
ICON = icns_path(REPO)

hiddenimports = collect_submodules("PyPedal", filter=keep_submodule)
hiddenimports += [
    "numpy",
    "pandas",
    "scipy",
    "networkx",
    "PySide6",
    "PySide6.QtCore",
    "PySide6.QtGui",
    "PySide6.QtWidgets",
]

datas = [
    item
    for item in collect_data_files("PyPedal")
    if "examples" not in Path(item[0]).parts and "griffon" not in item[0].lower()
]

a = Analysis(  # noqa: F821
    [str(REPO / "PyPedal" / "desktop" / "__main__.py")],
    pathex=[str(REPO)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "customtkinter",
        "matplotlib",
        "IPython",
        "pytest",
        "PyPedal.examples",
    ],
    noarchive=False,
)

a.binaries = [item for item in a.binaries if not drop_collected_path(str(item[0]))]
a.datas = [item for item in a.datas if not drop_collected_path(str(item[0]))]

pyz = PYZ(a.pure)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=BUNDLE_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=True,
)

coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name=BUNDLE_NAME,
)

app = BUNDLE(  # noqa: F821
    coll,
    name=f"{BUNDLE_NAME}.app",
    icon=str(ICON) if ICON is not None else None,
    bundle_identifier=BUNDLE_IDENTIFIER,
    info_plist=info_plist(VERSION),
)
