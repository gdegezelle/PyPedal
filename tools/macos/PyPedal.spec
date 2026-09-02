# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the 4.2-B macOS .app spike.

This is an engineering spike, not the 4.2-D distribution package.
Build with tools/macos/build_app.sh so output stays outside the repository.
Bundle identifier: org.pypedal.PyPedal

Analysis, PYZ, EXE, COLLECT, and BUNDLE are provided by PyInstaller when
it executes this spec.
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

SPECDIR = Path(SPECPATH).resolve()
REPO = SPECDIR.parents[1]

hiddenimports = collect_submodules(
    "PyPedal",
    filter=lambda name: "examples" not in name.split("."),
)
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
    excludes=["tkinter", "customtkinter", "matplotlib", "IPython"],
    noarchive=False,
)

pyz = PYZ(a.pure)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="PyPedal",
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
    name="PyPedal",
)

app = BUNDLE(  # noqa: F821
    coll,
    name="PyPedal.app",
    icon=None,
    bundle_identifier="org.pypedal.PyPedal",
    info_plist={
        "CFBundleName": "PyPedal",
        "CFBundleDisplayName": "PyPedal",
        "CFBundleIdentifier": "org.pypedal.PyPedal",
        "CFBundleShortVersionString": "4.1.0",
        "CFBundleVersion": "4.1.0",
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": "13.0",
    },
)
