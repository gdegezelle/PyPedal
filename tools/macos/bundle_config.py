"""Shared macOS bundle identity and content-filter rules.

PyInstaller executes ``PyPedal.spec`` in a separate process. This module
is the single place for bundle identifier, version lookup, icon discovery,
and exclusion filters so E can rebuild at 4.2.0 by bumping project metadata.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

BUNDLE_IDENTIFIER = "org.pypedal.PyPedal"
BUNDLE_NAME = "PyPedal"
VOLUME_NAME = "PyPedal"
MINIMUM_MACOS = "13.0"

# Qt / PySide pieces the desktop does not import. Dropping them is a size
# and hygiene choice, not a scientific one.
_QT_DROP_MARKERS = (
    "QtWebEngine",
    "QtWebEngineCore",
    "QtWebEngineWidgets",
    "QtWebEngineQuick",
    "QtQml",
    "QtQuick",
    "QtQuick3D",
    "QtMultimedia",
    "Qt3D",
    "QtCharts",
    "QtDataVisualization",
    "QtRemoteObjects",
    "QtSensors",
    "QtBluetooth",
    "QtNfc",
    "QtPositioning",
    "QtLocation",
    "QtSerialPort",
    "QtWebSockets",
    "QtPdf",
    "QtSql",
    "QtTest",
    "QtDesigner",
    "QtHelp",
    "QtUiTools",
    "QtHttpServer",
    "QtSpatialAudio",
    "QtTextToSpeech",
    "QtWebView",
    "QtXmlPatterns",
)

_PLUGIN_DROP_MARKERS = (
    "/qml/",
    "/sqldrivers/",
    "/multimedia/",
    "/geoservices/",
    "/canbus/",
    "/sensors/",
    "/webview/",
    "/designer/",
)

_FORBIDDEN_BUNDLE_MARKERS = (
    "griffon",
    "customtkinter",
    "tkinter",
    ".git/",
    "/tests/",
    "docs/engineering",
)


def macos_dir() -> Path:
    return Path(__file__).resolve().parent


def repository_root() -> Path:
    return macos_dir().parents[1]


def read_project_version(root: Path | None = None) -> str:
    """Return ``[project].version`` from pyproject.toml."""
    path = (root or repository_root()) / "pyproject.toml"
    with path.open("rb") as handle:
        return str(tomllib.load(handle)["project"]["version"])


def icns_path(root: Path | None = None) -> Path | None:
    """Return ``tools/macos/PyPedal.icns`` when a maintainer supplies it."""
    candidate = macos_dir() / "PyPedal.icns"
    if root is not None:
        candidate = Path(root) / "tools" / "macos" / "PyPedal.icns"
    return candidate if candidate.is_file() else None


def info_plist(version: str) -> dict[str, str | bool]:
    """CFBundle keys for the PyPedal.app Info.plist."""
    return {
        "CFBundleName": BUNDLE_NAME,
        "CFBundleDisplayName": BUNDLE_NAME,
        "CFBundleIdentifier": BUNDLE_IDENTIFIER,
        "CFBundleExecutable": BUNDLE_NAME,
        "CFBundlePackageType": "APPL",
        "CFBundleShortVersionString": version,
        "CFBundleVersion": version,
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": MINIMUM_MACOS,
    }


def keep_submodule(name: str) -> bool:
    """Hidden-import filter for ``collect_submodules('PyPedal')``."""
    parts = name.split(".")
    if "examples" in parts:
        return False
    if "pyp_app" in parts:
        return False
    if "pyp_tests" in parts:
        return False
    return True


def drop_collected_path(destination: str) -> bool:
    """True when a PyInstaller dest path should not enter the bundle."""
    text = destination.replace("\\", "/")
    lowered = text.lower()
    if "griffon" in lowered:
        return True
    if "customtkinter" in lowered or "/tkinter/" in lowered or lowered.endswith("/tkinter"):
        return True
    for marker in _QT_DROP_MARKERS:
        if marker in text:
            return True
    for marker in _PLUGIN_DROP_MARKERS:
        if marker in lowered:
            return True
    return False


def forbidden_bundle_hits(paths: list[str]) -> list[str]:
    """Return bundled paths that must not ship in PyPedal.app."""
    hits: list[str] = []
    for raw in paths:
        lowered = raw.replace("\\", "/").lower()
        for marker in _FORBIDDEN_BUNDLE_MARKERS:
            if marker in lowered:
                hits.append(raw)
                break
    return hits
