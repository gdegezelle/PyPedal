"""Architecture gates for the PySide6 desktop package.

Permanent rule: scientific modules must not import application or desktop.
Application must not import GUI toolkits.

The 4.2-A test that forbids application imports of analysis modules is an
A-scope restriction. 4.2-C will add application adapters that call
``pyp_nrm`` / ``pyp_metrics``. This file does not treat that A test as a
permanent invariant.
"""

from __future__ import annotations

import ast
from pathlib import Path

from _pedhelpers import REPO

PACKAGE = Path(REPO) / "PyPedal"
DESKTOP = PACKAGE / "desktop"
APPLICATION = PACKAGE / "application"

GUI_ROOTS = frozenset(
    {
        "customtkinter",
        "tkinter",
        "PySide6",
        "PySide2",
        "PyQt6",
        "PyQt5",
        "PyQt",
    }
)

ANALYSIS_ROOTS = frozenset(
    {
        "pyp_nrm",
        "pyp_metrics",
        "pyp_demog",
        "pyp_jbc",
        "pyp_network",
        "pyp_graphics",
        "pyp_reports",
        "pyp_snp",
        "pyp_db",
        "pyp_tests",
        "pyp_utils",
        "pyp_io",
        "pyp_chronology",
        "pyp_newclasses",
    }
)


def _iter_python(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*.py") if "__pycache__" not in path.parts)


def _imported_names(tree: ast.AST) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.name)
                names.add(alias.name.split(".", 1)[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
            names.add(node.module.split(".", 1)[0])
            for alias in node.names:
                names.add(f"{node.module}.{alias.name}")
    return names


def test_desktop_package_exists() -> None:
    assert (DESKTOP / "__init__.py").is_file()
    for name in ("main.py", "app.py", "main_window.py", "workers.py", "settings.py"):
        assert (DESKTOP / name).is_file(), name


def test_application_still_has_no_gui_imports() -> None:
    offenders: list[str] = []
    for path in _iter_python(APPLICATION):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = _imported_names(tree)
        hits = sorted(name for name in imported if name.split(".", 1)[0] in GUI_ROOTS)
        if hits:
            offenders.append(f"{path.relative_to(PACKAGE)}: {hits}")
    assert offenders == []


def test_desktop_does_not_import_analysis_modules() -> None:
    """Desktop reaches science through PyPedal.application, not pyp_* APIs."""
    offenders: list[str] = []
    for path in _iter_python(DESKTOP):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = _imported_names(tree)
        hits = sorted(
            name
            for name in imported
            if name in ANALYSIS_ROOTS or name.split(".")[-1] in ANALYSIS_ROOTS
        )
        if hits:
            offenders.append(f"{path.relative_to(PACKAGE)}: {hits}")
    assert offenders == []


def test_scientific_modules_do_not_import_desktop() -> None:
    offenders: list[str] = []
    for path in _iter_python(PACKAGE):
        rel = path.relative_to(PACKAGE)
        if rel.parts[0] in {"application", "desktop", "examples"}:
            continue
        if path.name in {"pyp_app.py", "__init__.py", "__main__.py", "__version__.py"}:
            continue
        if not path.name.startswith("pyp_"):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = _imported_names(tree)
        for name in imported:
            if name == "desktop" or name.startswith("desktop.") or "PyPedal.desktop" in name:
                offenders.append(f"{rel}: {name}")
    assert offenders == []


def test_desktop_init_does_not_import_pyside() -> None:
    source = (DESKTOP / "__init__.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported = _imported_names(tree)
    assert "PySide6" not in imported
    assert not any(name.startswith("PySide6") for name in imported)
