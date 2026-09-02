"""AST architecture gates for the Qt-free application layer."""

from __future__ import annotations

import ast
from pathlib import Path

from _pedhelpers import REPO

PACKAGE = Path(REPO) / "PyPedal"
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
        "qtpy",
        "QtCore",
        "QtWidgets",
        "QtGui",
    }
)

QT_TYPE_TOKENS = frozenset(
    {
        "QWidget",
        "QObject",
        "QModelIndex",
        "QVariant",
        "QApplication",
        "QMainWindow",
        "QAbstractTableModel",
        "QTableView",
        "QThread",
        "QSettings",
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


def _from_import_modules(tree: ast.AST) -> set[str]:
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.module:
                modules.add(node.module)
            for alias in node.names:
                modules.add(alias.name)
    return modules


def test_application_package_exists():
    assert (APPLICATION / "__init__.py").is_file()
    for name in ("session.py", "load.py", "tables.py", "errors.py", "jobs.py", "export.py"):
        assert (APPLICATION / name).is_file(), name


def test_application_does_not_import_gui_toolkits():
    offenders: list[str] = []
    for path in _iter_python(APPLICATION):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = _imported_names(tree)
        hits = sorted(name for name in imported if name.split(".", 1)[0] in GUI_ROOTS)
        if hits:
            offenders.append(f"{path.relative_to(PACKAGE)}: {hits}")
    assert offenders == []


def test_application_does_not_name_qt_types():
    offenders: list[str] = []
    for path in _iter_python(APPLICATION):
        source = path.read_text(encoding="utf-8")
        hits = sorted(token for token in QT_TYPE_TOKENS if token in source)
        if hits:
            offenders.append(f"{path.relative_to(PACKAGE)}: {hits}")
    assert offenders == []


def test_application_job_adapters_may_import_scientific_modules():
    """4.2-A forbade this because A had no job adapters.

    4.2-C introduces explicit adapters in ``jobs.py``. Application MAY
    import scientific modules. Scientific modules still MUST NOT import
    application or desktop.
    """
    jobs = APPLICATION / "jobs.py"
    tree = ast.parse(jobs.read_text(encoding="utf-8"), filename=str(jobs))
    imported = _imported_names(tree) | _from_import_modules(tree)
    scientific = {"pyp_nrm", "pyp_metrics"}
    hits = {name.split(".")[-1] for name in imported} & scientific
    assert hits == scientific


def test_scientific_modules_do_not_import_application_or_desktop():
    offenders: list[str] = []
    for path in _iter_python(PACKAGE):
        rel = path.relative_to(PACKAGE)
        if rel.parts[0] == "application":
            continue
        if rel.parts[0] == "examples":
            continue
        if path.name == "pyp_app.py":
            continue
        if path.name in {"__init__.py", "__main__.py", "__version__.py"}:
            continue
        if not path.name.startswith("pyp_"):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = _imported_names(tree) | _from_import_modules(tree)
        for name in imported:
            if name == "application" or name.startswith("application."):
                offenders.append(f"{rel}: {name}")
            if "PyPedal.application" in name or name.startswith("PyPedal.desktop"):
                offenders.append(f"{rel}: {name}")
            if name == "desktop" or name.startswith("desktop."):
                offenders.append(f"{rel}: {name}")
    assert offenders == []
