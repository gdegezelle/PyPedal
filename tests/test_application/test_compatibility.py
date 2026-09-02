"""GUI-toolkit isolation and pyp_app compatibility launcher."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

from _pedhelpers import REPO

from PyPedal.application import EXIT_STATUS as APP_EXIT_STATUS
from PyPedal.application import exit_status_for as app_exit_status_for
from PyPedal.application import normalize_sepchar as app_normalize_sepchar
from PyPedal.desktop.main import main as desktop_main
from PyPedal.pyp_app import EXIT_STATUS as SHIM_EXIT_STATUS
from PyPedal.pyp_app import exit_status_for as shim_exit_status_for
from PyPedal.pyp_app import main as pyp_app_main
from PyPedal.pyp_errors import PyPedalUsageError


def test_application_import_does_not_load_gui_toolkits():
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            """
import sys
before = set(sys.modules)
from PyPedal import application
after = set(sys.modules)
new = after - before
blocked = []
for name in new:
    root = name.split('.', 1)[0]
    if root in {'tkinter', 'customtkinter', 'PySide6', 'PySide2', 'PyQt6', 'PyQt5', 'PyQt'}:
        blocked.append(name)
if blocked:
    raise SystemExit('gui modules imported: ' + ', '.join(sorted(blocked)))
print('ok')
print(application.__name__)
""",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ok" in proc.stdout


def test_pyp_app_reexports_application_helpers():
    assert shim_exit_status_for is app_exit_status_for
    assert SHIM_EXIT_STATUS is APP_EXIT_STATUS
    assert shim_exit_status_for(PyPedalUsageError("x")) == 2
    assert app_normalize_sepchar(", ") == ","


def test_pyp_app_main_delegates_to_desktop():
    assert pyp_app_main(["--version"]) == 0
    assert desktop_main(["--version"]) == 0


def test_pyp_app_has_no_tk_or_customtkinter():
    source = (Path(REPO) / "PyPedal" / "pyp_app.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.name.split(".", 1)[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split(".", 1)[0])
    assert "customtkinter" not in names
    assert "tkinter" not in names
    assert "PySide6" not in names
    assert "CTk" not in source
    assert "class PyPedalApp" not in source


def test_plain_pypedal_import_does_not_require_gui_extra() -> None:
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            """
import sys
import PyPedal
import PyPedal.application as application
assert PyPedal.__version__
assert application.PedigreeSession is not None
blocked = [name for name in ('customtkinter', 'tkinter', 'PySide6') if name in sys.modules]
if blocked:
    raise SystemExit('gui modules imported: ' + ', '.join(blocked))
print('ok')
""",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ok" in proc.stdout
