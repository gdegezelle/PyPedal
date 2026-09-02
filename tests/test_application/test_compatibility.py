"""GUI-toolkit isolation and CustomTkinter compatibility wrappers."""

from __future__ import annotations

import subprocess
import sys

from PyPedal.application import EXIT_STATUS as APP_EXIT_STATUS
from PyPedal.application import exit_status_for as app_exit_status_for
from PyPedal.application import normalize_sepchar as app_normalize_sepchar
from PyPedal.pyp_app import EXIT_STATUS as CTK_EXIT_STATUS
from PyPedal.pyp_app import exit_status_for as ctk_exit_status_for
from PyPedal.pyp_app import normalize_sepchar as ctk_normalize_sepchar
from PyPedal.pyp_app import pedigree_open_options
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
    assert ctk_normalize_sepchar is app_normalize_sepchar
    assert ctk_exit_status_for is app_exit_status_for
    assert CTK_EXIT_STATUS is APP_EXIT_STATUS
    assert ctk_normalize_sepchar(", ") == ","
    assert ctk_exit_status_for(PyPedalUsageError("x")) == 2


def test_ctk_open_options_still_omit_pedigree_summary():
    opts = pedigree_open_options("/tmp/dogs.ped", "asdxbn", ", ", True)
    assert opts["sepchar"] == ","
    assert opts["messages"] == "quiet"
    assert "pedigree_summary" not in opts
    assert opts["pedformat"] == "asdxbn"
    assert opts["pedname"] == "dogs.ped"


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
