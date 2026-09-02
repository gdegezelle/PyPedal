"""Entry-point and import-laziness checks that do not need a QApplication."""

from __future__ import annotations

import subprocess
import sys

from _pedhelpers import REPO

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.desktop.main import main


def test_version_prints_and_exits_zero_without_qt() -> None:
    assert main(["--version"]) == 0


def test_module_version_subprocess_does_not_create_qapplication() -> None:
    script = (
        "import sys\n"
        "from PyPedal.desktop.main import main\n"
        "assert main(['--version']) == 0\n"
        "assert 'PySide6.QtWidgets' not in sys.modules\n"
        "try:\n"
        "    from PySide6.QtWidgets import QApplication\n"
        "except ImportError:\n"
        "    raise SystemExit(0)\n"
        "assert QApplication.instance() is None\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_python_m_desktop_version() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "PyPedal.desktop", "--version"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == PYPEDAL_VERSION


def test_self_test_runs_jobs_without_qapplication(tmp_path) -> None:
    from _pedhelpers import close_owned_pypedal_log_handlers

    pedigree = tmp_path / "mrode.ped"
    pedigree.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n", encoding="utf-8")
    try:
        assert main(["--self-test", str(pedigree)]) == 0
    finally:
        close_owned_pypedal_log_handlers()


def test_import_desktop_does_not_import_pyside() -> None:
    script = (
        "import sys\n"
        "import PyPedal.desktop\n"
        "assert not any(name == 'PySide6' or name.startswith('PySide6.') "
        "for name in sys.modules), sorted(n for n in sys.modules if 'PySide' in n)\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
