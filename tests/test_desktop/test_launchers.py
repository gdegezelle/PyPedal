"""Final 4.2 launcher cutover: one desktop bootstrap, no pypedal-qt."""

from __future__ import annotations

import subprocess
import sys
import tomllib
from importlib.metadata import entry_points
from pathlib import Path

from _pedhelpers import REPO

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.desktop.main import main as desktop_main
from PyPedal.pyp_app import main as pyp_app_main


def _project() -> dict:
    with (Path(REPO) / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]


def _named_entry_points(group: str) -> dict[str, str]:
    selected = entry_points().select(group=group)
    return {item.name: item.value for item in selected}


def test_console_and_gui_scripts_share_desktop_main() -> None:
    scripts = _named_entry_points("console_scripts")
    gui_scripts = _named_entry_points("gui_scripts")
    assert scripts.get("pypedal") == "PyPedal.desktop.main:main"
    assert gui_scripts.get("pypedal-gui") == "PyPedal.desktop.main:main"
    assert "pypedal-qt" not in scripts
    assert "pypedal-qt" not in gui_scripts
    project = _project()
    assert "pypedal-qt" not in project.get("scripts", {})
    assert project["scripts"]["pypedal"] == "PyPedal.desktop.main:main"
    assert project["gui-scripts"]["pypedal-gui"] == "PyPedal.desktop.main:main"


def test_python_m_pypedal_version_matches_desktop_bootstrap() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "PyPedal", "--version"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == PYPEDAL_VERSION
    assert desktop_main(["--version"]) == 0
    assert pyp_app_main(["--version"]) == 0


def test_python_m_pypedal_version_does_not_create_qapplication() -> None:
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


def test_gui_extra_is_pyside6_not_customtkinter() -> None:
    extras = _project()["optional-dependencies"]
    assert extras["gui"] == extras["desktop"]
    blob = "\n".join(item for group in extras.values() for item in group).lower()
    assert "customtkinter" not in blob
    assert "pyside6" in " ".join(extras["all"]).lower()
    assert "customtkinter" not in " ".join(extras["all"]).lower()
    assert "pyinstaller" not in " ".join(extras["all"]).lower()
