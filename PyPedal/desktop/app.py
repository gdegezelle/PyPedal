"""QApplication bootstrap for the PyPedal Qt desktop."""

from __future__ import annotations

import os
import sys
from pathlib import Path

from PySide6.QtCore import QCoreApplication
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QApplication

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.desktop.main_window import MainWindow
from PyPedal.desktop.settings import APPLICATION_NAME, ORGANIZATION_NAME, DesktopSettings


def user_data_directory() -> Path:
    """Directory for logs when the process cwd is not a project tree."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "PyPedal"
    return Path.home() / ".pypedal"


def prepare_frozen_cwd() -> Path | None:
    """Give a frozen app a writable cwd so pedigree logs are not cwd-relative.

    Finder launches often start with ``/`` or a useless directory. Scientific
    loads still use the absolute pedigree path chosen by the user.
    """
    if not getattr(sys, "frozen", False):
        return None
    target = user_data_directory()
    target.mkdir(parents=True, exist_ok=True)
    os.chdir(target)
    return target


def apply_application_identity() -> None:
    """Set Qt application metadata used by menus, About, and QSettings.

    Call this before constructing ``QApplication``. The constructor may
    still copy ``argv[0]`` (``__main__.py`` under ``python -m``) into the
    application name, so ``create_application`` reapplies the same values
    on the instance afterwards. ``applicationDisplayName`` is what macOS
    uses for About / Hide / Quit once Qt owns the menu.
    """
    QCoreApplication.setOrganizationName(ORGANIZATION_NAME)
    QCoreApplication.setApplicationName(APPLICATION_NAME)
    QCoreApplication.setApplicationVersion(PYPEDAL_VERSION)
    QGuiApplication.setApplicationDisplayName(APPLICATION_NAME)


def create_application(argv: list[str] | None = None) -> QApplication:
    """Return a QApplication using native platform style.

    Does not force Fusion or install a custom theme. Does not create
    widgets. Identity is established before construction and again after
    so a ``python -m PyPedal.desktop`` ``argv[0]`` cannot stick as the
    application name.
    """
    apply_application_identity()
    existing = QApplication.instance()
    if isinstance(existing, QApplication):
        app = existing
    elif existing is None:
        app = QApplication(argv if argv is not None else sys.argv)
    else:
        raise RuntimeError("A non-QApplication Qt application already exists")
    app.setOrganizationName(ORGANIZATION_NAME)
    app.setApplicationName(APPLICATION_NAME)
    app.setApplicationDisplayName(APPLICATION_NAME)
    app.setApplicationVersion(PYPEDAL_VERSION)
    return app


def run_desktop(pedigree: Path | None = None) -> int:
    """Show the main window and run the Qt event loop."""
    prepare_frozen_cwd()
    app = create_application()
    window = MainWindow(DesktopSettings())
    window.show()
    if pedigree is not None:
        window.open_path(Path(pedigree))
    return app.exec()
