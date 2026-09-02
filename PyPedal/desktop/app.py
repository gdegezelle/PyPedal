"""QApplication bootstrap for the PyPedal Qt desktop."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.desktop.main_window import MainWindow
from PyPedal.desktop.settings import ORGANIZATION_NAME, DesktopSettings

APPLICATION_NAME = "PyPedal"


def create_application(argv: list[str] | None = None) -> QApplication:
    """Return a QApplication using native platform style.

    Does not force Fusion or install a custom theme.
    """
    existing = QApplication.instance()
    if isinstance(existing, QApplication):
        app = existing
    elif existing is None:
        app = QApplication(argv if argv is not None else sys.argv)
    else:
        raise RuntimeError("A non-QApplication Qt application already exists")
    app.setOrganizationName(ORGANIZATION_NAME)
    app.setApplicationName(APPLICATION_NAME)
    app.setApplicationVersion(PYPEDAL_VERSION)
    return app


def run_desktop(pedigree: Path | None = None) -> int:
    """Show the main window and run the Qt event loop."""
    app = create_application()
    window = MainWindow(DesktopSettings())
    window.show()
    if pedigree is not None:
        window.open_path(Path(pedigree))
    return app.exec()
