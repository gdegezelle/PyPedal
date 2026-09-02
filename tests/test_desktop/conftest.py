"""Qt test isolation for the PySide6 desktop suite.

Redirect QSettings to a temporary directory and prefer the offscreen
platform plugin. Tests must not write the maintainer's real preferences.
"""

from __future__ import annotations

import os
import shutil
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def pytest_configure(config: object) -> None:
    try:
        from PySide6.QtCore import QCoreApplication, QSettings
    except ImportError:
        return
    root = tempfile.mkdtemp(prefix="pypedal-qsettings-")
    config._pypedal_qsettings_root = root  # type: ignore[attr-defined]
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, root)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.SystemScope, root)
    QCoreApplication.setOrganizationName("PyPedalTests")
    QCoreApplication.setApplicationName("PyPedalDesktopTests")


def pytest_unconfigure(config: object) -> None:
    root = getattr(config, "_pypedal_qsettings_root", None)
    if root:
        shutil.rmtree(root, ignore_errors=True)
