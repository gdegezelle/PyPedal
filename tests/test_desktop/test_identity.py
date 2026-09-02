"""Qt application identity for the PySide6 desktop."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QCoreApplication
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QApplication

from PyPedal.desktop.app import apply_application_identity, create_application, prepare_frozen_cwd

_TEST_ORG = "PyPedalTests"
_TEST_APP = "PyPedalDesktopTests"


def _restore_test_identity() -> None:
    QCoreApplication.setOrganizationName(_TEST_ORG)
    QCoreApplication.setApplicationName(_TEST_APP)
    QGuiApplication.setApplicationDisplayName(_TEST_APP)


def test_apply_application_identity_sets_qt_metadata() -> None:
    try:
        apply_application_identity()
        assert QCoreApplication.organizationName() == "PyPedal"
        assert QCoreApplication.applicationName() == "PyPedal"
        assert QGuiApplication.applicationDisplayName() == "PyPedal"
        assert QCoreApplication.applicationVersion() == "4.1.0"
    finally:
        _restore_test_identity()


def test_create_application_sets_native_identity() -> None:
    try:
        app = create_application()
        assert isinstance(app, QApplication)
        assert app.organizationName() == "PyPedal"
        assert app.applicationName() == "PyPedal"
        assert app.applicationDisplayName() == "PyPedal"
        assert app.applicationVersion() == "4.1.0"
    finally:
        _restore_test_identity()


def test_prepare_frozen_cwd_is_a_noop_when_not_frozen() -> None:
    assert prepare_frozen_cwd() is None
