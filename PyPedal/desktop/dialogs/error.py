"""Error presentation. Does not discard the original exception."""

from __future__ import annotations

from PySide6.QtWidgets import QMessageBox, QWidget

from PyPedal.application import describe_exception
from PyPedal.pyp_errors import PyPedalError


def show_application_error(
    parent: QWidget | None,
    exc: BaseException,
    details: str,
) -> None:
    """Show a typed PyPedal error, or a concise unexpected error with details."""
    info = describe_exception(exc)
    box = QMessageBox(parent)
    box.setWindowTitle(info.title)
    box.setText(info.text)
    if isinstance(exc, PyPedalError):
        box.setIcon(QMessageBox.Icon.Warning)
    else:
        box.setIcon(QMessageBox.Icon.Critical)
        box.setDetailedText(details)
    box.exec()


def show_load_error(
    parent: QWidget | None,
    exc: BaseException,
    details: str,
) -> None:
    """Compatibility alias used by the open-pedigree path."""
    show_application_error(parent, exc, details)
