"""Thin compatibility launcher for the PyPedal desktop.

``import PyPedal.pyp_app`` succeeds. ``main()`` delegates to
``PyPedal.desktop.main.main``. This module contains no widgets and does
not import CustomTkinter, tkinter, or PySide6.
"""

from __future__ import annotations

from PyPedal.application.errors import EXIT_STATUS, exit_status_for

__all__ = ["EXIT_STATUS", "exit_status_for", "main"]


def main(argv: list[str] | None = None) -> int:
    """Launch the PySide6 desktop. Returns an exit status."""
    from PyPedal.desktop.main import main as run_desktop_main

    return run_desktop_main(argv)
