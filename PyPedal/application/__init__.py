"""Qt-free application layer between a future desktop and the library.

``PyPedal.application`` owns session/load orchestration and
GUI-independent error presentation. It does not implement scientific
formulas and it does not import GUI toolkits.

The intended dependency direction is::

    PyPedal.desktop      (4.2-B, not present yet)
            ↓
    PyPedal.application
            ↓
    PyPedal scientific modules

Scientific ``pyp_*.py`` modules must not import this package.
"""

from __future__ import annotations

from PyPedal.application.errors import (
    EXIT_STATUS,
    ApplicationErrorInfo,
    describe_exception,
    exit_status_for,
)
from PyPedal.application.load import (
    PedigreeOpenOptions,
    load_into_session,
    normalize_sepchar,
    resolve_source_path,
)
from PyPedal.application.session import PedigreeSession

__all__ = [
    "EXIT_STATUS",
    "ApplicationErrorInfo",
    "PedigreeOpenOptions",
    "PedigreeSession",
    "describe_exception",
    "exit_status_for",
    "load_into_session",
    "normalize_sepchar",
    "resolve_source_path",
]
