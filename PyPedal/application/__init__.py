"""Qt-free application layer between a future desktop and the library.

``PyPedal.application`` owns session/load orchestration, browse access,
and GUI-independent error presentation. It does not implement scientific
formulas and it does not import GUI toolkits.

The intended dependency direction is::

    PyPedal.desktop
            ↓
    PyPedal.application
            ↓
    PyPedal scientific modules

Scientific ``pyp_*.py`` modules must not import this package or
``PyPedal.desktop``. Application adapters in 4.2-C may call scientific
APIs such as ``pyp_nrm``; that direction is allowed. The reverse is not.
"""

from __future__ import annotations

from PyPedal.application.errors import (
    EXIT_STATUS,
    ApplicationErrorInfo,
    describe_exception,
    exit_status_for,
)
from PyPedal.application.inbreeding import (
    UNKNOWN_YEAR_TOKENS,
    InbreedingByYearRow,
    group_inbreeding_by_year,
)
from PyPedal.application.load import (
    PedigreeOpenOptions,
    load_into_session,
    normalize_sepchar,
    resolve_source_path,
)
from PyPedal.application.session import PedigreeSession
from PyPedal.application.tables import (
    BROWSE_COLUMNS,
    PedigreeTableSource,
    TableColumn,
)

__all__ = [
    "BROWSE_COLUMNS",
    "EXIT_STATUS",
    "UNKNOWN_YEAR_TOKENS",
    "ApplicationErrorInfo",
    "InbreedingByYearRow",
    "PedigreeOpenOptions",
    "PedigreeSession",
    "PedigreeTableSource",
    "TableColumn",
    "describe_exception",
    "exit_status_for",
    "group_inbreeding_by_year",
    "load_into_session",
    "normalize_sepchar",
    "resolve_source_path",
]
