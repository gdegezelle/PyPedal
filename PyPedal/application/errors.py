"""Application-level error policy.

PyPedal's scientific modules raise typed exceptions and never choose a
process exit status. Mapping those exceptions onto a console/bootstrap
exit code, and onto a short human-readable title plus text, is an
application decision.

This module does not create GUI dialogs. Callers keep the original
exception; presentation never discards it.
"""

from __future__ import annotations

from dataclasses import dataclass

from PyPedal import pyp_errors

# Distinct codes let a wrapping shell script or CI job branch on the kind
# of failure instead of parsing stderr. Values below 64 follow ordinary
# convention (2 = usage/input); those at 64 and above follow sysexits.h,
# where 70 is EX_SOFTWARE. Anything unrecognised gets 1. Zero is never
# used for a failure.
EXIT_STATUS: tuple[tuple[type[BaseException], int], ...] = (
    (pyp_errors.PyPedalInputError, 2),
    (pyp_errors.PyPedalUsageError, 2),
    (pyp_errors.PyPedalConfigurationError, 3),
    (pyp_errors.PyPedalDependencyError, 4),
    (pyp_errors.PyPedalNotImplementedError, 5),
    (pyp_errors.PyPedalValidationError, 65),
    (pyp_errors.PyPedalInternalError, 70),
    # 73 is EX_CANTCREAT: the output file could not be produced.
    (pyp_errors.PyPedalExportFormatError, 73),
    (pyp_errors.PyPedalError, 1),
)


@dataclass(frozen=True, slots=True)
class ApplicationErrorInfo:
    """GUI-independent presentation of an exception.

    ``exception`` is the original object. Downstream UI may format
    ``title`` and ``text`` but must not need to parse them to recover
    the failure.
    """

    title: str
    text: str
    exception: BaseException


def exit_status_for(exc: BaseException) -> int:
    """Map a PyPedal exception to the process exit status the CLI should use."""
    for klass, status in EXIT_STATUS:
        if isinstance(exc, klass):
            return status
    return 1


def describe_exception(exc: BaseException) -> ApplicationErrorInfo:
    """Return a short title and human-readable text for ``exc``.

    Unexpected exceptions are not swallowed or rewritten as PyPedal
    errors. They are labelled with their type name so a desktop can
    show them without hiding the original object.
    """
    text = str(exc)
    if not text:
        text = type(exc).__name__
    return ApplicationErrorInfo(
        title=type(exc).__name__,
        text=text,
        exception=exc,
    )
