"""Application error presentation and process-exit mapping."""

from __future__ import annotations

from PyPedal.application import EXIT_STATUS, describe_exception, exit_status_for
from PyPedal.pyp_errors import (
    PyPedalDependencyError,
    PyPedalExportFormatError,
    PyPedalInternalError,
    PyPedalOptionError,
    PyPedalPedigreeFormatError,
    PyPedalUsageError,
    PyPedalValidationError,
)


def test_exit_status_mapping_matches_historical_policy():
    cases = {
        PyPedalPedigreeFormatError("x"): 2,
        PyPedalUsageError("x"): 2,
        PyPedalOptionError("x"): 3,
        PyPedalDependencyError("x"): 4,
        PyPedalValidationError("x"): 65,
        PyPedalInternalError("x"): 70,
        PyPedalExportFormatError("FORMAT", "FIELD", 3, -999, 4): 73,
    }
    for exc, expected in cases.items():
        assert exit_status_for(exc) == expected


def test_unrecognised_error_is_still_failure():
    assert exit_status_for(RuntimeError("x")) == 1
    for _klass, status in EXIT_STATUS:
        assert status != 0


def test_describe_exception_keeps_the_original():
    exc = PyPedalUsageError("bad argument")
    info = describe_exception(exc)
    assert info.title == "PyPedalUsageError"
    assert info.text == "bad argument"
    assert info.exception is exc


def test_describe_exception_does_not_hide_unexpected_errors():
    exc = RuntimeError("boom")
    info = describe_exception(exc)
    assert info.title == "RuntimeError"
    assert info.text == "boom"
    assert info.exception is exc
