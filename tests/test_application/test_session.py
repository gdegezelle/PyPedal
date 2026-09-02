"""PedigreeSession ownership, replacement, failure retention, and cache."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import PedigreeOpenOptions, PedigreeSession, load_into_session
from PyPedal.pyp_errors import PyPedalError
from PyPedal.pyp_results import InbreedingResult


def test_session_starts_empty():
    session = PedigreeSession()
    assert session.is_empty
    assert session.state == "empty"
    assert session.pedigree is None
    assert session.source_path is None
    assert session.load_options is None
    assert session.inbreeding_result is None


def test_successful_load_installs_pedigree_and_path(tmp_path: Path):
    session = PedigreeSession()
    source = tmp_path / "demo.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    options = PedigreeOpenOptions(pedformat="asd", separator=" ", renumber=True)
    try:
        pedigree = load_into_session(session, source, options)
    finally:
        close_owned_pypedal_log_handlers()
    assert session.pedigree is pedigree
    assert session.source_path == source.resolve()
    assert session.load_options is not None
    assert session.load_options.pedformat == "asd"
    assert session.load_options.separator == " "
    assert session.state == "loaded"
    assert len(pedigree.pedigree) == 3


def test_failed_load_keeps_previous_pedigree(tmp_path: Path):
    session = PedigreeSession()
    good = tmp_path / "good.ped"
    good.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    bad = tmp_path / "bad.ped"
    bad.write_text("1 0 0\n2 0 0\n3 1\n", encoding="utf-8")
    try:
        first = load_into_session(session, good, PedigreeOpenOptions())
        session.inbreeding_result = InbreedingResult({"fx": {1: 0.0}, "metadata": {}})
        cached = session.inbreeding_result
        with pytest.raises(PyPedalError):
            load_into_session(session, bad, PedigreeOpenOptions())
    finally:
        close_owned_pypedal_log_handlers()
    assert session.pedigree is first
    assert session.source_path == good.resolve()
    assert session.inbreeding_result is cached


def test_successful_replacement_clears_inbreeding_cache(tmp_path: Path):
    session = PedigreeSession()
    first_path = tmp_path / "a.ped"
    first_path.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    second_path = tmp_path / "b.ped"
    second_path.write_text("10 0 0\n20 0 0\n", encoding="utf-8")
    try:
        load_into_session(session, first_path, PedigreeOpenOptions())
        session.inbreeding_result = InbreedingResult({"fx": {1: 0.125}, "metadata": {}})
        second = load_into_session(session, second_path, PedigreeOpenOptions())
    finally:
        close_owned_pypedal_log_handlers()
    assert session.pedigree is second
    assert session.inbreeding_result is None
    assert session.source_path == second_path.resolve()


def test_clear_returns_to_empty(tmp_path: Path):
    session = PedigreeSession()
    source = tmp_path / "demo.ped"
    source.write_text("1 0 0\n2 0 0\n", encoding="utf-8")
    try:
        load_into_session(session, source, PedigreeOpenOptions())
        session.inbreeding_result = InbreedingResult({"fx": {}, "metadata": {}})
    finally:
        close_owned_pypedal_log_handlers()
    session.clear()
    assert session.is_empty
    assert session.inbreeding_result is None
    assert session.source_path is None
    assert session.load_options is None
