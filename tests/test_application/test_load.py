"""Open-option normalisation, pathlib paths, progress, and error propagation."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import (
    PedigreeOpenOptions,
    PedigreeSession,
    load_into_session,
    normalize_sepchar,
    resolve_source_path,
)
from PyPedal.pyp_errors import PyPedalError, PyPedalPedigreeFormatError


def test_normalize_sepchar_empty_and_spaces_mean_a_space():
    assert normalize_sepchar("") == " "
    assert normalize_sepchar(" ") == " "
    assert normalize_sepchar("  ") == " "
    assert normalize_sepchar(None) == " "


def test_normalize_sepchar_comma_with_leftover_spaces_is_a_comma():
    assert normalize_sepchar(",") == ","
    assert normalize_sepchar(", ") == ","
    assert normalize_sepchar(" ,") == ","
    assert normalize_sepchar(" , ") == ","


def test_normalize_sepchar_keeps_tab():
    assert normalize_sepchar("\t") == "\t"


def test_open_options_normalize_pedformat_and_separator():
    raw = PedigreeOpenOptions(pedformat="  asdxbn  ", separator=", ", renumber=True)
    options = raw.normalized()
    assert options.pedformat == "asdxbn"
    assert options.separator == ","
    assert options.renumber is True
    assert options.messages == "quiet"
    assert options.pedigree_summary == 0


def test_open_options_empty_pedformat_falls_back_to_asd():
    assert PedigreeOpenOptions(pedformat="").normalized().pedformat == "asd"


def test_application_load_uses_quiet_and_no_summary(tmp_path: Path):
    options = PedigreeOpenOptions()
    assert options.messages == "quiet"
    assert options.pedigree_summary == 0
    source = tmp_path / "dogs.ped"
    payload = options.to_library_options(source)
    assert payload["messages"] == "quiet"
    assert payload["pedigree_summary"] == 0
    assert payload["pedname"] == source.name
    assert Path(payload["pedfile"]) == source


def test_resolve_source_path_is_absolute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path)
    relative = Path("nested.ped")
    relative.write_text("1 0 0\n", encoding="utf-8")
    resolved = resolve_source_path(relative)
    assert resolved.is_absolute()
    assert resolved == (tmp_path / "nested.ped").resolve()


def test_load_accepts_pathlib_path(tmp_path: Path):
    source = tmp_path / "space.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(session, source, PedigreeOpenOptions(separator=" "))
    finally:
        close_owned_pypedal_log_handlers()
    assert len(pedigree.pedigree) == 3


def test_progress_is_forwarded(tmp_path: Path):
    source = tmp_path / "progress.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    calls: list[tuple[int, int | None]] = []

    def progress(done: int, total: int | None) -> None:
        calls.append((done, total))

    session = PedigreeSession()
    try:
        load_into_session(session, source, PedigreeOpenOptions(), progress=progress)
    finally:
        close_owned_pypedal_log_handlers()
    assert calls


def test_pypedal_errors_propagate_typed(tmp_path: Path):
    source = tmp_path / "bad.ped"
    source.write_text("1 0 0\n2 0 0\n3 1\n", encoding="utf-8")
    session = PedigreeSession()
    with pytest.raises(PyPedalPedigreeFormatError) as raised:
        try:
            load_into_session(session, source, PedigreeOpenOptions())
        finally:
            close_owned_pypedal_log_handlers()
    assert isinstance(raised.value, PyPedalError)
    assert session.is_empty


def test_comma_csv_loads_with_normalized_separator(tmp_path: Path):
    source = tmp_path / "csv.ped"
    source.write_text(
        "1,30497,52843,f,03132018,A Day Before Sunrise de Mar&Mar\n"
        "2,12401,68419,f,01111995,A Galaxie Mii Jimajo Nubegin\n",
        encoding="utf-8",
    )
    session = PedigreeSession()
    options = PedigreeOpenOptions(pedformat="asdxbn", separator=", ", renumber=True)
    try:
        pedigree = load_into_session(session, source, options)
    finally:
        close_owned_pypedal_log_handlers()
    named = [animal for animal in pedigree.pedigree if animal.originalID == 1]
    assert len(named) == 1
    assert named[0].name == "A Day Before Sunrise de Mar&Mar"
    assert session.load_options is not None
    assert session.load_options.separator == ","


def test_quiet_application_load_omits_per_record_debug(tmp_path: Path) -> None:
    source = tmp_path / "quiet.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(session, source, PedigreeOpenOptions(separator=" "))
        logfile = Path(pedigree.kw["logfile"])
        text = logfile.read_text(encoding="utf-8")
        assert "Raw line read" not in text
        assert "Processing line:" not in text
        assert len(text) < 50_000
    finally:
        close_owned_pypedal_log_handlers()
