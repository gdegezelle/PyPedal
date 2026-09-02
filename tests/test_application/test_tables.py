"""PedigreeTableSource row/column access without copying animals."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import (
    BROWSE_COLUMNS,
    PedigreeOpenOptions,
    PedigreeSession,
    PedigreeTableSource,
    load_into_session,
)


def _load(tmp_path: Path, text: str, pedformat: str = "asd") -> PedigreeTableSource:
    source = tmp_path / "table.ped"
    source.write_text(text, encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(
            session,
            source,
            PedigreeOpenOptions(pedformat=pedformat, separator=" ", renumber=True),
        )
    finally:
        close_owned_pypedal_log_handlers()
    return PedigreeTableSource(pedigree)


def test_table_source_shape_and_raw_values(tmp_path: Path):
    source = _load(tmp_path, "1 0 0\n2 0 0\n3 1 2\n")
    assert source.row_count() == 3
    assert source.column_count() == len(BROWSE_COLUMNS)
    keys = [column.key for column in source.columns()]
    assert keys == [
        "originalID",
        "animalID",
        "sire",
        "dam",
        "year",
        "sex",
        "name",
        "fa",
    ]
    # Renumber assigns 1-based animalID; original IDs stay as recorded.
    assert source.value(0, 0) == 1
    assert source.value(0, 1) == 1
    assert source.value(2, 1) == 3
    assert source.value(2, 2) == 1
    assert source.value(2, 3) == 2


def test_table_source_does_not_copy_the_animal_list(tmp_path: Path):
    session = PedigreeSession()
    path = tmp_path / "live.ped"
    path.write_text("1 0 0\n2 0 0\n", encoding="utf-8")
    try:
        pedigree = load_into_session(session, path, PedigreeOpenOptions())
    finally:
        close_owned_pypedal_log_handlers()
    table = PedigreeTableSource(pedigree)
    assert table._pedigree is pedigree
    assert table._pedigree.pedigree is pedigree.pedigree
    assert table.row_count() == 2


def test_missing_parent_and_year_are_raw_domain_values(tmp_path: Path):
    source = _load(tmp_path, "1 0 0\n")
    sire = source.value(0, 2)
    dam = source.value(0, 3)
    year = source.value(0, 4)
    sex = source.value(0, 5)
    fa = source.value(0, 7)
    assert sire == 0
    assert dam == 0
    assert year is None
    assert sex == "u"
    assert fa == 0.0
    assert "unknown" not in {sire, dam, year, fa}


def test_string_identities_appear_in_name(tmp_path: Path):
    source = _load(tmp_path, "ALPHA 0 0\nBETA 0 0\nGAMMA ALPHA BETA\n", pedformat="ASD")
    names = [source.value(row, 6) for row in range(source.row_count())]
    assert "ALPHA" in names
    assert "GAMMA" in names
    animal_id = source.value(0, 1)
    original_id = source.value(0, 0)
    assert isinstance(animal_id, int)
    assert original_id is not None


def test_out_of_range_access_raises(tmp_path: Path):
    source = _load(tmp_path, "1 0 0\n")
    with pytest.raises(IndexError):
        source.value(1, 0)
    with pytest.raises(IndexError):
        source.value(0, 99)
    with pytest.raises(IndexError):
        source.column(-1)
