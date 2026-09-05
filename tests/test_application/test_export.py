"""Analysis CSV exports preserve source identity and do not rewrite stored coefficients."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import (
    EXPORT_ZERO_ABS,
    InbreedingByYearRow,
    PedigreeOpenOptions,
    PedigreeSession,
    export_inbreeding_csv,
    export_mating_group_csv,
    export_mating_pair_csv,
    export_relationship_csv,
    export_year_inbreeding_csv,
    load_into_session,
    serialize_coefficient,
    serialize_percent,
    write_text,
)
from PyPedal.application.jobs import PairwiseResult
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult, MatingCoIGroupResult

NEAR_ZERO = 6.661338147750939e-16
NEAR_ZERO_NEGATIVE = -6.661338147750939e-16
SMALL_NONZERO = 7.450580596923828e-07
CANONICAL_MATING_F = 0.10095650884805218

NAMED_RENUMBER_PED = """\
53,0,0,f,01012019,A'Vigdors Berta Beautiful
64,0,0,f,01012020,Heart Breaker
75,53,64,m,01012021,O'Malley Junior
"""

COLETTE_PED = """\
20196,0,0,f,01012019,Colette
20209,0,0,f,01012020,Colette
30000,20196,20209,m,01012021,Offspring
"""

UNNAMED_PED = """\
10 0 0
20 0 0
30 10 20
"""


def _animal(current_id: int, original_id: int, name: str | None = "") -> SimpleNamespace:
    return SimpleNamespace(animalID=current_id, originalID=original_id, name=name)


def _pedigree(*animals: SimpleNamespace) -> SimpleNamespace:
    return SimpleNamespace(pedigree=list(animals))


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _load(tmp_path: Path, text: str, *, pedformat: str, separator: str) -> PedigreeSession:
    source = tmp_path / "demo.ped"
    source.write_text(text, encoding="utf-8")
    session = PedigreeSession()
    load_into_session(
        session,
        source,
        PedigreeOpenOptions(pedformat=pedformat, separator=separator, renumber=True),
    )
    return session


def test_serialize_coefficient_clamps_residue_only() -> None:
    assert serialize_coefficient(NEAR_ZERO) == 0.0
    assert serialize_coefficient(NEAR_ZERO_NEGATIVE) == 0.0
    assert serialize_coefficient(0.0) == 0.0
    assert serialize_coefficient(SMALL_NONZERO) == SMALL_NONZERO
    assert serialize_coefficient(0.125) == 0.125
    assert serialize_percent(NEAR_ZERO) == 0.0
    assert serialize_percent(0.125) == 12.5
    assert serialize_percent(CANONICAL_MATING_F) == 10.1
    assert EXPORT_ZERO_ABS == 1e-12


def test_export_inbreeding_csv_named_renumbered_identity(tmp_path: Path) -> None:
    try:
        session = _load(tmp_path, NAMED_RENUMBER_PED, pedformat="asdxbn", separator=",")
        pedigree = session.pedigree
        assert pedigree is not None
        by_original = {animal.originalID: animal for animal in pedigree.pedigree}
        heart = by_original[64]
        assert heart.animalID != heart.originalID
        fx = {
            animal.animalID: 0.125 if animal.originalID == 64 else 0.0
            for animal in pedigree.pedigree
        }
        result = InbreedingResult({"fx": fx, "metadata": {}})
        dest = tmp_path / "f.csv"
        export_inbreeding_csv(dest, result, pedigree)
        rows = {int(row["original_id"]): row for row in _read_rows(dest)}
        assert list(_read_rows(dest)[0].keys()) == [
            "original_id",
            "name",
            "current_id",
            "f",
            "f_percent",
        ]
        assert "animal_id" not in dest.read_text(encoding="utf-8")
        assert rows[53]["name"] == "A'Vigdors Berta Beautiful"
        assert rows[64]["name"] == "Heart Breaker"
        assert rows[75]["name"] == "O'Malley Junior"
        assert int(rows[64]["current_id"]) == heart.animalID
        assert int(rows[64]["original_id"]) != int(rows[64]["current_id"])
        assert rows[64]["f"] == "0.125"
        assert rows[64]["f_percent"] == "12.5"
        assert rows[64]["f"] != "125"
        with dest.open(encoding="utf-8", newline="") as handle:
            reread = next(row for row in csv.DictReader(handle) if row["original_id"] == "64")
        assert reread["name"] == "Heart Breaker"
        with pytest.raises(PyPedalUsageError, match="already exists"):
            export_inbreeding_csv(dest, result, pedigree)
    finally:
        close_owned_pypedal_log_handlers()


def test_export_inbreeding_csv_duplicate_colette(tmp_path: Path) -> None:
    try:
        session = _load(tmp_path, COLETTE_PED, pedformat="asdxbn", separator=",")
        pedigree = session.pedigree
        assert pedigree is not None
        colettes = [animal for animal in pedigree.pedigree if animal.name == "Colette"]
        assert len(colettes) == 2
        originals = {animal.originalID for animal in colettes}
        currents = {animal.animalID for animal in colettes}
        assert originals == {20196, 20209}
        assert len(currents) == 2
        result = InbreedingResult(
            {"fx": {animal.animalID: 0.0 for animal in pedigree.pedigree}, "metadata": {}}
        )
        dest = tmp_path / "colette.csv"
        export_inbreeding_csv(dest, result, pedigree)
        rows = [row for row in _read_rows(dest) if row["name"] == "Colette"]
        assert len(rows) == 2
        assert {int(row["original_id"]) for row in rows} == {20196, 20209}
        assert len({row["current_id"] for row in rows}) == 2
    finally:
        close_owned_pypedal_log_handlers()


def test_export_inbreeding_csv_unnamed_keeps_name_column(tmp_path: Path) -> None:
    try:
        session = _load(tmp_path, UNNAMED_PED, pedformat="asd", separator=" ")
        pedigree = session.pedigree
        assert pedigree is not None
        animal = next(item for item in pedigree.pedigree if item.originalID == 10)
        assert animal.animalID != animal.originalID
        result = InbreedingResult(
            {"fx": {item.animalID: 0.0 for item in pedigree.pedigree}, "metadata": {}}
        )
        dest = tmp_path / "unnamed.csv"
        export_inbreeding_csv(dest, result, pedigree)
        rows = {int(row["original_id"]): row for row in _read_rows(dest)}
        assert "name" in rows[10]
        assert rows[10]["name"] == str(animal.name)
        assert int(rows[10]["current_id"]) == animal.animalID
        assert int(rows[10]["original_id"]) != int(rows[10]["current_id"])
    finally:
        close_owned_pypedal_log_handlers()


def test_export_inbreeding_csv_near_zero_and_small_nonzero(tmp_path: Path) -> None:
    pedigree = _pedigree(
        _animal(1, 101, "Residue"),
        _animal(2, 102, "Negative"),
        _animal(3, 103, "Tiny"),
    )
    fx = {1: NEAR_ZERO, 2: NEAR_ZERO_NEGATIVE, 3: SMALL_NONZERO}
    result = InbreedingResult({"fx": fx, "metadata": {}})
    dest = tmp_path / "residue.csv"
    export_inbreeding_csv(dest, result, pedigree)
    rows = {int(row["current_id"]): row for row in _read_rows(dest)}
    assert rows[1]["f"] == "0.0"
    assert rows[1]["f_percent"] == "0.0"
    assert rows[2]["f"] == "0.0"
    assert rows[3]["f"] != "0.0"
    assert float(rows[3]["f"]) == SMALL_NONZERO
    assert result.fx[1] == NEAR_ZERO
    assert result.fx[2] == NEAR_ZERO_NEGATIVE
    assert result.fx[3] == SMALL_NONZERO


def test_export_inbreeding_csv_missing_name_is_empty(tmp_path: Path) -> None:
    pedigree = _pedigree(_animal(1, 9, None))
    result = InbreedingResult({"fx": {1: 0.0}, "metadata": {}})
    dest = tmp_path / "missing.csv"
    export_inbreeding_csv(dest, result, pedigree)
    assert _read_rows(dest)[0]["name"] == ""


def test_export_inbreeding_csv_unknown_id_raises(tmp_path: Path) -> None:
    pedigree = _pedigree(_animal(1, 9, "Only"))
    result = InbreedingResult({"fx": {99: 0.125}, "metadata": {}})
    with pytest.raises(PyPedalUsageError, match="not in the loaded pedigree"):
        export_inbreeding_csv(tmp_path / "missing-id.csv", result, pedigree)


def test_export_csv_dialect_is_locale_independent(tmp_path: Path) -> None:
    pedigree = _pedigree(_animal(1, 53, "O'Malley, Junior"))
    result = InbreedingResult({"fx": {1: 0.125}, "metadata": {}})
    dest = tmp_path / "dialect.csv"
    export_inbreeding_csv(dest, result, pedigree)
    raw = dest.read_bytes()
    assert raw.startswith(b"original_id,name,current_id,f,f_percent")
    assert b"\r\n" in raw[:80]
    assert not raw.startswith(b"\xef\xbb\xbf")
    text = dest.read_text(encoding="utf-8")
    assert "0.125" in text
    assert ",125\n" not in text.replace("\r\n", "\n")
    assert "%" not in text
    rows = _read_rows(dest)
    assert rows[0]["name"] == "O'Malley, Junior"
    assert rows[0]["f"] == "0.125"
    assert rows[0]["f_percent"] == "12.5"


def test_export_relationship_csv_identity_and_near_zero(tmp_path: Path) -> None:
    pedigree = _pedigree(
        _animal(11, 20196, "Colette"),
        _animal(22, 20209, "Colette"),
    )
    pair = PairwiseResult(animal_a=11, animal_b=22, coefficient=NEAR_ZERO)
    dest = tmp_path / "relationship.csv"
    export_relationship_csv(dest, pair, pedigree)
    row = _read_rows(dest)[0]
    assert list(row.keys()) == [
        "animal_a_original_id",
        "animal_a_name",
        "animal_a_current_id",
        "animal_b_original_id",
        "animal_b_name",
        "animal_b_current_id",
        "relationship",
    ]
    assert row["animal_a_original_id"] == "20196"
    assert row["animal_b_original_id"] == "20209"
    assert row["animal_a_name"] == "Colette"
    assert row["animal_b_name"] == "Colette"
    assert row["animal_a_current_id"] == "11"
    assert row["animal_b_current_id"] == "22"
    assert row["relationship"] == "0.0"
    assert pair.coefficient == NEAR_ZERO
    raw_ok = PairwiseResult(animal_a=11, animal_b=22, coefficient=0.25)
    dest2 = tmp_path / "relationship-raw.csv"
    export_relationship_csv(dest2, raw_ok, pedigree)
    assert _read_rows(dest2)[0]["relationship"] == "0.25"


def test_export_mating_pair_and_group_csv(tmp_path: Path) -> None:
    pedigree = _pedigree(
        _animal(98001, 98685, "Hierners Heartbreaker"),
        _animal(97984, 98667, "Hierners Honeybear"),
    )
    pair = PairwiseResult(
        animal_a=98001,
        animal_b=97984,
        coefficient=CANONICAL_MATING_F,
    )
    dest = tmp_path / "mating.csv"
    export_mating_pair_csv(dest, pair, pedigree)
    row = _read_rows(dest)[0]
    assert row["animal_a_original_id"] == "98685"
    assert row["animal_a_name"] == "Hierners Heartbreaker"
    assert row["animal_a_current_id"] == "98001"
    assert row["animal_b_original_id"] == "98667"
    assert row["animal_b_name"] == "Hierners Honeybear"
    assert row["animal_b_current_id"] == "97984"
    assert row["f"] == str(CANONICAL_MATING_F)
    assert row["f_percent"] == "10.1"
    assert "%" not in dest.read_text(encoding="utf-8")
    assert pair.coefficient == CANONICAL_MATING_F

    group = MatingCoIGroupResult(
        {"matings": {(98001, 97984): 0.125, (97984, 98001): NEAR_ZERO}, "metadata": {}}
    )
    group_path = tmp_path / "mating_group.csv"
    export_mating_group_csv(group_path, group, pedigree)
    group_rows = _read_rows(group_path)
    assert list(group_rows[0].keys())[-2:] == ["f", "f_percent"]
    by_a = {row["animal_a_original_id"]: row for row in group_rows}
    assert by_a["98685"]["f"] == "0.125"
    assert by_a["98685"]["f_percent"] == "12.5"
    assert by_a["98667"]["f"] == "0.0"
    assert group.matings[(97984, 98001)] == NEAR_ZERO


def test_export_year_and_text_aggregates(tmp_path: Path) -> None:
    years = (
        InbreedingByYearRow(year=1990, mean=0.05, n=2),
        InbreedingByYearRow(year=1991, mean=NEAR_ZERO, n=1),
    )
    year_path = tmp_path / "years.csv"
    export_year_inbreeding_csv(year_path, years)
    rows = _read_rows(year_path)
    assert list(rows[0].keys()) == ["year", "n", "mean_f"]
    assert "original_id" not in rows[0]
    assert rows[0]["year"] == "1990"
    assert rows[0]["n"] == "2"
    assert rows[0]["mean_f"] == "0.05"
    assert rows[1]["mean_f"] == "0.0"
    note = tmp_path / "ne.txt"
    write_text(note, "Theoretical Ne from metadata: 4.0\n")
    assert note.read_text(encoding="utf-8").startswith("Theoretical Ne")
