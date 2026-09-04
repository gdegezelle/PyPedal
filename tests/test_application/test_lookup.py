"""Qt-free animal lookup: names are display data, not identities."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from _pedhelpers import NAMED_DUPLICATE_PED, close_owned_pypedal_log_handlers

from PyPedal.application import PedigreeOpenOptions, PedigreeSession, load_into_session
from PyPedal.application.lookup import (
    DEFAULT_RESULT_LIMIT,
    AnimalLookupIndex,
    format_animal_label,
)


def _animal(
    animal_id: int,
    original_id: int,
    name: str = "",
    sex: str = "f",
    by: int | None = 2020,
) -> SimpleNamespace:
    return SimpleNamespace(
        animalID=animal_id,
        originalID=original_id,
        name=name,
        sex=sex,
        by=by,
    )


def _index(*animals: SimpleNamespace) -> AnimalLookupIndex:
    pedigree = SimpleNamespace(pedigree=list(animals), kw={"missing_name": "Unknown_Name"})
    return AnimalLookupIndex.from_pedigree(pedigree)


def test_exact_current_id_ranks_first() -> None:
    index = _index(
        _animal(12, 99, "Zebra"),
        _animal(20, 12, "Other"),
    )
    result = index.search("12")
    assert [hit.animal_id for hit in result.hits] == [12, 20]
    assert result.hits[0].original_id == 99


def test_exact_original_id() -> None:
    index = _index(_animal(7, 555, "Pip"))
    result = index.search("555")
    assert len(result.hits) == 1
    assert result.hits[0].animal_id == 7
    assert result.hits[0].original_id == 555


def test_exact_name_case_insensitive() -> None:
    index = _index(_animal(3, 30, "Morning Bell Virgine", "f", 2022))
    lower = index.search("morning bell virginie")
    assert lower.hits == ()
    matched = index.search("  MORNING   bell  Virgine ")
    assert len(matched.hits) == 1
    assert matched.hits[0].animal_id == 3


def test_prefix_ranks_above_substring() -> None:
    index = _index(
        _animal(1, 10, "Heartbreaker Sue"),
        _animal(2, 20, "Hierners Heartbreaker"),
    )
    result = index.search("heart")
    assert [hit.animal_id for hit in result.hits] == [1, 2]


def test_substring_match() -> None:
    index = _index(_animal(8, 80, "Hierners Heartbreaker"))
    result = index.search("heartbreak")
    assert len(result.hits) == 1
    assert result.hits[0].animal_id == 8


def test_no_result() -> None:
    index = _index(_animal(1, 1, "Max"))
    result = index.search("zzzz-not-a-name")
    assert result.hits == ()
    assert result.truncated is False
    assert result.total == 0


def test_duplicate_names_return_both_and_neither_is_the_identity() -> None:
    index = _index(
        _animal(1, 101, "Bella", "f", 2019),
        _animal(2, 103, "Bella", "f", 2021),
    )
    result = index.search("bella")
    assert len(result.hits) == 2
    assert {hit.original_id for hit in result.hits} == {101, 103}
    assert {hit.animal_id for hit in result.hits} == {1, 2}
    labels = [hit.label for hit in result.hits]
    assert all("101" in label or "103" in label for label in labels)
    assert all("Bella" in label for label in labels)
    assert result.hits[0].original_id == 101
    assert result.hits[1].original_id == 103


def test_unnamed_animal_is_found_by_ids_not_name() -> None:
    index = _index(
        _animal(5, 50, "50", "m", 2018),
        _animal(6, 60, "Unknown_Name", "m", 2019),
        _animal(7, 70, "", "f", 2020),
    )
    assert index.named_count == 0
    assert index.search("50").hits[0].animal_id == 5
    assert index.search("5").hits[0].animal_id == 5
    assert index.search("unknown_name").hits == ()
    assert index.search("50").hits[0].name == ""
    unnamed = index.hit_for_animal_id(7)
    assert unnamed is not None
    assert unnamed.label.startswith("70")
    assert "None" not in unnamed.label
    assert "Unknown" not in unnamed.label


def test_stable_ranking_within_equal_rank() -> None:
    index = _index(
        _animal(9, 30, "Prefix Z"),
        _animal(4, 10, "Prefix A"),
        _animal(8, 20, "Prefix M"),
    )
    result = index.search("prefix")
    assert [hit.original_id for hit in result.hits] == [10, 20, 30]
    assert [hit.animal_id for hit in result.hits] == [4, 8, 9]


def test_result_limit_reports_truncated() -> None:
    animals = [_animal(i + 1, 1000 + i, f"Shared {i:03d}") for i in range(60)]
    index = _index(*animals)
    result = index.search("shared", limit=DEFAULT_RESULT_LIMIT)
    assert len(result.hits) == DEFAULT_RESULT_LIMIT
    assert result.truncated is True
    assert result.total == 60
    assert result.hits[0].original_id == 1000
    assert result.hits[-1].original_id == 1000 + DEFAULT_RESULT_LIMIT - 1


def test_empty_query_returns_nothing() -> None:
    index = _index(_animal(1, 1, "Max"))
    assert index.search("").hits == ()
    assert index.search("   ").hits == ()


def test_label_omits_placeholders() -> None:
    label = format_animal_label(
        name="",
        original_id=12345,
        sex="unknown",
        birth_year=None,
        animal_id=81231,
    )
    assert label == "12345 — ID 81231"
    named = format_animal_label(
        name="Bella",
        original_id=12345,
        sex="f",
        birth_year=2021,
        animal_id=81231,
    )
    assert named == "Bella — 12345 — ♀ — 2021 — ID 81231"


def test_session_replacement_invalidates_lookup(tmp_path: Path) -> None:
    first = tmp_path / "named.ped"
    first.write_text(NAMED_DUPLICATE_PED, encoding="utf-8")
    second = tmp_path / "plain.ped"
    second.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        load_into_session(
            session,
            first,
            PedigreeOpenOptions(pedformat="asdxbn", separator=",").normalized(),
        )
        assert session.animal_lookup is not None
        bellas = session.animal_lookup.search("bella")
        assert len(bellas.hits) == 2
        load_into_session(
            session,
            second,
            PedigreeOpenOptions(pedformat="asd", separator=" ").normalized(),
        )
        assert session.animal_lookup is not None
        assert session.animal_lookup.search("bella").hits == ()
        current = session.animal_lookup.search("1")
        assert len(current.hits) == 1
        assert current.hits[0].original_id == 1
        session.clear()
        assert session.animal_lookup is None
    finally:
        close_owned_pypedal_log_handlers()


def test_ordinary_load_builds_lookup_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    builds: list[int] = []
    original = AnimalLookupIndex.from_pedigree.__func__

    def counted(cls: type[AnimalLookupIndex], pedigree: object) -> AnimalLookupIndex:
        builds.append(1)
        return original(cls, pedigree)

    monkeypatch.setattr(AnimalLookupIndex, "from_pedigree", classmethod(counted))
    source = tmp_path / "demo.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        load_into_session(session, source, PedigreeOpenOptions(separator=" "))
        assert len(builds) == 1
        assert session.animal_lookup is not None
        session.rebuild_animal_lookup()
        assert len(builds) == 2
    finally:
        close_owned_pypedal_log_handlers()
