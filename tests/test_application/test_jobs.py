"""Application analysis job adapters: cache, output=False, mutation, progress."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import (
    PedigreeOpenOptions,
    PedigreeSession,
    ensure_inbreeding,
    load_into_session,
    parse_animal_id,
    require_pedigree,
    run_effective_founders,
    run_inbreeding,
    run_inbreeding_by_year,
    run_mating_coi,
    run_mating_coi_group,
    run_relationship,
    run_theoretical_ne,
    save_pedigree,
)
from PyPedal.pyp_errors import PyPedalError, PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult

MRODE = """\
# Pedigree from Mrode (2005) Table 2.1
1 0 0
2 0 0
3 1 2
4 1 0
5 4 3
6 5 2
"""


def _load(tmp_path: Path, text: str = "1 0 0\n2 0 0\n3 1 2\n") -> PedigreeSession:
    source = tmp_path / "demo.ped"
    source.write_text(text, encoding="utf-8")
    session = PedigreeSession()
    load_into_session(session, source, PedigreeOpenOptions(separator=" "))
    return session


def test_require_pedigree_without_load_raises() -> None:
    with pytest.raises(PyPedalUsageError, match="must be open"):
        require_pedigree(PedigreeSession())


def test_parse_animal_id_rejects_blank_and_non_integer() -> None:
    with pytest.raises(PyPedalUsageError, match="required"):
        parse_animal_id("  ")
    with pytest.raises(PyPedalUsageError, match="integer"):
        parse_animal_id("1.5")
    with pytest.raises(PyPedalUsageError, match="integer"):
        parse_animal_id("abc")
    assert parse_animal_id(" 5 ") == 5


def test_inbreeding_requires_pedigree() -> None:
    with pytest.raises(PyPedalUsageError):
        run_inbreeding(PedigreeSession())


def test_inbreeding_mrode_caches_mutates_and_writes_no_dat(tmp_path: Path) -> None:
    session = _load(tmp_path, MRODE)
    progress: list[tuple[int, int | None]] = []

    def report(done: int, total: int | None) -> None:
        progress.append((done, total))

    try:
        result = run_inbreeding(session, progress=report)
        animal_five = next(a for a in session.pedigree.pedigree if a.originalID == 5)
        assert animal_five.fa == 0.125
        assert result.fx[animal_five.animalID] == 0.125
        assert session.inbreeding_result is result
        assert isinstance(result, InbreedingResult)
        dat_files = list(tmp_path.glob("*.dat"))
        assert dat_files == []
        assert progress  # Meuwissen-Luo reports progress
    finally:
        close_owned_pypedal_log_handlers()


def test_year_analysis_reuses_cached_inbreeding(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "years.ped"
    source.write_text("1 0 0 1990\n2 0 0 1990\n3 1 2 2000\n", encoding="utf-8")
    session = PedigreeSession()
    import PyPedal.application.jobs as jobs

    calls: list[int] = []
    original = jobs.inbreeding

    def wrapped(pedigree, **kwargs):
        calls.append(1)
        return original(pedigree, **kwargs)

    monkeypatch.setattr(jobs, "inbreeding", wrapped)
    try:
        load_into_session(session, source, PedigreeOpenOptions(pedformat="asdy"))
        first = run_inbreeding_by_year(session)
        assert first.computed_inbreeding is True
        second = run_inbreeding_by_year(session)
        assert second.computed_inbreeding is False
        assert second.rows == first.rows
        assert len(calls) == 1
        assert any(row.year == 1990 for row in first.rows)
    finally:
        close_owned_pypedal_log_handlers()


def test_ensure_inbreeding_returns_cached_object(tmp_path: Path) -> None:
    session = _load(tmp_path)
    try:
        first, computed = ensure_inbreeding(session)
        assert computed is True
        second, computed_again = ensure_inbreeding(session)
        assert computed_again is False
        assert second is first
    finally:
        close_owned_pypedal_log_handlers()


def test_lacy_founders_output_false_and_cache(tmp_path: Path) -> None:
    session = _load(tmp_path)
    try:
        outcome = run_effective_founders(session)
        assert outcome.implicit_renumber is False
        assert outcome.result.fa_effective_founders > 0
        assert session.effective_founders_result is outcome.result
        assert list(tmp_path.glob("*.dat")) == []
    finally:
        close_owned_pypedal_log_handlers()


def test_lacy_implicit_renumber_clears_inbreeding_cache(tmp_path: Path) -> None:
    source = tmp_path / "unnumbered.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        load_into_session(
            session,
            source,
            PedigreeOpenOptions(separator=" ", renumber=False),
        )
        session.inbreeding_result = InbreedingResult({"fx": {1: 0.0}, "metadata": {}})
        outcome = run_effective_founders(session)
        assert outcome.implicit_renumber is True
        assert session.inbreeding_result is None
        assert session.effective_founders_result is outcome.result
    finally:
        close_owned_pypedal_log_handlers()


def test_relationship_and_mating_on_mrode(tmp_path: Path) -> None:
    session = _load(tmp_path, MRODE)
    try:
        related = run_relationship(session, 4, 3)
        assert related.coefficient == 0.25
        unrelated = run_relationship(session, 1, 2)
        assert unrelated.coefficient == 0.0
        mating = run_mating_coi(session, 4, 3)
        assert mating.coefficient == 0.125
        group = run_mating_coi_group(session, [(4, 3), (1, 2)])
        assert group.matings[(4, 3)] == 0.125
        assert group.matings[(1, 2)] == 0.0
        with pytest.raises(PyPedalUsageError):
            run_mating_coi_group(session, [])
        with pytest.raises(PyPedalError):
            run_relationship(session, 4, 99)
    finally:
        close_owned_pypedal_log_handlers()


def test_theoretical_ne_from_metadata(tmp_path: Path) -> None:
    session = _load(tmp_path, MRODE)
    try:
        value = run_theoretical_ne(session)
        assert value > 0
        assert session.theoretical_ne == value
        assert list(tmp_path.glob("*.dat")) == []
    finally:
        close_owned_pypedal_log_handlers()


def test_save_pedigree_writes_and_refuses_overwrite(tmp_path: Path) -> None:
    session = _load(tmp_path)
    dest = tmp_path / "out.ped"
    try:
        saved = save_pedigree(session, dest)
        assert saved == dest.resolve()
        assert dest.is_file()
        with pytest.raises(PyPedalUsageError, match="already exists"):
            save_pedigree(session, dest)
        save_pedigree(session, dest, overwrite=True)
    finally:
        close_owned_pypedal_log_handlers()


def test_typed_inbreeding_failure_propagates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    session = _load(tmp_path)
    import PyPedal.application.jobs as jobs

    def boom(*_args, **_kwargs):
        raise PyPedalUsageError("forced")

    monkeypatch.setattr(jobs, "inbreeding", boom)
    try:
        with pytest.raises(PyPedalUsageError, match="forced"):
            run_inbreeding(session)
    finally:
        close_owned_pypedal_log_handlers()
