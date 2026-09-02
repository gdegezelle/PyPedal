"""Group existing inbreeding coefficients by year without recomputing them."""

from __future__ import annotations

from pathlib import Path

from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import (
    PedigreeOpenOptions,
    PedigreeSession,
    group_inbreeding_by_year,
    load_into_session,
)
from PyPedal.pyp_results import InbreedingResult


def test_group_inbreeding_by_year_uses_supplied_fx(tmp_path: Path):
    source = tmp_path / "years.ped"
    source.write_text("1 0 0 1990\n2 0 0 1990\n3 1 2 2000\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(
            session,
            source,
            PedigreeOpenOptions(pedformat="asdy"),
        )
    finally:
        close_owned_pypedal_log_handlers()
    fx = {animal.animalID: 0.0 for animal in pedigree.pedigree}
    child = next(animal for animal in pedigree.pedigree if animal.originalID == 3)
    fx[child.animalID] = 0.125
    rows = group_inbreeding_by_year(pedigree, fx)
    assert [(row.year, row.mean, row.n) for row in rows] == [
        (1990, 0.0, 2),
        (2000, 0.125, 1),
    ]


def test_unknown_years_are_omitted(tmp_path: Path):
    source = tmp_path / "mixed.ped"
    source.write_text("1 0 0 .\n2 0 0 1995\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(
            session,
            source,
            PedigreeOpenOptions(pedformat="asdy"),
        )
    finally:
        close_owned_pypedal_log_handlers()
    fx = {animal.animalID: 0.05 for animal in pedigree.pedigree}
    rows = group_inbreeding_by_year(pedigree, fx)
    assert len(rows) == 1
    assert rows[0].year == 1995
    assert rows[0].n == 1


def test_cached_inbreeding_result_can_be_grouped(tmp_path: Path):
    source = tmp_path / "cache.ped"
    source.write_text("1 0 0 2001\n2 0 0 2001\n", encoding="utf-8")
    session = PedigreeSession()
    try:
        pedigree = load_into_session(
            session,
            source,
            PedigreeOpenOptions(pedformat="asdy"),
        )
    finally:
        close_owned_pypedal_log_handlers()
    fx = {animal.animalID: 0.0 for animal in pedigree.pedigree}
    session.inbreeding_result = InbreedingResult({"fx": fx, "metadata": {}})
    assert session.inbreeding_result is not None
    rows = group_inbreeding_by_year(pedigree, session.inbreeding_result.fx)
    assert len(rows) == 1
    assert rows[0].mean == 0.0
    assert rows[0].n == 2
