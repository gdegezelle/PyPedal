"""Project existing inbreeding coefficients by recorded birth year.

This module does not compute inbreeding. It groups an already-computed
``fx`` mapping (typically ``InbreedingResult.fx``) using each animal's
stored ``by`` field. Meuwissen-Luo and every other recurrence stay in
``pyp_nrm``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from PyPedal.pyp_results import AnimalId

if TYPE_CHECKING:
    from PyPedal.pyp_newclasses import NewPedigree

# Tokens the CustomTkinter GUI skipped when grouping by year. Unknown
# chronology is ``None`` in PyPedal 4.x; ``0`` and ``-999`` are leftover
# missing-year sentinels. Grouping must not treat those as birth years.
UNKNOWN_YEAR_TOKENS: tuple[object, ...] = (None, "", 0, -999)


@dataclass(frozen=True, slots=True)
class InbreedingByYearRow:
    """Mean of existing *F* values for one recorded birth year."""

    year: int
    mean: float
    n: int


def group_inbreeding_by_year(
    pedigree: NewPedigree,
    fx: Mapping[AnimalId, float],
) -> tuple[InbreedingByYearRow, ...]:
    """Group ``fx`` by ``animal.by`` without recomputing coefficients.

    Animals whose year is an unknown-year token, or whose current
    ``animalID`` is absent from ``fx``, are omitted. The mean is the
    arithmetic average of the supplied coefficients.
    """
    by_year: dict[int, list[float]] = {}
    for animal in pedigree.pedigree:
        year = getattr(animal, "by", None)
        if year in UNKNOWN_YEAR_TOKENS:
            continue
        if not isinstance(year, int):
            continue
        coef = fx.get(animal.animalID)
        if coef is None:
            continue
        by_year.setdefault(year, []).append(float(coef))
    rows = [
        InbreedingByYearRow(year=year, mean=sum(values) / len(values), n=len(values))
        for year, values in by_year.items()
    ]
    rows.sort(key=lambda row: row.year)
    return tuple(rows)
