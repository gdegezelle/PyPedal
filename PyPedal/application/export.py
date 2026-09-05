"""Explicit user-requested export of analysis results.

Exports happen only when the caller asks. Analysis jobs never call this
module. Encoding is UTF-8. There is no Excel/OpenXML dependency.

Animal-level CSV identifies animals by source ``originalID``, stored
``name``, and current ``animalID``. Analysis APIs stay ID-based.
Coefficient serialization may clamp IEEE residue around zero; stored
result objects are not modified.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from PyPedal.application.inbreeding import InbreedingByYearRow
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult, MatingCoIGroupResult

if TYPE_CHECKING:
    from PyPedal.application.jobs import PairwiseResult
    from PyPedal.pyp_newclasses import NewPedigree

# Existing scientific comparison contract. Export-only: in-memory
# coefficients are not changed.
EXPORT_ZERO_ABS = 1e-12

INBREEDING_CSV_HEADERS = (
    "original_id",
    "name",
    "current_id",
    "f",
    "f_percent",
)
PAIR_IDENTITY_HEADERS = (
    "animal_a_original_id",
    "animal_a_name",
    "animal_a_current_id",
    "animal_b_original_id",
    "animal_b_name",
    "animal_b_current_id",
)
RELATIONSHIP_CSV_HEADERS = (*PAIR_IDENTITY_HEADERS, "relationship")
MATING_CSV_HEADERS = (*PAIR_IDENTITY_HEADERS, "f", "f_percent")
YEAR_INBREEDING_CSV_HEADERS = ("year", "n", "mean_f")


def serialize_coefficient(value: float) -> float:
    """Normalize IEEE residue around zero for user-facing export only."""
    x = float(value)
    return 0.0 if abs(x) < EXPORT_ZERO_ABS else x


def serialize_percent(value: float) -> float:
    """Numeric percentage of a serialized coefficient. No percent sign."""
    return round(100.0 * serialize_coefficient(value), 2)


def write_csv(
    destination: Path,
    headers: Sequence[str],
    rows: Iterable[Sequence[object]],
    *,
    overwrite: bool = False,
) -> Path:
    """Write ``headers`` and ``rows`` as locale-independent UTF-8 CSV."""
    path = destination.expanduser().resolve()
    if path.exists() and not overwrite:
        raise PyPedalUsageError(f"{path} already exists.")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(list(headers))
        for row in rows:
            writer.writerow(list(row))
    return path


def write_text(
    destination: Path,
    text: str,
    *,
    overwrite: bool = False,
) -> Path:
    """Write UTF-8 text. Used for scalar results."""
    path = destination.expanduser().resolve()
    if path.exists() and not overwrite:
        raise PyPedalUsageError(f"{path} already exists.")
    path.write_text(text, encoding="utf-8")
    return path


def _animals_by_current_id(pedigree: NewPedigree) -> dict[Any, Any]:
    animals = getattr(pedigree, "pedigree", None)
    if animals is None:
        raise PyPedalUsageError("A pedigree must be open before exporting.")
    return {animal.animalID: animal for animal in animals}


def _require_animal(animals_by_id: Mapping[Any, Any], current_id: object) -> Any:
    animal = animals_by_id.get(current_id)
    if animal is None:
        raise PyPedalUsageError(f"Export identity {current_id!r} is not in the loaded pedigree.")
    return animal


def _export_name(animal: object) -> str:
    raw = getattr(animal, "name", None)
    if raw is None:
        return ""
    return str(raw)


def _identity_fields(animal: Any) -> tuple[object, str, object]:
    return (animal.originalID, _export_name(animal), animal.animalID)


def _inbreeding_row(animal: object, coefficient: float) -> tuple[object, str, object, float, float]:
    original_id, name, current_id = _identity_fields(animal)
    return (
        original_id,
        name,
        current_id,
        serialize_coefficient(coefficient),
        serialize_percent(coefficient),
    )


def _pair_identity_fields(
    animals_by_id: Mapping[Any, Any],
    animal_a_id: object,
    animal_b_id: object,
) -> tuple[object, str, object, object, str, object]:
    animal_a = _require_animal(animals_by_id, animal_a_id)
    animal_b = _require_animal(animals_by_id, animal_b_id)
    return (*_identity_fields(animal_a), *_identity_fields(animal_b))


def export_inbreeding_csv(
    destination: Path,
    result: InbreedingResult,
    pedigree: NewPedigree,
    *,
    overwrite: bool = False,
) -> Path:
    """Export inbreeding with source identity, name, and raw *F*.

    ``original_id`` is the file identity. ``current_id`` is the internal
    PyPedal ID. Stored ``fx`` values are not modified.
    """
    animals_by_id = _animals_by_current_id(pedigree)
    rows = (
        _inbreeding_row(_require_animal(animals_by_id, animal_id), coefficient)
        for animal_id, coefficient in result.fx.items()
    )
    return write_csv(destination, INBREEDING_CSV_HEADERS, rows, overwrite=overwrite)


def export_year_inbreeding_csv(
    destination: Path,
    rows: Sequence[InbreedingByYearRow],
    *,
    overwrite: bool = False,
) -> Path:
    payload = ((row.year, row.n, serialize_coefficient(row.mean)) for row in rows)
    return write_csv(
        destination,
        YEAR_INBREEDING_CSV_HEADERS,
        payload,
        overwrite=overwrite,
    )


def export_mating_group_csv(
    destination: Path,
    result: MatingCoIGroupResult,
    pedigree: NewPedigree,
    *,
    overwrite: bool = False,
) -> Path:
    animals_by_id = _animals_by_current_id(pedigree)
    payload = (
        (
            *_pair_identity_fields(animals_by_id, animal_a, animal_b),
            serialize_coefficient(coefficient),
            serialize_percent(coefficient),
        )
        for (animal_a, animal_b), coefficient in result.matings.items()
    )
    return write_csv(destination, MATING_CSV_HEADERS, payload, overwrite=overwrite)


def export_mating_pair_csv(
    destination: Path,
    result: PairwiseResult,
    pedigree: NewPedigree,
    *,
    overwrite: bool = False,
) -> Path:
    animals_by_id = _animals_by_current_id(pedigree)
    row = (
        *_pair_identity_fields(animals_by_id, result.animal_a, result.animal_b),
        serialize_coefficient(result.coefficient),
        serialize_percent(result.coefficient),
    )
    return write_csv(destination, MATING_CSV_HEADERS, (row,), overwrite=overwrite)


def export_relationship_csv(
    destination: Path,
    result: PairwiseResult,
    pedigree: NewPedigree,
    *,
    overwrite: bool = False,
) -> Path:
    animals_by_id = _animals_by_current_id(pedigree)
    row = (
        *_pair_identity_fields(animals_by_id, result.animal_a, result.animal_b),
        serialize_coefficient(result.coefficient),
    )
    return write_csv(
        destination,
        RELATIONSHIP_CSV_HEADERS,
        (row,),
        overwrite=overwrite,
    )
