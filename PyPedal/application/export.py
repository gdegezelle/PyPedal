"""Explicit user-requested export of analysis results.

Exports happen only when the caller asks. Analysis jobs never call this
module. Encoding is UTF-8. There is no Excel/OpenXML dependency.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from pathlib import Path

from PyPedal.application.inbreeding import InbreedingByYearRow
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult, MatingCoIGroupResult


def write_csv(
    destination: Path,
    headers: Sequence[str],
    rows: Iterable[Sequence[object]],
    *,
    overwrite: bool = False,
) -> Path:
    """Write ``headers`` and ``rows`` as UTF-8 CSV."""
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


def export_inbreeding_csv(
    destination: Path,
    result: InbreedingResult,
    *,
    overwrite: bool = False,
) -> Path:
    """Export raw *F* coefficients. Display rounding is not applied."""
    rows = ((animal_id, coefficient) for animal_id, coefficient in result.fx.items())
    return write_csv(destination, ("animal_id", "f"), rows, overwrite=overwrite)


def export_year_inbreeding_csv(
    destination: Path,
    rows: Sequence[InbreedingByYearRow],
    *,
    overwrite: bool = False,
) -> Path:
    payload = ((row.year, row.n, row.mean) for row in rows)
    return write_csv(
        destination,
        ("year", "n", "mean_f"),
        payload,
        overwrite=overwrite,
    )


def export_mating_group_csv(
    destination: Path,
    result: MatingCoIGroupResult,
    *,
    overwrite: bool = False,
) -> Path:
    payload = ((a, b, coefficient) for (a, b), coefficient in result.matings.items())
    return write_csv(
        destination,
        ("animal_a", "animal_b", "f"),
        payload,
        overwrite=overwrite,
    )
