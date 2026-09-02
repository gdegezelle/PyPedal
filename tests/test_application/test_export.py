"""Explicit CSV/text export does not rewrite stored coefficients."""

from __future__ import annotations

from pathlib import Path

import pytest

from PyPedal.application import (
    InbreedingByYearRow,
    export_inbreeding_csv,
    export_mating_group_csv,
    export_year_inbreeding_csv,
    write_text,
)
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult, MatingCoIGroupResult


def test_export_inbreeding_csv_utf8_raw_values(tmp_path: Path) -> None:
    result = InbreedingResult({"fx": {1: 0.125, 2: 0.0}, "metadata": {}})
    dest = tmp_path / "f.csv"
    export_inbreeding_csv(dest, result)
    text = dest.read_text(encoding="utf-8")
    assert "animal_id,f" in text
    assert "0.125" in text
    with pytest.raises(PyPedalUsageError, match="already exists"):
        export_inbreeding_csv(dest, result)


def test_export_year_and_mating_and_text(tmp_path: Path) -> None:
    years = (InbreedingByYearRow(year=1990, mean=0.05, n=2),)
    year_path = tmp_path / "years.csv"
    export_year_inbreeding_csv(year_path, years)
    assert "1990,2,0.05" in year_path.read_text(encoding="utf-8")
    group = MatingCoIGroupResult(
        {
            "matings": {(1, 2): 0.125},
            "metadata": {},
        }
    )
    mating_path = tmp_path / "matings.csv"
    export_mating_group_csv(mating_path, group)
    assert "1,2,0.125" in mating_path.read_text(encoding="utf-8")
    note = tmp_path / "ne.txt"
    write_text(note, "Theoretical Ne from metadata: 4.0\n")
    assert note.read_text(encoding="utf-8").startswith("Theoretical Ne")
