"""QAbstractTableModel over PedigreeTableSource — no per-row item copies."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt

from PyPedal.application import PedigreeTableSource
from PyPedal.desktop.models.pedigree_table import (
    MISSING_DISPLAY,
    PedigreeFilterProxy,
    PedigreeTableModel,
    format_display_value,
)


class _LazyAnimals:
    def __init__(self, count: int) -> None:
        self.count = count
        self.gets = 0

    def __len__(self) -> int:
        return self.count

    def __getitem__(self, index: int) -> SimpleNamespace:
        if index < 0 or index >= self.count:
            raise IndexError(index)
        self.gets += 1
        return SimpleNamespace(
            originalID=index + 1,
            animalID=index + 1,
            sireID=0,
            damID=0,
            by=None,
            sex="u",
            name=f"a{index + 1}",
            fa=0.125 if index == 4 else 0.0,
        )


def _source(count: int) -> tuple[PedigreeTableSource, _LazyAnimals]:
    animals = _LazyAnimals(count)
    pedigree = SimpleNamespace(pedigree=animals)
    return PedigreeTableSource(pedigree), animals


def test_model_construction_does_not_touch_rows():
    source, animals = _source(100_000)
    model = PedigreeTableModel(source)
    assert model.rowCount() == 100_000
    assert model.columnCount() == 8
    assert animals.gets == 0


def test_headers_use_table_column_titles():
    source, _animals = _source(3)
    model = PedigreeTableModel(source)
    headers = [
        model.headerData(column, Qt.Orientation.Horizontal, Qt.ItemDataRole.DisplayRole)
        for column in range(model.columnCount())
    ]
    assert headers == [
        "Original ID",
        "Animal ID",
        "Sire",
        "Dam",
        "Year",
        "Sex",
        "Name",
        "F",
    ]


def test_display_and_raw_roles():
    source, animals = _source(6)
    model = PedigreeTableModel(source)
    year = model.index(0, 4)
    fa = model.index(4, 7)
    name = model.index(0, 6)
    assert model.data(year, Qt.ItemDataRole.DisplayRole) == MISSING_DISPLAY
    assert model.data(year, Qt.ItemDataRole.UserRole) is None
    assert model.data(fa, Qt.ItemDataRole.DisplayRole) == "0.125000"
    assert model.data(fa, Qt.ItemDataRole.UserRole) == 0.125
    assert model.data(name, Qt.ItemDataRole.DisplayRole) == "a1"
    assert animals.gets == 5


def test_ids_are_not_formatted_as_floats():
    source, _animals = _source(1)
    model = PedigreeTableModel(source)
    animal_id = model.data(model.index(0, 1), Qt.ItemDataRole.DisplayRole)
    assert animal_id == "1"
    assert "." not in animal_id


def test_set_source_resets_and_empty_model():
    first, _a = _source(4)
    second, _b = _source(2)
    model = PedigreeTableModel(first)
    assert model.rowCount() == 4
    model.set_source(second)
    assert model.rowCount() == 2
    model.set_source(None)
    assert model.rowCount() == 0
    assert model.columnCount() == 8


def test_filter_matches_ids_and_names():
    source, _animals = _source(8)
    model = PedigreeTableModel(source)
    proxy = PedigreeFilterProxy()
    proxy.setSourceModel(model)
    proxy.set_query("a5")
    assert proxy.rowCount() == 1
    proxy.set_query("5")
    assert proxy.rowCount() == 1
    proxy.set_query("")
    assert proxy.rowCount() == 8


def test_format_display_value_does_not_mutate_raw():
    raw = 0.125
    text = format_display_value("fa", raw)
    assert text == "0.125000"
    assert raw == 0.125
