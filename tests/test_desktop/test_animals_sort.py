"""Animals table must not sort 98k rows on pedigree load."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from PyPedal.application import PedigreeTableSource
from PyPedal.desktop.models.pedigree_table import PedigreeFilterProxy
from PyPedal.desktop.pages.animals import AnimalsPage

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])

_ANIMAL_ID_COLUMN = 1
_NAME_COLUMN = 6


def _source(rows: list[tuple[int, str]]) -> PedigreeTableSource:
    animals = [
        SimpleNamespace(
            originalID=original,
            animalID=animal_id,
            sireID=0,
            damID=0,
            by=None,
            sex="u",
            name=name,
            fa=0.0,
        )
        for animal_id, (original, name) in enumerate(rows, start=1)
    ]
    return PedigreeTableSource(SimpleNamespace(pedigree=animals))


def _displayed_animal_ids(page: AnimalsPage) -> list[int]:
    ids: list[int] = []
    proxy = page.proxy
    for row in range(proxy.rowCount()):
        raw = proxy.data(proxy.index(row, _ANIMAL_ID_COLUMN), Qt.ItemDataRole.UserRole)
        ids.append(int(raw))
    return ids


def test_animals_page_starts_unsorted() -> None:
    page = AnimalsPage()
    assert page.proxy.sortColumn() == -1
    assert page.view.isSortingEnabled() is True
    assert page.view.horizontalHeader().sortIndicatorSection() == -1


def test_set_source_does_not_call_less_than(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = {"n": 0}
    original = PedigreeFilterProxy.lessThan

    def counted(self: PedigreeFilterProxy, left: object, right: object) -> bool:
        calls["n"] += 1
        return original(self, left, right)

    monkeypatch.setattr(PedigreeFilterProxy, "lessThan", counted)
    page = AnimalsPage()
    source = _source([(10, "c"), (20, "a"), (30, "b")] * 200)
    page.set_source(source)
    assert page.proxy.sortColumn() == -1
    assert page.proxy.rowCount() == 600
    assert calls["n"] == 0
    assert _displayed_animal_ids(page) == list(range(1, 601))


def test_header_sort_orders_by_name() -> None:
    page = AnimalsPage()
    page.set_source(_source([(1, "c"), (2, "a"), (3, "b")]))
    assert _displayed_animal_ids(page) == [1, 2, 3]
    page.view.sortByColumn(_NAME_COLUMN, Qt.SortOrder.AscendingOrder)
    assert page.proxy.sortColumn() == _NAME_COLUMN
    names = [
        page.proxy.data(page.proxy.index(row, _NAME_COLUMN), Qt.ItemDataRole.DisplayRole)
        for row in range(3)
    ]
    assert names == ["a", "b", "c"]
    assert _displayed_animal_ids(page) == [2, 3, 1]


def test_new_source_restores_pedigree_order_after_user_sort() -> None:
    page = AnimalsPage()
    first = _source([(1, "c"), (2, "a"), (3, "b")])
    page.set_source(first)
    page.view.sortByColumn(_NAME_COLUMN, Qt.SortOrder.AscendingOrder)
    assert page.proxy.sortColumn() == _NAME_COLUMN
    second = _source([(11, "z"), (12, "m"), (13, "n"), (14, "a")])
    page.set_source(second)
    assert page.proxy.sortColumn() == -1
    assert _displayed_animal_ids(page) == [1, 2, 3, 4]
    assert page.view.horizontalHeader().sortIndicatorSection() == -1
