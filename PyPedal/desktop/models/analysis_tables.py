"""Presentation table models for analysis results. No scientific formulas."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, override

from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QObject,
    QPersistentModelIndex,
    Qt,
)

from PyPedal.application import InbreedingByYearRow, InbreedingResult, MatingCoIGroupResult
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value

ModelIndex = QModelIndex | QPersistentModelIndex


def format_inbreeding_percent(raw: object) -> str:
    """Breeder-facing inbreeding *F* as a percentage. Raw coefficients are unchanged."""
    if raw is None:
        return "—"
    if not isinstance(raw, int | float):
        return str(raw)
    percent = round(float(raw) * 100, 2)
    if percent == 0.0:
        percent = 0.0
    return f"{percent:.2f}%"


def format_effective_founders(raw: object) -> str:
    """Two-decimal effective-founder display. The scientific value is unchanged."""
    if raw is None:
        return "—"
    if not isinstance(raw, int | float):
        return str(raw)
    return f"{float(raw):.2f}"


def format_count(raw: object) -> str:
    """Integer count with locale-independent thousands separators."""
    if raw is None:
        return "—"
    if isinstance(raw, bool) or not isinstance(raw, int):
        return str(raw)
    return f"{raw:,}"


def _id_sort_key(animal_id: object) -> tuple[int, float | str]:
    if isinstance(animal_id, bool):
        return (2, str(animal_id))
    if isinstance(animal_id, int | float):
        return (0, float(animal_id))
    return (1, str(animal_id).casefold())


def _coef_sort_key(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        return float("-inf")
    return float(value)


class InbreedingResultTableModel(QAbstractTableModel):
    """Virtual table over ``InbreedingResult.fx`` (id list + mapping)."""

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._ids: list[int | str] = []
        self._fx: dict[int | str, float] = {}

    def set_result(self, result: InbreedingResult | None) -> None:
        self.beginResetModel()
        if result is None:
            self._ids = []
            self._fx = {}
        else:
            self._fx = result.fx
            self._ids = list(result.fx.keys())
        self.endResetModel()

    @override
    def rowCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(self._ids)

    @override
    def columnCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return 2

    @override
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return ("Animal ID", "F (%)")[section] if 0 <= section < 2 else None
        return str(section + 1)

    @override
    def data(self, index: ModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid() or not (0 <= index.row() < len(self._ids)):
            return None
        animal_id = self._ids[index.row()]
        coef = self._fx.get(animal_id)
        if role == Qt.ItemDataRole.UserRole:
            return animal_id if index.column() == 0 else coef
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if index.column() == 0:
            return str(animal_id)
        return format_inbreeding_percent(coef)

    @override
    def sort(
        self,
        column: int,
        order: Qt.SortOrder = Qt.SortOrder.AscendingOrder,
    ) -> None:
        reverse = order == Qt.SortOrder.DescendingOrder
        self.layoutAboutToBeChanged.emit()
        if column == 0:
            self._ids.sort(key=_id_sort_key, reverse=reverse)
        else:
            fx = self._fx
            self._ids.sort(key=lambda animal_id: _coef_sort_key(fx.get(animal_id)), reverse=reverse)
        self.layoutChanged.emit()


class YearInbreedingTableModel(QAbstractTableModel):
    """Table of year, count, and mean F from the application projection."""

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._rows: tuple[InbreedingByYearRow, ...] = ()

    def set_rows(self, rows: Sequence[InbreedingByYearRow] | None) -> None:
        self.beginResetModel()
        self._rows = tuple(rows or ())
        self.endResetModel()

    @override
    def rowCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(self._rows)

    @override
    def columnCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return 3

    @override
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            titles = ("Year", "Animals", "Mean F")
            return titles[section] if 0 <= section < 3 else None
        return str(section + 1)

    @override
    def data(self, index: ModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid() or not (0 <= index.row() < len(self._rows)):
            return None
        row = self._rows[index.row()]
        values: tuple[object, ...] = (row.year, row.n, row.mean)
        raw = values[index.column()]
        if role == Qt.ItemDataRole.UserRole:
            return raw
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if index.column() == 2:
            return format_display_value(FA_COLUMN, raw)
        return str(raw)


class MatingResultTableModel(QAbstractTableModel):
    """Explicit mating-group results. Not a Cartesian product."""

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._pairs: list[tuple[int | str, int | str]] = []
        self._values: dict[tuple[int | str, int | str], float] = {}

    def set_result(self, result: MatingCoIGroupResult | None) -> None:
        self.beginResetModel()
        if result is None:
            self._pairs = []
            self._values = {}
        else:
            self._values = result.matings
            self._pairs = list(result.matings.keys())
        self.endResetModel()

    @override
    def rowCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(self._pairs)

    @override
    def columnCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return 3

    @override
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            titles = ("Animal A", "Animal B", "F (%)")
            return titles[section] if 0 <= section < 3 else None
        return str(section + 1)

    @override
    def data(self, index: ModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid() or not (0 <= index.row() < len(self._pairs)):
            return None
        pair = self._pairs[index.row()]
        coef = self._values.get(pair)
        values: tuple[object, ...] = (pair[0], pair[1], coef)
        raw = values[index.column()]
        if role == Qt.ItemDataRole.UserRole:
            return raw
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if index.column() == 2:
            return format_inbreeding_percent(raw)
        return str(raw)
