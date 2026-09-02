"""Qt table model over ``PedigreeTableSource``.

The model holds a source reference. It does not copy animals into items,
dicts, or widgets at construction time.
"""

from __future__ import annotations

from typing import Any, override

from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QObject,
    QPersistentModelIndex,
    QSortFilterProxyModel,
    Qt,
)

from PyPedal.application import BROWSE_COLUMNS, PedigreeTableSource

MISSING_DISPLAY = "—"
FA_DISPLAY_DIGITS = 6
ID_COLUMNS = frozenset({"originalID", "animalID"})
NAME_COLUMN = "name"
YEAR_COLUMN = "year"
FA_COLUMN = "fa"
UNKNOWN_YEAR_TOKENS = frozenset({None, "", 0, -999})

ModelIndex = QModelIndex | QPersistentModelIndex


def format_display_value(column_key: str, raw: object) -> str:
    """Format a raw application value for ``DisplayRole`` only."""
    if column_key == YEAR_COLUMN and raw in UNKNOWN_YEAR_TOKENS:
        return MISSING_DISPLAY
    if raw is None:
        return MISSING_DISPLAY
    if column_key == FA_COLUMN:
        if not isinstance(raw, int | float):
            return str(raw)
        rounded = round(float(raw), FA_DISPLAY_DIGITS)
        if rounded == 0.0:
            rounded = 0.0
        return f"{rounded:.{FA_DISPLAY_DIGITS}f}"
    return str(raw)


def _sort_key(value: object) -> tuple[int, float | str]:
    if value is None:
        return (0, "")
    if isinstance(value, bool):
        return (1, float(value))
    if isinstance(value, int | float):
        return (1, float(value))
    return (2, str(value).casefold())


class PedigreeTableModel(QAbstractTableModel):
    """On-demand table model backed by ``PedigreeTableSource``."""

    def __init__(
        self,
        source: PedigreeTableSource | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._source = source

    @property
    def source(self) -> PedigreeTableSource | None:
        return self._source

    def set_source(self, source: PedigreeTableSource | None) -> None:
        self.beginResetModel()
        self._source = source
        self.endResetModel()

    def refresh_column(self, column: int) -> None:
        """Notify views that one column changed (for later analysis updates)."""
        if self.rowCount() == 0:
            return
        top_left = self.index(0, column)
        bottom_right = self.index(self.rowCount() - 1, column)
        self.dataChanged.emit(top_left, bottom_right, [Qt.ItemDataRole.DisplayRole])

    @override
    def rowCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid() or self._source is None:
            return 0
        return self._source.row_count()

    @override
    def columnCount(self, parent: ModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        if self._source is None:
            return len(BROWSE_COLUMNS)
        return self._source.column_count()

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
            if 0 <= section < len(BROWSE_COLUMNS):
                return BROWSE_COLUMNS[section].title
            return None
        return str(section + 1)

    @override
    def data(self, index: ModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid() or self._source is None:
            return None
        raw = self._source.value(index.row(), index.column())
        column_key = BROWSE_COLUMNS[index.column()].key
        if role == Qt.ItemDataRole.DisplayRole:
            return format_display_value(column_key, raw)
        if role == Qt.ItemDataRole.UserRole:
            return raw
        return None


class PedigreeFilterProxy(QSortFilterProxyModel):
    """Filter original ID, animal ID, and name without copying the pedigree."""

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._query = ""
        self.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.setSortRole(Qt.ItemDataRole.UserRole)

    def set_query(self, text: str) -> None:
        self.beginFilterChange()
        self._query = text.strip()
        self.endFilterChange(QSortFilterProxyModel.Direction.Rows)

    @override
    def filterAcceptsRow(self, source_row: int, source_parent: ModelIndex) -> bool:
        if not self._query:
            return True
        model = self.sourceModel()
        if model is None:
            return False
        needle = self._query.casefold()
        for column, field in enumerate(BROWSE_COLUMNS):
            if field.key not in ID_COLUMNS and field.key != NAME_COLUMN:
                continue
            index = model.index(source_row, column, source_parent)
            raw = model.data(index, Qt.ItemDataRole.UserRole)
            display = str(raw if raw is not None else "")
            if field.key in ID_COLUMNS and display == self._query:
                return True
            if needle in display.casefold():
                return True
        return False

    @override
    def lessThan(self, left: ModelIndex, right: ModelIndex) -> bool:
        model = self.sourceModel()
        if model is None:
            return False
        left_key = _sort_key(model.data(left, Qt.ItemDataRole.UserRole))
        right_key = _sort_key(model.data(right, Qt.ItemDataRole.UserRole))
        return left_key < right_key
