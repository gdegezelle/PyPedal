"""Animals browse page with sortable, filterable table."""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import QLineEdit, QTableView, QVBoxLayout, QWidget

from PyPedal.application import BROWSE_COLUMNS, PedigreeTableSource
from PyPedal.desktop.models.pedigree_table import PedigreeFilterProxy, PedigreeTableModel

FILTER_DEBOUNCE_MS = 250
_ANIMAL_ID_COLUMN = next(i for i, column in enumerate(BROWSE_COLUMNS) if column.key == "animalID")


class AnimalsPage(QWidget):
    """Full pedigree table. No 500-row cap."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.model = PedigreeTableModel(parent=self)
        self.proxy = PedigreeFilterProxy(self)
        self.proxy.setSourceModel(self.model)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Search original ID, animal ID, or name")
        self.search.setClearButtonEnabled(True)
        self.search.setObjectName("animal_search")

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(FILTER_DEBOUNCE_MS)
        self._debounce.timeout.connect(self._apply_filter)
        self.search.textChanged.connect(self._schedule_filter)

        self.view = QTableView()
        self.view.setObjectName("animal_table")
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.setAlternatingRowColors(True)
        self.view.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.view.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self.view.verticalHeader().setVisible(False)
        header = self.view.horizontalHeader()
        header.setStretchLastSection(True)
        # Sorting stays available on header click. An active sort column
        # makes QSortFilterProxyModel re-sort every row on model reset
        # (millions of Python lessThan calls on a 98k pedigree).
        self.show_source_order()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(self.search)
        layout.addWidget(self.view)

    def show_source_order(self) -> None:
        """Show pedigree/source order. Header-click sorting remains enabled."""
        self.proxy.sort(-1)
        self.view.horizontalHeader().setSortIndicator(-1, Qt.SortOrder.AscendingOrder)

    def set_source(self, source: PedigreeTableSource | None) -> None:
        """Attach ``source`` in pedigree order. Does not sort on reset."""
        self.show_source_order()
        self.model.set_source(source)

    def _schedule_filter(self, _text: str) -> None:
        self._debounce.start()

    def _apply_filter(self) -> None:
        self.proxy.set_query(self.search.text())

    def apply_filter_now(self) -> None:
        """Apply the current search without waiting for debounce (tests)."""
        self._debounce.stop()
        self._apply_filter()

    def selected_animal_id(self) -> int | None:
        """Return the current-row animal ID, or ``None`` if nothing is selected."""
        selection = self.view.selectionModel()
        if selection is None:
            return None
        rows = selection.selectedRows()
        if not rows:
            return None
        raw = self.proxy.data(
            self.proxy.index(rows[0].row(), _ANIMAL_ID_COLUMN),
            Qt.ItemDataRole.UserRole,
        )
        if isinstance(raw, int):
            return raw
        if isinstance(raw, str) and raw.isdigit():
            return int(raw)
        return None
