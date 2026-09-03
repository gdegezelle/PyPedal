"""Reusable animal selector: search field, bounded popup, explicit choice."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QModelIndex, QObject, Qt, QTimer, Signal
from PySide6.QtGui import QKeyEvent, QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListView,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from PyPedal.application.lookup import (
    DEFAULT_RESULT_LIMIT,
    AnimalLookupHit,
    AnimalLookupIndex,
)

SEARCH_DEBOUNCE_MS = 200
_MORE_MATCHES_TEXT = "More matches exist — refine the search"
_HIT_ROLE = Qt.ItemDataRole.UserRole
_NONE_SELECTED = "No animal selected"


class AnimalSelector(QWidget):
    """Search by name, original ID, or current ID; commit only on explicit choice.

    Changing the search text never changes the selected animal ID. A pedigree
    reload replaces the index and clears the selection.
    """

    selection_changed = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        search_object_name: str = "animal_search",
    ) -> None:
        super().__init__(parent)
        self._index: AnimalLookupIndex | None = None
        self._selected: AnimalLookupHit | None = None
        self._hits: tuple[AnimalLookupHit, ...] = ()

        self.search = QLineEdit()
        self.search.setObjectName(search_object_name)
        self.search.setPlaceholderText("Name, original ID, or current ID")
        self.search.setClearButtonEnabled(True)
        self.search.installEventFilter(self)

        self.summary = QLabel(_NONE_SELECTED)
        self.summary.setObjectName(f"{search_object_name}_summary")
        self.summary.setWordWrap(True)

        self.clear_button = QPushButton("Clear")
        self.clear_button.setObjectName(f"{search_object_name}_clear")
        self.clear_button.setEnabled(False)
        self.clear_button.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.clear_button.clicked.connect(self.clear_selection)

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(SEARCH_DEBOUNCE_MS)
        self._debounce.timeout.connect(self._search_now)
        self.search.textChanged.connect(self._schedule_search)

        self._model = QStandardItemModel(self)
        self._popup = QFrame(None, Qt.WindowType.Popup)
        self._popup.setObjectName(f"{search_object_name}_popup")
        self._list = QListView(self._popup)
        self._list.setObjectName(f"{search_object_name}_results")
        self._list.setModel(self._model)
        self._list.setUniformItemSizes(True)
        self._list.setEditTriggers(QListView.EditTrigger.NoEditTriggers)
        self._list.setSelectionMode(QListView.SelectionMode.SingleSelection)
        self._list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._list.clicked.connect(self._on_result_clicked)

        popup_layout = QVBoxLayout(self._popup)
        popup_layout.setContentsMargins(0, 0, 0, 0)
        popup_layout.addWidget(self._list)
        self.destroyed.connect(self._popup.deleteLater)

        selected_row = QHBoxLayout()
        selected_row.addWidget(self.summary, 1)
        selected_row.addWidget(self.clear_button, 0)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(self.search)
        layout.addLayout(selected_row)

    def set_index(self, index: AnimalLookupIndex | None) -> None:
        """Install a lookup index and drop any previous selection."""
        self._index = index
        self.search.blockSignals(True)
        self.search.clear()
        self.search.blockSignals(False)
        self._hide_popup()
        self._model.clear()
        self._hits = ()
        self.clear_selection()

    def selected_animal_id(self) -> int | None:
        if self._selected is None:
            return None
        return self._selected.animal_id

    def selected_hit(self) -> AnimalLookupHit | None:
        return self._selected

    def select_animal_id(self, animal_id: int) -> bool:
        """Commit ``animal_id`` when it exists in the current index (tests)."""
        if self._index is None:
            return False
        hit = self._index.hit_for_animal_id(animal_id)
        if hit is None:
            return False
        self._commit(hit)
        return True

    def apply_search_now(self) -> None:
        """Run the current query without waiting for debounce (tests)."""
        self._debounce.stop()
        self._search_now()

    def result_labels(self) -> list[str]:
        """Selectable popup labels, excluding the 'more matches' marker."""
        labels: list[str] = []
        for row in range(self._model.rowCount()):
            item = self._model.item(row)
            if item is None or not item.isEnabled():
                continue
            labels.append(item.text())
        return labels

    def popup_is_visible(self) -> bool:
        return self._popup.isVisible()

    def choose_result_row(self, row: int) -> bool:
        """Commit the selectable popup row ``row`` (tests)."""
        item = self._model.item(row)
        if item is None or not item.isEnabled():
            return False
        self._choose_index(self._model.index(row, 0))
        return self._selected is not None

    def clear_selection(self) -> None:
        changed = self._selected is not None
        self._selected = None
        self.summary.setText(_NONE_SELECTED)
        self.clear_button.setEnabled(False)
        if changed:
            self.selection_changed.emit()

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:  # noqa: N802
        if watched is not self.search or event.type() != QEvent.Type.KeyPress:
            return super().eventFilter(watched, event)
        if not isinstance(event, QKeyEvent):
            return super().eventFilter(watched, event)
        key = event.key()
        if key == Qt.Key.Key_Down:
            if not self._popup.isVisible():
                self.apply_search_now()
            self._move_current(1)
            return True
        if key == Qt.Key.Key_Up:
            if self._popup.isVisible():
                self._move_current(-1)
            return True
        if key in {Qt.Key.Key_Return, Qt.Key.Key_Enter}:
            if self._popup.isVisible():
                self._choose_current()
                return True
            return False
        if key == Qt.Key.Key_Escape:
            if self._popup.isVisible():
                self._hide_popup()
                return True
            return False
        if key == Qt.Key.Key_Tab:
            self._hide_popup()
            return False
        return super().eventFilter(watched, event)

    def _schedule_search(self, _text: str) -> None:
        self._debounce.start()

    def _search_now(self) -> None:
        query = self.search.text()
        if self._index is None or not query.strip():
            self._hits = ()
            self._model.clear()
            self._hide_popup()
            return
        result = self._index.search(query, limit=DEFAULT_RESULT_LIMIT)
        self._hits = result.hits
        self._model.clear()
        for hit in result.hits:
            item = QStandardItem(hit.label)
            item.setEditable(False)
            item.setData(hit.animal_id, _HIT_ROLE)
            self._model.appendRow(item)
        if result.truncated:
            more = QStandardItem(_MORE_MATCHES_TEXT)
            more.setEnabled(False)
            more.setSelectable(False)
            more.setEditable(False)
            self._model.appendRow(more)
        if self._model.rowCount() == 0:
            self._hide_popup()
            return
        self._list.setCurrentIndex(self._model.index(0, 0))
        self._show_popup()

    def _show_popup(self) -> None:
        width = max(self.search.width(), 360)
        row_height = self._list.sizeHintForRow(0)
        if row_height <= 0:
            row_height = 24
        rows = min(max(self._model.rowCount(), 1), 12)
        self._popup.setFixedWidth(width)
        self._popup.setFixedHeight(row_height * rows + 4)
        origin = self.search.mapToGlobal(self.search.rect().bottomLeft())
        self._popup.move(origin)
        self._popup.show()

    def _hide_popup(self) -> None:
        self._popup.hide()

    def _move_current(self, step: int) -> None:
        if not self._popup.isVisible() or self._model.rowCount() == 0:
            return
        current = self._list.currentIndex().row()
        if current < 0:
            current = 0
        target = current + step
        selectable = [
            row
            for row in range(self._model.rowCount())
            if self._model.item(row) is not None and self._model.item(row).isEnabled()
        ]
        if not selectable:
            return
        if target not in selectable:
            if step > 0:
                later = [row for row in selectable if row > current]
                target = later[0] if later else selectable[-1]
            else:
                earlier = [row for row in selectable if row < current]
                target = earlier[-1] if earlier else selectable[0]
        self._list.setCurrentIndex(self._model.index(target, 0))

    def _choose_current(self) -> None:
        self._choose_index(self._list.currentIndex())

    def _on_result_clicked(self, index: QModelIndex) -> None:
        self._choose_index(index)

    def _choose_index(self, index: QModelIndex) -> None:
        if not index.isValid():
            return
        item = self._model.item(index.row())
        if item is None or not item.isEnabled():
            return
        animal_id = item.data(_HIT_ROLE)
        if not isinstance(animal_id, int) or self._index is None:
            return
        hit = self._index.hit_for_animal_id(animal_id)
        if hit is None:
            return
        self._commit(hit)

    def _commit(self, hit: AnimalLookupHit) -> None:
        self._selected = hit
        self.summary.setText(hit.label)
        self.clear_button.setEnabled(True)
        self._hide_popup()
        self.selection_changed.emit()
