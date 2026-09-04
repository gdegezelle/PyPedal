"""Reusable animal selector: QLineEdit plus bounded QCompleter popup."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QModelIndex, QObject, Qt, QTimer, Signal
from PySide6.QtGui import QKeyEvent, QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (
    QCompleter,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListView,
    QPushButton,
    QSizePolicy,
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
_SEX_SYMBOLS = {"f": "♀", "female": "♀", "m": "♂", "male": "♂"}


def primary_display_text(hit: AnimalLookupHit) -> str:
    """Human-facing editor text for a committed hit. Not the full label."""
    name = hit.name.strip()
    if name:
        return name
    original = str(hit.original_id).strip()
    if original:
        return original
    return str(hit.animal_id)


def selection_detail_text(hit: AnimalLookupHit) -> str:
    """Committed metadata line. Does not repeat the editor's primary value."""
    primary = primary_display_text(hit)
    parts: list[str] = []
    original = str(hit.original_id).strip()
    if original and original != primary:
        parts.append(original)
    sex = _SEX_SYMBOLS.get(hit.sex.strip().casefold())
    if sex:
        parts.append(sex)
    if hit.birth_year is not None:
        parts.append(str(hit.birth_year))
    current = f"ID {hit.animal_id}"
    if str(hit.animal_id) != primary:
        parts.append(current)
    return " — ".join(parts)


class CompleterHitModel(QStandardItemModel):
    """At most ``DEFAULT_RESULT_LIMIT`` lookup hits, never the full pedigree."""

    def set_hits(self, hits: tuple[AnimalLookupHit, ...], *, truncated: bool) -> None:
        self.clear()
        for hit in hits:
            item = QStandardItem(hit.label)
            item.setEditable(False)
            item.setData(hit.animal_id, _HIT_ROLE)
            item.setData(primary_display_text(hit), Qt.ItemDataRole.EditRole)
            self.appendRow(item)
        if truncated:
            more = QStandardItem(_MORE_MATCHES_TEXT)
            more.setEnabled(False)
            more.setSelectable(False)
            more.setEditable(False)
            self.appendRow(more)


class AnimalSelector(QWidget):
    """Search by name, original ID, or current ID; commit only on explicit choice.

    Suggestions come from ``AnimalLookupIndex`` into a bounded QCompleter
    model (at most 50 hits). The completer popup is a non-focus view with
    this editor as its focus proxy, so typing continues in the QLineEdit.
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
        self.summary.setWordWrap(False)
        self.summary.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Preferred,
        )

        self.clear_button = QPushButton("Clear")
        self.clear_button.setObjectName(f"{search_object_name}_clear")
        self.clear_button.setEnabled(False)
        self.clear_button.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.clear_button.clicked.connect(self.clear_selection)

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(SEARCH_DEBOUNCE_MS)
        self._debounce.timeout.connect(self._search_now)
        self.search.textChanged.connect(self._on_text_changed)

        self._model = CompleterHitModel(self)
        self._completer = QCompleter(self._model, self)
        self._completer.setCompletionMode(QCompleter.CompletionMode.UnfilteredPopupCompletion)
        self._completer.setMaxVisibleItems(12)
        self._completer.setWidget(self.search)
        popup = QListView()
        popup.setObjectName(f"{search_object_name}_results")
        popup.setUniformItemSizes(True)
        popup.setEditTriggers(QListView.EditTrigger.NoEditTriggers)
        popup.setSelectionMode(QListView.SelectionMode.SingleSelection)
        popup.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        popup.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        popup.setFocusProxy(self.search)
        self._completer.setPopup(popup)
        popup.installEventFilter(self)
        self._completer.activated.connect(self._on_completer_activated)

        editor_row = QHBoxLayout()
        editor_row.addWidget(self.search, 1)
        editor_row.addWidget(self.clear_button, 0)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addLayout(editor_row)
        layout.addWidget(self.summary)

    def set_index(self, index: AnimalLookupIndex | None) -> None:
        """Install a lookup index and drop any previous selection."""
        self._index = index
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
        self._search_now(force=True)

    def _search_now(self, force: bool = False) -> None:
        query = self.search.text()
        if self._index is None or not query.strip():
            self._hits = ()
            self._model.set_hits((), truncated=False)
            self._hide_popup()
            return
        if (
            not force
            and self._selected is not None
            and query == primary_display_text(self._selected)
        ):
            self._hide_popup()
            return
        result = self._index.search(query, limit=DEFAULT_RESULT_LIMIT)
        self._hits = result.hits
        self._model.set_hits(result.hits, truncated=result.truncated)
        if not result.hits:
            self._hide_popup()
            return
        self._show_popup()

    def result_labels(self) -> list[str]:
        """Selectable popup labels, excluding the 'more matches' marker."""
        return [hit.label for hit in self._hits]

    def popup_is_visible(self) -> bool:
        popup = self._completer.popup()
        return popup is not None and popup.isVisible()

    def choose_result_row(self, row: int) -> bool:
        """Commit the selectable popup row ``row`` (tests)."""
        if row < 0 or row >= len(self._hits):
            return False
        self._commit(self._hits[row])
        return self._selected is not None

    def clear_selection(self) -> None:
        """Clear editor, popup, summary, and committed ID."""
        self._debounce.stop()
        self._hide_popup()
        self._hits = ()
        self._model.set_hits((), truncated=False)
        self.search.blockSignals(True)
        self.search.clear()
        self.search.blockSignals(False)
        changed = self._selected is not None
        self._selected = None
        self.summary.setText(_NONE_SELECTED)
        self.clear_button.setEnabled(False)
        if changed:
            self.selection_changed.emit()

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:  # noqa: N802
        popup = self._completer.popup()
        if event.type() != QEvent.Type.KeyPress or not isinstance(event, QKeyEvent):
            return super().eventFilter(watched, event)
        if watched is not self.search and watched is not popup:
            return super().eventFilter(watched, event)
        key = event.key()
        if key == Qt.Key.Key_Down:
            if not self.popup_is_visible():
                self.apply_search_now()
            self._move_current(1)
            return True
        if key == Qt.Key.Key_Up:
            if self.popup_is_visible():
                self._move_current(-1)
            return True
        if key in {Qt.Key.Key_Return, Qt.Key.Key_Enter}:
            if self.popup_is_visible():
                self._choose_current()
                return True
            return False
        if key == Qt.Key.Key_Escape:
            if self.popup_is_visible():
                self._hide_popup()
                return True
            return False
        if key == Qt.Key.Key_Tab:
            self._hide_popup()
            return False
        return super().eventFilter(watched, event)

    def _on_text_changed(self, text: str) -> None:
        if self._selected is not None and text != primary_display_text(self._selected):
            self._forget_selection()
        self._debounce.start()

    def _forget_selection(self) -> None:
        if self._selected is None:
            return
        self._selected = None
        self.summary.setText(_NONE_SELECTED)
        self.clear_button.setEnabled(False)
        self.selection_changed.emit()

    def _show_popup(self) -> None:
        popup = self._completer.popup()
        if popup is None:
            return
        popup.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating, True)
        popup.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        popup.setFocusProxy(self.search)
        self._completer.setCompletionPrefix("")
        width = max(self.search.width(), 360)
        popup.setMinimumWidth(width)
        self._completer.complete()
        popup.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        popup.setFocusProxy(self.search)
        model = popup.model()
        if model is not None and model.rowCount() > 0:
            popup.setCurrentIndex(model.index(0, 0))

    def _hide_popup(self) -> None:
        popup = self._completer.popup()
        if popup is None:
            return
        popup.hide()
        if popup.isVisible():
            popup.close()

    def _move_current(self, step: int) -> None:
        popup = self._completer.popup()
        if popup is None or not self.popup_is_visible() or self._model.rowCount() == 0:
            return
        selectable = [
            row
            for row in range(self._model.rowCount())
            if self._model.item(row) is not None and self._model.item(row).isEnabled()
        ]
        if not selectable:
            return
        model = popup.model()
        if model is None:
            return
        current = popup.currentIndex().row()
        if current < 0:
            popup.setCurrentIndex(model.index(selectable[0], 0))
            return
        target = current + step
        if target not in selectable:
            if step > 0:
                later = [row for row in selectable if row > current]
                target = later[0] if later else selectable[-1]
            else:
                earlier = [row for row in selectable if row < current]
                target = earlier[-1] if earlier else selectable[0]
        popup.setCurrentIndex(model.index(target, 0))

    def _choose_current(self) -> None:
        popup = self._completer.popup()
        if popup is None:
            return
        self._choose_index(popup.currentIndex())

    def _on_completer_activated(self, value: object) -> None:
        index = value if isinstance(value, QModelIndex) else None
        if index is None or not index.isValid():
            popup = self._completer.popup()
            index = popup.currentIndex() if popup is not None else QModelIndex()
        self._choose_index(index)

    def _choose_index(self, index: QModelIndex) -> None:
        if not index.isValid():
            return
        row = index.row()
        if 0 <= row < len(self._hits):
            self._commit(self._hits[row])
            return
        animal_id = index.data(_HIT_ROLE)
        if not isinstance(animal_id, int) or self._index is None:
            return
        hit = self._index.hit_for_animal_id(animal_id)
        if hit is None:
            return
        self._commit(hit)

    def _commit(self, hit: AnimalLookupHit) -> None:
        self._debounce.stop()
        self._hits = ()
        self._model.set_hits((), truncated=False)
        self._hide_popup()
        self._selected = hit
        self.search.blockSignals(True)
        self.search.setText(primary_display_text(hit))
        self.search.blockSignals(False)
        self.summary.setText(selection_detail_text(hit))
        self.clear_button.setEnabled(True)
        self._hide_popup()
        self.selection_changed.emit()
