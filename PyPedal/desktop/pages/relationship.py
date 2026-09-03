"""Pairwise relationship page."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QLabel, QVBoxLayout, QWidget

from PyPedal.application import AnimalLookupIndex, PairwiseResult
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)
from PyPedal.desktop.widgets.animal_selector import AnimalSelector

_ABSENT = "—"


class RelationshipPage(QWidget):
    """Two explicitly selected animals and one relationship coefficient."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.result: PairwiseResult | None = None
        self._armed = False
        self.selector_a = AnimalSelector(search_object_name="relationship_id_a")
        self.selector_b = AnimalSelector(search_object_name="relationship_id_b")
        self.id_a = self.selector_a.search
        self.id_b = self.selector_b.search
        self.selector_a.selection_changed.connect(self._update_run_enabled)
        self.selector_b.selection_changed.connect(self._update_run_enabled)
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("relationship_value")
        form = QFormLayout()
        form.addRow("Animal A", self.selector_a)
        form.addRow("Animal B", self.selector_b)
        form.addRow("Relationship", self.value_label)
        self.run_button = make_run_button("Compute relationship")
        self.run_button.setEnabled(False)
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Relationship",
            "Coefficient of relationship for two animals. Search by display "
            "name, original ID, or current animal ID, then choose a result. "
            "Duplicate names are not unique identities — pick the exact animal. "
            "This is not a pedigree-wide relationship matrix.",
        )
        layout.addLayout(form)
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addStretch(1)

    def set_lookup(self, index: AnimalLookupIndex | None) -> None:
        self.selector_a.set_index(index)
        self.selector_b.set_index(index)

    def set_armed(self, armed: bool) -> None:
        self._armed = armed
        self._update_run_enabled()

    def selected_animal_a(self) -> int | None:
        return self.selector_a.selected_animal_id()

    def selected_animal_b(self) -> int | None:
        return self.selector_b.selected_animal_id()

    def show_empty(self) -> None:
        self.result = None
        self.value_label.setText(_ABSENT)
        self.export_button.setEnabled(False)
        self.selector_a.clear_selection()
        self.selector_b.clear_selection()

    def show_result(self, result: PairwiseResult) -> None:
        self.result = result
        self.value_label.setText(format_display_value(FA_COLUMN, result.coefficient))
        self.export_button.setEnabled(True)

    def _update_run_enabled(self) -> None:
        ready = (
            self._armed
            and self.selector_a.selected_animal_id() is not None
            and self.selector_b.selected_animal_id() is not None
        )
        self.run_button.setEnabled(ready)
