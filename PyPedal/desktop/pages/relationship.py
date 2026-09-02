"""Pairwise relationship page."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QLabel, QLineEdit, QVBoxLayout, QWidget

from PyPedal.application import PairwiseResult
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)

_ABSENT = "—"


class RelationshipPage(QWidget):
    """Two current animal IDs and one relationship coefficient."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.result: PairwiseResult | None = None
        self.id_a = QLineEdit()
        self.id_a.setObjectName("relationship_id_a")
        self.id_b = QLineEdit()
        self.id_b.setObjectName("relationship_id_b")
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("relationship_value")
        form = QFormLayout()
        form.addRow("Animal A", self.id_a)
        form.addRow("Animal B", self.id_b)
        form.addRow("Relationship", self.value_label)
        self.run_button = make_run_button("Compute relationship")
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Relationship",
            "Coefficient of relationship for two current (renumbered) animal IDs. "
            "This is not a pedigree-wide relationship matrix.",
        )
        layout.addLayout(form)
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addStretch(1)

    def show_empty(self) -> None:
        self.result = None
        self.value_label.setText(_ABSENT)
        self.export_button.setEnabled(False)

    def show_result(self, result: PairwiseResult) -> None:
        self.result = result
        self.value_label.setText(format_display_value(FA_COLUMN, result.coefficient))
        self.export_button.setEnabled(True)
