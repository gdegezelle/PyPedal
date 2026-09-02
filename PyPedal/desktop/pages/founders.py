"""Lacy effective founders page."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QLabel, QVBoxLayout, QWidget

from PyPedal.application import EffectiveFoundersResult, FoundersOutcome
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)

_ABSENT = "—"


class FoundersPage(QWidget):
    """Concise Lacy effective-founder display."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.result: EffectiveFoundersResult | None = None
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("founders_value")
        self.animals_label = QLabel(_ABSENT)
        self.founder_count_label = QLabel(_ABSENT)
        self.descendant_label = QLabel(_ABSENT)
        form = QFormLayout()
        form.addRow("Effective founders", self.value_label)
        form.addRow("Animals", self.animals_label)
        form.addRow("Founders", self.founder_count_label)
        form.addRow("Descendants", self.descendant_label)
        self.run_button = make_run_button("Run effective founders")
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Effective Founders",
            "Lacy's effective number of founders. This does not form a dense relationship matrix.",
        )
        layout.addLayout(form)
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addStretch(1)

    def show_empty(self) -> None:
        self.result = None
        for label in (
            self.value_label,
            self.animals_label,
            self.founder_count_label,
            self.descendant_label,
        ):
            label.setText(_ABSENT)
        self.export_button.setEnabled(False)

    def show_outcome(self, outcome: FoundersOutcome) -> None:
        self.result = outcome.result
        result = outcome.result
        self.value_label.setText(format_display_value(FA_COLUMN, result.fa_effective_founders))
        self.animals_label.setText(str(result.fa_animal_count))
        self.founder_count_label.setText(str(result.fa_founder_count))
        self.descendant_label.setText(str(result.fa_descendant_count))
        self.export_button.setEnabled(True)
