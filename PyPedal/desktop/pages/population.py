"""Theoretical Ne from pedigree metadata."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QLabel, QVBoxLayout, QWidget

from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)

_ABSENT = "—"


class PopulationPage(QWidget):
    """Theoretical Ne from sire and dam counts. Not a census Ne."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.value: float | None = None
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("ne_value")
        form = QFormLayout()
        form.addRow("Theoretical Ne from metadata", self.value_label)
        self.run_button = make_run_button("Compute theoretical Ne")
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Population",
            "Theoretical effective population size from the numbers of unique "
            "sires and dams in the pedigree metadata. This is not a universally "
            "estimated Ne.",
        )
        layout.addLayout(form)
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addStretch(1)

    def show_empty(self) -> None:
        self.value = None
        self.value_label.setText(_ABSENT)
        self.export_button.setEnabled(False)

    def show_value(self, value: float) -> None:
        self.value = value
        self.value_label.setText(format_display_value(FA_COLUMN, value))
        self.export_button.setEnabled(True)
