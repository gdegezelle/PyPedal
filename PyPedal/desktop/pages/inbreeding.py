"""Meuwissen-Luo inbreeding page."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QLabel, QTableView, QVBoxLayout, QWidget

from PyPedal.application import InbreedingResult
from PyPedal.desktop.models.analysis_tables import InbreedingResultTableModel
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)

_ABSENT = "—"


class InbreedingPage(QWidget):
    """Summary plus full Fx table. Run lives here; progress is global."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.result: InbreedingResult | None = None
        self.model = InbreedingResultTableModel(self)

        self.count_label = QLabel(_ABSENT)
        self.mean_label = QLabel(_ABSENT)
        self.min_label = QLabel(_ABSENT)
        self.max_label = QLabel(_ABSENT)
        self.positive_label = QLabel(_ABSENT)
        for name, widget in (
            ("inbreeding_count", self.count_label),
            ("inbreeding_mean", self.mean_label),
            ("inbreeding_min", self.min_label),
            ("inbreeding_max", self.max_label),
            ("inbreeding_positive", self.positive_label),
        ):
            widget.setObjectName(name)

        form = QFormLayout()
        form.addRow("Animals with results", self.count_label)
        form.addRow("Mean F", self.mean_label)
        form.addRow("Min F", self.min_label)
        form.addRow("Max F", self.max_label)
        form.addRow("Count F > 0", self.positive_label)

        self.view = QTableView()
        self.view.setObjectName("inbreeding_table")
        self.view.setModel(self.model)
        self.view.setSortingEnabled(True)
        self.view.setAlternatingRowColors(True)
        self.view.verticalHeader().setVisible(False)

        self.run_button = make_run_button("Run inbreeding")
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Inbreeding",
            "Meuwissen–Luo coefficients of inbreeding. This updates each animal's F.",
        )
        layout.addLayout(form)
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addWidget(self.view)

    def show_empty(self) -> None:
        self.result = None
        self.model.set_result(None)
        for label in (
            self.count_label,
            self.mean_label,
            self.min_label,
            self.max_label,
            self.positive_label,
        ):
            label.setText(_ABSENT)
        self.export_button.setEnabled(False)

    def show_result(self, result: InbreedingResult) -> None:
        self.result = result
        self.model.set_result(result)
        all_stats = result.metadata.get("all", {})
        nonzero = result.metadata.get("nonzero", {})
        self.count_label.setText(str(all_stats.get("f_count", _ABSENT)))
        self.mean_label.setText(format_display_value(FA_COLUMN, all_stats.get("f_avg")))
        self.min_label.setText(format_display_value(FA_COLUMN, all_stats.get("f_min")))
        self.max_label.setText(format_display_value(FA_COLUMN, all_stats.get("f_max")))
        self.positive_label.setText(str(nonzero.get("f_count", _ABSENT)))
        self.export_button.setEnabled(True)
