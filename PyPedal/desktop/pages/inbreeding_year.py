"""Inbreeding by birth year, grouped from cached Meuwissen-Luo results."""

from __future__ import annotations

from collections.abc import Sequence

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QLabel, QTableView, QVBoxLayout, QWidget

from PyPedal.application import InbreedingByYearRow
from PyPedal.desktop.models.analysis_tables import YearInbreedingTableModel
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)


class InbreedingYearPage(QWidget):
    """Year, animal count, and mean F. Reuses the inbreeding cache."""

    run_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.rows: tuple[InbreedingByYearRow, ...] = ()
        self.model = YearInbreedingTableModel(self)
        self.status = QLabel("No year summary yet.")
        self.status.setObjectName("year_status")
        self.view = QTableView()
        self.view.setObjectName("year_table")
        self.view.setModel(self.model)
        self.view.setSortingEnabled(True)
        self.view.setAlternatingRowColors(True)
        self.view.verticalHeader().setVisible(False)
        self.run_button = make_run_button("Summarize by year")
        self.run_button.clicked.connect(self.run_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Inbreeding by Year",
            "Groups cached Meuwissen–Luo coefficients by recorded birth year. "
            "If inbreeding has not been run yet, it is computed once.",
        )
        layout.addWidget(self.run_button)
        layout.addWidget(self.export_button)
        layout.addWidget(self.status)
        layout.addWidget(self.view)

    def show_empty(self) -> None:
        self.rows = ()
        self.model.set_rows(())
        self.status.setText("No year summary yet.")
        self.export_button.setEnabled(False)

    def show_rows(
        self,
        rows: Sequence[InbreedingByYearRow],
        *,
        computed_inbreeding: bool,
    ) -> None:
        self.rows = tuple(rows)
        self.model.set_rows(self.rows)
        if computed_inbreeding:
            self.status.setText("Meuwissen–Luo was run once, then grouped by year.")
        else:
            self.status.setText("Grouped from the cached inbreeding result.")
        self.export_button.setEnabled(bool(self.rows))
