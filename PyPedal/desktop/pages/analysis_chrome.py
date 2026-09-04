"""Shared chrome for analysis pages. Native Qt widgets, no QSS theme."""

from __future__ import annotations

from PySide6.QtWidgets import QFormLayout, QHeaderView, QLabel, QPushButton, QTableView, QVBoxLayout


def add_analysis_header(layout: QVBoxLayout, title: str, explanation: str) -> None:
    heading = QLabel(title)
    heading.setObjectName("analysis_title")
    note = QLabel(explanation)
    note.setWordWrap(True)
    note.setObjectName("analysis_explanation")
    layout.addWidget(heading)
    layout.addWidget(note)


def configure_growing_form(form: QFormLayout) -> None:
    """Let expanding fields (animal selectors) use leftover window width."""
    form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)


def configure_result_table(view: QTableView) -> None:
    """Native table settings so result columns remain usable when resized."""
    view.setSortingEnabled(True)
    view.setAlternatingRowColors(True)
    view.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
    view.verticalHeader().setVisible(False)
    view.horizontalHeader().setStretchLastSection(True)
    view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)


def make_run_button(text: str = "Run") -> QPushButton:
    button = QPushButton(text)
    button.setObjectName("analysis_run")
    return button


def make_export_button(text: str = "Export…") -> QPushButton:
    button = QPushButton(text)
    button.setObjectName("analysis_export")
    button.setEnabled(False)
    return button
