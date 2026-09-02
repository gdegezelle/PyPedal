"""Mating CoI page: one pair, plus an explicit small group list."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from PyPedal.application import MatingCoIGroupResult, PairwiseResult
from PyPedal.desktop.models.analysis_tables import MatingResultTableModel
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    make_export_button,
    make_run_button,
)

_ABSENT = "—"


class MatingPage(QWidget):
    """Single-pair mating CoI and optional explicit-pair group mode."""

    run_pair_requested = Signal()
    run_group_requested = Signal()
    export_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.pair_result: PairwiseResult | None = None
        self.group_result: MatingCoIGroupResult | None = None
        self.id_a = QLineEdit()
        self.id_a.setObjectName("mating_id_a")
        self.id_b = QLineEdit()
        self.id_b.setObjectName("mating_id_b")
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("mating_value")
        form = QFormLayout()
        form.addRow("Animal A", self.id_a)
        form.addRow("Animal B", self.id_b)
        form.addRow("Offspring F", self.value_label)

        self.pair_list = QListWidget()
        self.pair_list.setObjectName("mating_pairs")
        self.add_pair_button = QPushButton("Add pair to group")
        self.add_pair_button.setObjectName("mating_add_pair")
        self.add_pair_button.clicked.connect(self._add_pair)
        self.clear_pairs_button = QPushButton("Clear group")
        self.clear_pairs_button.clicked.connect(self.pair_list.clear)

        self.group_model = MatingResultTableModel(self)
        self.group_view = QTableView()
        self.group_view.setObjectName("mating_group_table")
        self.group_view.setModel(self.group_model)
        self.group_view.verticalHeader().setVisible(False)

        self.run_button = make_run_button("Run pair")
        self.run_button.clicked.connect(self.run_pair_requested.emit)
        self.run_group_button = make_run_button("Run group")
        self.run_group_button.setObjectName("mating_run_group")
        self.run_group_button.clicked.connect(self.run_group_requested.emit)
        self.export_button = make_export_button()
        self.export_button.clicked.connect(self.export_requested.emit)

        buttons = QHBoxLayout()
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.run_group_button)
        buttons.addWidget(self.export_button)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        add_analysis_header(
            layout,
            "Mating",
            "Prospective offspring inbreeding for explicit pairs. Group mode "
            "evaluates only the pairs you add; it does not mate every animal "
            "with every other animal.",
        )
        layout.addLayout(form)
        layout.addWidget(self.add_pair_button)
        layout.addWidget(self.clear_pairs_button)
        layout.addWidget(self.pair_list)
        layout.addLayout(buttons)
        layout.addWidget(self.group_view)

    def group_pairs(self) -> list[tuple[str, str]]:
        pairs: list[tuple[str, str]] = []
        for index in range(self.pair_list.count()):
            item = self.pair_list.item(index)
            if item is None:
                continue
            left, right = item.text().split(",", 1)
            pairs.append((left, right))
        return pairs

    def _add_pair(self) -> None:
        text_a = self.id_a.text().strip()
        text_b = self.id_b.text().strip()
        if not text_a or not text_b:
            return
        self.pair_list.addItem(f"{text_a},{text_b}")

    def show_empty(self) -> None:
        self.pair_result = None
        self.group_result = None
        self.value_label.setText(_ABSENT)
        self.group_model.set_result(None)
        self.pair_list.clear()
        self.export_button.setEnabled(False)

    def show_pair(self, result: PairwiseResult) -> None:
        self.pair_result = result
        self.value_label.setText(format_display_value(FA_COLUMN, result.coefficient))
        self.export_button.setEnabled(True)

    def show_group(self, result: MatingCoIGroupResult) -> None:
        self.group_result = result
        self.group_model.set_result(result)
        self.export_button.setEnabled(True)
