"""Mating CoI page: one pair, plus an explicit small group list."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from PyPedal.application import AnimalLookupIndex, MatingCoIGroupResult, PairwiseResult
from PyPedal.desktop.models.analysis_tables import MatingResultTableModel, format_inbreeding_percent
from PyPedal.desktop.pages.analysis_chrome import (
    add_analysis_header,
    configure_result_table,
    make_export_button,
    make_run_button,
)
from PyPedal.desktop.widgets.animal_selector import AnimalSelector

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
        self._armed = False
        self.selector_a = AnimalSelector(search_object_name="mating_id_a")
        self.selector_b = AnimalSelector(search_object_name="mating_id_b")
        self.id_a = self.selector_a.search
        self.id_b = self.selector_b.search
        self.selector_a.selection_changed.connect(self._on_selection_changed)
        self.selector_b.selection_changed.connect(self._on_selection_changed)
        self.value_label = QLabel(_ABSENT)
        self.value_label.setObjectName("mating_value")
        form = QFormLayout()
        form.addRow("Animal A", self.selector_a)
        form.addRow("Animal B", self.selector_b)
        form.addRow("Offspring inbreeding", self.value_label)

        self.pair_list = QListWidget()
        self.pair_list.setObjectName("mating_pairs")
        self.add_pair_button = QPushButton("Add pair to group")
        self.add_pair_button.setObjectName("mating_add_pair")
        self.add_pair_button.setEnabled(False)
        self.add_pair_button.clicked.connect(self._add_pair)
        self.clear_pairs_button = QPushButton("Clear group")
        self.clear_pairs_button.clicked.connect(self._clear_group)

        self.group_model = MatingResultTableModel(self)
        self.group_view = QTableView()
        self.group_view.setObjectName("mating_group_table")
        self.group_view.setModel(self.group_model)
        configure_result_table(self.group_view)

        self.run_button = make_run_button("Run pair")
        self.run_button.setEnabled(False)
        self.run_button.clicked.connect(self.run_pair_requested.emit)
        self.run_group_button = make_run_button("Run group")
        self.run_group_button.setObjectName("mating_run_group")
        self.run_group_button.setEnabled(False)
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
            "Prospective offspring inbreeding for explicit pairs. Search by "
            "display name, original ID, or current animal ID, then choose a "
            "result. Sex is shown so a mistaken selection is visible; animals "
            "are not swapped. Group mode evaluates only the pairs you add; it "
            "does not mate every animal with every other animal.",
        )
        layout.addLayout(form)
        layout.addWidget(self.add_pair_button)
        layout.addWidget(self.clear_pairs_button)
        layout.addWidget(self.pair_list)
        layout.addLayout(buttons)
        layout.addWidget(self.group_view)

    def set_lookup(self, index: AnimalLookupIndex | None) -> None:
        self.selector_a.set_index(index)
        self.selector_b.set_index(index)

    def set_armed(self, armed: bool) -> None:
        self._armed = armed
        self._update_actions()

    def selected_animal_a(self) -> int | None:
        return self.selector_a.selected_animal_id()

    def selected_animal_b(self) -> int | None:
        return self.selector_b.selected_animal_id()

    def group_pairs(self) -> list[tuple[str, str]]:
        pairs: list[tuple[str, str]] = []
        for index in range(self.pair_list.count()):
            item = self.pair_list.item(index)
            if item is None:
                continue
            payload = item.data(Qt.ItemDataRole.UserRole)
            if isinstance(payload, tuple) and len(payload) == 2:
                pairs.append((str(payload[0]), str(payload[1])))
                continue
            text = item.text()
            if "," in text:
                left, right = text.split(",", 1)
                pairs.append((left.strip(), right.strip()))
        return pairs

    def _add_pair(self) -> None:
        animal_a = self.selector_a.selected_animal_id()
        animal_b = self.selector_b.selected_animal_id()
        hit_a = self.selector_a.selected_hit()
        hit_b = self.selector_b.selected_hit()
        if animal_a is None or animal_b is None or hit_a is None or hit_b is None:
            return
        item = QListWidgetItem(f"{hit_a.label}  ×  {hit_b.label}")
        item.setData(Qt.ItemDataRole.UserRole, (animal_a, animal_b))
        self.pair_list.addItem(item)
        self._update_actions()

    def _clear_group(self) -> None:
        self.pair_list.clear()
        self._update_actions()

    def show_empty(self) -> None:
        self.pair_result = None
        self.group_result = None
        self.value_label.setText(_ABSENT)
        self.group_model.set_result(None)
        self.pair_list.clear()
        self.export_button.setEnabled(False)
        self.selector_a.clear_selection()
        self.selector_b.clear_selection()
        self._update_actions()

    def show_pair(self, result: PairwiseResult) -> None:
        self.pair_result = result
        self.value_label.setText(format_inbreeding_percent(result.coefficient))
        self._update_export()

    def show_group(self, result: MatingCoIGroupResult) -> None:
        self.group_result = result
        self.group_model.set_result(result)
        self._update_export()

    def _on_selection_changed(self) -> None:
        self._invalidate_pair_result()
        self._update_actions()

    def _invalidate_pair_result(self) -> None:
        self.pair_result = None
        self.value_label.setText(_ABSENT)
        self._update_export()

    def _update_export(self) -> None:
        self.export_button.setEnabled(self.pair_result is not None or self.group_result is not None)

    def _update_actions(self) -> None:
        both = (
            self.selector_a.selected_animal_id() is not None
            and self.selector_b.selected_animal_id() is not None
        )
        self.add_pair_button.setEnabled(self._armed and both)
        self.run_button.setEnabled(self._armed and both)
        self.run_group_button.setEnabled(self._armed and self.pair_list.count() > 0)
        self._update_export()
