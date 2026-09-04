"""Desktop analysis navigation, busy state, workers, cache, and F refresh."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QApplication, QFormLayout, QMessageBox

from PyPedal.application import FoundersOutcome, PedigreeOpenOptions
from PyPedal.application.tables import BROWSE_COLUMNS
from PyPedal.desktop.main_window import (
    PAGE_INBREEDING,
    PAGE_YEAR,
    MainWindow,
)
from PyPedal.desktop.models.analysis_tables import (
    InbreedingResultTableModel,
    MatingResultTableModel,
    format_inbreeding_percent,
)
from PyPedal.desktop.models.pedigree_table import FA_COLUMN, format_display_value
from PyPedal.desktop.pages.founders import FoundersPage
from PyPedal.desktop.pages.inbreeding import InbreedingPage
from PyPedal.desktop.settings import DesktopSettings
from PyPedal.pyp_results import EffectiveFoundersResult, InbreedingResult, MatingCoIGroupResult

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])

LOAD_TIMEOUT_MS = 30_000
MRODE = """\
1 0 0
2 0 0
3 1 2
4 1 0
5 4 3
6 5 2
"""
FA_INDEX = next(i for i, column in enumerate(BROWSE_COLUMNS) if column.key == "fa")


def _settings(tmp_path: Path) -> DesktopSettings:
    ini = tmp_path / "desktop.ini"
    return DesktopSettings(QSettings(str(ini), QSettings.Format.IniFormat))


def _wait_idle(qtbot: object, window: MainWindow) -> None:
    qtbot.waitUntil(
        lambda: window._job is None and not window._busy,
        timeout=LOAD_TIMEOUT_MS,
    )


def _load_mrode(qtbot: object, tmp_path: Path) -> MainWindow:
    source = tmp_path / "mrode.ped"
    source.write_text(MRODE, encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    window.open_path(source, PedigreeOpenOptions(separator=" ").normalized())
    _wait_idle(qtbot, window)
    return window


def test_analysis_navigation_and_run_disabled_when_empty(qtbot: object, tmp_path: Path) -> None:
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    assert window.nav.count() == 8
    assert window.inbreeding_page.run_button.isEnabled() is False
    window.analysis_inbreeding_action.trigger()
    assert window.nav.currentRow() == PAGE_INBREEDING
    window.analysis_year_action.trigger()
    assert window.nav.currentRow() == PAGE_YEAR


def test_inbreeding_run_refreshes_f_column_and_enables_export(
    qtbot: object, tmp_path: Path
) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        assert window.inbreeding_page.run_button.isEnabled() is True
        assert window.inbreeding_page.export_button.isEnabled() is False
        window.run_inbreeding_analysis()
        assert window._busy is True
        status = window.status_operation.text()
        assert "Calculating inbreeding" in status
        assert "6" in status
        assert window.progress.isHidden() is False
        assert window.inbreeding_page.run_button.isEnabled() is False
        assert window.about_action.isEnabled() is True
        _wait_idle(qtbot, window)
        assert window.session.inbreeding_result is not None
        assert window.inbreeding_page.export_button.isEnabled() is True
        assert window.inbreeding_page.model.rowCount() == 6
        animal_five = next(a for a in window.session.pedigree.pedigree if a.originalID == 5)
        row = next(
            index
            for index in range(window.animals_page.model.rowCount())
            if window.animals_page.model.source is not None
            and window.animals_page.model.source.value(index, 1) == animal_five.animalID
        )
        displayed = window.animals_page.model.data(window.animals_page.model.index(row, FA_INDEX))
        assert displayed == format_display_value(FA_COLUMN, 0.125)
        assert window.inbreeding_page.count_label.text() == "6"
    finally:
        close_owned_pypedal_log_handlers()


def test_year_analysis_reuses_cache_without_second_ml(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "years.ped"
    source.write_text("1 0 0 1990\n2 0 0 1990\n3 1 2 2000\n", encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    import PyPedal.application.jobs as jobs

    calls: list[int] = []
    original = jobs.inbreeding

    def wrapped(pedigree, **kwargs):
        calls.append(1)
        return original(pedigree, **kwargs)

    monkeypatch.setattr(jobs, "inbreeding", wrapped)
    try:
        window.open_path(source, PedigreeOpenOptions(pedformat="asdy").normalized())
        _wait_idle(qtbot, window)
        window.run_year_analysis()
        _wait_idle(qtbot, window)
        assert window.year_page.rows
        assert "Meuwissen" in window.year_page.status.text()
        window.run_year_analysis()
        _wait_idle(qtbot, window)
        assert len(calls) == 1
        assert "cached" in window.year_page.status.text()
    finally:
        close_owned_pypedal_log_handlers()


def test_relationship_and_mating_and_ne(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.relationship_page.selector_a.select_animal_id(4)
        window.relationship_page.selector_b.select_animal_id(3)
        assert window.relationship_page.run_button.isEnabled() is True
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert window.relationship_page.result is not None
        assert window.relationship_page.result.coefficient == 0.25
        window.mating_page.selector_a.select_animal_id(4)
        window.mating_page.selector_b.select_animal_id(3)
        assert window.mating_page.run_button.isEnabled() is True
        window.run_mating_pair()
        _wait_idle(qtbot, window)
        assert window.mating_page.pair_result is not None
        assert window.mating_page.pair_result.coefficient == 0.125
        assert window.mating_page.value_label.text() == "12.50%"
        window.run_theoretical_ne_analysis()
        assert window.population_page.value is not None
        assert window.population_page.export_button.isEnabled() is True
    finally:
        close_owned_pypedal_log_handlers()


def test_relationship_result_clears_when_selection_changes(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        page = window.relationship_page
        page.selector_a.select_animal_id(4)
        page.selector_b.select_animal_id(3)
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert page.result is not None
        assert page.value_label.text() != "—"
        assert page.export_button.isEnabled() is True

        page.selector_a.clear_selection()
        assert page.result is None
        assert page.value_label.text() == "—"
        assert page.run_button.isEnabled() is False
        assert page.export_button.isEnabled() is False

        page.selector_a.select_animal_id(4)
        page.selector_b.select_animal_id(3)
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert page.result is not None
        page.selector_b.clear_selection()
        assert page.result is None
        assert page.value_label.text() == "—"

        page.selector_a.select_animal_id(4)
        page.selector_b.select_animal_id(3)
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        page.selector_a.search.setText("5")
        assert page.selected_animal_a() is None
        assert page.result is None
        assert page.value_label.text() == "—"

        page.selector_a.select_animal_id(5)
        page.selector_b.select_animal_id(3)
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert page.result is not None
        page.selector_a.select_animal_id(4)
        assert page.result is None
        assert page.value_label.text() == "—"
        assert page.selected_animal_a() == 4
    finally:
        close_owned_pypedal_log_handlers()


def test_mating_pair_result_clears_but_group_is_kept(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        page = window.mating_page
        page.selector_a.select_animal_id(4)
        page.selector_b.select_animal_id(3)
        page.add_pair_button.click()
        assert page.pair_list.count() == 1
        window.run_mating_pair()
        _wait_idle(qtbot, window)
        assert page.pair_result is not None
        assert page.pair_result.coefficient == 0.125
        assert page.value_label.text() == "12.50%"

        page.selector_a.clear_selection()
        assert page.pair_result is None
        assert page.value_label.text() == "—"
        assert page.run_button.isEnabled() is False
        assert page.add_pair_button.isEnabled() is False
        assert page.pair_list.count() == 1

        page.selector_a.select_animal_id(4)
        window.run_mating_pair()
        _wait_idle(qtbot, window)
        assert page.pair_result is not None
        page.selector_b.search.setText("2")
        assert page.selected_animal_b() is None
        assert page.pair_result is None
        assert page.value_label.text() == "—"
        assert page.pair_list.count() == 1
    finally:
        close_owned_pypedal_log_handlers()


def test_mating_pair_displays_percent_not_fraction(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        page = window.mating_page
        assert page.value_label.text() == "—"
        page.selector_a.select_animal_id(4)
        page.selector_b.select_animal_id(3)
        window.run_mating_pair()
        _wait_idle(qtbot, window)
        assert page.pair_result is not None
        assert page.pair_result.coefficient == 0.125
        assert page.value_label.text() == "12.50%"
        page.selector_a.clear_selection()
        assert page.value_label.text() == "—"
        assert "0%" not in page.value_label.text()
    finally:
        close_owned_pypedal_log_handlers()


def test_mating_group_table_displays_percent() -> None:
    raw = 0.10095650884805218
    assert format_inbreeding_percent(raw) == "10.10%"
    assert format_inbreeding_percent(None) == "—"
    model = MatingResultTableModel()
    model.set_result(MatingCoIGroupResult({"matings": {(98001, 97984): raw}, "metadata": {}}))
    header = model.headerData(2, Qt.Orientation.Horizontal, Qt.ItemDataRole.DisplayRole)
    assert header == "F (%)"
    index = model.index(0, 2)
    assert model.data(index) == "10.10%"
    assert model.data(index, Qt.ItemDataRole.UserRole) == raw


def test_compute_stays_disabled_without_explicit_selection(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.relationship_page.id_a.setText("4")
        window.relationship_page.id_b.setText("99")
        window.relationship_page.selector_a.apply_search_now()
        window.relationship_page.selector_b.apply_search_now()
        assert window.relationship_page.selected_animal_a() is None
        assert window.relationship_page.selected_animal_b() is None
        assert window.relationship_page.run_button.isEnabled() is False
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert window.relationship_page.result is None
    finally:
        close_owned_pypedal_log_handlers()


def test_unexpected_error_shows_traceback_details(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: list[tuple[str, str]] = []

    def fake_exec(self: QMessageBox) -> int:
        captured.append((self.text(), self.detailedText()))
        return 0

    monkeypatch.setattr(QMessageBox, "exec", fake_exec)

    import PyPedal.application.jobs as jobs

    def boom(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("unexpected boom")

    monkeypatch.setattr(jobs, "inbreeding", boom)
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.run_inbreeding_analysis()
        _wait_idle(qtbot, window)
        assert captured
        _text, details = captured[0]
        assert "RuntimeError" in details
        assert "unexpected boom" in details
        assert "Traceback" in details
    finally:
        close_owned_pypedal_log_handlers()


def test_failed_load_retains_inbreeding_cache(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "PyPedal.desktop.main_window.show_application_error",
        lambda *_args, **_kwargs: None,
    )
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.run_inbreeding_analysis()
        _wait_idle(qtbot, window)
        cached = window.session.inbreeding_result
        assert cached is not None
        bad = tmp_path / "bad.ped"
        bad.write_text("1 0 0\n2 0 0\n3 1\n", encoding="utf-8")
        window.open_path(bad, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.session.inbreeding_result is cached
        assert window.inbreeding_page.export_button.isEnabled() is True
    finally:
        close_owned_pypedal_log_handlers()


def test_new_pedigree_clears_inbreeding_cache(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.run_inbreeding_analysis()
        _wait_idle(qtbot, window)
        assert window.session.inbreeding_result is not None
        other = tmp_path / "other.ped"
        other.write_text("10 0 0\n20 0 0\n", encoding="utf-8")
        window.open_path(other, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.session.inbreeding_result is None
        assert window.inbreeding_page.export_button.isEnabled() is False
    finally:
        close_owned_pypedal_log_handlers()


def test_close_clears_analysis_cache_and_pages(qtbot: object, tmp_path: Path) -> None:
    window = _load_mrode(qtbot, tmp_path)
    try:
        window.run_inbreeding_analysis()
        _wait_idle(qtbot, window)
        assert window.session.inbreeding_result is not None
        window.close_pedigree()
        assert window.session.inbreeding_result is None
        assert window.inbreeding_page.export_button.isEnabled() is False
        assert window.inbreeding_page.run_button.isEnabled() is False
    finally:
        close_owned_pypedal_log_handlers()


def _form_labels(form: QFormLayout) -> list[str]:
    labels: list[str] = []
    for row in range(form.rowCount()):
        item = form.itemAt(row, QFormLayout.ItemRole.LabelRole)
        widget = None if item is None else item.widget()
        if widget is not None:
            labels.append(widget.text())
    return labels


def test_founders_page_displays_two_decimals_not_raw() -> None:
    raw = 193.31434658869796
    result = EffectiveFoundersResult(
        {
            "fa_animal_count": 98001,
            "fa_founder_count": 7604,
            "fa_descendant_count": 91312,
            "fa_effective_founders": raw,
        }
    )
    page = FoundersPage()
    page.show_outcome(FoundersOutcome(result=result, implicit_renumber=False))
    assert page.result is result
    assert page.result.fa_effective_founders == raw
    assert page.value_label.text() == "193.31"
    assert page.animals_label.text() == "98,001"
    assert page.founder_count_label.text() == "7,604"
    assert page.descendant_label.text() == "91,312"


def test_inbreeding_summary_uses_percentages_and_breeder_labels() -> None:
    mean = 0.09313044278029989
    minimum = 0.0
    maximum = 0.546875
    result = InbreedingResult(
        {
            "fx": {1: 0.125, 2: 0.0},
            "metadata": {
                "all": {
                    "f_count": 98001,
                    "f_avg": mean,
                    "f_min": minimum,
                    "f_max": maximum,
                },
                "nonzero": {"f_count": 84442},
            },
        }
    )
    page = InbreedingPage()
    page.show_result(result)
    assert page.result is result
    assert page.result.metadata["all"]["f_avg"] == mean
    assert _form_labels(page.form) == [
        "Animals with results",
        "Mean inbreeding",
        "Minimum inbreeding",
        "Maximum inbreeding",
        "Animals with F > 0",
    ]
    assert page.count_label.text() == "98,001"
    assert page.mean_label.text() == "9.31%"
    assert page.min_label.text() == "0.00%"
    assert page.max_label.text() == "54.69%"
    assert page.positive_label.text() == "84,442"


def test_inbreeding_table_displays_percent_and_sorts_numerically() -> None:
    result = InbreedingResult(
        {
            "fx": {1: 0.09, 2: 0.125, 3: 0.0},
            "metadata": {},
        }
    )
    model = InbreedingResultTableModel()
    model.set_result(result)
    header = model.headerData(1, Qt.Orientation.Horizontal, Qt.ItemDataRole.DisplayRole)
    assert header == "F (%)"
    sample = model.index(1, 1)
    assert model.data(sample) == "12.50%"
    assert model.data(sample, Qt.ItemDataRole.UserRole) == 0.125
    model.sort(1, Qt.SortOrder.AscendingOrder)
    displayed = [model.data(model.index(row, 1)) for row in range(3)]
    raw = [model.data(model.index(row, 1), Qt.ItemDataRole.UserRole) for row in range(3)]
    assert displayed == ["0.00%", "9.00%", "12.50%"]
    assert raw == [0.0, 0.09, 0.125]
    assert "12.50%" < "9.00%"
