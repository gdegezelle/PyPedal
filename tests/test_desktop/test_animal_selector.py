"""Desktop animal selector: explicit choice, duplicates, keyboard, reload."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import NAMED_DUPLICATE_PED, close_owned_pypedal_log_handlers

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QApplication

from PyPedal.application import PedigreeOpenOptions
from PyPedal.desktop.main_window import MainWindow
from PyPedal.desktop.settings import DesktopSettings

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])

LOAD_TIMEOUT_MS = 30_000


def _settings(tmp_path: Path) -> DesktopSettings:
    ini = tmp_path / "desktop.ini"
    return DesktopSettings(QSettings(str(ini), QSettings.Format.IniFormat))


def _wait_idle(qtbot: object, window: MainWindow) -> None:
    qtbot.waitUntil(
        lambda: window._job is None and not window._busy,
        timeout=LOAD_TIMEOUT_MS,
    )


def _current_id(window: MainWindow, original_id: int) -> int:
    pedigree = window.session.pedigree
    assert pedigree is not None
    for animal in pedigree.pedigree:
        if int(animal.originalID) == original_id:
            return int(animal.animalID)
    raise AssertionError(f"missing originalID {original_id}")


def _load_named(qtbot: object, tmp_path: Path) -> MainWindow:
    source = tmp_path / "named.ped"
    source.write_text(NAMED_DUPLICATE_PED, encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    window.open_path(
        source,
        PedigreeOpenOptions(pedformat="asdxbn", separator=",").normalized(),
    )
    _wait_idle(qtbot, window)
    return window


def test_relationship_select_by_name_enables_compute(qtbot: object, tmp_path: Path) -> None:
    window = _load_named(qtbot, tmp_path)
    try:
        bella_old = _current_id(window, 101)
        bella_young = _current_id(window, 103)
        max_id = _current_id(window, 102)
        page = window.relationship_page
        assert page.run_button.isEnabled() is False
        page.selector_a.search.setText("bella")
        page.selector_a.apply_search_now()
        labels = page.selector_a.result_labels()
        assert len(labels) == 2
        assert all("Bella" in label for label in labels)
        assert any("101" in label for label in labels)
        assert any("103" in label for label in labels)
        assert page.selector_a.choose_result_row(0) is True
        assert page.selected_animal_a() == bella_old
        page.selector_b.search.setText("max")
        page.selector_b.apply_search_now()
        assert page.selector_b.choose_result_row(0) is True
        assert page.selected_animal_b() == max_id
        assert page.run_button.isEnabled() is True
        window.run_relationship_analysis()
        _wait_idle(qtbot, window)
        assert page.result is not None
        assert page.result.animal_a == bella_old
        assert page.result.animal_b == max_id
        page.selector_a.search.setText("spot")
        page.selector_a.apply_search_now()
        assert page.selected_animal_a() == bella_old
        page.selector_a.search.setText("bella")
        page.selector_a.apply_search_now()
        assert page.selector_a.choose_result_row(1) is True
        assert page.selected_animal_a() == bella_young
    finally:
        close_owned_pypedal_log_handlers()


def test_mating_select_by_name_and_duplicates(qtbot: object, tmp_path: Path) -> None:
    window = _load_named(qtbot, tmp_path)
    try:
        bella_young = _current_id(window, 103)
        max_id = _current_id(window, 102)
        page = window.mating_page
        page.selector_a.search.setText("Bella")
        page.selector_a.apply_search_now()
        assert len(page.selector_a.result_labels()) == 2
        assert page.selector_a.choose_result_row(1) is True
        assert page.selected_animal_a() == bella_young
        page.selector_b.search.setText("102")
        page.selector_b.apply_search_now()
        assert page.selector_b.choose_result_row(0) is True
        assert page.selected_animal_b() == max_id
        assert page.run_button.isEnabled() is True
        window.run_mating_pair()
        _wait_idle(qtbot, window)
        assert page.pair_result is not None
        assert page.pair_result.animal_a == bella_young
        assert page.pair_result.animal_b == max_id
    finally:
        close_owned_pypedal_log_handlers()


def test_enter_selects_highlighted_result(qtbot: object, tmp_path: Path) -> None:
    window = _load_named(qtbot, tmp_path)
    try:
        max_id = _current_id(window, 102)
        selector = window.relationship_page.selector_a
        window.show()
        qtbot.waitExposed(window)
        selector.search.setFocus(Qt.FocusReason.OtherFocusReason)
        selector.search.setText("Max")
        selector.apply_search_now()
        qtbot.waitUntil(selector.popup_is_visible, timeout=2000)
        qtbot.keyClick(selector.search, Qt.Key.Key_Return)
        assert selector.selected_animal_id() == max_id
        assert selector.popup_is_visible() is False
    finally:
        close_owned_pypedal_log_handlers()


def test_unnamed_pedigree_selects_by_current_id(qtbot: object, tmp_path: Path) -> None:
    source = tmp_path / "mrode.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n", encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    window.open_path(source, PedigreeOpenOptions(separator=" ").normalized())
    _wait_idle(qtbot, window)
    try:
        page = window.relationship_page
        page.selector_a.search.setText("4")
        page.selector_a.apply_search_now()
        assert page.selector_a.choose_result_row(0) is True
        assert page.selected_animal_a() == 4
        page.selector_b.search.setText("3")
        page.selector_b.apply_search_now()
        assert page.selector_b.choose_result_row(0) is True
        assert page.run_button.isEnabled() is True
        assert window.session.animal_lookup is not None
        assert window.session.animal_lookup.search("max").hits == ()
    finally:
        close_owned_pypedal_log_handlers()


def test_reload_clears_animal_selections(qtbot: object, tmp_path: Path) -> None:
    window = _load_named(qtbot, tmp_path)
    try:
        max_id = _current_id(window, 102)
        window.relationship_page.selector_a.select_animal_id(max_id)
        window.mating_page.selector_b.select_animal_id(max_id)
        assert window.relationship_page.selected_animal_a() == max_id
        other = tmp_path / "other.ped"
        other.write_text("10 0 0\n20 0 0\n", encoding="utf-8")
        window.open_path(other, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.relationship_page.selected_animal_a() is None
        assert window.mating_page.selected_animal_b() is None
        assert window.relationship_page.run_button.isEnabled() is False
    finally:
        close_owned_pypedal_log_handlers()
