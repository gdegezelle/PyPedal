"""Named Griffon desktop acceptance: select by name, then compute."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers, named_griffon_path

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QApplication, QToolButton

from PyPedal.application import PedigreeOpenOptions
from PyPedal.desktop.main_window import (
    PAGE_MATING,
    PAGE_RELATIONSHIP,
    MainWindow,
)
from PyPedal.desktop.settings import DesktopSettings

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])

LOAD_TIMEOUT_MS = 180_000
ANALYSIS_TIMEOUT_MS = 60_000
A_EXPECTED = 0.20191301769610437
F_EXPECTED = 0.10095650884805218
NAME_A = "Hierners Heartbreaker"
NAME_B = "Morning Bell Virgine"
CURRENT_A = 98001
CURRENT_B = 97984


def _settings(tmp_path: Path) -> DesktopSettings:
    ini = tmp_path / "desktop.ini"
    return DesktopSettings(QSettings(str(ini), QSettings.Format.IniFormat))


def _focus_editor(qtbot: object, selector: object) -> None:
    editor = selector.search
    window = editor.window()
    if window is not None:
        window.activateWindow()
    if QApplication.focusWidget() is editor:
        return
    qtbot.mouseClick(editor, Qt.MouseButton.LeftButton)
    if QApplication.focusWidget() is not editor:
        editor.setFocus(Qt.FocusReason.MouseFocusReason)
    qtbot.waitUntil(lambda: QApplication.focusWidget() is editor, timeout=5000)


def _type_query(qtbot: object, selector: object, text: str) -> None:
    editor = selector.search
    _focus_editor(qtbot, selector)
    qtbot.keyClicks(editor, text[0])
    qtbot.waitUntil(selector.popup_is_visible, timeout=5000)
    focused = QApplication.focusWidget()
    assert focused is editor
    if len(text) > 1:
        qtbot.keyClicks(focused, text[1:])
    assert editor.text() == text


def _commit_query(qtbot: object, selector: object, text: str, original_id: int) -> None:
    _type_query(qtbot, selector, text)
    selector.apply_search_now()
    qtbot.waitUntil(selector.popup_is_visible, timeout=5000)
    labels = selector.result_labels()
    row = next(i for i, label in enumerate(labels) if str(original_id) in label)
    assert selector.choose_result_row(row) is True
    qtbot.waitUntil(lambda: not selector.popup_is_visible(), timeout=2000)


def _clear_both(page: object) -> None:
    page.selector_a.clear_button.click()
    page.selector_b.clear_button.click()
    assert page.selector_a.selected_animal_id() is None
    assert page.selector_b.selected_animal_id() is None
    assert page.selector_a.search.text() == ""
    assert page.selector_b.search.text() == ""
    assert page.selector_a.search.isEnabled() is True
    assert page.selector_b.search.isEnabled() is True


@pytest.mark.integration
def test_named_griffon_relationship_and_mating_by_name(qtbot: object, tmp_path: Path) -> None:
    source = tmp_path / "griffonbruxellois_2026_named_pyp.ped"
    shutil.copy(named_griffon_path(), source)
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    try:
        window.open_path(
            source,
            PedigreeOpenOptions(pedformat="asdxbn", separator=",").normalized(),
        )
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and not window.session.is_empty,
            timeout=LOAD_TIMEOUT_MS,
        )
        pedigree = window.session.pedigree
        assert pedigree is not None
        assert len(pedigree.pedigree) == 98_001

        rel = window.relationship_page
        window.show()
        qtbot.waitExposed(window)
        window.nav.setCurrentRow(PAGE_RELATIONSHIP)
        qtbot.waitUntil(
            lambda: window.stack.currentWidget() is rel and rel.selector_a.search.isVisible(),
            timeout=5000,
        )
        _type_query(qtbot, rel.selector_a, "Hierners Heart")
        rel.selector_a.apply_search_now()
        qtbot.waitUntil(rel.selector_a.popup_is_visible, timeout=5000)
        qtbot.keyClick(rel.selector_a.search, Qt.Key.Key_Return)
        assert rel.selected_animal_a() == CURRENT_A
        assert rel.selector_a.search.text() == NAME_A
        assert rel.selector_a.summary.text() == "98685 — ♂ — 2024 — ID 98001"
        assert NAME_A not in rel.selector_a.summary.text()
        qtbot.waitUntil(lambda: not rel.selector_a.popup_is_visible(), timeout=2000)
        qtbot.keyClick(rel.selector_a.search, Qt.Key.Key_Tab)
        _type_query(qtbot, rel.selector_b, "Morning Bell Virg")
        rel.selector_b.apply_search_now()
        qtbot.waitUntil(rel.selector_b.popup_is_visible, timeout=5000)
        qtbot.keyClick(rel.selector_b.search, Qt.Key.Key_Return)
        assert rel.selected_animal_b() == CURRENT_B
        assert rel.selector_b.search.text() == NAME_B
        assert rel.selector_b.summary.text() == "98667 — ♀ — 2022 — ID 97984"
        assert NAME_B not in rel.selector_b.summary.text()
        assert rel.run_button.isEnabled() is True
        window.run_relationship_analysis()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and rel.result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert rel.result is not None
        assert abs(rel.result.coefficient - A_EXPECTED) < 1e-12
        _clear_both(rel)
        assert rel.result is None
        assert rel.value_label.text() == "—"
        assert rel.run_button.isEnabled() is False
        _commit_query(qtbot, rel.selector_a, "Colette", 20196)
        _commit_query(qtbot, rel.selector_b, "Colette", 20209)
        assert rel.selector_a.search.text() == "Colette"
        assert rel.selector_b.search.text() == "Colette"
        assert rel.selected_animal_a() != rel.selected_animal_b()
        window.run_relationship_analysis()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and rel.result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert rel.result is not None
        assert rel.value_label.text() != "—"

        mating = window.mating_page
        window.nav.setCurrentRow(PAGE_MATING)
        qtbot.waitUntil(
            lambda: window.stack.currentWidget() is mating and mating.selector_a.search.isVisible(),
            timeout=5000,
        )
        _commit_query(qtbot, mating.selector_a, "Hierners Heart", 98685)
        _commit_query(qtbot, mating.selector_b, "Morning Bell Virg", 98667)
        assert mating.selected_animal_a() == CURRENT_A
        assert mating.selected_animal_b() == CURRENT_B
        assert mating.run_button.isEnabled() is True
        window.run_mating_pair()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and mating.pair_result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert mating.pair_result is not None
        assert abs(mating.pair_result.coefficient - F_EXPECTED) < 1e-12
        assert mating.value_label.text() == "10.10%"

        _clear_both(mating)
        assert mating.pair_result is None
        assert mating.value_label.text() == "—"
        assert mating.run_button.isEnabled() is False
        _commit_query(qtbot, mating.selector_a, "Colette", 20196)
        _commit_query(qtbot, mating.selector_b, "Colette", 20209)
        assert mating.run_button.isEnabled() is True
        window.run_mating_pair()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and mating.pair_result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert mating.pair_result is not None
        assert mating.value_label.text() != "—"

        _clear_both(mating)
        assert mating.pair_result is None
        _commit_query(qtbot, mating.selector_a, "Hierners Heart", 98685)
        _commit_query(qtbot, mating.selector_b, "Morning Bell Virg", 98667)
        window.run_mating_pair()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and mating.pair_result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert mating.pair_result is not None
        assert abs(mating.pair_result.coefficient - F_EXPECTED) < 1e-12
        assert mating.value_label.text() == "10.10%"
        assert mating.selector_a.search.findChildren(QToolButton) == []
        assert mating.selector_b.search.findChildren(QToolButton) == []
        assert mating.selector_a.search.isClearButtonEnabled() is False
    finally:
        close_owned_pypedal_log_handlers()
