"""Named Griffon desktop acceptance: select by name, then compute."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers, named_griffon_path

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from PyPedal.application import PedigreeOpenOptions
from PyPedal.desktop.main_window import MainWindow
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
        rel.selector_a.search.setText(NAME_A)
        rel.selector_a.apply_search_now()
        assert rel.selector_a.choose_result_row(0) is True
        assert rel.selected_animal_a() == CURRENT_A
        rel.selector_b.search.setText(NAME_B)
        rel.selector_b.apply_search_now()
        assert rel.selector_b.choose_result_row(0) is True
        assert rel.selected_animal_b() == CURRENT_B
        assert rel.run_button.isEnabled() is True
        window.run_relationship_analysis()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and rel.result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert rel.result is not None
        assert abs(rel.result.coefficient - A_EXPECTED) < 1e-12

        mating = window.mating_page
        mating.selector_a.search.setText(NAME_A)
        mating.selector_a.apply_search_now()
        assert mating.selector_a.choose_result_row(0) is True
        mating.selector_b.search.setText(NAME_B)
        mating.selector_b.apply_search_now()
        assert mating.selector_b.choose_result_row(0) is True
        assert mating.run_button.isEnabled() is True
        window.run_mating_pair()
        qtbot.waitUntil(
            lambda: window._job is None and not window._busy and mating.pair_result is not None,
            timeout=ANALYSIS_TIMEOUT_MS,
        )
        assert mating.pair_result is not None
        assert abs(mating.pair_result.coefficient - F_EXPECTED) < 1e-12
    finally:
        close_owned_pypedal_log_handlers()
