"""L4: window smoke, worker load, failed-load retention."""

from __future__ import annotations

from pathlib import Path

import pytest
from _pedhelpers import close_owned_pypedal_log_handlers

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QApplication

from PyPedal.application import PedigreeOpenOptions
from PyPedal.application.lookup import AnimalLookupIndex
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


def test_main_window_smoke(qtbot: object, tmp_path: Path) -> None:
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    window.show()
    assert window.windowTitle().startswith("PyPedal")
    assert window.nav.count() == 8
    assert window.stack.count() == 8
    assert window.open_action.text() == "Open…"


def test_worker_loads_small_pedigree(qtbot: object, tmp_path: Path) -> None:
    source = tmp_path / "three.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    try:
        window.open_path(source, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.session.is_empty is False
        assert window.session.source_path == source.resolve()
        assert window.animals_page.model.rowCount() == 3
        assert window.close_action.isEnabled() is True
        assert "3" in window.status_count.text()
        assert window.status_file.text() == source.name
        assert window.metadata_page._banner.isHidden()
        assert window.metadata_page._form_host.isHidden() is False
        assert window.metadata_page._hint.isHidden()
        assert window.settings.recent_files()[0] == source.resolve()
        assert window.animals_page.proxy.sortColumn() == -1
        ids = [
            window.animals_page.proxy.data(
                window.animals_page.proxy.index(row, 1),
                Qt.ItemDataRole.UserRole,
            )
            for row in range(3)
        ]
        assert ids == [1, 2, 3]
    finally:
        close_owned_pypedal_log_handlers()


def test_failed_load_retains_previous_pedigree(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "PyPedal.desktop.main_window.show_application_error",
        lambda *args, **kwargs: None,
    )
    good = tmp_path / "good.ped"
    good.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    bad = tmp_path / "bad.ped"
    bad.write_text("1 0 0\n2 0 0\n3 1\n", encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    try:
        window.open_path(good, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.animals_page.model.rowCount() == 3
        installed = window.session.pedigree
        window.open_path(bad, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert window.session.pedigree is installed
        assert window.session.source_path == good.resolve()
        assert window.animals_page.model.rowCount() == 3
        assert window.settings.recent_files() == [good.resolve()]
    finally:
        close_owned_pypedal_log_handlers()


def test_desktop_display_does_not_rebuild_lookup(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    builds: list[int] = []
    original = AnimalLookupIndex.from_pedigree.__func__

    def counted(cls: type[AnimalLookupIndex], pedigree: object) -> AnimalLookupIndex:
        builds.append(1)
        return original(cls, pedigree)

    monkeypatch.setattr(AnimalLookupIndex, "from_pedigree", classmethod(counted))
    source = tmp_path / "three.ped"
    source.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    try:
        window.open_path(source, PedigreeOpenOptions(separator=" ").normalized())
        _wait_idle(qtbot, window)
        assert len(builds) == 1
        assert window.session.animal_lookup is not None
        assert len(window.session.animal_lookup) == 3
        assert window.relationship_page.selector_a.select_animal_id(1) is True
    finally:
        close_owned_pypedal_log_handlers()
