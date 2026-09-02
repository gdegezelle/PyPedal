"""L3: open dialog, pages, menus, recent files, settings, About."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QSettings
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QApplication, QMessageBox

from PyPedal.application import PedigreeTableSource
from PyPedal.desktop.dialogs.open_pedigree import OpenPedigreeDialog
from PyPedal.desktop.main_window import MainWindow
from PyPedal.desktop.pages.animals import AnimalsPage
from PyPedal.desktop.pages.metadata import MetadataPage
from PyPedal.desktop.settings import MAX_RECENT_FILES, DesktopSettings

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])


def _settings(tmp_path: Path) -> DesktopSettings:
    ini = tmp_path / "desktop.ini"
    return DesktopSettings(QSettings(str(ini), QSettings.Format.IniFormat))


def _window(qtbot: object, tmp_path: Path) -> MainWindow:
    window = MainWindow(_settings(tmp_path))
    qtbot.addWidget(window)
    return window


def test_open_dialog_uses_application_options_and_normalize(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    dialog = OpenPedigreeDialog(settings)
    pedigree = tmp_path / "dogs.ped"
    pedigree.write_text("1 0 0\n", encoding="utf-8")
    dialog.path_edit.setText(str(pedigree))
    dialog.format_edit.setText(" asdxb ")
    dialog.separator_edit.setText(", ")
    dialog.renumber_box.setChecked(True)
    options = dialog.selected_options()
    assert dialog.selected_path() == pedigree
    assert options.pedformat == "asdxb"
    assert options.separator == ","
    assert options.renumber is True
    assert options.messages == "quiet"
    assert options.pedigree_summary == 0


def test_open_dialog_empty_separator_is_a_space(tmp_path: Path) -> None:
    dialog = OpenPedigreeDialog(_settings(tmp_path))
    dialog.separator_edit.setText("")
    assert dialog.selected_options().separator == " "


def test_metadata_empty_state() -> None:
    page = MetadataPage()
    assert "No pedigree is open" in page._banner.text()
    assert page._banner.isHidden() is False


def test_animals_search_does_not_copy_rows() -> None:
    animals = [
        SimpleNamespace(
            originalID=i,
            animalID=i,
            sireID=0,
            damID=0,
            by=None,
            sex="u",
            name=f"dog{i}",
            fa=0.0,
        )
        for i in range(1, 6)
    ]
    page = AnimalsPage()
    page.model.set_source(PedigreeTableSource(SimpleNamespace(pedigree=animals)))
    page.search.setText("dog3")
    assert page.proxy.rowCount() == 5
    page.apply_filter_now()
    assert page.proxy.rowCount() == 1


def test_empty_window_disables_close_and_shows_empty_pages(qtbot: object, tmp_path: Path) -> None:
    window = _window(qtbot, tmp_path)
    assert window.close_action.isEnabled() is False
    assert window.open_action.isEnabled() is True
    assert "No pedigree is open" in window.metadata_page._banner.text()
    assert window.animals_page.model.rowCount() == 0
    assert window.status_file.text() == "No pedigree"
    assert window.status_operation.text() == "Ready"


def test_about_shows_version_author_maintainer_license(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: list[tuple[str, str]] = []

    def fake_about(_parent: object, title: str, text: str) -> None:
        captured.append((title, text))

    monkeypatch.setattr(QMessageBox, "about", fake_about)
    window = _window(qtbot, tmp_path)
    window.show_about()
    assert captured
    title, text = captured[0]
    assert title == "About PyPedal"
    assert "4.1.0" in text
    assert "John B. Cole" in text
    assert "Geert Degezelle" in text
    assert "LGPL-2.1-or-later" in text
    assert "GNU General Public License" not in text


def test_settings_persist_open_choices_in_isolated_ini(tmp_path: Path) -> None:
    ini = tmp_path / "desktop.ini"
    settings = DesktopSettings(QSettings(str(ini), QSettings.Format.IniFormat))
    path = tmp_path / "dogs.ped"
    path.write_text("1 0 0\n", encoding="utf-8")
    settings.remember_successful_open(path, "asdxb", ",", True)
    assert settings.last_pedformat() == "asdxb"
    assert settings.last_separator() == ","
    assert settings.last_renumber() is True
    assert settings.last_directory() == tmp_path
    assert settings.recent_files() == [path]
    text = ini.read_text(encoding="utf-8")
    assert "asdxb" in text
    home = str(Path.home() / "Library" / "Preferences")
    assert home not in str(ini)


def test_recent_files_policy_successful_only_most_recent_first(
    tmp_path: Path,
) -> None:
    """Successful loads only; most recent first; duplicates dropped; cap 10.

    Failed loads are never recorded. Options remembered are the last
    successful global pedformat/separator/renumber, not per-file memory.
    """
    settings = _settings(tmp_path)
    paths = []
    for index in range(12):
        path = tmp_path / f"dogs{index}.ped"
        path.write_text("1 0 0\n", encoding="utf-8")
        paths.append(path)
        settings.remember_successful_open(path, "asd", " ", True)
    settings.remember_successful_open(paths[11], "asdxb", ",", False)
    recent = settings.recent_files()
    assert len(recent) == MAX_RECENT_FILES
    assert recent[0] == paths[11]
    assert paths[0] not in recent
    assert paths[1] not in recent
    assert settings.last_pedformat() == "asdxb"
    assert settings.last_separator() == ","
    assert settings.last_renumber() is False


def test_missing_recent_file_is_reported_and_removed(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    notices: list[str] = []

    def fake_info(_parent: object, title: str, text: str) -> None:
        notices.append(f"{title}: {text}")

    monkeypatch.setattr(QMessageBox, "information", fake_info)
    settings = _settings(tmp_path)
    missing = tmp_path / "gone.ped"
    settings.remember_successful_open(missing, "asd", " ", True)
    window = MainWindow(settings)
    qtbot.addWidget(window)
    window.open_recent_path(missing)
    assert notices
    assert "no longer available" in notices[0]
    assert settings.recent_files() == []


def test_busy_disables_open_and_close_event_is_ignored(
    qtbot: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(QMessageBox, "information", lambda *args, **kwargs: None)
    window = _window(qtbot, tmp_path)
    window._set_busy(True, "Loading dogs.ped")
    assert window.open_action.isEnabled() is False
    assert window.recent_menu.isEnabled() is False
    assert window.close_action.isEnabled() is False
    event = QCloseEvent()
    window.closeEvent(event)
    assert event.isAccepted() is False
    assert window._busy is True
