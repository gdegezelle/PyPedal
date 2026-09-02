"""Native QMainWindow shell for the PyPedal Qt desktop."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QSplitter,
    QStackedWidget,
    QWidget,
)

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.application import PedigreeOpenOptions, PedigreeSession, PedigreeTableSource
from PyPedal.desktop.dialogs.error import show_load_error
from PyPedal.desktop.dialogs.open_pedigree import OpenPedigreeDialog
from PyPedal.desktop.pages.animals import AnimalsPage
from PyPedal.desktop.pages.metadata import MetadataPage
from PyPedal.desktop.settings import DesktopSettings
from PyPedal.desktop.workers import LoadJob

PAGE_METADATA = 0
PAGE_ANIMALS = 1


class MainWindow(QMainWindow):
    """Restrained scientific desktop shell."""

    def __init__(
        self,
        settings: DesktopSettings | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"PyPedal {PYPEDAL_VERSION}")
        self.resize(960, 640)
        self.settings = settings or DesktopSettings()
        self.session = PedigreeSession()
        self._job: LoadJob | None = None
        self._busy = False

        self._build_pages()
        self._build_menus()
        self._build_status()
        self._restore_geometry()
        self._refresh_recent_menu()
        self._sync_empty_state()

    def _build_pages(self) -> None:
        self.nav = QListWidget()
        self.nav.setObjectName("navigation")
        self.nav.addItem(QListWidgetItem("Metadata"))
        self.nav.addItem(QListWidgetItem("Animals"))
        self.nav.setCurrentRow(PAGE_METADATA)
        self.nav.setMaximumWidth(160)

        self.metadata_page = MetadataPage()
        self.animals_page = AnimalsPage()
        self.stack = QStackedWidget()
        self.stack.addWidget(self.metadata_page)
        self.stack.addWidget(self.animals_page)
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)

        splitter = QSplitter()
        splitter.addWidget(self.nav)
        splitter.addWidget(self.stack)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)

    def _build_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        self.open_action = QAction("Open…", self)
        self.open_action.setObjectName("action_open")
        self.open_action.setShortcut(QKeySequence.StandardKey.Open)
        self.open_action.triggered.connect(self.open_pedigree)
        file_menu.addAction(self.open_action)

        self.recent_menu = file_menu.addMenu("Open Recent")
        self.recent_menu.setObjectName("menu_recent")

        self.close_action = QAction("Close Pedigree", self)
        self.close_action.setObjectName("action_close")
        self.close_action.triggered.connect(self.close_pedigree)
        file_menu.addAction(self.close_action)
        file_menu.addSeparator()

        self.quit_action = QAction("Quit", self)
        self.quit_action.setObjectName("action_quit")
        self.quit_action.setMenuRole(QAction.MenuRole.QuitRole)
        self.quit_action.setShortcut(QKeySequence.StandardKey.Quit)
        self.quit_action.triggered.connect(self.close)
        file_menu.addAction(self.quit_action)

        help_menu = self.menuBar().addMenu("&Help")
        self.about_action = QAction("About PyPedal", self)
        self.about_action.setObjectName("action_about")
        self.about_action.setMenuRole(QAction.MenuRole.AboutRole)
        self.about_action.triggered.connect(self.show_about)
        help_menu.addAction(self.about_action)

    def _build_status(self) -> None:
        self.status_file = QLabel("No pedigree")
        self.status_file.setObjectName("status_file")
        self.status_count = QLabel("")
        self.status_count.setObjectName("status_count")
        self.status_operation = QLabel("Ready")
        self.status_operation.setObjectName("status_operation")
        self.progress = QProgressBar()
        self.progress.setObjectName("status_progress")
        self.progress.setMaximumWidth(160)
        self.progress.setTextVisible(False)
        self.progress.hide()
        status = self.statusBar()
        status.addWidget(self.status_file, 1)
        status.addWidget(self.status_count)
        status.addPermanentWidget(self.status_operation)
        status.addPermanentWidget(self.progress)

    def _restore_geometry(self) -> None:
        geometry = self.settings.window_geometry()
        if geometry is not None:
            self.restoreGeometry(geometry)
        state = self.settings.window_state()
        if state is not None:
            self.restoreState(state)

    def _refresh_recent_menu(self) -> None:
        self.recent_menu.clear()
        recent = self.settings.recent_files()
        if not recent:
            empty = QAction("(None)", self)
            empty.setEnabled(False)
            self.recent_menu.addAction(empty)
            return
        for path in recent:
            action = QAction(path.name, self)
            action.setData(str(path))
            action.setToolTip(str(path))
            action.triggered.connect(self._open_recent_action)
            self.recent_menu.addAction(action)

    def _open_recent_action(self) -> None:
        action = self.sender()
        if not isinstance(action, QAction):
            return
        raw = action.data()
        if not raw:
            return
        self.open_recent_path(Path(str(raw)))

    def open_recent_path(self, path: Path) -> None:
        """Open a recent file using the last successful load options.

        Failed paths are never stored. A missing file is reported and
        removed from the recent list. Options are the last successful
        pedformat/separator/renumber, not a per-file memory.
        """
        if not path.is_file():
            QMessageBox.information(
                self,
                "Recent file missing",
                f"{path.name} is no longer available and will be removed from Open Recent.",
            )
            self.settings.remove_recent_file(path)
            self._refresh_recent_menu()
            return
        options = PedigreeOpenOptions(
            pedformat=self.settings.last_pedformat(),
            separator=self.settings.last_separator(),
            renumber=self.settings.last_renumber(),
            messages="quiet",
            pedigree_summary=0,
        ).normalized()
        self._start_load(path, options)

    def open_pedigree(self) -> None:
        if self._busy:
            return
        dialog = OpenPedigreeDialog(self.settings, self)
        if dialog.exec() != OpenPedigreeDialog.DialogCode.Accepted:
            return
        path = dialog.selected_path()
        if path is None:
            return
        self._start_load(path, dialog.selected_options())

    def open_path(
        self,
        path: Path,
        options: PedigreeOpenOptions | None = None,
    ) -> None:
        """Open ``path`` with stored or supplied options (startup / tests)."""
        chosen = (
            options
            or PedigreeOpenOptions(
                pedformat=self.settings.last_pedformat(),
                separator=self.settings.last_separator(),
                renumber=self.settings.last_renumber(),
                messages="quiet",
                pedigree_summary=0,
            ).normalized()
        )
        self._start_load(path, chosen)

    def _start_load(self, path: Path, options: PedigreeOpenOptions) -> None:
        if self._busy:
            return
        self._set_busy(True, f"Loading {path.name}")
        job = LoadJob.start(path, options)
        job.worker.progress.connect(self._on_progress)
        job.worker.succeeded.connect(self._on_load_succeeded)
        job.worker.failed.connect(self._on_load_failed)
        job.thread.finished.connect(self._on_job_finished)
        self._job = job

    def _on_progress(self, done: int, total: object) -> None:
        if not isinstance(total, int) or total <= 0:
            self.progress.setRange(0, 0)
            return
        self.progress.setRange(0, total)
        self.progress.setValue(done)

    def _on_load_succeeded(self, session: object) -> None:
        if not isinstance(session, PedigreeSession):
            return
        self.session = session
        pedigree = session.pedigree
        source = session.source_path
        options = session.load_options
        if pedigree is None or source is None or options is None:
            return
        separator = options.separator if options.separator is not None else " "
        self.settings.remember_successful_open(
            source,
            options.pedformat,
            separator,
            options.renumber,
        )
        self._refresh_recent_menu()
        self.animals_page.model.set_source(PedigreeTableSource(pedigree))
        self.metadata_page.show_session(self.session)
        self._sync_empty_state()

    def _on_load_failed(self, exc: object, details: str) -> None:
        if isinstance(exc, BaseException):
            show_load_error(self, exc, details)
        self._sync_empty_state()

    def _on_job_finished(self) -> None:
        self._job = None
        self._set_busy(False, "Ready")

    def _set_busy(self, busy: bool, operation: str) -> None:
        self._busy = busy
        self.open_action.setEnabled(not busy)
        self.recent_menu.setEnabled(not busy)
        self.close_action.setEnabled(not busy and not self.session.is_empty)
        self.status_operation.setText(operation)
        if busy:
            self.progress.setRange(0, 0)
            self.progress.show()
        else:
            self.progress.hide()
            self.progress.setRange(0, 1)
            self.progress.setValue(0)

    def close_pedigree(self) -> None:
        if self._busy or self.session.is_empty:
            return
        self.session.clear()
        self.animals_page.model.set_source(None)
        self.animals_page.search.clear()
        self.animals_page.apply_filter_now()
        self.metadata_page.show_empty()
        self._sync_empty_state()

    def _sync_empty_state(self) -> None:
        empty = self.session.is_empty
        self.close_action.setEnabled(not empty and not self._busy)
        if empty:
            self.status_file.setText("No pedigree")
            self.status_count.setText("")
            return
        source = self.session.source_path
        name = source.name if source is not None else "pedigree"
        self.status_file.setText(name)
        pedigree = self.session.pedigree
        count = len(pedigree.pedigree) if pedigree is not None else 0
        self.status_count.setText(f"{count:,} animals")

    def show_about(self) -> None:
        QMessageBox.about(
            self,
            "About PyPedal",
            (
                f"<b>PyPedal</b><br>"
                f"Version {PYPEDAL_VERSION}<br><br>"
                "Original author: John B. Cole, PhD<br>"
                "Maintainer: Geert Degezelle<br><br>"
                "License: GNU LGPL-2.1-or-later"
            ),
        )

    def closeEvent(self, event: QCloseEvent) -> None:
        if self._busy:
            QMessageBox.information(
                self,
                "Load in progress",
                "A pedigree is still loading. Wait for it to finish before quitting.",
            )
            event.ignore()
            return
        self.settings.set_window_geometry(self.saveGeometry())
        self.settings.set_window_state(self.saveState())
        event.accept()
