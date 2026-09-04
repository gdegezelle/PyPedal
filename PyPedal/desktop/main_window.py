"""Native QMainWindow shell for the PyPedal Qt desktop."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QFileDialog,
    QInputDialog,
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
from PyPedal.application import (
    BROWSE_COLUMNS,
    FoundersOutcome,
    InbreedingResult,
    MatingCoIGroupResult,
    PairwiseResult,
    PedigreeOpenOptions,
    PedigreeSession,
    PedigreeTableSource,
    YearInbreedingOutcome,
    export_inbreeding_csv,
    export_mating_group_csv,
    export_metadata_pdf,
    export_three_gen_pdf,
    export_year_inbreeding_csv,
    parse_animal_id,
    run_effective_founders,
    run_inbreeding,
    run_inbreeding_by_year,
    run_mating_coi,
    run_mating_coi_group,
    run_relationship,
    run_theoretical_ne,
    save_pedigree,
    write_text,
)
from PyPedal.desktop.dialogs.error import show_application_error
from PyPedal.desktop.dialogs.open_pedigree import OpenPedigreeDialog
from PyPedal.desktop.pages.animals import AnimalsPage
from PyPedal.desktop.pages.founders import FoundersPage
from PyPedal.desktop.pages.inbreeding import InbreedingPage
from PyPedal.desktop.pages.inbreeding_year import InbreedingYearPage
from PyPedal.desktop.pages.mating import MatingPage
from PyPedal.desktop.pages.metadata import MetadataPage
from PyPedal.desktop.pages.population import PopulationPage
from PyPedal.desktop.pages.relationship import RelationshipPage
from PyPedal.desktop.settings import DesktopSettings
from PyPedal.desktop.workers import AnalysisJob, LoadJob, WorkFn

PAGE_METADATA = 0
PAGE_ANIMALS = 1
PAGE_INBREEDING = 2
PAGE_YEAR = 3
PAGE_FOUNDERS = 4
PAGE_RELATIONSHIP = 5
PAGE_MATING = 6
PAGE_POPULATION = 7

_FA_COLUMN = next(index for index, column in enumerate(BROWSE_COLUMNS) if column.key == "fa")
_NAV_LABELS = (
    "Metadata",
    "Animals",
    "Inbreeding",
    "Inbreeding by Year",
    "Effective Founders",
    "Relationship",
    "Mating",
    "Population",
)


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
        self._job: LoadJob | AnalysisJob | None = None
        self._busy = False
        self._analysis_success: Callable[[object], None] | None = None

        self._build_pages()
        self._build_menus()
        self._build_status()
        self._connect_analysis()
        self._restore_geometry()
        self._refresh_recent_menu()
        self._sync_empty_state()

    def _build_pages(self) -> None:
        self.nav = QListWidget()
        self.nav.setObjectName("navigation")
        for label in _NAV_LABELS:
            self.nav.addItem(QListWidgetItem(label))
        self.nav.setCurrentRow(PAGE_METADATA)
        self.nav.setMaximumWidth(180)

        self.metadata_page = MetadataPage()
        self.animals_page = AnimalsPage()
        self.inbreeding_page = InbreedingPage()
        self.year_page = InbreedingYearPage()
        self.founders_page = FoundersPage()
        self.relationship_page = RelationshipPage()
        self.mating_page = MatingPage()
        self.population_page = PopulationPage()
        self.stack = QStackedWidget()
        for page in (
            self.metadata_page,
            self.animals_page,
            self.inbreeding_page,
            self.year_page,
            self.founders_page,
            self.relationship_page,
            self.mating_page,
            self.population_page,
        ):
            self.stack.addWidget(page)
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

        self.save_action = QAction("Save Pedigree As…", self)
        self.save_action.setObjectName("action_save")
        self.save_action.setShortcut(QKeySequence.StandardKey.SaveAs)
        self.save_action.triggered.connect(self.save_pedigree_as)
        file_menu.addAction(self.save_action)

        self.export_menu = file_menu.addMenu("Export")
        self.export_menu.setObjectName("menu_export")
        self.pdf_metadata_action = QAction("Metadata Report as PDF…", self)
        self.pdf_metadata_action.setObjectName("action_pdf_metadata")
        self.pdf_metadata_action.triggered.connect(self.export_metadata_pdf_report)
        self.export_menu.addAction(self.pdf_metadata_action)
        self.pdf_three_gen_action = QAction("Three-Generation Pedigree as PDF…", self)
        self.pdf_three_gen_action.setObjectName("action_pdf_three_gen")
        self.pdf_three_gen_action.triggered.connect(self.export_three_gen_pdf_report)
        self.export_menu.addAction(self.pdf_three_gen_action)

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

        analysis_menu = self.menuBar().addMenu("&Analysis")
        self.analysis_inbreeding_action = QAction("Inbreeding", self)
        self.analysis_inbreeding_action.triggered.connect(
            lambda: self.nav.setCurrentRow(PAGE_INBREEDING)
        )
        analysis_menu.addAction(self.analysis_inbreeding_action)
        self.analysis_year_action = QAction("Inbreeding by Year", self)
        self.analysis_year_action.triggered.connect(lambda: self.nav.setCurrentRow(PAGE_YEAR))
        analysis_menu.addAction(self.analysis_year_action)
        self.analysis_founders_action = QAction("Effective Founders", self)
        self.analysis_founders_action.triggered.connect(
            lambda: self.nav.setCurrentRow(PAGE_FOUNDERS)
        )
        analysis_menu.addAction(self.analysis_founders_action)
        self.analysis_relationship_action = QAction("Relationship", self)
        self.analysis_relationship_action.triggered.connect(
            lambda: self.nav.setCurrentRow(PAGE_RELATIONSHIP)
        )
        analysis_menu.addAction(self.analysis_relationship_action)
        self.analysis_mating_action = QAction("Mating", self)
        self.analysis_mating_action.triggered.connect(lambda: self.nav.setCurrentRow(PAGE_MATING))
        analysis_menu.addAction(self.analysis_mating_action)
        self.analysis_population_action = QAction("Population", self)
        self.analysis_population_action.triggered.connect(
            lambda: self.nav.setCurrentRow(PAGE_POPULATION)
        )
        analysis_menu.addAction(self.analysis_population_action)

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

    def _connect_analysis(self) -> None:
        self.inbreeding_page.run_requested.connect(self.run_inbreeding_analysis)
        self.inbreeding_page.export_requested.connect(self.export_inbreeding)
        self.year_page.run_requested.connect(self.run_year_analysis)
        self.year_page.export_requested.connect(self.export_year)
        self.founders_page.run_requested.connect(self.run_founders_analysis)
        self.founders_page.export_requested.connect(self.export_founders)
        self.relationship_page.run_requested.connect(self.run_relationship_analysis)
        self.relationship_page.export_requested.connect(self.export_relationship)
        self.mating_page.run_pair_requested.connect(self.run_mating_pair)
        self.mating_page.run_group_requested.connect(self.run_mating_group)
        self.mating_page.export_requested.connect(self.export_mating)
        self.population_page.run_requested.connect(self.run_theoretical_ne_analysis)
        self.population_page.export_requested.connect(self.export_ne)

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
        """Open a recent file using the last successful load options."""
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
        job.worker.failed.connect(self._on_job_failed)
        job.thread.finished.connect(self._on_job_finished)
        self._job = job

    def _start_analysis(
        self,
        label: str,
        work: WorkFn,
        on_success: Callable[[object], None],
    ) -> None:
        if self._busy:
            return
        self._analysis_success = on_success
        self._set_busy(True, label)
        job = AnalysisJob.start(work)
        job.worker.progress.connect(self._on_progress)
        job.worker.succeeded.connect(self._on_analysis_succeeded)
        job.worker.failed.connect(self._on_job_failed)
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
        self._show_loaded_pedigree()

    def _on_analysis_succeeded(self, result: object) -> None:
        handler = self._analysis_success
        self._analysis_success = None
        if handler is not None:
            handler(result)

    def _on_job_failed(self, exc: object, details: str) -> None:
        if isinstance(exc, BaseException):
            show_application_error(self, exc, details)
        self._sync_empty_state()

    def _on_job_finished(self) -> None:
        self._job = None
        self._analysis_success = None
        self._set_busy(False, "Ready")

    def _analysis_buttons(self) -> tuple:
        return (
            self.inbreeding_page.run_button,
            self.year_page.run_button,
            self.founders_page.run_button,
            self.population_page.run_button,
        )

    def _sync_selector_pages(self, *, armed: bool) -> None:
        self.relationship_page.set_armed(armed)
        self.mating_page.set_armed(armed)

    def _set_busy(self, busy: bool, operation: str) -> None:
        self._busy = busy
        empty = self.session.is_empty
        self.open_action.setEnabled(not busy)
        self.recent_menu.setEnabled(not busy)
        self.save_action.setEnabled(not busy and not empty)
        self.pdf_metadata_action.setEnabled(not busy and not empty)
        self.pdf_three_gen_action.setEnabled(not busy and not empty)
        self.close_action.setEnabled(not busy and not empty)
        for button in self._analysis_buttons():
            button.setEnabled(not busy and not empty)
        self._sync_selector_pages(armed=not busy and not empty)
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
        self.animals_page.set_source(None)
        self.animals_page.search.clear()
        self.animals_page.apply_filter_now()
        self.metadata_page.show_empty()
        self.relationship_page.set_lookup(None)
        self.mating_page.set_lookup(None)
        self._clear_analysis_pages()
        self._sync_empty_state()

    def _clear_analysis_pages(self) -> None:
        self.inbreeding_page.show_empty()
        self.year_page.show_empty()
        self.founders_page.show_empty()
        self.relationship_page.show_empty()
        self.mating_page.show_empty()
        self.population_page.show_empty()

    def _show_loaded_pedigree(self) -> None:
        pedigree = self.session.pedigree
        if pedigree is None:
            return
        self.animals_page.set_source(PedigreeTableSource(pedigree))
        self.metadata_page.show_session(self.session)
        self.relationship_page.set_lookup(self.session.animal_lookup)
        self.mating_page.set_lookup(self.session.animal_lookup)
        self._clear_analysis_pages()
        self._sync_empty_state()

    def _refresh_f_column(self) -> None:
        self.animals_page.model.refresh_column(_FA_COLUMN)

    def _sync_empty_state(self) -> None:
        empty = self.session.is_empty
        self.close_action.setEnabled(not empty and not self._busy)
        self.save_action.setEnabled(not empty and not self._busy)
        self.pdf_metadata_action.setEnabled(not empty and not self._busy)
        self.pdf_three_gen_action.setEnabled(not empty and not self._busy)
        for button in self._analysis_buttons():
            button.setEnabled(not empty and not self._busy)
        self._sync_selector_pages(armed=not empty and not self._busy)
        if empty:
            self.status_file.setText("No pedigree")
            self.status_count.setText("")
            self.setWindowFilePath("")
            return
        source = self.session.source_path
        name = source.name if source is not None else "pedigree"
        self.status_file.setText(name)
        self.setWindowFilePath(str(source) if source is not None else "")
        pedigree = self.session.pedigree
        count = len(pedigree.pedigree) if pedigree is not None else 0
        self.status_count.setText(f"{count:,} animals")

    def run_inbreeding_analysis(self) -> None:
        session = self.session
        self._start_analysis(
            "Calculating inbreeding",
            lambda progress: run_inbreeding(session, progress=progress),
            self._apply_inbreeding_result,
        )

    def _apply_inbreeding_result(self, result: object) -> None:
        if not isinstance(result, InbreedingResult):
            return
        self.inbreeding_page.show_result(result)
        self._refresh_f_column()

    def run_year_analysis(self) -> None:
        session = self.session
        if session.inbreeding_result is not None:
            try:
                outcome = run_inbreeding_by_year(session)
            except Exception as exc:
                show_application_error(self, exc, "")
                return
            self._apply_year_outcome(outcome)
            return
        self._start_analysis(
            "Calculating inbreeding",
            lambda progress: run_inbreeding_by_year(session, progress=progress),
            self._apply_year_outcome,
        )

    def _apply_year_outcome(self, outcome: object) -> None:
        if not isinstance(outcome, YearInbreedingOutcome):
            return
        if outcome.computed_inbreeding and self.session.inbreeding_result is not None:
            self.inbreeding_page.show_result(self.session.inbreeding_result)
            self._refresh_f_column()
        self.year_page.show_rows(outcome.rows, computed_inbreeding=outcome.computed_inbreeding)

    def run_founders_analysis(self) -> None:
        session = self.session
        self._start_analysis(
            "Calculating effective founders",
            lambda _progress: run_effective_founders(session),
            self._apply_founders_outcome,
        )

    def _apply_founders_outcome(self, outcome: object) -> None:
        if not isinstance(outcome, FoundersOutcome):
            return
        self.founders_page.show_outcome(outcome)
        if outcome.implicit_renumber:
            QMessageBox.information(
                self,
                "Pedigree renumbered",
                "Lacy effective founders automatically renumbered this pedigree. "
                "Cached inbreeding and pairwise results were cleared.",
            )
            self._show_loaded_pedigree()
            self.founders_page.show_outcome(outcome)

    def run_relationship_analysis(self) -> None:
        animal_a = self.relationship_page.selected_animal_a()
        animal_b = self.relationship_page.selected_animal_b()
        if animal_a is None or animal_b is None:
            return
        session = self.session
        self._start_analysis(
            "Computing relationship",
            lambda _progress: run_relationship(session, animal_a, animal_b),
            self._apply_relationship_result,
        )

    def _apply_relationship_result(self, result: object) -> None:
        if isinstance(result, PairwiseResult):
            self.relationship_page.show_result(result)

    def run_mating_pair(self) -> None:
        animal_a = self.mating_page.selected_animal_a()
        animal_b = self.mating_page.selected_animal_b()
        if animal_a is None or animal_b is None:
            return
        session = self.session
        self._start_analysis(
            "Computing mating CoI",
            lambda _progress: run_mating_coi(session, animal_a, animal_b),
            self._apply_mating_pair_result,
        )

    def _apply_mating_pair_result(self, result: object) -> None:
        if isinstance(result, PairwiseResult):
            self.mating_page.show_pair(result)

    def run_mating_group(self) -> None:
        try:
            pairs = [
                (
                    parse_animal_id(str(left), label="Animal A"),
                    parse_animal_id(str(right), label="Animal B"),
                )
                for left, right in self.mating_page.group_pairs()
            ]
        except Exception as exc:
            show_application_error(self, exc, "")
            return
        session = self.session
        self._start_analysis(
            "Computing mating group",
            lambda _progress: run_mating_coi_group(session, pairs),
            self._apply_mating_group_result,
        )

    def _apply_mating_group_result(self, result: object) -> None:
        if isinstance(result, MatingCoIGroupResult):
            self.mating_page.show_group(result)

    def run_theoretical_ne_analysis(self) -> None:
        try:
            value = run_theoretical_ne(self.session)
        except Exception as exc:
            show_application_error(self, exc, "")
            return
        self.population_page.show_value(value)

    def save_pedigree_as(self) -> None:
        if self._busy or self.session.is_empty:
            return
        path, _selected = QFileDialog.getSaveFileName(
            self,
            "Save Pedigree As",
            str(self.session.source_path or Path.home()),
            "Pedigree files (*.ped);;All files (*)",
        )
        if not path:
            return
        destination = Path(path)
        session = self.session
        self._start_analysis(
            f"Saving {destination.name}",
            lambda _progress: save_pedigree(session, destination, overwrite=True),
            lambda _result: None,
        )

    def _choose_pdf_path(self, title: str, suggested: str) -> Path | None:
        path, _selected = QFileDialog.getSaveFileName(
            self,
            title,
            suggested,
            "PDF files (*.pdf);;All files (*)",
        )
        if not path:
            return None
        return Path(path)

    def export_metadata_pdf_report(self) -> None:
        if self._busy or self.session.is_empty:
            return
        suggested = "metadata.pdf"
        source = self.session.source_path
        if source is not None:
            suggested = str(source.with_name(f"{source.stem}_metadata.pdf"))
        path = self._choose_pdf_path("Export metadata report", suggested)
        if path is None:
            return
        session = self.session
        self._start_analysis(
            f"Writing {path.name}",
            lambda _progress: export_metadata_pdf(session, path, overwrite=True),
            lambda _result: None,
        )

    def export_three_gen_pdf_report(self) -> None:
        if self._busy or self.session.is_empty:
            return
        suggested_id = self.animals_page.selected_animal_id()
        text, accepted = QInputDialog.getText(
            self,
            "Three-generation pedigree",
            "Current animal ID:",
            text="" if suggested_id is None else str(suggested_id),
        )
        if not accepted:
            return
        try:
            animal_id = parse_animal_id(text, label="Animal ID")
        except Exception as exc:
            show_application_error(self, exc, "")
            return
        suggested = f"three_generation_{animal_id}.pdf"
        path = self._choose_pdf_path("Export three-generation pedigree", suggested)
        if path is None:
            return
        session = self.session
        self._start_analysis(
            f"Writing {path.name}",
            lambda _progress: export_three_gen_pdf(session, animal_id, path, overwrite=True),
            lambda _result: None,
        )

    def _choose_export_path(self, title: str, suggested: str) -> Path | None:
        path, _selected = QFileDialog.getSaveFileName(
            self,
            title,
            suggested,
            "CSV files (*.csv);;Text files (*.txt);;All files (*)",
        )
        if not path:
            return None
        return Path(path)

    def export_inbreeding(self) -> None:
        result = self.inbreeding_page.result
        if result is None:
            return
        path = self._choose_export_path("Export inbreeding", "inbreeding.csv")
        if path is None:
            return
        try:
            export_inbreeding_csv(path, result, overwrite=True)
        except Exception as exc:
            show_application_error(self, exc, "")

    def export_year(self) -> None:
        if not self.year_page.rows:
            return
        path = self._choose_export_path("Export inbreeding by year", "inbreeding_by_year.csv")
        if path is None:
            return
        try:
            export_year_inbreeding_csv(path, self.year_page.rows, overwrite=True)
        except Exception as exc:
            show_application_error(self, exc, "")

    def export_founders(self) -> None:
        result = self.founders_page.result
        if result is None:
            return
        path = self._choose_export_path("Export effective founders", "effective_founders.txt")
        if path is None:
            return
        text = (
            f"Effective founders: {result.fa_effective_founders}\n"
            f"Animals: {result.fa_animal_count}\n"
            f"Founders: {result.fa_founder_count}\n"
            f"Descendants: {result.fa_descendant_count}\n"
        )
        try:
            write_text(path, text, overwrite=True)
        except Exception as exc:
            show_application_error(self, exc, "")

    def export_relationship(self) -> None:
        result = self.relationship_page.result
        if result is None:
            return
        path = self._choose_export_path("Export relationship", "relationship.txt")
        if path is None:
            return
        text = (
            f"Animal A: {result.animal_a}\n"
            f"Animal B: {result.animal_b}\n"
            f"Relationship: {result.coefficient}\n"
        )
        try:
            write_text(path, text, overwrite=True)
        except Exception as exc:
            show_application_error(self, exc, "")

    def export_mating(self) -> None:
        group = self.mating_page.group_result
        pair = self.mating_page.pair_result
        if group is not None:
            path = self._choose_export_path("Export mating group", "mating_group.csv")
            if path is None:
                return
            try:
                export_mating_group_csv(path, group, overwrite=True)
            except Exception as exc:
                show_application_error(self, exc, "")
            return
        if pair is None:
            return
        path = self._choose_export_path("Export mating", "mating.txt")
        if path is None:
            return
        text = (
            f"Animal A: {pair.animal_a}\n"
            f"Animal B: {pair.animal_b}\n"
            f"Offspring F: {pair.coefficient}\n"
        )
        try:
            write_text(path, text, overwrite=True)
        except Exception as exc:
            show_application_error(self, exc, "")

    def export_ne(self) -> None:
        value = self.population_page.value
        if value is None:
            return
        path = self._choose_export_path("Export theoretical Ne", "theoretical_ne.txt")
        if path is None:
            return
        try:
            write_text(
                path,
                f"Theoretical Ne from metadata: {value}\n",
                overwrite=True,
            )
        except Exception as exc:
            show_application_error(self, exc, "")

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
                "Operation in progress",
                "Wait for the current operation to finish before quitting.",
            )
            event.ignore()
            return
        self.settings.set_window_geometry(self.saveGeometry())
        self.settings.set_window_state(self.saveState())
        event.accept()
