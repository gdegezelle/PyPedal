"""Background workers for pedigree load and scientific jobs.

Workers run on a ``QThread``. They never touch widgets. Progress is a
queued signal so the GUI thread can update the bar.

Load uses a scratch ``PedigreeSession`` so a failed open cannot replace
the displayed pedigree. Analysis jobs mutate the live session on the
worker thread (Meuwissen-Luo writes ``animal.fa`` in place). One
long-running job at a time.

Worker jobs: pedigree load, inbreeding, Lacy founders, relationship,
mating, pedigree save.

UI-thread jobs: inbreeding-by-year when a cache already exists,
theoretical Ne, CSV/text export, table refresh.
"""

from __future__ import annotations

import traceback
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal, Slot

from PyPedal.application import PedigreeOpenOptions, PedigreeSession, load_into_session

ProgressFn = Callable[[int, int | None], None]
WorkFn = Callable[[ProgressFn], object]


class LoadWorker(QObject):
    """Load a pedigree into a *scratch* session on a worker thread.

    The window's displayed session is not mutated here. On success the GUI
    thread adopts ``session``. A failed load never replaces the displayed
    pedigree because the scratch session is discarded.
    """

    progress = Signal(int, object)
    succeeded = Signal(object)
    failed = Signal(object, str)
    finished = Signal()

    def __init__(self, source: Path, options: PedigreeOpenOptions) -> None:
        super().__init__()
        self._source = source
        self._options = options
        self.session = PedigreeSession()

    def _emit_progress(self, done: int, total: int | None) -> None:
        self.progress.emit(done, total)

    @Slot()
    def run(self) -> None:
        try:
            load_into_session(
                self.session,
                self._source,
                self._options,
                progress=self._emit_progress,
            )
        except Exception as exc:
            self.failed.emit(exc, traceback.format_exc())
        else:
            self.succeeded.emit(self.session)
        finally:
            self.finished.emit()


class LoadJob:
    """Own one worker/thread pair and tear it down deterministically."""

    def __init__(self, worker: LoadWorker, thread: QThread) -> None:
        self.worker = worker
        self.thread = thread

    @classmethod
    def start(cls, source: Path, options: PedigreeOpenOptions) -> LoadJob:
        worker = LoadWorker(source, options)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        job = cls(worker, thread)
        thread.start()
        return job


class AnalysisWorker(QObject):
    """Run an application job on a worker thread."""

    progress = Signal(int, object)
    succeeded = Signal(object)
    failed = Signal(object, str)
    finished = Signal()

    def __init__(self, work: WorkFn) -> None:
        super().__init__()
        self._work = work

    def _emit_progress(self, done: int, total: int | None) -> None:
        self.progress.emit(done, total)

    @Slot()
    def run(self) -> None:
        try:
            result = self._work(self._emit_progress)
        except Exception as exc:
            self.failed.emit(exc, traceback.format_exc())
        else:
            self.succeeded.emit(result)
        finally:
            self.finished.emit()


class AnalysisJob:
    """Own one analysis worker/thread pair and tear it down deterministically."""

    def __init__(self, worker: AnalysisWorker, thread: QThread) -> None:
        self.worker = worker
        self.thread = thread

    @classmethod
    def start(cls, work: WorkFn) -> AnalysisJob:
        worker = AnalysisWorker(work)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        job = cls(worker, thread)
        thread.start()
        return job
