"""Background pedigree load worker.

The worker runs ``load_into_session`` on a ``QThread``. It never touches
widgets. Progress is a queued signal so the GUI thread can update the bar.
"""

from __future__ import annotations

import traceback
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal, Slot

from PyPedal.application import PedigreeOpenOptions, PedigreeSession, load_into_session


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
