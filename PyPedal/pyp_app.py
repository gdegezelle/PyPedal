#!/usr/bin/env python3
"""
Cross-platform desktop GUI for PyPedal.

This replaces the old wxPython/Wax interface. CustomTkinter runs on
macOS, Windows, and Linux and does not require a browser or a web stack.

Written for PyPedal 4.0. SPDX-License-Identifier: LGPL-2.1-or-later.
"""

from __future__ import annotations

import os
import sys
import threading
import traceback
from typing import Any, Optional

from PyPedal import pyp_errors, pyp_io, pyp_metrics, pyp_newclasses, pyp_nrm
from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.application.errors import EXIT_STATUS as EXIT_STATUS
from PyPedal.application.errors import exit_status_for as exit_status_for
from PyPedal.application.load import normalize_sepchar as normalize_sepchar


def _require_customtkinter():
    try:
        import customtkinter as ctk
        from tkinter import filedialog, messagebox
    except ImportError as exc:
        raise SystemExit(
            "The PyPedal GUI needs CustomTkinter.\n"
            "Install it with:  pip install 'PyPedal[gui]'\n"
            f"Original error: {exc}"
        ) from exc
    return ctk, filedialog, messagebox


def pedigree_open_options(
    path: str,
    pedformat: str,
    sepchar: str,
    renumber: bool,
) -> dict:
    """Options dict the desktop app passes to ``NewPedigree``.

    Separator normalisation lives in ``PyPedal.application``. This wrapper
    keeps the CustomTkinter dict shape, including omitting
    ``pedigree_summary`` so the library default (1) is unchanged. The
    application-layer ``PedigreeOpenOptions`` default is
    ``pedigree_summary=0`` for future PySide6 loads.
    """
    return {
        "pedfile": path,
        "pedformat": (pedformat or "").strip() or "asd",
        "sepchar": normalize_sepchar(sepchar),
        "renumber": bool(renumber),
        "messages": "quiet",
        "pedname": os.path.basename(path),
    }


GUI_PREVIEW_ROWS = 500
GUI_PROGRESS_POLL_MS = 80


def gui_progress_mode(done: int, total: int | None) -> tuple[str, float | None]:
    """Map a scientific progress event to a progress-bar mode.

    Known ``total`` becomes determinate. Unknown ``total`` stays
    indeterminate. This helper does not touch Tk widgets.
    """
    if total is None or total <= 0:
        return "indeterminate", None
    fraction = float(done) / float(total)
    if fraction < 0.0:
        fraction = 0.0
    elif fraction > 1.0:
        fraction = 1.0
    return "determinate", fraction


class GuiProgressBridge:
    """Publish the latest ``(done, total)`` from a worker thread.

    The scientific callback only stores state. Tk widgets must be updated
    on the UI thread by reading ``latest``. Rapid events coalesce: the UI
    renders the most recent snapshot, not every callback.
    """

    def __init__(self) -> None:
        self.latest: tuple[int, int | None] | None = None

    def __call__(self, done: int, total: int | None) -> None:
        self.latest = (done, total)


def gui_control_states(busy: bool) -> dict:
    """Enabled/disabled intent for GUI controls while work is running.

    About stays available. Open and analysis actions do not.
    """
    return {
        "open": not busy,
        "analyses": not busy,
        "about": True,
    }


def format_preview_caption(shown: int, total: int) -> str:
    return f"Showing {shown:,} of {total:,}"


def _format_inbreeding(
    result: Any,
    max_rows: int = GUI_PREVIEW_ROWS,
    result_file: Optional[str] = None,
) -> str:
    if isinstance(result, tuple):
        result = result[0]
    if not isinstance(result, dict):
        return str(result)

    metadata = result.get("metadata")
    if isinstance(metadata, dict) and "all" not in metadata and "f_count" in metadata:
        metadata = {"all": metadata, "nonzero": {}}
    lines = []
    if isinstance(metadata, dict):
        lines.append(pyp_io.summary_inbreeding(metadata))
    else:
        lines.append(str(result))

    fx = result.get("fx") or {}
    if fx:
        items = sorted(fx.items(), key=lambda item: item[0])
        total = len(items)
        shown = items[:max_rows]
        lines.append("Coefficients by animal")
        lines.append("-" * 40)
        if total > max_rows:
            lines.append(format_preview_caption(len(shown), total))
            if result_file:
                lines.append(f"Full coefficients are in {os.path.basename(result_file)}")
        for animal_id, coef in shown:
            lines.append(f"  {animal_id}: {pyp_io.format_display_coefficient(coef)}")
        if all(float(coef) == 0.0 for coef in fx.values()):
            lines.append("")
            lines.append(
                "No inbreeding in this pedigree (every coefficient is 0)."
            )
    return "\n".join(lines)


def _list_animals(pedigree, max_rows: int = GUI_PREVIEW_ROWS) -> str:
    animals = list(pedigree.pedigree)
    total = len(animals)
    shown = animals[:max_rows]
    lines = [
        f"{'ID':>8}  {'Sire':>8}  {'Dam':>8}  {'Year':>6}  {'Sex':<4}  Name",
        "-" * 64,
    ]
    if total > max_rows:
        lines.insert(0, format_preview_caption(len(shown), total))
    for animal in shown:
        year = getattr(animal, "by", None)
        year_label = "" if year is None else year
        sex = getattr(animal, "sex", "")
        name = getattr(animal, "name", "")
        lines.append(
            f"{animal.animalID:>8}  {animal.sireID:>8}  {animal.damID:>8}  "
            f"{year_label:>6}  {str(sex):<4}  {name}"
        )
    return "\n".join(lines)


def apply_pedigree_load(session, attempted_path, pedigree=None, error=None):
    """Install a loaded pedigree, or leave the previous one active.

    Call this on the UI thread. A worker must not assign ``session.pedigree``
    or ``session.filename``. On failure the previous pedigree is left in
    place; if there was none, the session stays empty.
    """
    attempted = os.path.basename(attempted_path) or attempted_path
    if error is None and pedigree is not None:
        session.pedigree = pedigree
        session.filename = attempted_path
        return {
            "ok": True,
            "output": None,
            "status": f"Loaded {attempted}",
        }
    if session.filename:
        active = os.path.basename(session.filename)
        note = f"{active} remains the active pedigree."
        status = f"Load of {attempted} failed — {active} remains active"
    else:
        note = "No pedigree is active."
        status = f"Load of {attempted} failed — no pedigree is active"
    body = error if error else "Unknown error"
    return {
        "ok": False,
        "output": f"Load of {attempted} failed.\n{note}\n\n{body}",
        "status": status,
    }


class PyPedalApp:
    """Main window: load a pedigree and run the analyses the old GUI offered."""

    def __init__(self) -> None:
        ctk, filedialog, messagebox = _require_customtkinter()
        self.ctk = ctk
        self.filedialog = filedialog
        self.messagebox = messagebox

        ctk.set_appearance_mode("system")
        ctk.set_default_color_theme("blue")

        self.root = ctk.CTk()
        self.root.title(f"PyPedal {PYPEDAL_VERSION}")
        self.root.geometry("960x640")
        self.root.minsize(720, 480)

        self.pedigree = None
        self.filename = ""
        self._busy = False
        self._gui_progress = GuiProgressBridge()
        self._progress_poll_id = None
        self._progress_bar_mode = "indeterminate"

        self._build_layout()

    def _build_layout(self) -> None:
        ctk = self.ctk

        header = ctk.CTkFrame(self.root, fg_color="transparent")
        header.pack(fill="x", padx=16, pady=(16, 8))

        ctk.CTkLabel(
            header,
            text="PyPedal",
            font=ctk.CTkFont(size=22, weight="bold"),
        ).pack(side="left")
        ctk.CTkLabel(
            header,
            text="Pedigree analysis",
            font=ctk.CTkFont(size=14),
        ).pack(side="left", padx=(10, 0))

        controls = ctk.CTkFrame(self.root)
        controls.pack(fill="x", padx=16, pady=8)

        self.open_button = ctk.CTkButton(
            controls, text="Open pedigree…", command=self.open_pedigree, width=140
        )
        self.open_button.pack(side="left", padx=8, pady=10)

        ctk.CTkLabel(controls, text="Format").pack(side="left", padx=(8, 4))
        self.format_var = ctk.StringVar(value="asdxb")
        ctk.CTkEntry(controls, textvariable=self.format_var, width=90).pack(side="left")

        ctk.CTkLabel(controls, text="Separator").pack(side="left", padx=(12, 4))
        self.sep_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            controls,
            textvariable=self.sep_var,
            width=48,
            placeholder_text="space",
        ).pack(side="left")

        self.renumber_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(controls, text="Renumber", variable=self.renumber_var).pack(
            side="left", padx=12
        )

        body = ctk.CTkFrame(self.root, fg_color="transparent")
        body.pack(fill="both", expand=True, padx=16, pady=8)

        sidebar = ctk.CTkFrame(body, width=200)
        sidebar.pack(side="left", fill="y", padx=(0, 12))
        sidebar.pack_propagate(False)

        actions = [
            ("Metadata", self.show_metadata),
            ("List animals", self.list_animals),
            ("Inbreeding", self.calc_inbreeding),
            ("Effective founders", self.calc_founders),
            ("Inbreeding by year", self.inbreeding_by_year),
        ]
        self.action_buttons = []
        for label, command in actions:
            button = ctk.CTkButton(sidebar, text=label, command=command, anchor="w")
            button.pack(fill="x", padx=10, pady=6)
            self.action_buttons.append(button)

        self.about_button = ctk.CTkButton(
            sidebar, text="About", command=self.show_about, fg_color="gray"
        )
        self.about_button.pack(fill="x", padx=10, pady=(18, 6))

        self.output = ctk.CTkTextbox(body, font=ctk.CTkFont(family="Menlo", size=13))
        self.output.pack(side="left", fill="both", expand=True)
        self._write(
            "Open a .ped file to get started.\n\n"
            "Format codes describe pedigree columns, for example:\n"
            "  asd    animal, sire, dam\n"
            "  asdxb  animal, sire, dam, sex, birth year\n"
            "  asdgy  animal, sire, dam, generation, birth year\n"
        )

        self.status = ctk.CTkLabel(self.root, text="No pedigree loaded", anchor="w")
        self.status.pack(fill="x", padx=16, pady=(0, 12))
        self.progress = ctk.CTkProgressBar(self.root, mode="indeterminate")

    def _write(self, text: str) -> None:
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)

    def _set_status(self, text: str) -> None:
        self.status.configure(text=text)

    def _set_busy(self, busy: bool) -> None:
        self._busy = busy
        states = gui_control_states(busy)
        open_state = "normal" if states["open"] else "disabled"
        analysis_state = "normal" if states["analyses"] else "disabled"
        about_state = "normal" if states["about"] else "disabled"
        self.open_button.configure(state=open_state)
        for button in self.action_buttons:
            button.configure(state=analysis_state)
        self.about_button.configure(state=about_state)
        if busy:
            self._progress_bar_mode = "indeterminate"
            self.progress.configure(mode="indeterminate")
            self.progress.set(0)
            self.progress.pack(fill="x", padx=16, pady=(0, 8), before=self.status)
            self.progress.start()
            self._schedule_progress_poll()
        else:
            self._cancel_progress_poll()
            self._apply_progress_snapshot()
            self.progress.stop()
            self._progress_bar_mode = "indeterminate"
            self.progress.configure(mode="indeterminate")
            self.progress.pack_forget()

    def _schedule_progress_poll(self) -> None:
        if self._progress_poll_id is not None:
            return
        self._progress_poll_id = self.root.after(
            GUI_PROGRESS_POLL_MS, self._poll_progress
        )

    def _cancel_progress_poll(self) -> None:
        poll_id = self._progress_poll_id
        if poll_id is not None:
            self.root.after_cancel(poll_id)
            self._progress_poll_id = None

    def _poll_progress(self) -> None:
        self._progress_poll_id = None
        self._apply_progress_snapshot()
        if self._busy:
            self._schedule_progress_poll()

    def _apply_progress_snapshot(self) -> None:
        snapshot = self._gui_progress.latest
        if snapshot is None:
            return
        done, total = snapshot
        mode, fraction = gui_progress_mode(done, total)
        if mode == "determinate":
            if self._progress_bar_mode != "determinate":
                self.progress.stop()
                self.progress.configure(mode="determinate")
                self._progress_bar_mode = "determinate"
            self.progress.set(fraction)
        elif self._progress_bar_mode != "indeterminate":
            self.progress.configure(mode="indeterminate")
            self._progress_bar_mode = "indeterminate"
            self.progress.start()

    def _need_pedigree(self) -> bool:
        if self.pedigree is None:
            self.messagebox.showinfo("PyPedal", "Open a pedigree file first.")
            return False
        return True

    def _run_background(self, work, title: str, finish=None) -> None:
        if self._busy:
            return
        self._gui_progress = GuiProgressBridge()
        self._set_busy(True)
        self._set_status(f"{title}…")
        complete = finish or self._finish_background

        def runner():
            result = None
            error = None
            try:
                result = work()
            except pyp_errors.PyPedalError as exc:
                # A PyPedal exception describes something the user can act on --
                # a malformed record, a missing option, an absent dependency --
                # so show the message, not a traceback. These paths must not
                # call sys.exit(0) and take the GUI down while reporting success.
                error = f"{type(exc).__name__}\n\n{exc}"
            except Exception:
                # Anything else is unexpected and the traceback is the useful
                # part; it is still caught, so the GUI survives.
                error = traceback.format_exc()
            self.root.after(0, lambda: complete(title, result, error))

        threading.Thread(target=runner, daemon=True).start()

    def _finish_background(self, title: str, result, error=None) -> None:
        self._set_busy(False)
        text = error if error is not None else result
        self._write(text)
        self._set_status(f"{title} finished — {os.path.basename(self.filename)}" if self.filename else title)

    def _finish_load(self, attempted_path: str, title: str, result, error=None) -> None:
        self._set_busy(False)
        outcome = apply_pedigree_load(
            self,
            attempted_path,
            pedigree=None if error else result,
            error=error,
        )
        if outcome["ok"]:
            self._write(result.metadata.stringme())
            self._set_status(outcome["status"])
            return
        self._write(outcome["output"])
        self._set_status(outcome["status"])

    def open_pedigree(self) -> None:
        path = self.filedialog.askopenfilename(
            title="Open pedigree",
            filetypes=[("Pedigree files", "*.ped"), ("All files", "*.*")],
        )
        if not path:
            return

        options = pedigree_open_options(
            path,
            self.format_var.get(),
            self.sep_var.get(),
            self.renumber_var.get(),
        )

        def work():
            ped = pyp_newclasses.NewPedigree(options)
            ped.load(progress=self._gui_progress)
            return ped

        self._run_background(
            work,
            "Loading pedigree",
            finish=lambda title, result, error: self._finish_load(path, title, result, error),
        )

    def show_metadata(self) -> None:
        if not self._need_pedigree():
            return
        self._write(self.pedigree.metadata.stringme())
        self._set_status(f"Metadata — {os.path.basename(self.filename)}")

    def list_animals(self) -> None:
        if not self._need_pedigree():
            return
        self._write(_list_animals(self.pedigree))
        self._set_status(f"{len(self.pedigree.pedigree)} animals")

    def calc_inbreeding(self) -> None:
        if not self._need_pedigree():
            return

        def work():
            result = pyp_nrm.inbreeding(
                self.pedigree, method="meu_luo", progress=self._gui_progress
            )
            result_file = f"{self.pedigree.kw.get('filetag')}_inbreeding.dat"
            return _format_inbreeding(result, result_file=result_file)

        self._run_background(work, "Calculating inbreeding")

    def calc_founders(self) -> None:
        if not self._need_pedigree():
            return

        def work():
            result = pyp_metrics.effective_founders_lacy(self.pedigree)
            fe = result.get("fa_effective_founders", result)
            lines = [
                "Effective founders (Lacy)",
                "=" * 40,
                f"f_e = {fe}",
            ]
            for key, value in result.items():
                if key != "fa_effective_founders":
                    lines.append(f"{key}: {value}")
            return "\n".join(lines)

        self._run_background(work, "Calculating effective founders")

    def inbreeding_by_year(self) -> None:
        if not self._need_pedigree():
            return

        def work():
            result = pyp_nrm.inbreeding(
                self.pedigree, method="meu_luo", progress=self._gui_progress
            )
            if isinstance(result, tuple):
                result = result[0]
            fx = result.get("fx", {}) if isinstance(result, dict) else {}
            by_year = {}
            for animal in self.pedigree.pedigree:
                year = getattr(animal, "by", None)
                coef = fx.get(animal.animalID)
                if year in (None, "", 0, -999) or coef is None:
                    continue
                by_year.setdefault(year, []).append(float(coef))
            if not by_year:
                return "No birth-year and inbreeding pairs were available."
            lines = ["Average inbreeding by birth year", "=" * 40]
            for year in sorted(by_year):
                values = by_year[year]
                mean = sum(values) / len(values)
                lines.append(f"{year}: {mean:.4f}  (n={len(values)})")
            return "\n".join(lines)

        self._run_background(work, "Summarizing inbreeding by year")

    def show_about(self) -> None:
        self.messagebox.showinfo(
            "About PyPedal",
            f"PyPedal {PYPEDAL_VERSION}\n"
            "Tools for pedigree analysis.\n\n"
            "This desktop app replaces the old wxPython GUI.\n"
            "It uses CustomTkinter so it runs on macOS, Windows, and Linux.\n\n"
            "Original author: John B. Cole\n"
            "License: GNU LGPL",
        )

    def run(self) -> None:
        self.root.mainloop()


def main(argv: Optional[list] = None) -> int:
    """
    Application entry point. Returns the process exit status.

    This is one of only three places in PyPedal permitted to end the process,
    and the only one that chooses a status from a failure. A PyPedal exception
    reaching here is reported as a message rather than a traceback -- it
    describes something the user can act on -- and mapped to a distinct non-zero
    status by ``exit_status_for``.
    """
    _ = argv or sys.argv[1:]
    try:
        app = PyPedalApp()
        app.run()
    except pyp_errors.PyPedalError as exc:
        print(f"[ERROR]: {exc}", file=sys.stderr)
        return exit_status_for(exc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
