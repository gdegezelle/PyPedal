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


# Exception class -> process exit status.
#
# PyPedal is a library: it raises, and never decides that the program is over.
# Choosing a status is an APPLICATION decision, so it happens here and nowhere
# else. Distinct codes let a wrapping shell script or CI job branch on the kind
# of failure instead of parsing stderr.
#
# The values below 64 follow ordinary convention (2 = usage/input); those at 64
# and above follow sysexits.h, where 70 is EX_SOFTWARE, an internal error.
# Anything unrecognised gets 1. Zero is never used for a failure -- the whole
# point is that these paths must never exit 0 and tell the shell they succeeded.
EXIT_STATUS = (
    (pyp_errors.PyPedalInputError, 2),
    (pyp_errors.PyPedalUsageError, 2),
    (pyp_errors.PyPedalConfigurationError, 3),
    (pyp_errors.PyPedalDependencyError, 4),
    (pyp_errors.PyPedalNotImplementedError, 5),
    (pyp_errors.PyPedalValidationError, 65),
    (pyp_errors.PyPedalInternalError, 70),
    # 73 is EX_CANTCREAT: the output file could not be produced. Distinct from
    # 70 because it is not a PyPedal defect in the general case -- the value
    # genuinely does not fit the declared format -- and distinct from 2 because
    # the caller's pedigree file is usually not what is wrong.
    (pyp_errors.PyPedalExportFormatError, 73),
    (pyp_errors.PyPedalError, 1),
)


def exit_status_for(exc: BaseException) -> int:
    """Map a PyPedal exception to the process exit status the CLI should use."""
    for klass, status in EXIT_STATUS:
        if isinstance(exc, klass):
            return status
    return 1


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


def normalize_sepchar(raw: Optional[str]) -> str:
    """Map the GUI separator field to a ``sepchar``.

    The field defaults to empty (meaning a space). A comma with leftover
    spaces from that default -- ``', '`` or ``' ,'`` -- must still be a
    comma, or a CSV with no space after the delimiter is read as one
    column. A tab is kept. Only spaces are trimmed, so a tab is not
    stripped to empty and then replaced by a space.
    """
    if raw is None:
        return " "
    text = str(raw)
    if text == "\t":
        return "\t"
    stripped = text.strip(" ")
    if stripped == "":
        return " "
    return stripped


def pedigree_open_options(
    path: str,
    pedformat: str,
    sepchar: str,
    renumber: bool,
) -> dict:
    """Options dict the desktop app passes to ``NewPedigree``."""
    return {
        "pedfile": path,
        "pedformat": (pedformat or "").strip() or "asd",
        "sepchar": normalize_sepchar(sepchar),
        "renumber": bool(renumber),
        "messages": "quiet",
        "pedname": os.path.basename(path),
    }


def _format_inbreeding(result: Any) -> str:
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
        lines.append("Coefficients by animal")
        lines.append("-" * 40)
        for animal_id, coef in sorted(fx.items(), key=lambda item: item[0]):
            lines.append(f"  {animal_id}: {pyp_io.format_display_coefficient(coef)}")
        if all(float(coef) == 0.0 for coef in fx.values()):
            lines.append("")
            lines.append(
                "No inbreeding in this pedigree (every coefficient is 0).\n"
                "That is expected for new_lacy.ped. Try mrode.ped (format asd)\n"
                "or hartlandclark.ped (format asdb) to see non-zero values."
            )
    return "\n".join(lines)


def _list_animals(pedigree) -> str:
    lines = [
        f"{'ID':>8}  {'Sire':>8}  {'Dam':>8}  {'Year':>6}  {'Sex':<4}  Name",
        "-" * 64,
    ]
    for animal in pedigree.pedigree:
        year = getattr(animal, "by", None)
        year_label = "" if year is None else year
        sex = getattr(animal, "sex", "")
        name = getattr(animal, "name", "")
        lines.append(
            f"{animal.animalID:>8}  {animal.sireID:>8}  {animal.damID:>8}  "
            f"{year_label:>6}  {str(sex):<4}  {name}"
        )
    return "\n".join(lines)


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

        ctk.CTkButton(controls, text="Open pedigree…", command=self.open_pedigree, width=140).pack(
            side="left", padx=8, pady=10
        )

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
        for label, command in actions:
            ctk.CTkButton(sidebar, text=label, command=command, anchor="w").pack(
                fill="x", padx=10, pady=6
            )

        ctk.CTkButton(sidebar, text="About", command=self.show_about, fg_color="gray").pack(
            fill="x", padx=10, pady=(18, 6)
        )

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

    def _write(self, text: str) -> None:
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)

    def _set_status(self, text: str) -> None:
        self.status.configure(text=text)

    def _need_pedigree(self) -> bool:
        if self.pedigree is None:
            self.messagebox.showinfo("PyPedal", "Open a pedigree file first.")
            return False
        return True

    def _run_background(self, work, title: str) -> None:
        if self._busy:
            return
        self._busy = True
        self._set_status(f"{title}…")

        def runner():
            try:
                result = work()
            except pyp_errors.PyPedalError as exc:
                # A PyPedal exception describes something the user can act on --
                # a malformed record, a missing option, an absent dependency --
                # so show the message, not a traceback. These paths must not
                # call sys.exit(0) and take the GUI down while reporting success.
                result = f"{type(exc).__name__}\n\n{exc}"
            except Exception:
                # Anything else is unexpected and the traceback is the useful
                # part; it is still caught, so the GUI survives.
                result = traceback.format_exc()
            self.root.after(0, lambda: self._finish_background(title, result))

        threading.Thread(target=runner, daemon=True).start()

    def _finish_background(self, title: str, result: str) -> None:
        self._busy = False
        self._write(result)
        self._set_status(f"{title} finished — {os.path.basename(self.filename)}" if self.filename else title)

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
            ped.load()
            self.pedigree = ped
            self.filename = path
            return ped.metadata.stringme()

        self._run_background(work, "Loading pedigree")

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
            result = pyp_nrm.inbreeding(self.pedigree, method="meu_luo")
            return _format_inbreeding(result)

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
            result = pyp_nrm.inbreeding(self.pedigree, method="meu_luo")
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
