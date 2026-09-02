"""Canonical PyPedal desktop process entry.

``pypedal``, ``pypedal-gui``, and ``python -m PyPedal`` all call ``main``.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

from PyPedal.__version__ import version as PYPEDAL_VERSION
from PyPedal.application import exit_status_for


def _missing_pyside6_message(exc: BaseException) -> str:
    return (
        "The PyPedal desktop needs PySide6.\n"
        "Install it with:  pip install 'PyPedal[gui]'\n"
        f"Original error: {exc}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pypedal",
        description="PyPedal desktop application.",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the PyPedal version and exit.",
    )
    parser.add_argument(
        "--self-test",
        type=Path,
        metavar="PEDIGREE",
        help="Engineering smoke: load PEDIGREE and run representative jobs, then exit.",
    )
    parser.add_argument(
        "pedigree",
        nargs="?",
        type=Path,
        help="Optional pedigree file to open after startup.",
    )
    return parser


def _run_self_test(pedigree: Path) -> int:
    """Load a tiny pedigree and exercise analysis adapters. No GUI."""
    from PyPedal.application import (
        PedigreeOpenOptions,
        PedigreeSession,
        load_into_session,
        run_effective_founders,
        run_inbreeding,
        run_mating_coi,
        run_relationship,
        run_theoretical_ne,
    )

    source = pedigree.expanduser().resolve()
    previous = Path.cwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        try:
            session = PedigreeSession()
            load_into_session(session, source, PedigreeOpenOptions(separator=" ").normalized())
            run_inbreeding(session)
            run_effective_founders(session)
            loaded = session.pedigree
            if loaded is None:
                raise RuntimeError("self-test load produced no pedigree")
            ids = [int(animal.animalID) for animal in loaded.pedigree[:2]]
            if len(ids) >= 2:
                run_relationship(session, ids[0], ids[1])
                run_mating_coi(session, ids[0], ids[1])
            run_theoretical_ne(session)
        finally:
            os.chdir(previous)
    print("SELF-TEST OK")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Process-entry for the desktop. Returns an exit status."""
    args = build_parser().parse_args(argv if argv is not None else sys.argv[1:])
    if args.version:
        print(PYPEDAL_VERSION)
        return 0
    if args.self_test is not None:
        try:
            return _run_self_test(args.self_test)
        except Exception as exc:
            print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
            return exit_status_for(exc)
    try:
        from PyPedal.desktop.app import run_desktop
    except ImportError as exc:
        print(_missing_pyside6_message(exc), file=sys.stderr)
        return 4
    return run_desktop(args.pedigree)
