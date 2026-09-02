"""Temporary Qt desktop entry point.

``pypedal`` and ``pypedal-gui`` still launch CustomTkinter during 4.2-B.
This module is the development PySide6 path.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PyPedal.__version__ import version as PYPEDAL_VERSION


def _missing_pyside6_message(exc: BaseException) -> str:
    return (
        "The PyPedal Qt desktop needs PySide6.\n"
        "Install it with:  pip install 'PyPedal[desktop]'\n"
        f"Original error: {exc}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pypedal-qt",
        description="PyPedal Qt desktop (development entry point).",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the PyPedal version and exit.",
    )
    parser.add_argument(
        "pedigree",
        nargs="?",
        type=Path,
        help="Optional pedigree file to open after startup.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Process-entry for the Qt desktop. Returns an exit status."""
    args = build_parser().parse_args(argv if argv is not None else sys.argv[1:])
    if args.version:
        print(PYPEDAL_VERSION)
        return 0
    try:
        from PyPedal.desktop.app import run_desktop
    except ImportError as exc:
        print(_missing_pyside6_message(exc), file=sys.stderr)
        return 4
    return run_desktop(args.pedigree)
