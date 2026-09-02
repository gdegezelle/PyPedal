#!/usr/bin/env python3
"""Build PyPedal.app outside the repository.

Usage (from the repository root, with the macos-app extra installed):

    python tools/macos/build_app.py
    python tools/macos/build_app.py /tmp/pypedal-macos-app
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from bundle_config import (
    BUNDLE_NAME,
    forbidden_bundle_hits,
    icns_path,
    macos_dir,
    read_project_version,
    repository_root,
)

_COCOA_MARKERS = ("libqcocoa", "qcocoa")
_OFFSCREEN_MARKERS = ("libqoffscreen", "qoffscreen")


def _bundle_files(app: Path) -> list[Path]:
    return [path for path in app.rglob("*") if path.is_file()]


def validate_app(app: Path) -> None:
    """Fail if the bundle is missing required pieces or includes Griffon/CTk."""
    if not app.is_dir():
        raise SystemExit(f"expected application bundle at {app}")
    files = _bundle_files(app)
    names = [str(path) for path in files]
    blob = "\n".join(names)
    lowered = blob.lower()
    hits = forbidden_bundle_hits(names)
    if hits:
        raise SystemExit("forbidden paths in bundle:\n  " + "\n  ".join(hits[:20]))
    if not any(marker in lowered for marker in _COCOA_MARKERS):
        raise SystemExit("Cocoa Qt platform plugin not found in bundle")
    executable = app / "Contents" / "MacOS" / BUNDLE_NAME
    if not executable.is_file():
        raise SystemExit(f"missing executable {executable}")
    plist = app / "Contents" / "Info.plist"
    if not plist.is_file():
        raise SystemExit("missing Info.plist")
    print("app:", app)
    print("files:", len(files))
    print("cocoa plugin: present")
    print(
        "offscreen plugin:",
        "present" if any(marker in lowered for marker in _OFFSCREEN_MARKERS) else "absent",
    )
    print("icon:", "present" if icns_path() is not None else "not supplied")


def build_app(out: Path) -> Path:
    repo = repository_root()
    spec = macos_dir() / "PyPedal.spec"
    out.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--distpath",
        str(out / "dist"),
        "--workpath",
        str(out / "work"),
        str(spec),
    ]
    subprocess.run(command, cwd=repo, check=True)
    app = out / "dist" / f"{BUNDLE_NAME}.app"
    validate_app(app)
    return app


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the PyPedal macOS application bundle.")
    parser.add_argument(
        "out",
        nargs="?",
        default=str(Path(os.environ.get("TMPDIR", "/tmp")) / "pypedal-macos-app"),
        help="Output directory (default: $TMPDIR/pypedal-macos-app)",
    )
    args = parser.parse_args(argv)
    out = Path(args.out).expanduser().resolve()
    app = build_app(out)
    version = read_project_version()
    print(f"version: {version}")
    print(f"MACOS APP OK: {app}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
