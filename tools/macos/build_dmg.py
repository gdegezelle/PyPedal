#!/usr/bin/env python3
"""Build an engineering PyPedal.dmg outside the repository.

Usage:

    python tools/macos/build_dmg.py --app /tmp/pypedal-macos-app/dist/PyPedal.app
    python tools/macos/build_dmg.py --app ... --out /tmp/pypedal-macos-dmg
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from bundle_config import BUNDLE_NAME, VOLUME_NAME


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, capture_output=True, text=True)


def build_dmg(app: Path, destination: Path, *, volume_name: str = VOLUME_NAME) -> Path:
    """Copy ``app`` into a staging folder with an Applications link and image it."""
    if not app.is_dir() or app.suffix != ".app":
        raise SystemExit(f"expected a .app bundle, got {app}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.parent / f"{destination.stem}-staging"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    shutil.copytree(app, staging / f"{BUNDLE_NAME}.app", symlinks=True)
    (staging / "Applications").symlink_to("/Applications")
    if destination.exists():
        destination.unlink()
    _run(
        [
            "hdiutil",
            "create",
            "-volname",
            volume_name,
            "-srcfolder",
            str(staging),
            "-ov",
            "-format",
            "UDZO",
            str(destination),
        ]
    )
    shutil.rmtree(staging)
    return destination


def validate_dmg(dmg: Path, *, volume_name: str = VOLUME_NAME) -> None:
    """Mount, check layout, and detach."""
    attach = _run(["hdiutil", "attach", "-nobrowse", "-readonly", str(dmg)])
    mount = None
    for line in attach.stdout.splitlines():
        if "/Volumes/" in line:
            mount = Path(line.split("\t")[-1].strip())
    if mount is None:
        raise SystemExit(f"could not determine mount point:\n{attach.stdout}")
    try:
        app = mount / f"{BUNDLE_NAME}.app"
        applications = mount / "Applications"
        if not app.is_dir():
            raise SystemExit(f"mounted image missing {app}")
        if not applications.is_symlink() and not applications.exists():
            raise SystemExit("mounted image missing Applications link")
        print("dmg:", dmg)
        print("volume:", volume_name)
        print("mount:", mount)
        print("app on image: present")
        print("Applications link: present")
    finally:
        subprocess.run(["hdiutil", "detach", str(mount), "-quiet"], check=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a PyPedal engineering DMG.")
    parser.add_argument("--app", required=True, type=Path, help="Path to PyPedal.app")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Directory for PyPedal.dmg (default: $TMPDIR/pypedal-macos-dmg)",
    )
    args = parser.parse_args(argv)
    out_dir = args.out
    if out_dir is None:
        out_dir = Path(os.environ.get("TMPDIR", "/tmp")) / "pypedal-macos-dmg"
    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    dmg = out_dir / f"{BUNDLE_NAME}.dmg"
    app = args.app.expanduser().resolve()
    build_dmg(app, dmg)
    validate_dmg(dmg)
    print("size:", dmg.stat().st_size)
    print(f"MACOS DMG OK: {dmg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
