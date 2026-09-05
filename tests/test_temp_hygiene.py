"""Regression tests for test-owned temporary directories and logfile handles."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path

from _pedhelpers import (
    REPO,
    cleanup_owned_temp_paths,
    close_owned_pypedal_log_handlers,
    load_corpus,
    owned_temp_dir,
    owned_temp_snapshot,
)

from PyPedal.pyp_newclasses import (
    _PYPEDAL_OWNED_HANDLER,
    PYPEDAL_LOGGER_NAME,
    load_pedigree,
)


def test_owned_temp_dir_is_removed_by_cleanup_helper():
    start = owned_temp_snapshot()
    path = Path(owned_temp_dir(prefix="pypedal_hygiene_"))
    marker = path / "marker.txt"
    marker.write_text("owned", encoding="utf-8")
    assert path.is_dir()
    assert marker.is_file()
    close_owned_pypedal_log_handlers()
    cleanup_owned_temp_paths(start=start)
    assert not path.exists()


def test_cleanup_does_not_remove_paths_registered_before_snapshot():
    start = owned_temp_snapshot()
    earlier = Path(owned_temp_dir(prefix="pypedal_hygiene_keep_"))
    (earlier / "keep.txt").write_text("keep", encoding="utf-8")
    later_start = owned_temp_snapshot()
    later = Path(owned_temp_dir(prefix="pypedal_hygiene_drop_"))
    (later / "drop.txt").write_text("drop", encoding="utf-8")
    close_owned_pypedal_log_handlers()
    cleanup_owned_temp_paths(start=later_start)
    assert earlier.is_dir()
    assert not later.exists()
    cleanup_owned_temp_paths(start=start)
    assert not earlier.exists()


def test_load_corpus_registers_an_owned_directory():
    start = owned_temp_snapshot()
    ped = load_corpus("mrode.ped")
    assert len(ped.pedigree) == 6
    pedfile = Path(ped.kw["pedfile"])
    tmp = pedfile.parent
    assert tmp.is_dir()
    close_owned_pypedal_log_handlers()
    cleanup_owned_temp_paths(start=start)
    assert not tmp.exists()


def test_closed_logfile_can_be_deleted_immediately(tmp_path):
    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    load_pedigree(
        options={
            "pedfile": str(pedfile),
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }
    )
    package = logging.getLogger(PYPEDAL_LOGGER_NAME)
    owned = [
        handler for handler in package.handlers if getattr(handler, _PYPEDAL_OWNED_HANDLER, False)
    ]
    assert owned
    logfile = Path(owned[0].baseFilename)
    assert logfile.is_file()
    close_owned_pypedal_log_handlers()
    logfile.unlink()
    assert not logfile.exists()


def test_failure_path_teardown_removes_owned_temp_dir():
    probe = Path(__file__).with_name("_hygiene_fail_probe.py")
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            str(probe),
            "-o",
            "python_files=_hygiene_fail_probe.py",
            "--no-repo-guard",
            "-q",
            "-s",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    owned = None
    for stream in (proc.stdout, proc.stderr):
        for line in stream.splitlines():
            if line.startswith("PYPEDAL_FAILPATH_OWNED="):
                owned = line.split("=", 1)[1]
    assert owned, proc.stdout + proc.stderr
    assert not os.path.exists(owned)
