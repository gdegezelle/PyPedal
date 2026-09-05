"""Intentional failing test used by test_temp_hygiene.py. Not collected by default."""

from pathlib import Path

from _pedhelpers import owned_temp_dir


def test_intentional_failure_still_owns_a_temp_dir():
    path = owned_temp_dir(prefix="pypedal_failpath_")
    (Path(path) / "sentinel.txt").write_text(path, encoding="utf-8")
    print(f"PYPEDAL_FAILPATH_OWNED={path}")
    assert False, "intentional failure to prove teardown still runs"
