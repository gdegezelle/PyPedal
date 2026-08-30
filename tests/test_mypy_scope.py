"""Mypy checks production library modules and is reporting-only for 4.0."""
import tomllib
from pathlib import Path

from _pedhelpers import REPO


def _mypy_tool():
    with (Path(REPO) / "pyproject.toml").open("rb") as fh:
        return tomllib.load(fh)["tool"]["mypy"]


def test_mypy_scope_is_production_modules_not_examples():
    cfg = _mypy_tool()
    files = cfg["files"]
    if isinstance(files, str):
        files = [files]
    assert files == ["PyPedal"]
    exclude = cfg["exclude"]
    if isinstance(exclude, str):
        exclude = [exclude]
    blob = " ".join(exclude)
    assert "PyPedal/examples/" in blob


def test_mypy_does_not_pin_a_zero_error_release_requirement():
    """4.0 records annotation debt; it does not require a clean mypy run."""
    text = (Path(REPO) / "pyproject.toml").read_text(encoding="utf-8")
    section = text.split("[tool.mypy]", 1)[1].split("[", 1)[0]
    compact = " ".join(line.lstrip("# ").strip() for line in section.splitlines())
    compact = compact.lower()
    assert "reporting-only" in compact
    assert "does not fail" in compact
