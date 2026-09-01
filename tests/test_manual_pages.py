"""Helpers and existence checks for the public MkDocs manual."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest
import tomllib

from _pedhelpers import REPO
from test_user_manual import REQUIRED_PAGES

pytestmark = pytest.mark.docs

MANUAL = Path(REPO) / "docs" / "manual"
USER = MANUAL
MKDOCS = Path(REPO) / "mkdocs.yml"


def _user_markdown():
    return sorted(USER.glob("*.md"))


def _mkdocs_nav_files():
    text = MKDOCS.read_text(encoding="utf-8")
    files = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.endswith(".md") and ":" in stripped:
            files.append(stripped.rsplit(":", 1)[-1].strip())
    return files


def test_doc1_required_living_manual_files_exist():
    assert USER.is_dir()
    for name in REQUIRED_PAGES:
        path = USER / name
        assert path.is_file(), path
        assert path.stat().st_size > 0


def test_doc2_mkdocs_nav_files_exist():
    assert MKDOCS.is_file()
    nav = _mkdocs_nav_files()
    assert nav, "mkdocs.yml nav listed no markdown files"
    for name in nav:
        assert (MANUAL / name).is_file(), name
    for required in REQUIRED_PAGES:
        assert (USER / required).is_file(), required


def test_doc3_living_docs_are_not_python2_or_setup_py():
    """Those strings may appear only as explicit exclusions, not as advice."""
    tokens = (
        "Python 2.7",
        "python 2.7",
        "Python 3.6+",
        "python 3.6+",
        "setup.py install",
        "SourceForge",
    )
    hits = []
    for path in _user_markdown():
        text = path.read_text(encoding="utf-8")
        lower = text.lower()
        for token in tokens:
            start = 0
            needle = token.lower()
            while True:
                idx = lower.find(needle, start)
                if idx < 0:
                    break
                window = lower[max(0, idx - 80): idx + len(needle) + 80]
                if "not" not in window and "do not" not in window and "no " not in window:
                    hits.append(f"{path.name}: {token!r}")
                start = idx + 1
    assert hits == []


def test_doc4_genes_is_not_advertised_as_supported():
    for path in _user_markdown():
        text = path.read_text(encoding="utf-8")
        if "GENES" not in text:
            continue
        lowered = text.replace("*", "").lower()
        assert (
            "not supported" in lowered
            or "not a supported" in lowered
            or "not part of" in lowered
            or "removed" in lowered
        ), path.name


def test_doc5_mating_coi_is_advertised_as_supported():
    index = (USER / "index.md").read_text(encoding="utf-8")
    mating = (USER / "mating.md").read_text(encoding="utf-8")
    limitations = (USER / "limitations.md").read_text(encoding="utf-8")
    assert "mating_coi" in index or "test mating" in index.lower() or "prospective" in index.lower()
    assert "mating_coi" in mating
    assert "test mating" in mating.lower() or "prospective" in mating.lower()
    compact = mating.lower().replace(" ", "").replace("_", "")
    assert "aij" in compact or "a<sub>ij</sub>" in mating.lower()
    assert "read-only" in mating.lower()
    lim = limitations.replace("*", "").lower()
    assert "mating_coi" in lim


def test_doc6_site_branding_is_release_line_not_rc3():
    text = MKDOCS.read_text(encoding="utf-8")
    assert "site_name: PyPedal" in text
    assert "4.0.0-rc3" not in text
    assert "docs_dir: docs/manual" in text


def test_doc7_format_codes_follow_the_loader():
    text = (USER / "pedigree-format-codes.md").read_text(encoding="utf-8")
    assert "`y`" in text
    assert "`b`" in text
    assert "dam" in text.lower()
    assert "`p`" in text
    assert "PEDIGREE_FORMAT_CODES.txt" in text


def test_docs_extra_is_not_a_runtime_dependency():
    with open(Path(REPO) / "pyproject.toml", "rb") as fh:
        project = tomllib.load(fh)
    extras = project["project"]["optional-dependencies"]
    assert "docs" in extras
    assert "mkdocs" in " ".join(extras["docs"]).lower()
    runtime = " ".join(project["project"]["dependencies"]).lower()
    assert "mkdocs" not in runtime
    assert "mkdocs" not in " ".join(extras.get("all", [])).lower()


@pytest.mark.skipif(
    importlib.util.find_spec("mkdocs") is None,
    reason="mkdocs extra is not installed in this interpreter",
)
def test_mkdocs_strict_build_outside_the_repository(tmp_path):
    site = tmp_path / "site"
    result = subprocess.run(
        [sys.executable, "-m", "mkdocs", "build", "--strict", "-d", str(site)],
        cwd=REPO,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": ""},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (site / "index.html").is_file()
    assert not (Path(REPO) / "site").exists()
