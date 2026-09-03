"""Public user-manual facts (docs/manual/)."""
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import REPO

pytestmark = pytest.mark.docs

MANUAL = Path(REPO) / "docs" / "manual"
MKDOCS = Path(REPO) / "mkdocs.yml"

REQUIRED_PAGES = (
    "index.md",
    "installation.md",
    "first-pedigree.md",
    "first-analysis.md",
    "ids-and-missing-parents.md",
    "pedigree-input.md",
    "pedigree-formats.md",
    "animal-identities-and-names.md",
    "birth-dates-and-chronology.md",
    "validation.md",
    "saving-and-exporting.md",
    "inbreeding.md",
    "relationships.md",
    "mating.md",
    "effective-founders.md",
    "effective-ancestors.md",
    "founder-genome-equivalents.md",
    "lacy-and-boichard.md",
    "generation-intervals.md",
    "gene-dropping.md",
    "pdf-reports.md",
    "graphics.md",
    "sqlite.md",
    "desktop-application.md",
    "genomics.md",
    "limitations.md",
    "large-pedigrees.md",
    "recipes.md",
    "configuration.md",
    "pedigree-format-codes.md",
    "object-model.md",
    "api-overview.md",
    "glossary.md",
    "references.md",
    "notices.md",
)

FORBIDDEN_LABELS = (
    r"\bBL1\b",
    r"\bBL2\b",
    r"\bBL4B\b",
    r"Finding 31",
    r"Finding 36",
    r"Finding 37",
    r"Finding 39",
    r"MEDIUM-",
    r"AGE-3a",
    r"RC1-S1",
    r"holdout protocol",
    r"adjudication finding",
    r"Phase 4",
)

MRODE = "1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n"


def _manual_text():
    return "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(MANUAL.glob("*.md"))
    )


def _load_mrode(tmp_path):
    pedfile = tmp_path / "mrode.ped"
    pedfile.write_text(MRODE, encoding="utf-8")
    return load_pedigree(
        options={
            "pedfile": str(pedfile),
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )


def test_manual_pages_exist():
    assert MANUAL.is_dir()
    for name in REQUIRED_PAGES:
        path = MANUAL / name
        assert path.is_file(), name
        assert path.stat().st_size > 0


def test_mkdocs_nav_lists_manual_pages():
    text = MKDOCS.read_text(encoding="utf-8")
    assert "docs_dir: docs/manual" in text
    listed = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.endswith(".md") and ":" in stripped:
            listed.append(stripped.rsplit(":", 1)[-1].strip())
    assert set(listed) == set(REQUIRED_PAGES)


def test_manual_has_no_internal_engineering_labels():
    blob = _manual_text()
    hits = []
    for pattern in FORBIDDEN_LABELS:
        if re.search(pattern, blob, flags=re.IGNORECASE):
            hits.append(pattern)
    assert hits == [], hits
    lowered = blob.lower()
    assert "rc4 contract" not in lowered
    assert "the rc4" not in lowered


def test_python_support_and_no_setup_py_install():
    text = (MANUAL / "installation.md").read_text(encoding="utf-8")
    assert "Python 3.12" in text
    lowered = text.lower()
    assert "python 2.7" not in lowered
    assert "this chapter does not cover" in lowered
    assert "sourceforge" in lowered
    assert "numarray" in lowered
    assert "pypi" in lowered
    assert "no published package" in lowered or "not been published" in lowered


def test_mating_is_documented():
    text = (MANUAL / "mating.md").read_text(encoding="utf-8").lower()
    assert "mating_coi" in text
    assert "mating_coi_group" in text
    assert "a_ij" in text.replace(" ", "") or "a<sub>ij</sub>" in (MANUAL / "mating.md").read_text(encoding="utf-8").lower()
    assert "read-only" in text
    assert "[(4, 3), (1, 2)]" in (MANUAL / "mating.md").read_text(encoding="utf-8") or "[(1, 2)" in (MANUAL / "recipes.md").read_text(encoding="utf-8")


def test_pdf_reports_are_documented():
    text = (MANUAL / "pdf-reports.md").read_text(encoding="utf-8")
    assert "pdf_pedigree_metadata" in text
    assert "pdf_three_gen_ped" in text
    assert "15-slot" in text or "15 slot" in text.lower()
    assert "ReportLab" in text
    assert "headless" in text.lower()


def test_unknown_chronology_is_none():
    text = (MANUAL / "birth-dates-and-chronology.md").read_text(encoding="utf-8")
    assert "`None`" in text or "None" in text
    assert "1800" in text and "1900" in text
    assert "legacy_missing_byear_token" in text


def test_griffon_sample_and_dense_nrm_warning():
    text = (MANUAL / "large-pedigrees.md").read_text(encoding="utf-8")
    assert "griffonbruxellois_2026_named_pyp.ped" in text
    assert "asdxbn" in text
    assert "98,001" in text
    assert "6,689" in text
    assert "915" in text
    assert "3,997" in text
    assert "11.018378975785259" in text
    assert "asdxb" in text
    assert '"sepchar": ","' in text
    assert "deterministic dataset regression" in text.lower()
    assert "80 GB" in text or "80 gb" in text.lower()
    assert "meu_luo" in text
    assert "no automatic switch" in text.lower()
    assert "curated project" in text.lower()


def test_cole_attribution_and_notices():
    notices = (MANUAL / "notices.md").read_text(encoding="utf-8")
    assert "John B. Cole" in notices
    assert "Geert Degezelle" in notices
    assert "LGPL" in notices
    assert "relicensed copy" in notices.lower() or "not a relicensed" in notices.lower()
    assert "Rick Muller" in notices
    intro = (MANUAL / "index.md").read_text(encoding="utf-8")
    assert "John B. Cole" in intro
    assert "Copyright (C) 2025-2026 Geert Degezelle" not in _manual_text()


def test_first_pedigree_and_inbreeding_snippets(tmp_path):
    ped = _load_mrode(tmp_path)
    assert len(ped.pedigree) == 6
    result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    assert result["fx"][5] == 0.125


def test_relationship_and_mating_snippets(tmp_path):
    ped = _load_mrode(tmp_path)
    assert pyp_metrics.relationship(4, 3, ped) == 0.25
    assert pyp_metrics.relationship(1, 2, ped) == 0.0
    assert pyp_metrics.mating_coi(4, 3, ped) == 0.125
    assert pyp_metrics.mating_coi(1, 2, ped) == 0.0
    assert pyp_metrics.mating_coi(5, 5, ped) == 0.5625
    got = pyp_metrics.mating_coi_group([(4, 3), (1, 2)], ped)
    assert got["matings"][(4, 3)] == 0.125
    assert got["matings"][(1, 2)] == 0.0
    with pytest.raises(PyPedalUsageError):
        pyp_metrics.relationship(99999, 3, ped)


@pytest.mark.skipif(
    importlib.util.find_spec("reportlab") is None,
    reason="reports extra is not installed in this interpreter",
)
def test_pdf_metadata_snippet(tmp_path):
    from PyPedal import pyp_reports

    ped = _load_mrode(tmp_path)
    out = tmp_path / "metadata.pdf"
    path = pyp_reports.pdf_pedigree_metadata(ped, reportfile=str(out))
    assert Path(path).is_file()
    assert Path(path).stat().st_size > 0


@pytest.mark.skipif(
    importlib.util.find_spec("mkdocs") is None,
    reason="mkdocs extra is not installed in this interpreter",
)
def test_mkdocs_strict_build_for_manual(tmp_path):
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
