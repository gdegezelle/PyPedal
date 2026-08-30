"""RC4 an earlier revision: headless PDF pedigree reports.

Artifact checks use pypdf. Production generation uses ReportLab.
Tests write only under tmp_path / isolated tmpdirs.
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

from PyPedal import pyp_reports
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import REPO, chdir_tmp, load_corpus, load_corpus_from_path

needs_reportlab = pytest.mark.skipif(
    importlib.util.find_spec("reportlab") is None,
    reason="ReportLab (reports extra) is not installed",
)
needs_pypdf = pytest.mark.skipif(
    importlib.util.find_spec("pypdf") is None,
    reason="pypdf (test extra) is not installed",
)


def _pdf_text(path):
    from pypdf import PdfReader
    reader = PdfReader(str(path))
    return "\n".join((page.extract_text() or "") for page in reader.pages)


def _page_count(path):
    from pypdf import PdfReader
    return len(PdfReader(str(path)).pages)


def _assert_valid_pdf(path, min_pages=1):
    from pypdf import PdfReader
    path = Path(path)
    assert path.is_file()
    data = path.read_bytes()
    assert len(data) > 64
    assert data.startswith(b"%PDF")
    reader = PdfReader(str(path))
    assert len(reader.pages) >= min_pages
    return reader


def _write_rows(tmp_path, rows, name="fixture.ped"):
    path = tmp_path / name
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


FIFTEEN_SLOT_ROWS = [
    "1 0 0 GGP_SSS",
    "2 0 0 GGP_SSD",
    "3 0 0 GGP_SDS",
    "4 0 0 GGP_SDD",
    "5 0 0 GGP_DSS",
    "6 0 0 GGP_DSD",
    "7 0 0 GGP_DDS",
    "8 0 0 GGP_DDD",
    "9 1 2 GP_SS",
    "10 3 4 GP_SD",
    "11 5 6 GP_DS",
    "12 7 8 GP_DD",
    "13 9 10 Parent_Sire",
    "14 11 12 Parent_Dam",
    "15 13 14 Subject_Proband",
]


def _load_rows(tmp_path, rows, pedformat="asd", **overrides):
    path = _write_rows(tmp_path, rows)
    return load_corpus_from_path(str(path), pedformat, **overrides)


@needs_reportlab
@needs_pypdf
class TestHelpers:
    def test_initialize_restores_frame_geometry(self):
        ped = load_corpus("mrode.ped")
        settings = pyp_reports._pdfInitialize(ped)
        calcs = settings["_pdfCalcs"]
        assert "_frame_width" in calcs
        assert "_frame_height" in calcs
        assert calcs["_frame_width"] == calcs["_right_margin"] - calcs["_left_margin"]
        assert calcs["_frame_height"] == calcs["_top_margin"] - calcs["_bottom_margin"]
        assert calcs["_frame_width"] > 0
        assert calcs["_frame_height"] > 0
        assert calcs["_page_width"] > 0
        assert calcs["_page_height"] > 0

    def test_letter_and_a4_page_sizes(self):
        letter_ped = load_corpus("mrode.ped", paper_size="letter")
        a4_ped = load_corpus("mrode.ped", paper_size="A4")
        letter = pyp_reports._pdfInitialize(letter_ped)["_pdfCalcs"]
        a4 = pyp_reports._pdfInitialize(a4_ped)["_pdfCalcs"]
        assert letter["_page"] != a4["_page"]
        assert a4["_frame_width"] > 0
        assert letter["_frame_width"] > 0

    def test_title_page_wraps_long_text(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "title_wrap.pdf"
        path = pyp_reports.pdf_pedigree_metadata(
            ped,
            titlepage=1,
            reporttitle="A very long title that must wrap across several centred lines without overflowing",
            reportauthor="A very long author line that should also wrap rather than run off the page",
            reportfile=str(out),
        )
        reader = _assert_valid_pdf(path, min_pages=2)
        text = _pdf_text(path)
        assert "very long title" in text
        assert "very long author" in text
        assert len(reader.pages) >= 2


@needs_reportlab
@needs_pypdf
class TestMeanMetricHelper:
    def test_pdf_mean_metric_by_still_writes(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "mean_metric.pdf"
        rc = pyp_reports.pdfMeanMetricBy(
            ped,
            {1990: 0.1, 1991: 0.25},
            titlepage=0,
            reportfile=str(out),
        )
        assert rc == 1
        _assert_valid_pdf(out)
        text = _pdf_text(out)
        assert "1990" in text
        assert "0.1" in text

    def test_snake_case_metric_alias(self):
        assert pyp_reports.pdf_mean_metric_by is pyp_reports.pdfMeanMetricBy
        assert pyp_reports.mean_metric_by is pyp_reports.meanMetricBy


@needs_reportlab
@needs_pypdf
class TestMetadataPdf:
    def test_artifact_and_content(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "mrode_metadata.pdf"
        path = pyp_reports.pdf_pedigree_metadata(
            ped,
            titlepage=0,
            reporttitle="Mrode metadata",
            reportfile=str(out),
        )
        assert path == str(out)
        reader = _assert_valid_pdf(path)
        text = _pdf_text(path)
        assert "Records" in text
        assert "6" in text
        assert "Founders" in text
        assert "Unknown birth years" in text
        assert "1800" not in text
        assert "1900" not in text
        assert reader.pages

    def test_title_page_optional(self, tmp_path):
        ped = load_corpus("mrode.ped")
        plain = tmp_path / "plain.pdf"
        titled = tmp_path / "titled.pdf"
        pyp_reports.pdf_pedigree_metadata(ped, titlepage=0, reportfile=str(plain))
        pyp_reports.pdf_pedigree_metadata(
            ped,
            titlepage=1,
            reporttitle="Custom Metadata Title",
            reportauthor="Test Author",
            reportfile=str(titled),
        )
        assert _page_count(plain) >= 1
        assert _page_count(titled) >= 2
        titled_text = _pdf_text(titled)
        assert "Custom Metadata Title" in titled_text
        assert "Test Author" in titled_text

    def test_unknown_chronology_not_fake_year(self, tmp_path):
        ped = load_corpus("mrode.ped")
        assert all(animal.by is None and animal.bd is None for animal in ped.pedigree)
        out = tmp_path / "unknown_years.pdf"
        pyp_reports.pdf_pedigree_metadata(ped, reportfile=str(out))
        text = _pdf_text(out)
        assert "1800" not in text
        assert "1900" not in text
        assert "01011800" not in text
        assert "None" not in text

    def test_real_1800_and_1900_years_are_shown(self, tmp_path):
        rows = ["1 0 0 1800", "2 0 0 1900"]
        ped = _load_rows(tmp_path, rows, pedformat="asdy")
        assert ped.pedigree[0].by == 1800
        assert ped.pedigree[1].by == 1900
        out = tmp_path / "real_years.pdf"
        pyp_reports.pdf_pedigree_metadata(ped, reportfile=str(out))
        text = _pdf_text(out)
        assert "1800" in text
        assert "1900" in text

    def test_nina_title_does_not_crash(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "nina.pdf"
        path = pyp_reports.pdf_pedigree_metadata(
            ped,
            titlepage=1,
            reporttitle="Niña",
            reportfile=str(out),
        )
        _assert_valid_pdf(path, min_pages=2)
        text = _pdf_text(path)
        assert "Ni" in text

    def test_overwrite(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "overwrite.pdf"
        pyp_reports.pdf_pedigree_metadata(
            ped, titlepage=1, reporttitle="First Title", reportfile=str(out)
        )
        first_size = out.stat().st_size
        pyp_reports.pdf_pedigree_metadata(
            ped, titlepage=1, reporttitle="Second Title", reportfile=str(out)
        )
        assert out.stat().st_size > 64
        text = _pdf_text(out)
        assert "Second Title" in text
        assert "First Title" not in text
        assert first_size > 64

    def test_missing_parent_directory_raises(self, tmp_path):
        ped = load_corpus("mrode.ped")
        missing = tmp_path / "no_such_dir" / "out.pdf"
        with pytest.raises(OSError):
            pyp_reports.pdf_pedigree_metadata(ped, reportfile=str(missing))
        assert not missing.exists()

    def test_pathlike_and_default_name(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "pathlike.pdf"
        path = pyp_reports.pdf_pedigree_metadata(ped, reportfile=out)
        assert path == os.fspath(out)
        with chdir_tmp() as work:
            default = pyp_reports.pdf_pedigree_metadata(ped, reportfile="")
            assert default.endswith("_metadata.pdf")
            assert Path(work, default).is_file() or Path(default).is_file()


@needs_reportlab
@needs_pypdf
class TestThreeGenerationPdf:
    def test_mrode_half_founder_and_unknown_parent(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "mrode_3gen.pdf"
        path = pyp_reports.pdf_three_gen_ped(5, ped, reportfile=str(out))
        assert path == str(out)
        _assert_valid_pdf(path)
        text = _pdf_text(path)
        assert "Pedigree for" in text
        assert "Unknown Parent" in text
        assert "1800" not in text
        subject = next(a for a in ped.pedigree if a.animalID == 5)
        assert str(subject.originalID) in text
        sire = next(a for a in ped.pedigree if a.animalID == subject.sireID)
        dam = next(a for a in ped.pedigree if a.animalID == subject.damID)
        assert str(sire.originalID) in text
        assert str(dam.originalID) in text

    def test_fifteen_slot_ancestry(self, tmp_path):
        ped = _load_rows(tmp_path, FIFTEEN_SLOT_ROWS, pedformat="asdn")
        subject = next(a for a in ped.pedigree if a.name == "Subject_Proband")
        out = tmp_path / "fifteen.pdf"
        pyp_reports.pdf_three_gen_ped(subject.animalID, ped, reportfile=str(out))
        text = _pdf_text(out)
        for token in (
            "Subject_Proband",
            "Parent_Sire",
            "Parent_Dam",
            "GP_SS",
            "GP_SD",
            "GP_DS",
            "GP_DD",
            "GGP_SSS",
            "GGP_SSD",
            "GGP_SDS",
            "GGP_SDD",
            "GGP_DSS",
            "GGP_DSD",
            "GGP_DDS",
            "GGP_DDD",
        ):
            assert token in text, token
        assert "Unknown Parent" not in text

    def test_half_founder(self, tmp_path):
        rows = ["10 0 0 KnownSire", "20 10 0 HalfFounder"]
        ped = _load_rows(tmp_path, rows, pedformat="asdn")
        subject = next(a for a in ped.pedigree if a.name == "HalfFounder")
        out = tmp_path / "half.pdf"
        pyp_reports.pdf_three_gen_ped(subject.animalID, ped, reportfile=str(out))
        text = _pdf_text(out)
        assert "HalfFounder" in text
        assert "KnownSire" in text
        assert "Unknown Parent" in text
        assert "-1" not in text

    def test_multi_subject_pages_and_atomic_invalid(self, tmp_path):
        ped = load_corpus("mrode.ped")
        good = tmp_path / "two_subjects.pdf"
        path = pyp_reports.pdf_three_gen_ped([5, 6], ped, reportfile=str(good))
        assert _page_count(path) >= 2

        target = tmp_path / "must_not_write.pdf"
        target.write_text("sentinel", encoding="utf-8")
        with pytest.raises(PyPedalUsageError):
            pyp_reports.pdf_three_gen_ped([5, 999], ped, reportfile=str(target))
        assert target.read_text(encoding="utf-8") == "sentinel"

        missing = tmp_path / "never_created.pdf"
        with pytest.raises(PyPedalUsageError):
            pyp_reports.pdf_three_gen_ped(999, ped, reportfile=str(missing))
        assert not missing.exists()

    def test_rejects_call_name_strings(self, tmp_path):
        ped = _load_rows(tmp_path, FIFTEEN_SLOT_ROWS, pedformat="asdn")
        out = tmp_path / "names.pdf"
        with pytest.raises(PyPedalUsageError):
            pyp_reports.pdf_three_gen_ped("Subject_Proband", ped, reportfile=str(out))
        assert not out.exists()

    def test_historical_aliases(self, tmp_path):
        ped = load_corpus("mrode.ped")
        assert pyp_reports.pdfPedigreeMetadata is pyp_reports.pdf_pedigree_metadata
        assert pyp_reports.pdf3GenPed is pyp_reports.pdf_three_gen_ped
        assert pyp_reports.pdf_3_gen_ped is pyp_reports.pdf_three_gen_ped
        meta = tmp_path / "alias_meta.pdf"
        three = tmp_path / "alias_3gen.pdf"
        assert pyp_reports.pdfPedigreeMetadata(ped, 0, "", "", str(meta))
        assert pyp_reports.pdf3GenPed(5, ped, 0, "", "", str(three))
        _assert_valid_pdf(meta)
        _assert_valid_pdf(three)

    def test_none_chronology_three_gen(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "none_chrono.pdf"
        pyp_reports.pdf_three_gen_ped(5, ped, reportfile=str(out))
        text = _pdf_text(out)
        assert "1800" not in text
        assert "1900" not in text
        assert "01011800" not in text


class TestMissingReportLab:
    def test_calling_pdf_raises_typed_dependency_error(self, tmp_path):
        ped = load_corpus("mrode.ped")
        out = tmp_path / "blocked.pdf"
        code = r"""
import sys
class BlockReportLab:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "reportlab" or fullname.startswith("reportlab."):
            raise ModuleNotFoundError("blocked for test")
        return None
sys.meta_path.insert(0, BlockReportLab())
for name in list(sys.modules):
    if name == "reportlab" or name.startswith("reportlab."):
        del sys.modules[name]
from PyPedal.pyp_errors import PyPedalDependencyError
from PyPedal import pyp_reports
try:
    pyp_reports.pdf_pedigree_metadata(None, reportfile="x.pdf")
except PyPedalDependencyError:
    print("dependency-ok")
else:
    raise SystemExit("expected PyPedalDependencyError")
"""
        env = dict(os.environ)
        env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(tmp_path),
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "dependency-ok" in proc.stdout
        assert not out.exists()

    def test_import_pypedal_without_reportlab(self, tmp_path):
        code = r"""
import sys
class BlockReportLab:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "reportlab" or fullname.startswith("reportlab."):
            raise ModuleNotFoundError("blocked for test")
        return None
sys.meta_path.insert(0, BlockReportLab())
for name in list(sys.modules):
    if name == "reportlab" or name.startswith("reportlab."):
        del sys.modules[name]
import PyPedal
from PyPedal import pyp_newclasses, pyp_nrm, pyp_reports
print("imported", PyPedal.__version__)
print("reports-module", pyp_reports.__name__)
"""
        env = dict(os.environ)
        env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(tmp_path),
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "imported" in proc.stdout
