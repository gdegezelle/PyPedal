"""PH4-DOC-1..13: an earlier revision living-manual product/API contract.

Historical manuals and RC1/RC2 notes are out of scope.
"""
import os
import shutil
import tempfile
from pathlib import Path

import pytest

from PyPedal import pyp_app, pyp_graphics, pyp_snp, pyp_utils
from PyPedal.pyp_errors import (
    PyPedalNotImplementedError,
    PyPedalPedigreeStructureError,
)
from PyPedal.pyp_newclasses import load_pedigree

from test_manual_pages import REQUIRED_PAGES, USER, _mkdocs_nav_files, _user_markdown

pytestmark = pytest.mark.docs

PHASE4_PAGES = (
    "object-model.md",
    "ids-and-missing-parents.md",
    "genomics.md",
    "graphics.md",
    "desktop-application.md",
    "limitations.md",
    "references.md",
)


def _text(name):
    return (USER / name).read_text(encoding="utf-8")


def _plain(name):
    return _text(name).replace("*", "").lower()


def _load_mrode():
    examples = Path(__import__("PyPedal").__file__).resolve().parent / "examples"
    work = Path(tempfile.mkdtemp())
    dest = work / "mrode.ped"
    shutil.copy(examples / "mrode.ped", dest)
    ped = load_pedigree(
        options={
            "pedfile": str(dest),
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )
    return ped, work


def test_ph4_doc1_phase4_pages_exist():
    for name in PHASE4_PAGES:
        path = USER / name
        assert path.is_file(), path
        assert path.stat().st_size > 0


def test_ph4_doc2_phase4_pages_exist_on_disk():
    for name in PHASE4_PAGES:
        path = USER / name
        assert path.is_file(), name
    for name in REQUIRED_PAGES:
        assert (USER / name).is_file(), name


def test_ph4_doc3_limitations_has_mandatory_categories():
    text = _plain("limitations.md")
    for token in (
        "unsupported",
        "removed",
        "limited",
        "bounded",
        "external",
        "ambiguous",
    ):
        assert token in text, token


def test_ph4_doc4_genes_removed_or_not_supported():
    lim = _plain("limitations.md")
    assert "genes" in lim
    assert "removed" in lim
    assert "not supported" in lim


def test_ph4_doc5_mating_coi_supported_bounds():
    lim = _plain("limitations.md")
    assert "mating_coi" in lim
    assert "supported" in lim
    assert "gens" in lim
    assert "renumber" in lim
    mating_block = lim.split("test mating")[1].split("tabular inbreeding")[0]
    assert "pypedalnotimplementederror" not in mating_block
    assert "may be reconsidered later" not in mating_block


def test_ph4_doc6_gen_coeff_pattie_unsupported():
    lim = _plain("limitations.md")
    assert "gen_coeff" in lim
    assert "pattie" in lim
    assert "not supported" in lim


def test_ph4_doc7_unknown_chronology_is_none():
    lim = _plain("limitations.md")
    gen = _plain("birth-dates-and-chronology.md")
    assert "none" in lim
    assert "legacy_missing_byear_token" in lim or "legacy" in lim
    assert "none" in gen
    assert "igen" in gen
    assert "age" in gen


def test_ph4_doc8_tabular_large_pedigree():
    lim = _plain("limitations.md")
    assert "tabular" in lim
    assert "automatic switch" in lim
    assert "meu_luo" in lim


def test_ph4_doc9_sex_chrometype_refusal():
    lim = _plain("limitations.md")
    gene = _plain("gene-dropping.md")
    assert "chrometype" in lim and "sex" in lim
    assert "pypedalusagerror" in lim or "usageerror" in lim
    assert "chrometype" in gene and "sex" in gene


def test_ph4_doc10_vanraden_1992_vs_2008():
    genomic = _plain("genomics.md")
    inbreeding = _plain("inbreeding.md")
    assert "1992" in genomic and "2008" in genomic
    assert 'method="vanraden"' in genomic or "method='vanraden'" in genomic
    assert "pedigree" in genomic and "genomic" in genomic
    assert "1992" in inbreeding and "2008" in inbreeding
    assert "grm=" in genomic
    assert "g_matrix" in genomic
    assert "no supported" in genomic or "not" in genomic


def test_ph4_doc11_originalid_vs_animalid():
    ids = _plain("ids-and-missing-parents.md")
    assert "originalid" in ids
    assert "animalid" in ids
    assert "idmap" in ids
    assert "half-founder" in ids
    assert "animalid - 1" in ids or "animalid-1" in ids.replace(" ", "")


def test_object_model_documents_factual_derived_cached_state():
    text = _plain("object-model.md")
    assert "factual" in text
    assert "derived" in text
    assert "cached" in text or "computed" in text
    assert "animal.fa" in text or "fa" in text
    assert "f_computed" in text
    assert "overwritten" in text
    glossary = _plain("glossary.md")
    assert "fa" in glossary


def test_ph4_doc12_wxpython_is_historical():
    hits = []
    for path in _user_markdown():
        text = path.read_text(encoding="utf-8")
        lower = text.lower()
        start = 0
        while True:
            idx = lower.find("wxpython", start)
            if idx < 0:
                break
            window = lower[max(0, idx - 80): idx + 40]
            if not any(
                token in window
                for token in ("old", "historical", "replaces", "replaced")
            ):
                hits.append(f"{path.name}: wxPython presented as current")
            start = idx + 1
    assert hits == []


def test_ph4_doc13_references_page_and_core_citations():
    text = _text("references.md")
    plain = text.replace("*", "").lower()
    assert (USER / "references.md").stat().st_size > 0
    for token in (
        "mrode",
        "meuwissen",
        "vanraden",
        "lacy",
        "boichard",
        "suwanlee",
        "baumung",
        "pattie",
        "cole, j. b. 2007",
    ):
        assert token in plain, token
    assert "1992" in text and "2008" in text
    nav = _mkdocs_nav_files()
    assert "references.md" in nav


def test_ph4_pdf_reports_documented_as_supported():
    """PDF pedigree reports are a supported headless ReportLab capability."""
    graphics = _plain("graphics.md")
    lim = _plain("limitations.md")
    reports = _plain("pdf-reports.md")
    assert "pdf_pedigree_metadata" in reports or "pdfpedigreemetadata" in reports
    assert "pdf_three_gen_ped" in reports or "pdf3genped" in reports
    assert "pdf-reports.md" in _text("graphics.md")
    pdf_block = "pdf pedigree reports".join(
        lim.split("pdf pedigree reports")[1:]
    )
    assert "supported" in pdf_block
    assert "removed from the 4.0" not in pdf_block
    assert "not part of the supported" not in pdf_block
    assert "pdf reports](pdf-reports.md)" in graphics or "pdf-reports.md" in _text("graphics.md")


def test_ph4_stale_product_claims_absent():
    banned = (
        "setup.py install",
        "sourceforge",
        "production/stable",
        "v4.0.0 final",
    )
    hits = []
    for path in _user_markdown():
        lower = path.read_text(encoding="utf-8").lower()
        for token in banned:
            start = 0
            while True:
                idx = lower.find(token, start)
                if idx < 0:
                    break
                window = lower[max(0, idx - 80): idx + len(token) + 80]
                if not any(
                    marker in window
                    for marker in ("not", "do not", "no ", "absent")
                ):
                    hits.append(f"{path.name}: {token}")
                start = idx + 1
    assert hits == []


def test_ph4_smoke_objects_and_ids():
    ped, _work = _load_mrode()
    animal = ped.pedigree[4]
    assert animal.animalID == 5
    assert animal.originalID == 5
    assert animal.sireID == 4
    assert animal.damID == 3
    assert animal.founder == "n"
    assert ped.pedigree[animal.animalID - 1] is animal
    assert ped.idmap[animal.originalID] == animal.animalID
    half = [a for a in ped.pedigree if str(a.sireID) != "0" and str(a.damID) == "0"]
    assert any(a.originalID == 4 for a in half)
    assert all(a.founder == "n" for a in half)


def test_ph4_smoke_reorder_renumber():
    ped, _work = _load_mrode()
    assert ped.kw["renumber"] is True
    assert ped.kw["pedigree_is_renumbered"] is True
    ids = [a.animalID for a in ped.pedigree]
    assert ids == list(range(1, 7))
    snapshot = list(ped.pedigree)
    ordered = pyp_utils.fast_reorder(snapshot, missingparent=0)
    assert [a.animalID for a in ordered] == ids
    with pytest.raises(PyPedalPedigreeStructureError):
        pyp_utils.reorder([snapshot[0], snapshot[0]], missingparent=0)


def test_ph4_smoke_genomic_grm():
    work = Path(tempfile.mkdtemp())
    (work / "ped.ped").write_text("1 0 0\n2 0 0\n3 1 2\n4 1 2\n")
    (work / "geno.txt").write_text(
        "1 chip1 4 0120\n2 chip1 4 1201\n3 chip1 4 2012\n4 chip1 4 0000\n"
    )
    cwd = os.getcwd()
    try:
        os.chdir(work)
        ped = load_pedigree(
            options={
                "pedfile": str(work / "ped.ped"),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
                "snpfile": str(work / "geno.txt"),
            }
        )
        grm = pyp_snp.form_grm_from_snp(ped)
        fg = pyp_snp.compute_genomic_inbreeding_from_grm(ped, grm=grm)
        hom = pyp_snp.compute_genomic_homozygosity_from_snp(ped)
    finally:
        os.chdir(cwd)
    assert round(fg[1], 4) == 0.4667
    assert hom[4] == 1.0


def test_ph4_smoke_graphics_draw_pedigree():
    pytest.importorskip("pydot")
    if shutil.which("dot") is None:
        pytest.skip("Graphviz dot is not on PATH")
    ped, work = _load_mrode()
    stem = work / "mrode_graph"
    ok = pyp_graphics.draw_pedigree(
        ped, gfilename=str(stem), gformat="png", gdot=0, gtitle="mrode"
    )
    png = Path(str(stem) + ".png")
    assert ok == 1
    assert png.is_file()
    assert png.stat().st_size > 0


def test_ph4_smoke_desktop_app_import():
    assert callable(pyp_app.main)
    assert callable(pyp_app.exit_status_for)
    assert pyp_app.exit_status_for(PyPedalNotImplementedError("x", "y")) == 5
