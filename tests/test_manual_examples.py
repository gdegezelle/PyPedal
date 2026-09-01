"""PH5-1..12: an earlier revision product-hygiene living-manual / example contract."""
import inspect
import re
from pathlib import Path

import pytest

from PyPedal.pyp_newclasses import NewPedigree, loadPedigree

from _pedhelpers import REPO, EXAMPLES, chdir_tmp
from test_examples_integration import KNOWN_FAILING, NOT_EXAMPLES, SCRIPTS
from test_manual_pages import USER, _mkdocs_nav_files, _user_markdown

pytestmark = pytest.mark.docs

FORMAT_CODES = Path(REPO) / "PyPedal" / "PEDIGREE_FORMAT_CODES.txt"
INI = Path(REPO) / "PyPedal" / "PyPedal.ini"
README = Path(REPO) / "README.md"
NOTES = Path(REPO) / "notes"

# Loader attribute map that used to drift from PEDIGREE_FORMAT_CODES.txt.
FORMAT_ATTR = {
    "a": "animalID",
    "s": "sireID",
    "d": "damID",
    "y": "by",
    "b": "bd",
    "D": "damName",
    "S": "sireName",
    "A": "name",
    "p": "gencoeff",
    "h": "herd",
    "H": "originalHerd",
    "G": "genomicInbreeding",
    "Y": "genomicHomozygosity",
    "T": "traits",
    "P": "SNPgenotype",
}


def test_ph5_1_examples_page_exists():
    assert (USER / "recipes.md").is_file()
    assert (USER / "recipes.md").stat().st_size > 0


def test_ph5_2_options_page_exists():
    assert (USER / "configuration.md").is_file()
    assert (USER / "configuration.md").stat().st_size > 0


def test_ph5_3_format_codes_txt_is_not_version_2():
    text = FORMAT_CODES.read_text(encoding="utf-8")
    assert "VERSION: 2." not in text
    assert "2.0.0a20" not in text
    assert "PyPedal 4.0 pedigree format codes" in text


def test_ph5_4_format_codes_agree_with_loader():
    with chdir_tmp():
        ped = NewPedigree(
            {
                "pedfile": "probe.ped",
                "messages": "quiet",
                "pedigree_summary": 0,
            }
        )
    for code, attr in FORMAT_ATTR.items():
        assert code in ped.pedformat_codes, code
        assert ped.new_animal_attr[code] == attr, (code, attr)
    src = Path(REPO) / "PyPedal" / "pyp_newclasses.py"
    parse_src = Path(REPO) / "PyPedal" / "_pyp_parse.py"
    loader = src.read_text(encoding="utf-8") + "\n" + parse_src.read_text(encoding="utf-8")
    assert '("birthyear", "y"' in loader
    assert '("birthdate", "b"' in loader
    assert "index('T')" not in loader
    assert "pedformat_locations['traits']" not in loader
    txt = FORMAT_CODES.read_text(encoding="utf-8")
    assert "This is DAM" in txt
    assert "birth YEAR" in txt
    assert "birth DATE" in txt
    assert "CALCULATION is not supported" in txt
    assert "Does NOT currently populate" in txt
    assert "not a documented load recipe" in txt


def test_ph5_5_pypedal_ini_is_current_template():
    text = INI.read_text(encoding="utf-8")
    assert "sourceforge" not in text.lower()
    assert "VERSION: 2." not in text
    assert "2.0.0a9" not in text
    assert "PyPedal 4.0 example configuration" in text
    assert "renumber = 1" in text
    assert "reorder = 0" in text
    assert "legacy_missing_byear_token" in text
    assert "missing_byear = 1800" not in text
    assert "missing_parent = 0" in text


def test_ph5_8_supported_examples_not_in_known_failing():
    must_pass = {
        "new_snp.py",
        "new_snp2.py",
        "new_methods.py",
        "new_networkx.py",
        "new_sqlite.py",
        "new_graphics.py",
        "new_ids.py",
        "new_jbc.py",
        "new_inbreeding2.py",
        "new_options.py",
        "duplicates.py",
        "new_reporting.py",
        "new_lacy.py",
        "new_db.py",
    }
    overlap = must_pass & set(KNOWN_FAILING)
    assert not overlap, overlap
    for name in must_pass:
        assert name in SCRIPTS, name
        assert name not in NOT_EXAMPLES


def test_ph5_9_known_failing_has_current_reasons():
    # an earlier revision listed three simulate()-lifecycle failures. an earlier revision repaired
    # that lifecycle. an earlier revision removed the remaining dyad_census call from
    # new_simulate.py. Active KNOWN_FAILING must stay empty.
    assert KNOWN_FAILING == {}


def test_ph5_10_user_docs_do_not_link_archived_notes_as_current():
    hits = []
    for path in list(_user_markdown()) + [README]:
        text = path.read_text(encoding="utf-8")
        if "notes/QUICKSTART" in text or "notes/REPOSITORY_COMPARISON" in text:
            hits.append(path.name)
    assert hits == []


def test_ph5_11_no_stale_g_matrix_keyword_in_active_examples():
    hits = []
    for name in SCRIPTS:
        path = Path(EXAMPLES) / name
        if "g_matrix" in path.read_text(encoding="utf-8"):
            hits.append(name)
    assert hits == []


def test_ph5_12_active_examples_do_not_advertise_removed_paths():
    needles = (
        "getCursorSQA",
        "import wx",
        "from wx",
    )
    hits = []
    for name in SCRIPTS:
        text = (Path(EXAMPLES) / name).read_text(encoding="utf-8")
        if "GENES" in text and "not supported" not in text.lower():
            hits.append(f"{name}: GENES")
        for needle in needles:
            if needle in text:
                hits.append(f"{name}: {needle}")
    assert hits == []


def test_ph5_pdf_reports_are_supported_in_docs_and_example():
    """PDF pedigree reports are a supported 4.0 capability after an earlier revision."""
    index = (USER / "index.md").read_text(encoding="utf-8")
    reports = (USER / "pdf-reports.md").read_text(encoding="utf-8")
    limitations = (USER / "limitations.md").read_text(encoding="utf-8")
    graphics = (USER / "graphics.md").read_text(encoding="utf-8")
    example = (Path(EXAMPLES) / "new_reporting.py").read_text(encoding="utf-8")
    assert "pdf-reports.md" in _mkdocs_nav_files()
    assert "pdf_pedigree_metadata" in reports
    assert "pdf_three_gen_ped" in reports
    assert "15-slot" in reports or "15 slot" in reports.lower() or "15-slot" in reports
    assert ".[reports]" in reports or "'.[reports]'" in reports
    assert "pdfPedigreeMetadata" in reports
    assert "pdf3GenPed" in reports
    intro = index.lower()
    assert "pdf" in intro or "reportlab" in intro
    assert "what this manual does not treat as supported" not in intro
    lim = limitations.replace("*", "").lower()
    pdf_block = "pdf pedigree reports".join(
        lim.split("pdf pedigree reports")[1:]
    ).split("tabular inbreeding")[0]
    assert "not part of the supported" not in pdf_block
    assert "are not restored" not in pdf_block
    assert "unsupported pdf" not in pdf_block
    assert "removed from the 4.0" not in graphics.lower()
    assert "not part of the 4.0" not in graphics.lower()
    assert "pdf_pedigree_metadata" in example
    assert "pdf_three_gen_ped" in example
    assert "wx" not in example.lower()
    assert "customtkinter" not in example.lower()
    assert "pdf pedigree reports are removed" not in example.lower()


def test_ph5_option_defaults_match_newpedigree():
    with chdir_tmp():
        ped = NewPedigree(
            {
                "pedfile": "probe.ped",
                "messages": "quiet",
                "pedigree_summary": 0,
            }
        )
    assert ped.kw["renumber"] is True
    assert ped.kw["reorder"] is False
    assert "missing_byear" not in ped.kw
    assert ped.kw["legacy_missing_byear_token"] is None
    assert ped.kw["missing_parent"] == 0
    assert ped.kw["pedformat"] == "asd"
    assert ped.kw["sepchar"] == " "
    assert ped.kw["form_nrm"] is False
    assert ped.kw["matrix_type"] == "sparse"
    assert ped.kw["snpfile"] is False
    assert ped.kw["simulate_pedigree"] is False


def test_ph5_load_pedigree_default_pedsource_is_file():
    assert inspect.signature(loadPedigree).parameters["pedsource"].default == "file"


def test_ph5_format_codes_txt_lists_required_codes():
    text = FORMAT_CODES.read_text(encoding="utf-8")
    for code in (
        "a", "s", "d", "A", "S", "D", "y", "b", "T", "P", "G", "Y", "p", "h", "H",
    ):
        assert re.search(rf"^{code} =", text, re.MULTILINE), code
