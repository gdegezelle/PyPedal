"""SCI-DOC-1..10: an earlier revision scientific living-manual contracts.

Historical manuals and RC1/RC2 notes are out of scope.
"""
import shutil
from pathlib import Path

import pytest

from PyPedal import pyp_metrics, pyp_nrm, pyp_utils
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_newclasses import load_pedigree

from test_manual_pages import REQUIRED_PAGES, USER, _user_markdown
from _pedhelpers import owned_temp_dir

pytestmark = pytest.mark.docs

PHASE3_PAGES = (
    "inbreeding.md",
    "relationships.md",
    "effective-founders.md",
    "lacy-and-boichard.md",
    "gene-dropping.md",
    "generation-intervals.md",
)


def _load_example(name, pedformat):
    examples = Path(__import__("PyPedal").__file__).resolve().parent / "examples"
    work = Path(owned_temp_dir())
    dest = work / name
    shutil.copy(examples / name, dest)
    ped = load_pedigree(
        options={
            "pedfile": str(dest),
            "pedformat": pedformat,
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )
    return ped, work


def test_sci_doc1_phase3_pages_exist():
    for name in PHASE3_PAGES:
        path = USER / name
        assert path.is_file(), path
        assert path.stat().st_size > 0


def test_sci_doc2_phase3_pages_exist_on_disk():
    for name in PHASE3_PAGES:
        path = USER / name
        assert path.is_file(), name
    for name in REQUIRED_PAGES:
        assert (USER / name).is_file(), name


def test_sci_doc4_gen_coeff_pattie_not_supported():
    for path in _user_markdown():
        text = path.read_text(encoding="utf-8")
        if "gen_coeff" not in text and "Pattie" not in text:
            continue
        lowered = " ".join(text.replace("*", "").lower().split())
        assert (
            "not supported" in lowered
            or "does not compute" in lowered
            or "refused" in lowered
            or "leave false" in lowered
        ), path.name


def test_sci_doc6_large_pedigree_tabular_warning():
    text = (USER / "large-pedigrees.md").read_text(encoding="utf-8")
    lowered = text.replace("*", "").lower()
    assert "automatic switch" in lowered
    assert "meu_luo" in lowered
    assert "tabular" in lowered
    assert "griffon" in lowered


def test_sci_doc7_ng_is_not_raw_founder_count():
    text = (USER / "gene-dropping.md").read_text(encoding="utf-8")
    lowered = text.replace("*", "").lower()
    compact = text.replace(" ", "")
    assert "1/(2" in compact or "1 / (2" in text
    assert "historical founders" in lowered or "raw founder" in lowered
    assert "98,001" in text or "98001" in text


def test_sci_doc8_unknown_chronology_is_none():
    text = (USER / "birth-dates-and-chronology.md").read_text(encoding="utf-8")
    lowered = text.replace("*", "").lower()
    assert "none" in lowered
    assert "1800" in text
    assert "igen" in lowered


def test_sci_doc9_relationship_unresolved_ids_named():
    text = (USER / "relationships.md").read_text(encoding="utf-8")
    assert "PyPedalUsageError" in text
    assert "0.0" in text
    assert "unresolved" in text.lower()


def test_sci_doc10_inbreeding_mrode():
    ped, _ = _load_example("mrode.ped", "asd")
    result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    assert result["fx"][5] == 0.125


def test_sci_doc10_relationship_mrode():
    ped, _ = _load_example("mrode.ped", "asd")
    assert pyp_metrics.relationship(4, 3, ped) == 0.25
    assert pyp_metrics.relationship(1, 2, ped) == 0.0
    with pytest.raises(PyPedalUsageError):
        pyp_metrics.relationship(99999, 3, ped)


def test_sci_doc10_lacy_new_lacy():
    ped, _ = _load_example("new_lacy.ped", "asd")
    got = pyp_metrics.effective_founders_lacy(ped)
    assert got["fa_founder_count"] == 3
    assert got["fa_effective_founders"] == pytest.approx(2.909090909, abs=1e-9)


def test_sci_doc10_boichard_figure2():
    ped, _ = _load_example("boichard2.ped", "asdg")
    fe = pyp_metrics.a_effective_founders_boichard(ped)
    fa = pyp_metrics.a_effective_ancestors_definite(ped)
    assert fe == pytest.approx(5.56, abs=0.05)
    assert fa == pytest.approx(2.94, abs=0.05)


def test_sci_doc10_gene_drop_boichard2():
    ped, _ = _load_example("boichard2.ped", "asdg")
    ng = pyp_metrics.effective_founder_genomes(
        ped, rounds=200, seed=31, output=False
    )
    assert ng == pytest.approx(2.56, abs=0.2)
    with pytest.raises(PyPedalUsageError):
        pyp_metrics.effective_founder_genomes(
            ped, chrometype="sex", output=False
        )


def test_sci_doc10_set_generation_mrode():
    ped, _ = _load_example("mrode.ped", "asd")
    assert [a.igen for a in ped.pedigree] == [-999.0] * 6
    assert pyp_utils.set_generation(ped) is True
    assert [a.igen for a in ped.pedigree] == [1, 1, 2, 2, 3, 4]
    assert [a.gen for a in ped.pedigree] == [-999.0] * 6


def test_sci_doc10_generation_intervals_example():
    ped, _ = _load_example("generations.ped", "asdbx")
    mean = pyp_metrics.generation_intervals(ped)["mean"]
    assert 13.0 < mean < 15.0
