"""
an earlier revision workstream A: current 4.0rc3 caller/callee contracts.

The December-2025 pypedal3 audit recorded signature drift among
NewPedigree.renumber, pyp_utils.renumber, and inbreeding_tabular ->
fast_a_matrix. These tests ask only whether *this* tree's callers match
*this* tree's callees. They do not import pypedal3 signatures.
"""
import os
import tempfile

from PyPedal import pyp_nrm, pyp_utils
from PyPedal.pyp_newclasses import NewPedigree

from _pedhelpers import chdir_tmp, load_corpus, load_corpus_from_path

# Mrode (2005) Table 2.1: animal 5 is the offspring of a half-sib mating.
MRODE_F5 = 0.125


def _gappy_rows_path():
    tmp = tempfile.mkdtemp(prefix="pypedal_sig_")
    path = os.path.join(tmp, "gappy.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("10 0 0\n20 10 0\n30 20 10\n")
    return path


def test_sig1_newpedigree_renumber_zero_arg_establishes_invariants():
    """NewPedigree.renumber() takes no extra arguments and finalizes IDs."""
    ped = load_corpus_from_path(_gappy_rows_path(), "asd", renumber=False)
    originals = [a.originalID for a in ped.pedigree]
    assert originals == [10, 20, 30]

    ok = ped.renumber()
    assert ok is True
    ids = [a.animalID for a in ped.pedigree]
    assert ids == [1, 2, 3]
    assert all(ped.pedigree[i].animalID == i + 1 for i in range(3))
    assert ped.kw["pedigree_is_renumbered"] is True
    assert ped.idmap[10] == 1
    assert ped.idmap[20] == 2
    assert ped.idmap[30] == 3
    assert ped.backmap[1] == 10
    assert ped.backmap[2] == 20
    assert ped.backmap[3] == 30
    missing = str(ped.kw["missing_parent"])
    present = {a.animalID for a in ped.pedigree}
    for animal in ped.pedigree:
        for parent in (animal.sireID, animal.damID):
            if str(parent) != missing:
                assert parent in present


def test_sig2_pyp_utils_renumber_accepts_animal_list():
    """pyp_utils.renumber's first argument is a list of animal objects."""
    ped = load_corpus_from_path(_gappy_rows_path(), "asd", renumber=False)
    out = pyp_utils.renumber(
        ped.pedigree,
        missingparent=ped.kw["missing_parent"],
        animaltype="new",
        io="no",
    )
    assert isinstance(out, list)
    assert [a.animalID for a in out] == [1, 2, 3]
    assert all(out[i].animalID == i + 1 for i in range(3))


def test_sig2_pyp_utils_renumber_also_accepts_newpedigree():
    """Existing dual contract: a NewPedigree is iterable and is returned."""
    with chdir_tmp():
        ped = load_corpus_from_path(_gappy_rows_path(), "asd", renumber=False)
        out = pyp_utils.renumber(
            ped,
            missingparent=ped.kw["missing_parent"],
            animaltype="new",
            io="no",
        )
    assert isinstance(out, NewPedigree)
    assert out is ped
    assert [a.animalID for a in ped.pedigree] == [1, 2, 3]


def test_sig3_inbreeding_tabular_calls_fast_a_matrix_with_list_and_kw(monkeypatch):
    captured = {}
    real = pyp_nrm.fast_a_matrix

    def wrapper(pedigree, pedopts, save=False, method="sparse", debug=False, fill=1):
        captured["pedigree_type"] = type(pedigree)
        captured["is_list"] = isinstance(pedigree, list)
        captured["opts_is_dict"] = isinstance(pedopts, dict)
        captured["has_missing_parent"] = "missing_parent" in pedopts
        captured["method"] = method
        return real(pedigree, pedopts, save=save, method=method, debug=debug, fill=fill)

    monkeypatch.setattr(pyp_nrm, "fast_a_matrix", wrapper)
    ped = load_corpus("mrode.ped")
    result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    assert captured["is_list"] is True
    assert captured["opts_is_dict"] is True
    assert captured["has_missing_parent"] is True
    fx = {int(k): float(v) for k, v in result["fx"].items()}
    assert abs(fx[5] - MRODE_F5) < 1e-3


def test_sig4_mrode_tabular_animal_5_unchanged():
    ped = load_corpus("mrode.ped")
    result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    fx = {int(k): float(v) for k, v in result["fx"].items()}
    assert abs(fx[5] - MRODE_F5) < 1e-12
    assert fx[6] > 0.0
