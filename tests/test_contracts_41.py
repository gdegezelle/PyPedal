"""4.1-B contracts: refuse plausible-wrong scientific sentinels and invalid args."""
import os
import tempfile
import warnings

import pytest

from _pedhelpers import chdir_tmp, load_corpus, load_corpus_from_path

from PyPedal import pyp_io, pyp_metrics, pyp_network, pyp_nrm
from PyPedal.pyp_errors import PyPedalError, PyPedalUsageError
from PyPedal.pyp_results import InbreedingResult

LACY_APPENDIX_A_FE = 2.909090909090909
MRODE_NE = 4.8
MEMMAP_NAMES = (
    "fast_a_matrix_mmap.bin",
    "fast_a_matrix_r_mmap.bin",
    "fast_partial_a_matrix_mmap.bin",
)


def _rows_ped(rows):
    tmp = tempfile.mkdtemp(prefix="pypedal_contracts_")
    path = os.path.join(tmp, "tiny.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, "asd")


def _memory_error(*_args, **_kwargs):
    raise MemoryError("simulated allocation failure")


def _assert_no_memmap_residue(directory):
    leftover = [
        name for name in MEMMAP_NAMES
        if os.path.exists(os.path.join(directory, name))
    ]
    assert leftover == [], leftover


def test_relationship_unrelated_pair_is_zero():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        assert pyp_metrics.relationship(1, 2, ped) == 0.0


def test_relationship_allocation_failure_raises(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm, "lil_matrix", _memory_error)
        with pytest.raises(PyPedalError, match="fast_a_matrix"):
            pyp_metrics.relationship(1, 2, ped)
        _assert_no_memmap_residue(tmp)


def test_inbreeding_success_is_inbreeding_result():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
        assert isinstance(result, InbreedingResult)
        assert isinstance(result, dict)
        assert result["fx"][5] == pytest.approx(0.125)


def test_inbreeding_allocation_failure_raises(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")

        def boom(_n):
            raise MemoryError("simulated allocation failure")

        monkeypatch.setattr(pyp_nrm, "_alloc_float_list", boom)
        with pytest.raises(PyPedalError, match="inbreeding_meuwissen_luo"):
            pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
        _assert_no_memmap_residue(tmp)


def test_a_coefficients_no_inbred_is_empty_dict():
    with chdir_tmp():
        ped = _rows_ped(["1 0 0", "2 0 0"])
        result = pyp_metrics.a_coefficients(ped, output=False)
        assert result == {}


def test_a_coefficients_allocation_failure_raises(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm.np, "zeros", _memory_error)
        monkeypatch.setattr(pyp_nrm, "lil_matrix", _memory_error)
        with pytest.raises(PyPedalError, match="fast_a_matrix"):
            pyp_metrics.a_coefficients(ped, output=False)
        _assert_no_memmap_residue(tmp)


def test_a_matrix_one_animal_is_identity():
    with chdir_tmp():
        ped = _rows_ped(["1 0 0"])
        matrix = pyp_nrm.a_matrix(ped)
        assert matrix.shape == (1, 1)
        assert float(matrix[0, 0]) == pytest.approx(1.0)


def test_a_matrix_allocation_failure_raises(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm.np, "zeros", _memory_error)
        with pytest.raises(PyPedalError, match="a_matrix"):
            pyp_nrm.a_matrix(ped)
        _assert_no_memmap_residue(tmp)


def test_a_inverse_from_file_missing_raises():
    with chdir_tmp() as tmp:
        missing = os.path.join(tmp, "no-such-inverse.pkl")
        with pytest.raises(PyPedalUsageError, match="a_inverse_from_file"):
            pyp_io.a_inverse_from_file(missing)


def test_fast_a_matrix_memoryerror_no_float32_memmap(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm.np, "zeros", _memory_error)
        with pytest.raises(PyPedalError) as caught:
            pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="dense")
        message = str(caught.value)
        assert "fast_a_matrix" in message
        assert "6-by-6" in message
        assert "float64" in message
        assert "meu_luo" in message
        _assert_no_memmap_residue(tmp)
        assert not any(name.endswith(".bin") for name in os.listdir(tmp))


def test_fast_a_matrix_r_memoryerror_no_float32_memmap(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm.np, "zeros", _memory_error)
        with pytest.raises(PyPedalError, match="fast_a_matrix_r"):
            pyp_nrm.fast_a_matrix_r(ped.pedigree, ped.kw, method="dense")
        _assert_no_memmap_residue(tmp)


def test_fast_partial_a_matrix_memoryerror_no_float32_memmap(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm.np, "zeros", _memory_error)
        with pytest.raises(PyPedalError, match="fast_partial_a_matrix"):
            pyp_nrm.fast_partial_a_matrix(
                ped.pedigree, 1, [1, 2], ped.kw, method="dense")
        _assert_no_memmap_residue(tmp)


def test_sparse_fast_a_matrix_memoryerror_raises(monkeypatch):
    with chdir_tmp() as tmp:
        ped = load_corpus("mrode.ped")
        monkeypatch.setattr(pyp_nrm, "lil_matrix", _memory_error)
        with pytest.raises(PyPedalError, match="fast_a_matrix"):
            pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="sparse")
        _assert_no_memmap_residue(tmp)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"method": "not-a-method"},
        {"gens": -1},
        {"gens": True},
    ],
)
def test_inbreeding_invalid_arguments_raise(kwargs):
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        with pytest.raises(PyPedalUsageError, match="inbreeding"):
            pyp_nrm.inbreeding(ped, output=False, **kwargs)


def test_foundercoi_invalid_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        opts = dict(ped.kw)
        opts["foundercoi"] = 2
        with pytest.raises(PyPedalUsageError, match="foundercoi"):
            pyp_nrm.fast_a_matrix(ped.pedigree, opts, method="dense")
        for accepted in (0, 1, False, True):
            opts["foundercoi"] = accepted
            matrix = pyp_nrm.fast_a_matrix(ped.pedigree, opts, method="dense")
            assert matrix is not None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"rounds": 0},
        {"loci": 0},
        {"seed": None},
        {"seed": "not-an-int"},
        {"rounds": True},
    ],
)
def test_dropped_ancestral_invalid_arguments_raise(kwargs):
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        with pytest.raises(PyPedalUsageError, match="dropped_ancestral_inbreeding"):
            pyp_metrics.dropped_ancestral_inbreeding(ped, **kwargs)


def test_fast_a_matrix_invalid_method_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        with pytest.raises(PyPedalUsageError, match="method"):
            pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="huge")


def test_inbreeding_aguilar_invalid_amethod_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        with pytest.raises(PyPedalUsageError, match="amethod"):
            pyp_nrm.inbreeding_aguilar(ped, amethod=9)
        with pytest.raises(PyPedalUsageError, match="amethod"):
            pyp_nrm.inbreeding_aguilar(ped, amethod=True)


def test_find_ancestors_g_invalid_gens_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        graph = pyp_network.ped_to_graph(ped)
        with pytest.raises(PyPedalUsageError, match="find_ancestors_g"):
            pyp_network.find_ancestors_g(graph, 5, {}, "three")
        with pytest.raises(PyPedalUsageError, match="find_ancestors_g"):
            pyp_network.find_ancestors_g(graph, 5, {}, -1)
        limited = pyp_network.find_ancestors_g(graph, 5, {}, 1)
        assert isinstance(limited, dict)
        assert limited


def test_lacy_already_renumbered_emits_no_deprecation():
    with chdir_tmp():
        ped = load_corpus("new_lacy.ped")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", DeprecationWarning)
            result = pyp_metrics.effective_founders_lacy(ped, output=False)
        matching = [
            item for item in caught
            if issubclass(item.category, DeprecationWarning)
            and "renumber" in str(item.message).lower()
        ]
        assert matching == []
        assert result["fa_effective_founders"] == pytest.approx(LACY_APPENDIX_A_FE)


def test_lacy_auto_renumber_warns_and_keeps_result():
    with chdir_tmp():
        numbered = load_corpus("new_lacy.ped")
        expected = pyp_metrics.effective_founders_lacy(
            numbered, output=False)["fa_effective_founders"]
        ped = load_corpus(
            "new_lacy.ped",
            renumber=False,
            pedigree_is_renumbered=False,
            reorder=False,
        )
        with pytest.warns(DeprecationWarning, match="renumber"):
            result = pyp_metrics.effective_founders_lacy(ped, output=False)
        assert result["fa_effective_founders"] == pytest.approx(expected)
        assert ped.kw.get("pedigree_is_renumbered") is True


def test_theoretical_ne_returns_float_and_respects_output():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        ne_file = ped.kw["filetag"] + "_ne_from_metadata_.dat"
        value = pyp_metrics.theoretical_ne_from_metadata(ped, output=False)
        assert value == pytest.approx(MRODE_NE)
        assert isinstance(value, float)
        assert not os.path.exists(ne_file)
        written = pyp_metrics.theoretical_ne_from_metadata(ped, output=True)
        assert written == pytest.approx(value)
        assert os.path.exists(ne_file)


def test_theoretical_ne_failure_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        ped.metadata.num_unique_sires = 0
        with pytest.raises(PyPedalError, match="theoretical_ne_from_metadata"):
            pyp_metrics.theoretical_ne_from_metadata(ped, output=False)


def test_summary_inbreeding_rejects_non_dict():
    with pytest.raises(PyPedalUsageError, match="summary_inbreeding"):
        pyp_io.summary_inbreeding("0")


def test_min_max_f_invalid_forma_raises():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        with pytest.raises(PyPedalUsageError, match="forma"):
            pyp_metrics.min_max_f(ped, forma="huge")
