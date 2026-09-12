"""Deprecation of legacy ``pedigree_completeness``."""
import warnings

import pytest

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import load_corpus, load_corpus_from_path, write_temp_pedigree


def test_pedigree_completeness_emits_deprecation_warning():
    ped = load_corpus("new_lacy.ped")
    with pytest.warns(DeprecationWarning, match="equivalent_complete_generations"):
        pyp_metrics.pedigree_completeness(ped, 4)


def test_invalid_gens_raises_without_deprecation_warning():
    ped = load_corpus("new_lacy.ped")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        with pytest.raises(PyPedalUsageError):
            pyp_metrics.pedigree_completeness(ped, 0)
    completeness = [
        item for item in caught
        if issubclass(item.category, DeprecationWarning)
        and "pedigree_completeness" in str(item.message)
    ]
    assert completeness == []


def test_default_load_does_not_invoke_completeness():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        ped = load_corpus("mrode.ped")
    completeness = [
        item for item in caught
        if issubclass(item.category, DeprecationWarning)
        and "pedigree_completeness" in str(item.message)
    ]
    assert completeness == []
    assert ped.kw.get("pedcomp") is False
    assert all(float(animal.pedcomp) == float(ped.kw["missing_pedcomp"]) for animal in ped.pedigree)


def test_opt_in_pedcomp_load_emits_deprecation_warning():
    path = write_temp_pedigree(["1 0 0", "2 0 0", "3 1 2"])
    with pytest.warns(DeprecationWarning, match="equivalent_complete_generations"):
        load_corpus_from_path(path, "asd", pedcomp=True, pedcomp_gens=3)
