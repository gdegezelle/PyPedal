"""Lacy founder-row recurrence: independent phantom controls.

Permanent tests compare production ``effective_founders_lacy`` with
``a_effective_founders_lacy`` and the independent Lacy oracle. Public
strict/absorb modes are refused.
"""
import os
import unittest
import warnings

from _pedhelpers import chdir_tmp, corpus, load_corpus, load_corpus_from_path, write_temp_pedigree
from oracles import lacy_f_e

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError


def _write_rows(rows):
    return write_temp_pedigree(rows)


def _load(rows):
    path = _write_rows(rows)
    try:
        return load_corpus_from_path(path, "asd")
    finally:
        os.remove(path)


def _f_e(fn, ped, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        with chdir_tmp():
            return fn(ped, **kw)["fa_effective_founders"]


class TestLacyAppendixAUnchanged(unittest.TestCase):

    def test_new_lacy_appendix_a(self):
        ped = load_corpus("new_lacy.ped")
        for fn in (pyp_metrics.effective_founders_lacy, pyp_metrics.a_effective_founders_lacy):
            self.assertAlmostEqual(2.909090909090909, _f_e(fn, ped), places=9)
        want, gate = lacy_f_e(corpus("new_lacy.ped"))
        self.assertTrue(gate)
        self.assertAlmostEqual(2.909090909090909, want, places=9)


class TestLacyRecurrenceMatchesOracleAndPopulation(unittest.TestCase):

    CASES = (
        ("two_equal_founders", ["1 0 0", "2 0 0", "3 1 2"]),
        ("skewed_founders", ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 1 3"]),
        ("one_parent_branch", ["1 0 0", "2 1 0", "3 2 0"]),
        ("half_founder", ["1 0 0", "2 0 0", "3 1 0", "4 3 2"]),
        ("inbred_descendant", ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"]),
        (
            "multiple_founder_paths",
            ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"],
        ),
    )

    def test_synthetics_agree_with_oracle_and_a_effective(self):
        for name, rows in self.CASES:
            path = _write_rows(rows)
            try:
                ped = load_corpus_from_path(path, "asd")
                family = _f_e(pyp_metrics.effective_founders_lacy, ped)
                population = _f_e(pyp_metrics.a_effective_founders_lacy, ped)
                want, gate = lacy_f_e(path)
                with self.subTest(name=name):
                    self.assertTrue(gate, "%s: oracle q must sum to one" % name)
                    self.assertAlmostEqual(want, family, places=9)
                    self.assertAlmostEqual(want, population, places=9)
            finally:
                os.remove(path)


class TestLacyModesAndPhantomSources(unittest.TestCase):

    def _half_founder_ped(self):
        return _load(["1 0 0", "2 0 0", "3 1 0", "4 3 2"])

    def test_phantom_includes_unknown_parent_as_source(self):
        ped = self._half_founder_ped()
        phantom = _f_e(pyp_metrics.effective_founders_lacy, ped, mode="phantom")
        default = _f_e(pyp_metrics.a_effective_founders_lacy, ped)
        self.assertAlmostEqual(phantom, default, places=9)
        slots = pyp_metrics.lacy_phantom_slots(ped)
        self.assertEqual(1, len(slots))
        self.assertEqual("dam", slots[0][1])

    def test_strict_and_absorb_are_refused_on_both_entry_points(self):
        ped = self._half_founder_ped()
        for fn in (pyp_metrics.effective_founders_lacy,
                   pyp_metrics.a_effective_founders_lacy):
            for mode in ("strict", "absorb"):
                with self.subTest(routine=fn.__name__, mode=mode):
                    with self.assertRaises(PyPedalUsageError):
                        _f_e(fn, ped, mode=mode)

    def test_phantom_matches_independent_oracle(self):
        rows = ["1 0 0", "2 0 0", "3 1 0", "4 3 2"]
        path = _write_rows(rows)
        try:
            ped = load_corpus_from_path(path, "asd")
            family = _f_e(pyp_metrics.effective_founders_lacy, ped)
            want, gate = lacy_f_e(path)
            self.assertTrue(gate)
            self.assertAlmostEqual(want, family, places=9)
        finally:
            os.remove(path)

    def test_corpus_modes_agree_with_existing_controls(self):
        ped = load_corpus("mrode.ped")
        self.assertAlmostEqual(
            2.797814208, _f_e(pyp_metrics.effective_founders_lacy, ped), places=6
        )
        ped = load_corpus("hartlandclark.ped")
        self.assertAlmostEqual(
            5.831988609, _f_e(pyp_metrics.effective_founders_lacy, ped), places=6
        )
        ped = load_corpus("generations.ped", "asdbx")
        self.assertAlmostEqual(
            4.612612613, _f_e(pyp_metrics.effective_founders_lacy, ped), places=6
        )

    def test_inbred_descendant_does_not_drop_founder_paths(self):
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"]
        path = _write_rows(rows)
        try:
            ped = load_corpus_from_path(path, "asd")
            family = _f_e(pyp_metrics.effective_founders_lacy, ped)
            want, gate = lacy_f_e(path)
            self.assertTrue(gate)
            self.assertAlmostEqual(want, family, places=9)
        finally:
            os.remove(path)


if __name__ == "__main__":
    unittest.main()
