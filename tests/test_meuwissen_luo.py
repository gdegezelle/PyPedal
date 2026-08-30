"""Meuwissen and Luo (1992) inbreeding: independent controls and numbering safety.

Permanent tests compare the production implementation with the independent
tabular / Meuwissen-Luo oracle. They do not embed a previous-release snapshot.
"""
import os
import random
import tempfile
import unittest
from types import SimpleNamespace

from _pedhelpers import chdir_tmp, corpus, load_corpus, load_corpus_from_path
from oracles import (
    inbreeding_meuwissen_luo as oracle_ml,
    inbreeding_tabular as oracle_tabular,
    oracle_inbreeding,
    read_pedigree,
    renumber,
)

from PyPedal import pyp_errors, pyp_nrm


def _fx(ped):
    result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
    by_current = {int(k): float(v) for k, v in result["fx"].items()}
    by_original = {}
    for animal in ped.pedigree:
        by_original[int(animal.originalID)] = by_current[int(animal.animalID)]
    return by_original


def _write_rows(rows):
    tmp = tempfile.mkstemp(suffix=".ped")[1]
    with open(tmp, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return tmp


def _load_asd(rows, **overrides):
    path = _write_rows(rows)
    try:
        return load_corpus_from_path(path, "asd", **overrides)
    finally:
        os.remove(path)


def _assert_matches_oracle(test, ped, path, places=12):
    expected = oracle_inbreeding(path)
    got = _fx(ped)
    test.assertEqual(set(expected), set(got))
    for animal_id, want in expected.items():
        test.assertAlmostEqual(want, got[animal_id], places=places, msg="animal %s" % animal_id)


class TestMeuwissenLuoAgainstIndependentOracle(unittest.TestCase):

    def test_mrode_every_coefficient(self):
        ped = load_corpus("mrode.ped")
        _assert_matches_oracle(self, ped, corpus("mrode.ped"))
        self.assertAlmostEqual(0.125, _fx(ped)[5], places=12)

    def test_unrelated_founders(self):
        rows = ["1 0 0", "2 0 0", "3 0 0"]
        ped = _load_asd(rows)
        got = _fx(ped)
        self.assertEqual({1: 0.0, 2: 0.0, 3: 0.0}, got)

    def test_parent_offspring(self):
        rows = ["1 0 0", "2 0 0", "3 1 2"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)
        self.assertEqual(0.0, _fx(ped)[3])

    def test_half_siblings(self):
        rows = ["1 0 0", "2 0 0", "3 0 0", "4 1 2", "5 1 3"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)

    def test_full_siblings(self):
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)

    def test_inbred_individual(self):
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)
        self.assertGreater(_fx(ped)[5], 0.0)

    def test_one_known_parent(self):
        rows = ["1 0 0", "2 1 0", "3 0 1"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)
        got = _fx(ped)
        self.assertEqual(0.0, got[2])
        self.assertEqual(0.0, got[3])

    def test_repeated_common_ancestors(self):
        rows = [
            "1 0 0",
            "2 0 0",
            "3 1 2",
            "4 1 2",
            "5 3 4",
            "6 3 4",
            "7 5 6",
        ]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)
        got = _fx(ped)
        self.assertGreater(got[5], 0.0)
        self.assertGreater(got[7], got[5])

    def test_deep_ancestry(self):
        rows = ["1 0 0", "2 0 0"]
        for i in range(3, 21):
            rows.append("%d %d %d" % (i, i - 2, i - 1))
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)

    def test_multiple_paths_to_the_same_ancestor(self):
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 1", "6 5 4"]
        ped = _load_asd(rows)
        path = _write_rows(rows)
        try:
            _assert_matches_oracle(self, ped, path)
        finally:
            os.remove(path)

    def test_empty_pedigree_returns_empty_dict(self):
        pedobj = SimpleNamespace(
            pedigree=[],
            kw={"missing_parent": 0, "debug_messages": False},
        )
        self.assertEqual({}, pyp_nrm.inbreeding_meuwissen_luo(pedobj))

    def test_one_animal_founder(self):
        ped = _load_asd(["1 0 0"])
        self.assertEqual({1: 0.0}, _fx(ped))

    def test_seeded_random_acyclic_pedigrees_match_tabular_and_oracle(self):
        for seed in (1, 7, 31, 99):
            rng = random.Random(seed)
            n = 40
            rows = ["1 0 0", "2 0 0"]
            for i in range(3, n + 1):
                sire = rng.choice([0] + list(range(1, i)))
                dam = rng.choice([0] + list(range(1, i)))
                rows.append("%d %d %d" % (i, sire, dam))
            path = _write_rows(rows)
            try:
                ped = load_corpus_from_path(path, "asd")
                got = _fx(ped)
                expected = oracle_inbreeding(path)
                self.assertEqual(set(expected), set(got))
                for animal_id, want in expected.items():
                    self.assertAlmostEqual(
                        want, got[animal_id], places=12,
                        msg="production seed=%s original=%s" % (seed, animal_id),
                    )
                numbered, _back = renumber(read_pedigree(path))
                tabular = oracle_tabular(numbered)
                ml = oracle_ml(numbered)
                for i in range(1, n + 1):
                    self.assertAlmostEqual(
                        tabular[i], ml[i], places=12,
                        msg="oracle seed=%s i=%s" % (seed, i),
                    )
            finally:
                os.remove(path)


class TestMeuwissenLuoNumberingSafety(unittest.TestCase):
    """Unnumbered or parent-after-child lists must not compute silent wrong F."""

    def test_unnumbered_original_ids_raise_usage_error(self):
        ped = _load_asd(
            ["10 0 0", "20 0 0", "30 10 20"],
            renumber=False,
            pedigree_is_renumbered=False,
        )
        with self.assertRaises(pyp_errors.PyPedalUsageError) as raised:
            pyp_nrm.inbreeding_meuwissen_luo(ped)
        self.assertIn("renumber", str(raised.exception).lower())
        self.assertFalse(any(getattr(animal, "fa", 0.0) for animal in ped.pedigree))

    def test_public_inbreeding_meu_luo_raises_on_unnumbered(self):
        ped = _load_asd(
            ["10 0 0", "20 0 0", "30 10 20"],
            renumber=False,
            pedigree_is_renumbered=False,
        )
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            pyp_nrm.inbreeding(ped, method="meu_luo", output=False)

    def test_parent_after_child_raises_usage_error(self):
        animal_child = SimpleNamespace(animalID=1, sireID=2, damID=0)
        animal_parent = SimpleNamespace(animalID=2, sireID=0, damID=0)
        pedobj = SimpleNamespace(
            pedigree=[animal_child, animal_parent],
            kw={"missing_parent": 0, "debug_messages": False},
        )
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            pyp_nrm.inbreeding_meuwissen_luo(pedobj)

    def test_does_not_renumber_the_pedigree(self):
        ped = _load_asd(
            ["10 0 0", "20 0 0", "30 10 20"],
            renumber=False,
            pedigree_is_renumbered=False,
        )
        before = [animal.animalID for animal in ped.pedigree]
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            pyp_nrm.inbreeding_meuwissen_luo(ped)
        self.assertEqual(before, [animal.animalID for animal in ped.pedigree])

    def test_renumbered_mrode_still_computes(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_nrm.inbreeding_meuwissen_luo(ped)
        self.assertAlmostEqual(0.125, got[5], places=12)


if __name__ == "__main__":
    unittest.main()
