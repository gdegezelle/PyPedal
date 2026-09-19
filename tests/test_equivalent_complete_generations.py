"""Equivalent Complete Generations (Maignel, Boichard & Verrier 1996, p.50).

Expected values are hand-checked slot sums ``Σ known_slots(g) / 2^g``.
They are not taken from the production recurrence.
"""
import hashlib
import os
import statistics
import struct
import unittest
import warnings
from collections import deque

import pytest

from PyPedal import pyp_errors, pyp_metrics

from _pedhelpers import (
    load_canonical_griffon,
    load_corpus_from_path,
    write_temp_pedigree,
)


def _load(rows, **overrides):
    return load_corpus_from_path(write_temp_pedigree(rows), "asd", **overrides)


def _ecg(ped):
    return pyp_metrics.equivalent_complete_generations(ped, output=False)


def _by_original(ped, original_id):
    for animal in ped.pedigree:
        if int(animal.originalID) == int(original_id):
            return animal
    raise AssertionError("original ID %s is not in the pedigree" % original_id)


def _oracle_ecg(ped, animal_id):
    """TEST-ONLY Maignel slot enumeration. Does not call production code."""
    missing = ped.kw["missing_parent"]
    by_id = {int(animal.animalID): animal for animal in ped.pedigree}

    def _known(parent_id):
        return not (parent_id == missing or str(parent_id) == str(missing))

    total = 0.0
    start = by_id[int(animal_id)]
    queue = deque()
    if _known(start.sireID):
        queue.append((int(start.sireID), 1))
    if _known(start.damID):
        queue.append((int(start.damID), 1))
    while queue:
        parent_id, depth = queue.popleft()
        total += 1.0 / (2 ** depth)
        parent = by_id[parent_id]
        if _known(parent.sireID):
            queue.append((int(parent.sireID), depth + 1))
        if _known(parent.damID):
            queue.append((int(parent.damID), depth + 1))
    return total


class TestHandCheckableECG(unittest.TestCase):
    def test_founder_is_zero(self):
        ped = _load(["1 0 0"])
        self.assertEqual(0.0, _ecg(ped)[int(_by_original(ped, 1).animalID)])

    def test_half_founder_with_known_founder_parent_is_half(self):
        ped = _load(["1 0 0", "2 1 0"])
        self.assertEqual(0.5, _ecg(ped)[int(_by_original(ped, 2).animalID)])

    def test_two_known_founder_parents_is_one(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2"])
        self.assertEqual(1.0, _ecg(ped)[int(_by_original(ped, 3).animalID)])

    def test_two_complete_generations_is_two(self):
        ped = _load([
            "1 0 0", "2 0 0", "3 0 0", "4 0 0",
            "5 1 2", "6 3 4",
            "7 5 6",
        ])
        self.assertEqual(2.0, _ecg(ped)[int(_by_original(ped, 7).animalID)])

    def test_half_founder_whose_known_parent_has_both_parents_is_one(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2", "4 3 0"])
        self.assertEqual(1.0, _ecg(ped)[int(_by_original(ped, 4).animalID)])

    def test_full_sib_collapse_counts_repeated_slots(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"])
        self.assertEqual(2.0, _ecg(ped)[int(_by_original(ped, 5).animalID)])

    def test_three_complete_generations_is_three(self):
        ped = _load([
            "1 0 0", "2 0 0", "3 0 0", "4 0 0",
            "5 0 0", "6 0 0", "7 0 0", "8 0 0",
            "9 1 2", "10 3 4", "11 5 6", "12 7 8",
            "13 9 10", "14 11 12",
            "15 13 14",
        ])
        self.assertEqual(3.0, _ecg(ped)[int(_by_original(ped, 15).animalID)])

    def test_deeper_incomplete_case_is_one_and_three_quarters(self):
        ped = _load(["1 0 0", "2 0 0", "3 0 0", "4 1 2", "5 3 0", "6 4 5"])
        self.assertEqual(1.75, _ecg(ped)[int(_by_original(ped, 6).animalID)])


class TestSlotVersusIdentity(unittest.TestCase):
    def test_full_sib_ecg_uses_slots_and_differs_from_legacy_completeness(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"])
        mapping = _ecg(ped)
        animal = _by_original(ped, 5)
        self.assertEqual(2.0, mapping[int(animal.animalID)])
        self.assertEqual(2.0, float(animal.ecg))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            pyp_metrics.pedigree_completeness(ped, 3)
        self.assertNotEqual(2.0, float(animal.pedcomp))
        self.assertLess(float(animal.pedcomp), 1.0)


class TestIndependentSlotOracle(unittest.TestCase):
    CASES = (
        (["1 0 0"], 1),
        (["1 0 0", "2 1 0"], 2),
        (["1 0 0", "2 0 0", "3 1 2"], 3),
        (["1 0 0", "2 0 0", "3 0 0", "4 0 0", "5 1 2", "6 3 4", "7 5 6"], 7),
        (["1 0 0", "2 0 0", "3 1 2", "4 3 0"], 4),
        (["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"], 5),
        (
            [
                "1 0 0", "2 0 0", "3 0 0", "4 0 0",
                "5 0 0", "6 0 0", "7 0 0", "8 0 0",
                "9 1 2", "10 3 4", "11 5 6", "12 7 8",
                "13 9 10", "14 11 12",
                "15 13 14",
            ],
            15,
        ),
        (["1 0 0", "2 0 0", "3 0 0", "4 1 2", "5 3 0", "6 4 5"], 6),
    )

    def test_oracle_agrees_with_production(self):
        for rows, original_id in self.CASES:
            with self.subTest(original_id=original_id, rows=rows):
                ped = _load(rows)
                animal = _by_original(ped, original_id)
                got = _ecg(ped)[int(animal.animalID)]
                want = _oracle_ecg(ped, animal.animalID)
                self.assertEqual(want, got)


class TestECGApi(unittest.TestCase):
    def test_keys_are_current_animal_ids_and_values_are_stored(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2"])
        mapping = _ecg(ped)
        self.assertEqual({int(animal.animalID) for animal in ped.pedigree}, set(mapping))
        for animal in ped.pedigree:
            self.assertEqual(mapping[int(animal.animalID)], float(animal.ecg))

    def test_output_false_does_not_write_a_dat_file(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2"])
        path = ped.kw["filetag"] + "_ecg_.dat"
        _ecg(ped)
        self.assertFalse(os.path.exists(path))

    def test_output_true_writes_dat(self):
        ped = _load(["1 0 0", "2 0 0", "3 1 2"])
        pyp_metrics.equivalent_complete_generations(ped, output=True)
        self.assertTrue(os.path.exists(ped.kw["filetag"] + "_ecg_.dat"))

    def test_non_topological_pedigree_raises(self):
        ped = _load(
            ["3 1 2", "1 0 0", "2 0 0"],
            reorder=False,
            renumber=False,
            pedigree_is_renumbered=False,
        )
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            _ecg(ped)
        self.assertIn("does not appear earlier", str(caught.exception))


CANONICAL_ECG_N = 97002
CANONICAL_ECG_MEAN = 15.49766453668286
CANONICAL_ECG_MEDIAN = 17.069160033259667
CANONICAL_ECG_MIN = 0.0
CANONICAL_ECG_MAX = 26.542230867556363
CANONICAL_ECG_ZERO_COUNT = 6659
CANONICAL_ECG_VECTOR_SHA256 = (
    "ec6d177dcdf3eeef398052a1d2337cdf11d0bf2b1fa710cf88073fea5dca51ee"
)
CANONICAL_ECG_SELECTED = {
    98685: 25.12890448849202,
    98667: 25.151318398479063,
    37482: 23.166416324607525,
    54587: 22.167251976010363,
    20196: 0.0,
    20209: 11.512431582735644,
}


def _ecg_vector_sha256(values):
    digest = hashlib.sha256()
    digest.update(struct.pack(">Q", len(values)))
    for value in values:
        digest.update(struct.pack(">d", value))
    return digest.hexdigest()


def _neumaier_mean(values):
    total = 0.0
    correction = 0.0
    for value in values:
        t = total + value
        if abs(total) >= abs(value):
            correction += (total - t) + value
        else:
            correction += (value - t) + total
        total = t
    return (total + correction) / float(len(values))


@pytest.mark.integration
def test_canonical_griffon_ecg_summary():
    ped = load_canonical_griffon()
    mapping = _ecg(ped)
    assert len(ped.pedigree) == CANONICAL_ECG_N
    assert len(mapping) == CANONICAL_ECG_N
    values = [mapping[int(animal.animalID)] for animal in ped.pedigree]
    missing = ped.kw["missing_parent"]
    founders = sum(
        1
        for animal in ped.pedigree
        if (animal.sireID == missing or str(animal.sireID) == str(missing))
        and (animal.damID == missing or str(animal.damID) == str(missing))
    )
    zeros = sum(1 for value in values if value == 0.0)
    assert zeros == CANONICAL_ECG_ZERO_COUNT
    assert founders == CANONICAL_ECG_ZERO_COUNT
    assert min(values) == CANONICAL_ECG_MIN
    assert max(values) == CANONICAL_ECG_MAX
    assert _neumaier_mean(values) == CANONICAL_ECG_MEAN
    assert statistics.median(values) == CANONICAL_ECG_MEDIAN
    positional = [mapping[i] for i in range(1, CANONICAL_ECG_N + 1)]
    assert _ecg_vector_sha256(positional) == CANONICAL_ECG_VECTOR_SHA256
    by_oid = {
        int(animal.originalID): mapping[int(animal.animalID)]
        for animal in ped.pedigree
    }
    for original_id, expected in CANONICAL_ECG_SELECTED.items():
        assert by_oid[original_id] == expected


if __name__ == "__main__":
    unittest.main()
