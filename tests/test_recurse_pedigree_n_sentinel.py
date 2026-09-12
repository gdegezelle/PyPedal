"""Missing-parent sentinel handling in ``recurse_pedigree_n``.

``pedigree_completeness`` traces each parental side by calling
``recurse_pedigree_n`` with that parent's ID. For a half-founder one of
those IDs is the missing-parent sentinel. The append was already guarded,
but the subsequent ``pedigree[anid - 1]`` lookup was not, so sentinel 0
read ``pedigree[-1]`` (the last real animal) and contaminated completeness.
"""
import unittest
import warnings

from PyPedal import pyp_metrics, pyp_nrm

from _pedhelpers import load_corpus, load_corpus_from_path, write_temp_pedigree

HALF_FOUNDER_ROWS = [
    "1 0 0",
    "2 0 0",
    "3 1 2",
    "4 3 0",
    "5 3 2",
]

# Same genealogy, different source order. After a normal load the half-founder
# is original ID 4; the last source row is no longer animal 5.
HALF_FOUNDER_REORDERED_ROWS = [
    "5 3 2",
    "4 3 0",
    "3 1 2",
    "2 0 0",
    "1 0 0",
]


def _legacy_completeness(ped, gens=4):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        return pyp_metrics.pedigree_completeness(ped, gens)


def _load_rows(rows, **overrides):
    return load_corpus_from_path(write_temp_pedigree(rows), "asd", **overrides)


def _by_original(ped, original_id):
    for animal in ped.pedigree:
        if int(animal.originalID) == int(original_id):
            return animal
    raise AssertionError("original ID %s is not in the pedigree" % original_id)


class TestSentinelCallReturnsEmpty(unittest.TestCase):
    def test_missing_parent_contributes_nothing(self):
        ped = _load_rows(HALF_FOUNDER_ROWS)
        missing = ped.kw["missing_parent"]
        got = pyp_nrm.recurse_pedigree_n(ped, missing, [], 4)
        self.assertEqual([], got)


class TestHalfFounderCompleteness(unittest.TestCase):
    def test_animal_4_at_gens_4_is_one_tenth(self):
        ped = _load_rows(HALF_FOUNDER_ROWS)
        _legacy_completeness(ped, 4)
        animal = _by_original(ped, 4)
        self.assertEqual(0.1, float(animal.pedcomp))
        self.assertEqual(3 / 30, float(animal.pedcomp))

    def test_both_parents_known_siblings_are_unchanged(self):
        ped = _load_rows(HALF_FOUNDER_ROWS)
        _legacy_completeness(ped, 4)
        self.assertEqual(0.0, float(_by_original(ped, 1).pedcomp))
        self.assertEqual(0.0, float(_by_original(ped, 2).pedcomp))
        self.assertEqual(2 / 30, float(_by_original(ped, 3).pedcomp))
        self.assertEqual(4 / 30, float(_by_original(ped, 5).pedcomp))


class TestOrderIndependence(unittest.TestCase):
    def test_source_row_order_does_not_give_the_half_founder_last_item_ancestry(self):
        forward = _load_rows(HALF_FOUNDER_ROWS)
        reversed_rows = _load_rows(HALF_FOUNDER_REORDERED_ROWS)
        _legacy_completeness(forward, 4)
        _legacy_completeness(reversed_rows, 4)
        self.assertEqual(
            float(_by_original(forward, 4).pedcomp),
            float(_by_original(reversed_rows, 4).pedcomp),
        )
        self.assertEqual(0.1, float(_by_original(forward, 4).pedcomp))
        self.assertEqual(0.1, float(_by_original(reversed_rows, 4).pedcomp))
        self.assertEqual(
            float(_by_original(forward, 5).pedcomp),
            float(_by_original(reversed_rows, 5).pedcomp),
        )


class TestHalfFounderFreeRegression(unittest.TestCase):
    """new_lacy.ped has no half-founders; gens=4 values are pinned from 4.2.2."""

    LACY_PEDCOMP = (0.0, 0.0, 0.0, 2 / 30, 2 / 30, 4 / 30, 4 / 30)

    def test_new_lacy_completeness_is_bit_identical_to_4_2_2(self):
        ped = load_corpus("new_lacy.ped")
        summary = _legacy_completeness(ped, 4)
        self.assertEqual(list(self.LACY_PEDCOMP), [float(a.pedcomp) for a in ped.pedigree])
        self.assertEqual(0.4, summary["sum"])
        self.assertEqual(7, summary["n"])
        self.assertEqual(0.0, summary["min"])
        self.assertEqual(4 / 30, summary["max"])
        self.assertEqual(4 / 30, summary["range"])
        self.assertEqual(0.4 / 7, summary["average"])
        self.assertEqual(0.4, summary["nonfounder_sum"])
        self.assertEqual(4, summary["nonfounder_n"])
        self.assertEqual(2 / 30, summary["nonfounder_min"])
        self.assertEqual(4 / 30, summary["nonfounder_max"])
        self.assertEqual(2 / 30, summary["nonfounder_range"])
        self.assertEqual(0.1, summary["nonfounder_average"])


if __name__ == "__main__":
    unittest.main()
