"""Invalid-ID handling for ``related_animals`` and ``common_ancestors``.

The missing-parent sentinel is not an animal. Passing it used to walk
``pedigree[-1]`` via ``recurse_pedigree_idonly`` and return real ancestry.
Non-existent IDs were swallowed into ``[]``, which is also the legitimate
no-ancestor result.
"""
import unittest

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import load_corpus_from_path, write_temp_pedigree

ROWS = [
    "1 0 0",
    "2 0 0",
    "3 1 2",
    "4 3 0",
    "5 3 2",
]


def _load(**overrides):
    return load_corpus_from_path(write_temp_pedigree(ROWS), "asd", **overrides)


class TestRelatedAnimalsInvalidIds(unittest.TestCase):
    def test_missing_parent_raises(self):
        ped = _load()
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.related_animals(ped.kw["missing_parent"], ped)

    def test_nonexistent_id_raises(self):
        ped = _load()
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.related_animals(99999, ped)

    def test_valid_call_is_unchanged(self):
        ped = _load()
        self.assertEqual([4, 3, 1, 2], pyp_metrics.related_animals(4, ped))
        self.assertEqual([5, 3, 1, 2], pyp_metrics.related_animals(5, ped))
        self.assertEqual([1], pyp_metrics.related_animals(1, ped))


class TestCommonAncestorsInvalidIds(unittest.TestCase):
    def test_missing_parent_raises(self):
        ped = _load()
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.common_ancestors(4, ped.kw["missing_parent"], ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.common_ancestors(ped.kw["missing_parent"], 4, ped)

    def test_nonexistent_id_raises(self):
        ped = _load()
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.common_ancestors(4, 99999, ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.common_ancestors(99999, 4, ped)

    def test_valid_shared_ancestors_are_unchanged(self):
        ped = _load()
        self.assertEqual([1, 2, 3], pyp_metrics.common_ancestors(4, 5, ped))

    def test_no_common_ancestors_still_returns_empty(self):
        ped = _load()
        self.assertEqual([], pyp_metrics.common_ancestors(1, 2, ped))


if __name__ == "__main__":
    unittest.main()
