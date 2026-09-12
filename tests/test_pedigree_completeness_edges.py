"""Edge-case validation for legacy ``pedigree_completeness``."""
import unittest
import warnings

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import load_corpus_from_path, write_temp_pedigree


def _legacy_completeness(ped, gens=4):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        return pyp_metrics.pedigree_completeness(ped, gens)


def _load_rows(rows):
    return load_corpus_from_path(write_temp_pedigree(rows), "asd")


class TestCompletenessGensValidation(unittest.TestCase):
    def test_gens_zero_raises(self):
        ped = _load_rows(["1 0 0", "2 0 0", "3 1 2"])
        with self.assertRaises(PyPedalUsageError):
            _legacy_completeness(ped, 0)

    def test_negative_gens_raises(self):
        ped = _load_rows(["1 0 0", "2 0 0", "3 1 2"])
        with self.assertRaises(PyPedalUsageError):
            _legacy_completeness(ped, -1)


class TestFoundersOnlySummary(unittest.TestCase):
    def test_empty_nonfounder_statistics_are_zero_not_sentinels(self):
        ped = _load_rows(["1 0 0", "2 0 0"])
        summary = _legacy_completeness(ped, 4)
        self.assertEqual(2, summary["n"])
        self.assertEqual(0, summary["nonfounder_n"])
        self.assertEqual(0.0, summary["nonfounder_sum"])
        self.assertEqual(0.0, summary["nonfounder_min"])
        self.assertEqual(0.0, summary["nonfounder_max"])
        self.assertEqual(0.0, summary["nonfounder_range"])
        self.assertEqual(0.0, summary["nonfounder_average"])
        self.assertGreaterEqual(summary["nonfounder_range"], 0.0)


if __name__ == "__main__":
    unittest.main()
