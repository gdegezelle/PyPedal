"""RC4 an earlier revision -- ``mating_coi_group`` explicit pair evaluation."""
import unittest
import warnings

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import chdir_tmp, load_corpus
from test_mating_coi import snapshot, studbook, write_pedigree


class TestModernPairInput(unittest.TestCase):
    def test_mrode_modern_pairs(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group([(1, 2), (1, 6)], ped)
        self.assertEqual([(1, 2), (1, 6)], list(got["matings"]))
        self.assertEqual(0.0, got["matings"][(1, 2)])
        self.assertAlmostEqual(0.125, got["matings"][(1, 6)], places=12)
        meta = got["metadata"]["all"]
        self.assertEqual(2, meta["count"])
        self.assertAlmostEqual(0.125, meta["sum"], places=12)
        self.assertEqual(0.0, meta["min"])
        self.assertAlmostEqual(0.125, meta["max"], places=12)
        self.assertAlmostEqual(0.125, meta["range"], places=12)
        self.assertAlmostEqual(0.0625, meta["mean"], places=12)
        nz = got["metadata"]["nonzero"]
        self.assertEqual(1, nz["count"])
        self.assertAlmostEqual(0.125, nz["mean"], places=12)

    def test_list_pairs_are_accepted(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group([[1, 3]], ped)
        self.assertAlmostEqual(0.25, got["matings"][(1, 3)], places=12)


class TestLegacyStringInput(unittest.TestCase):
    def test_legacy_numeric_strings(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group(["1_2", "1_6"], ped)
        self.assertEqual([(1, 2), (1, 6)], list(got["matings"]))
        self.assertEqual(0.0, got["matings"][(1, 2)])
        self.assertAlmostEqual(0.125, got["matings"][(1, 6)], places=12)


class TestDuplicateAndReversed(unittest.TestCase):
    def test_duplicate_exact_pair_is_evaluated_once(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group([(1, 6), (1, 2), (1, 6)], ped)
        self.assertEqual([(1, 6), (1, 2)], list(got["matings"]))
        self.assertEqual(2, got["metadata"]["all"]["count"])
        self.assertAlmostEqual(0.125, got["metadata"]["all"]["sum"], places=12)

    def test_reversed_pairs_are_distinct_keys_with_equal_coi(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group([(1, 6), (6, 1)], ped)
        self.assertEqual([(1, 6), (6, 1)], list(got["matings"]))
        self.assertAlmostEqual(
            got["matings"][(1, 6)], got["matings"][(6, 1)], places=12
        )
        self.assertEqual(2, got["metadata"]["all"]["count"])
        self.assertAlmostEqual(0.25, got["metadata"]["all"]["sum"], places=12)


class TestEmptyAndNonzeroMetadata(unittest.TestCase):
    def test_empty_group_is_structurally_valid(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        got = pyp_metrics.mating_coi_group([], ped)
        self.assertEqual({}, got["matings"])
        for block in (got["metadata"]["all"], got["metadata"]["nonzero"]):
            self.assertEqual(0, block["count"])
            self.assertEqual(0.0, block["sum"])
            self.assertIsNone(block["min"])
            self.assertIsNone(block["max"])
            self.assertIsNone(block["range"])
            self.assertIsNone(block["mean"])
        self.assertEqual(before, snapshot(ped))

    def test_all_zero_matings_leave_nonzero_extrema_none(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group([(1, 2), (2, 4)], ped)
        self.assertEqual(2, got["metadata"]["all"]["count"])
        self.assertEqual(0.0, got["metadata"]["all"]["sum"])
        nz = got["metadata"]["nonzero"]
        self.assertEqual(0, nz["count"])
        self.assertEqual(0.0, nz["sum"])
        self.assertIsNone(nz["min"])
        self.assertIsNone(nz["max"])
        self.assertIsNone(nz["range"])
        self.assertIsNone(nz["mean"])


class TestMalformedInput(unittest.TestCase):
    def test_short_tuple_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group([(1,)], ped)
        self.assertEqual(before, snapshot(ped))

    def test_long_tuple_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group([(1, 2, 3)], ped)

    def test_legacy_single_token_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group(["1"], ped)

    def test_legacy_three_token_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group(["1_2_3"], ped)

    def test_invalid_ids_raise(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group([(1, 99999)], ped)
        self.assertEqual(before, snapshot(ped))

    def test_unsupported_gens_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group([(1, 2)], ped, gens=1)
        self.assertEqual(before, snapshot(ped))


class TestNamesCompatibility(unittest.TestCase):
    def test_names_one_resolves_unique_string_ids(self):
        ped = write_pedigree(
            "alpha 0 0\nbeta 0 0\ngamma alpha beta\n", pedformat="ASD"
        )
        self.assertIn("alpha", ped.namemap)
        with chdir_tmp():
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                got = pyp_metrics.mating_coi_group(
                    [("alpha", "beta")], ped, names=1
                )
        self.assertTrue(any(item.category is DeprecationWarning for item in caught))
        pair = next(iter(got["matings"]))
        self.assertEqual(0.0, got["matings"][pair])
        self.assertEqual(
            (
                int(ped.idmap[ped.namemap["alpha"]]),
                int(ped.idmap[ped.namemap["beta"]]),
            ),
            pair,
        )

    def test_names_one_does_not_search_call_name(self):
        ped = write_pedigree(
            "1 0 0 Freddy\n2 0 0 Sally\n3 1 2 Kid\n", pedformat="asdn"
        )
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi_group([("Freddy", "Sally")], ped, names=1)
        self.assertEqual(before, snapshot(ped))


class TestGroupZeroMutation(unittest.TestCase):
    def test_group_and_self_leave_state(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            pyp_metrics.mating_coi_group([(1, 2), (5, 5), (1, 6)], ped)
        self.assertEqual(before, snapshot(ped))

    def test_studbook_group_does_not_grow(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group(["1_2", "1_3"], ped)
        self.assertEqual(0.0, got["matings"][(1, 2)])
        self.assertAlmostEqual(0.25, got["matings"][(1, 3)], places=12)
        self.assertEqual(before, snapshot(ped))
