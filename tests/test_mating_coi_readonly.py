"""
 successor -- ``mating_coi`` is a read-only float; ``addanimal``
must not leave a failed append behind.

WHAT THIS FILE IS ABOUT
-----------------------
The historical mating-COI path appended a hypothetical offspring and
could return ``-999.9`` while leaking that animal. RC1 refused both
mating APIs. RC4 an earlier revision implements ``mating_coi`` as relationship
mathematics: a float, or ``PyPedalUsageError``, and zero pedigree
mutation. ``addanimal`` still restores call-local state on historical
failure.

The scientific pairwise controls live in ``test_mating_coi.py``. This
file keeps the RC1 regression intent: no phantom leak, no ``-999.9``,
and ``addanimal`` rollback.
"""
import copy
import os
import unittest

from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import owned_temp_dir, chdir_tmp, load_corpus, load_example, load_corpus_from_path, load_griffon_1871_1890

BASELINE = "bb2ee3a"

# Distinctive originals so animalID 1..n cannot be mistaken for identity.
STUD_ROWS = [
    "100 0 0",
    "200 0 0",
    "300 100 200",
]


def studbook(**overrides):
    tmp = owned_temp_dir(prefix="rc1s1_")
    path = os.path.join(tmp, "stud.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(STUD_ROWS) + "\n")
    return load_corpus_from_path(path, "asd", **overrides)


def snapshot(ped):
    return {
        "n": len(ped.pedigree),
        "objects": [id(animal) for animal in ped.pedigree],
        "animal_ids": [int(animal.animalID) for animal in ped.pedigree],
        "original_ids": [animal.originalID for animal in ped.pedigree],
        "sire_ids": [animal.sireID for animal in ped.pedigree],
        "dam_ids": [animal.damID for animal in ped.pedigree],
        "renumbered_ids": [animal.renumberedID for animal in ped.pedigree],
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "namemap": dict(ped.namemap),
        "namebackmap": dict(ped.namebackmap),
        "offspring": [
            (
                animal.originalID,
                dict(animal.sons),
                dict(animal.daus),
                dict(animal.unks),
            )
            for animal in ped.pedigree
        ],
        "meta": ped.metadata.num_records,
        "kw": copy.deepcopy(dict(ped.kw)),
        "has_nrm": getattr(ped, "nrm", None) is not None,
    }


def assert_same_state(test, before, ped):
    test.assertEqual(before, snapshot(ped))


# ===========================================================================
# Anti-vacuity and current-master characterisation (green on bb2ee3a)
# ===========================================================================
class TestFixtureIsNotVacuous(unittest.TestCase):
    def test_studbook_original_and_animal_ids_differ(self):
        ped = studbook()
        oids = [int(a.originalID) for a in ped.pedigree]
        aids = [int(a.animalID) for a in ped.pedigree]
        self.assertEqual([100, 200, 300], oids)
        self.assertEqual([1, 2, 3], aids)
        self.assertNotEqual(oids, aids)
        self.assertEqual({100: 1, 200: 2, 300: 3}, ped.idmap)

    def test_griffons_are_a_normal_renumbered_studbook(self):
        ped = load_griffon_1871_1890()
        self.assertGreaterEqual(len(ped.pedigree), 3)
        self.assertTrue(
            any(int(a.originalID) != int(a.animalID) for a in ped.pedigree)
        )


class TestCurrentMasterLeakInverted(unittest.TestCase):
    """The bb2ee3a sentinel leak is gone. Characterisation inverted."""

    def test_mrode_compatible_ids_return_float_without_append(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped, 0)
        self.assertEqual(0.0, got)
        self.assertIsInstance(got, float)
        self.assertNotEqual(-999.9, got)
        assert_same_state(self, before, ped)

    def test_studbook_valid_ids_compute_without_append(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped, 0)
        self.assertEqual(0.0, got)
        assert_same_state(self, before, ped)

    def test_invalid_ids_raise_usage_error_without_append(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(99999, 88888, ped, 0)
        assert_same_state(self, before, ped)

    def test_group_now_computes_without_append(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi_group(["1_2"], ped, names=0, gens=0)
        self.assertEqual(0.0, got["matings"][(1, 2)])
        assert_same_state(self, before, ped)

    def test_direct_addanimal_renumbered_parents_no_longer_leaks(self):
        """ADD-1 inverted — historical failure no longer leaves the append."""
        ped = studbook()
        before = snapshot(ped)
        new_id = max(ped.idmap.keys()) + 1
        ok = ped.addanimal(new_id, 1, 2)
        self.assertFalse(ok)
        assert_same_state(self, before, ped)


class TestSuccessfulAddanimalUnchanged(unittest.TestCase):
    """ADD-3 — green on the baseline and after the repair."""

    def test_original_id_parents_still_append_and_map(self):
        ped = studbook()
        before = snapshot(ped)
        ok = ped.addanimal(400, 100, 200)
        self.assertTrue(ok)
        self.assertEqual(before["n"] + 1, len(ped.pedigree))
        added = ped.pedigree[-1]
        self.assertEqual(400, added.originalID)
        self.assertEqual(4, added.animalID)
        self.assertEqual(1, added.sireID)
        self.assertEqual(2, added.damID)
        self.assertEqual(4, ped.idmap[400])
        self.assertEqual(400, ped.backmap[4])
        self.assertEqual(before["objects"], [id(a) for a in ped.pedigree[: before["n"]]])
        self.assertEqual(before["meta"], ped.metadata.num_records)
        self.assertEqual("loader", ped.kw["newanimal_caller"])


# ===========================================================================
# Desired 4.0.0-rc1 contract. Was xfail(strict=True) on bb2ee3a.
# ===========================================================================
class TestMC1ValidIdsReturnFloat(unittest.TestCase):
    def test_mating_coi_valid_mrode_ids_return_float(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped)
        self.assertEqual(0.0, got)
        self.assertIsInstance(got, float)
        assert_same_state(self, before, ped)


class TestMC2InvalidIdsRefuse(unittest.TestCase):
    def test_mating_coi_invalid_ids_raise_not_sentinel(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(99999, 88888, ped)
        assert_same_state(self, before, ped)


class TestMC3GensZeroSucceeds(unittest.TestCase):
    def test_explicit_gens_zero_computes(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 6, ped, 0)
        self.assertAlmostEqual(0.125, got, places=12)


class TestMC4GensMinusOneSucceeds(unittest.TestCase):
    def test_gens_minus_one_matches_gens_zero(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped, -1)
        self.assertEqual(0.0, got)
        assert_same_state(self, before, ped)


class TestMC5GroupComputesViaMatingCoi(unittest.TestCase):
    def test_group_calls_mating_coi_for_supplied_pairs(self):
        ped = studbook()
        before = snapshot(ped)
        called = []
        original = pyp_metrics.mating_coi

        def wrapped(*args, **kwargs):
            called.append((args, kwargs))
            return original(*args, **kwargs)

        pyp_metrics.mating_coi = wrapped
        try:
            with chdir_tmp():
                got = pyp_metrics.mating_coi_group(["1_2"], ped, names=0, gens=0)
        finally:
            pyp_metrics.mating_coi = original
        self.assertEqual(1, len(called))
        self.assertEqual(0.0, got["matings"][(1, 2)])
        assert_same_state(self, before, ped)


class TestMC6SuccessLeavesState(unittest.TestCase):
    def test_full_snapshot_unchanged_on_success(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            pyp_metrics.mating_coi(1, 2, ped, 0)
        assert_same_state(self, before, ped)


class TestMC7StudbookComputes(unittest.TestCase):
    def test_griffons_compute_and_do_not_grow(self):
        ped = load_griffon_1871_1890()
        before = snapshot(ped)
        ids = [int(a.animalID) for a in ped.pedigree]
        with chdir_tmp():
            got = pyp_metrics.mating_coi(ids[0], ids[1], ped, 0)
        self.assertIsInstance(got, float)
        self.assertNotEqual(-999.9, got)
        assert_same_state(self, before, ped)


class TestMC8RepeatedSuccess(unittest.TestCase):
    def test_three_successes_leave_constant_state(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            for _ in range(3):
                pyp_metrics.mating_coi(1, 2, ped, 0)
        assert_same_state(self, before, ped)


class TestMC9RepeatedGroupSuccess(unittest.TestCase):
    def test_three_group_calls_leave_constant_state(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            for _ in range(3):
                pyp_metrics.mating_coi_group(["1_2", "1_3"], ped)
        assert_same_state(self, before, ped)


class TestMC10NoSentinel(unittest.TestCase):
    def test_neither_entry_returns_the_numeric_sentinel(self):
        ped = studbook()
        with chdir_tmp():
            self.assertNotEqual(-999.9, pyp_metrics.mating_coi(1, 2, ped, 0))
            self.assertNotEqual(-999.9, pyp_metrics.mating_coi(1, 2, ped, -1))
            group = pyp_metrics.mating_coi_group(["1_2"], ped)
        self.assertNotEqual(-999.9, group["matings"][(1, 2)])
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(99999, 88888, ped, 0)


class TestADD2FailedAddanimalRestored(unittest.TestCase):
    def test_historical_failure_leaves_pre_call_state(self):
        ped = studbook()
        before = snapshot(ped)
        ok = ped.addanimal(max(ped.idmap.keys()) + 1, 1, 2)
        self.assertFalse(ok)
        assert_same_state(self, before, ped)


class TestADD4NoLeakIntoInbreeding(unittest.TestCase):
    def test_failed_addanimal_then_inbreeding_matches_fresh_load(self):
        dirty = studbook()
        clean = studbook()
        self.assertFalse(dirty.addanimal(max(dirty.idmap.keys()) + 1, 1, 2))
        with chdir_tmp():
            got = pyp_nrm.inbreeding(dirty, method="tabular", output=False)
            expect = pyp_nrm.inbreeding(clean, method="tabular", output=False)
        self.assertEqual(expect["fx"], got["fx"])
        self.assertEqual(len(clean.pedigree), len(dirty.pedigree))


class TestADD5RollbackIncludesMaps(unittest.TestCase):
    def test_maps_kw_metadata_offspring_match_pre_call(self):
        ped = studbook()
        before = snapshot(ped)
        self.assertFalse(ped.addanimal(max(ped.idmap.keys()) + 1, 1, 2))
        after = snapshot(ped)
        self.assertEqual(before["idmap"], after["idmap"])
        self.assertEqual(before["backmap"], after["backmap"])
        self.assertEqual(before["namemap"], after["namemap"])
        self.assertEqual(before["namebackmap"], after["namebackmap"])
        self.assertEqual(before["offspring"], after["offspring"])
        self.assertEqual(before["meta"], after["meta"])
        self.assertEqual(before["kw"], after["kw"])
        self.assertEqual(before["n"], after["n"])
        self.assertEqual(before["objects"], after["objects"])
