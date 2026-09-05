"""
 -- ``NewPedigree.delete_animals()`` must delete the requested
stable identities, not stale list positions.

WHAT THIS FILE IS ABOUT
-----------------------
After a normal load, ``idmap`` maps ``originalID -> animalID`` and
``animalID - 1`` is the current list index. ``delete_animals`` uses that
arithmetic for every ID in the request:

    for animalID in animal_list:
        anidx = self.idmap[animalID] - 1
        del maps...
        del self.pedigree[anidx]

The first deletion shrinks the list. Later ``idmap`` values are still the
pre-deletion animalIDs, so ``anidx`` now names a different animal -- or
falls off the end of the list.

Request identity is the ``idmap`` key, which after load+renumber is
``originalID``. Distinctive originals (10, 20, ..., 999) keep Finding-37
renumbered ``animalID`` values from being mistaken for identity.

The deletion loop does not rewrite parent slots.  now
refuses a referenced-parent deletion before any mutation unless every
referencing animal is also explicitly requested. Dangling success is
no longer a supported outcome. This file still proves 
identity: requested originals are removed, not stale list slots.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with the repaired
batch-delete contract marked ``xfail(strict=True)``. Those markers are
removed here. Characterisation of the stale-index failure is inverted,
not deleted.

``self.subTest`` was kept out of the tests that were xfailed.

WHAT IS *NOT* CLAIMED HERE
--------------------------
's reorder/renumber engine is closed. 's
``missingparent=self.kw['missing_parent']`` line is closed. 
oldsave is closed. This file does not ask ``delete_animals`` to cascade
descendants, to rewrite orphan slots to the sentinel, or to rebuild
offspring lists when ``set_offspring`` is off.
"""
import inspect
import os
import unittest

from PyPedal.pyp_errors import PyPedalPedigreeStructureError, PyPedalUsageError
from PyPedal.pyp_newclasses import NewPedigree

from _pedhelpers import owned_temp_dir, load_corpus_from_path

BASELINE = "f7963ff"

FOUNDER = 10
HALF_SIRE = 20
HALF_DAM = 30
ORDINARY = 40
ORDINARY2 = 50
CHILD = 60
LEAF = 999

ALL = (FOUNDER, HALF_SIRE, HALF_DAM, ORDINARY, ORDINARY2, CHILD, LEAF)
PRE_AID = {
    FOUNDER: 1,
    HALF_SIRE: 2,
    HALF_DAM: 3,
    ORDINARY: 4,
    ORDINARY2: 5,
    CHILD: 6,
    LEAF: 7,
}
PRE_GRAPH = {
    FOUNDER: ("MISSING", "MISSING"),
    HALF_SIRE: (FOUNDER, "MISSING"),
    HALF_DAM: ("MISSING", FOUNDER),
    ORDINARY: (FOUNDER, HALF_SIRE),
    ORDINARY2: (HALF_SIRE, HALF_DAM),
    CHILD: (ORDINARY, HALF_DAM),
    LEAF: (ORDINARY2, CHILD),
}

ROWS = [
    "%s 0 0" % FOUNDER,
    "%s %s 0" % (HALF_SIRE, FOUNDER),
    "%s 0 %s" % (HALF_DAM, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (ORDINARY2, HALF_SIRE, HALF_DAM),
    "%s %s %s" % (CHILD, ORDINARY, HALF_DAM),
    "%s %s %s" % (LEAF, ORDINARY2, CHILD),
]
CUSTOM = -999
CUSTOM_ROWS = [
    "%s %s %s" % (FOUNDER, CUSTOM, CUSTOM),
    "%s %s %s" % (HALF_SIRE, FOUNDER, CUSTOM),
    "%s %s %s" % (HALF_DAM, CUSTOM, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (ORDINARY2, HALF_SIRE, HALF_DAM),
    "%s %s %s" % (CHILD, ORDINARY, HALF_DAM),
    "%s %s %s" % (LEAF, ORDINARY2, CHILD),
]


def rows_to_ped(rows=None, pedformat="asd", **overrides):
    rows = ROWS if rows is None else rows
    tmp = owned_temp_dir(prefix="m6_")
    path = os.path.join(tmp, "m6.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, **overrides)


def isolate(ped):
    """
    Keep load-time renumbering (so idmap values are 1-based animalIDs)
    and disable the post-delete reorder/renumber/offspring path so the
    deletion loop can be scored on its own.
    """
    ped.kw["renumber"] = False
    ped.kw["reorder"] = 0
    ped.kw["set_generations"] = False
    ped.kw["set_ancestors"] = False
    ped.kw["set_sexes"] = False
    ped.kw["assign_sexes"] = False
    ped.kw["set_alleles"] = False
    ped.kw["form_nrm"] = False
    ped.kw["set_offspring"] = False
    return ped


def originals(ped):
    return [int(a.originalID) for a in ped.pedigree]


def stable_graph(ped):
    missing = ped.kw["missing_parent"]
    by_aid = {a.animalID: int(a.originalID) for a in ped.pedigree}
    edges = {}
    for animal in ped.pedigree:
        def parent(pid):
            if pid == missing or str(pid) == str(missing):
                return "MISSING"
            if pid in by_aid:
                return by_aid[pid]
            return ("DANGLING", pid)
        edges[int(animal.originalID)] = (
            parent(animal.sireID), parent(animal.damID))
    return edges


def expected_after(deleted):
    remaining = [oid for oid in ALL if oid not in deleted]
    edges = {}
    for child in remaining:
        sire, dam = PRE_GRAPH[child]

        def slot(parent):
            if parent == "MISSING":
                return "MISSING"
            if parent in deleted:
                return ("DANGLING", PRE_AID[parent])
            return parent

        edges[child] = (slot(sire), slot(dam))
    return set(remaining), edges


def snapshot_ids(ped):
    return {
        "originals": originals(ped),
        "objects": [id(a) for a in ped.pedigree],
        "slots": {
            int(a.originalID): (a.sireID, a.damID, a.animalID, a.renumberedID)
            for a in ped.pedigree
        },
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "n": len(ped.pedigree),
        "meta": ped.metadata.num_records,
    }


def assert_refused_unchanged(test, ped, requested, before, exc=PyPedalPedigreeStructureError):
    with test.assertRaises(exc):
        ped.delete_animals(requested)
    after = snapshot_ids(ped)
    test.assertEqual(before, after)
    for oid in requested:
        if oid in ped.idmap:
            test.assertIn(oid, set(originals(ped)))


def assert_deleted_exactly(test, ped, requested, ok=True):
    remaining, edges = expected_after(set(requested))
    test.assertTrue(ok)
    test.assertEqual(remaining, set(originals(ped)))
    test.assertEqual(edges, stable_graph(ped))
    test.assertEqual(len(remaining), ped.metadata.num_records)
    for oid in requested:
        test.assertNotIn(oid, ped.idmap)
        test.assertNotIn(oid, {int(a.originalID) for a in ped.pedigree})
    for oid in remaining:
        test.assertIn(oid, ped.idmap)


# ===========================================================================
# Anti-vacuity. Green on the broken tree and the repaired tree.
# ===========================================================================
class TestFixtureIsNotVacuous(unittest.TestCase):
    def test_load_renumbers_to_sequential_animal_ids(self):
        ped = rows_to_ped()
        self.assertEqual(list(ALL), originals(ped))
        self.assertEqual(
            [PRE_AID[oid] for oid in ALL],
            [int(a.animalID) for a in ped.pedigree])
        self.assertEqual(PRE_AID, {int(k): int(v) for k, v in ped.idmap.items()})
        self.assertNotEqual(
            [int(a.animalID) for a in ped.pedigree],
            [int(a.originalID) for a in ped.pedigree])

    def test_requested_ids_are_idmap_keys_not_list_indices(self):
        ped = rows_to_ped()
        self.assertIn(LEAF, ped.idmap)
        self.assertNotEqual(LEAF, ped.idmap[LEAF])
        self.assertEqual(7, ped.idmap[LEAF])

    def test_an_earlier_delete_shifts_a_later_list_position(self):
        """
        After removing animalID 2 (original 20) the animal that used to
        sit at index 5 (original 60, animalID 6) sits at index 4, but
        idmap[60] is still 6. That is the arithmetic the loop reuses.
        """
        ped = rows_to_ped()
        isolate(ped)
        anidx = ped.idmap[HALF_SIRE] - 1
        self.assertEqual(1, anidx)
        self.assertEqual(HALF_SIRE, int(ped.pedigree[anidx].originalID))
        del ped.pedigree[anidx]
        self.assertEqual(CHILD, int(ped.pedigree[4].originalID))
        self.assertEqual(LEAF, int(ped.pedigree[5].originalID))
        stale = ped.idmap[CHILD] - 1
        self.assertEqual(5, stale)
        self.assertEqual(LEAF, int(ped.pedigree[stale].originalID))
        self.assertNotEqual(CHILD, int(ped.pedigree[stale].originalID))


class TestSignatureAndIdentity(unittest.TestCase):
    def test_public_signature_is_animal_list(self):
        params = list(inspect.signature(NewPedigree.delete_animals).parameters)
        self.assertEqual(["self", "animal_list"], params)

    def test_pre_delete_graph_matches_the_oracle(self):
        ped = rows_to_ped()
        self.assertEqual(PRE_GRAPH, stable_graph(ped))


# ===========================================================================
# The f7963ff reproducer, inverted. 999 is no longer deleted for [20, 60].
# ===========================================================================
class TestTheReproducerNoLongerReproduces(unittest.TestCase):
    def test_da_repro_b_referenced_parent_pair_is_refused_not_stale_shifted(self):
        """
        Was: request [20, 60] removed 20 and 999, left 60. 
        now refuses because 40 and 50 survive and still name 20.
        999 is not deleted as a stale-index victim.
        """
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [HALF_SIRE, CHILD], before)
        self.assertIn(LEAF, originals(ped))
        self.assertIn(CHILD, originals(ped))

    def test_da_repro_e_open_three_set_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(
            self, ped, [HALF_SIRE, ORDINARY, CHILD], before)
        self.assertIn(ORDINARY2, originals(ped))

    def test_descending_open_parent_set_is_also_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [CHILD, HALF_SIRE], before)

    def test_later_index_after_an_earlier_delete_no_longer_raises(self):
        """Was: [60, 999] IndexError'd on idmap[999]-1 after the list shrank."""
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([CHILD, LEAF])
        assert_deleted_exactly(self, ped, [CHILD, LEAF], ok=ok)


# ===========================================================================
# DA-1, DA-9, DA-10 -- single deletion already works.
# ===========================================================================
class TestSingleDeletionControl(unittest.TestCase):
    def test_da1_single_leaf_removes_exactly_the_requested_identity(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF])
        assert_deleted_exactly(self, ped, [LEAF], ok=ok)

    def test_da1_single_leaf_default_post_path_also_works(self):
        ped = rows_to_ped()
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertEqual(
            set(ALL) - {LEAF},
            {int(a.originalID) for a in ped.pedigree})
        self.assertEqual(
            {k: v for k, v in PRE_GRAPH.items() if k != LEAF},
            stable_graph(ped))

    def test_da9_founder_deletion_is_refused_while_children_survive(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [FOUNDER], before)

    def test_da10_half_founder_deletion_is_refused_while_children_survive(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [HALF_SIRE], before)


# ===========================================================================
# Desired batch contract. Was xfail(strict=True) on f7963ff.
# ===========================================================================
class TestDesiredBatchContract(unittest.TestCase):
    def test_da2_two_nonadjacent_closed_set_removes_exactly_the_requested(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([ORDINARY2, LEAF])
        assert_deleted_exactly(self, ped, [ORDINARY2, LEAF], ok=ok)

    def test_da2_open_nonadjacent_parent_set_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [HALF_SIRE, CHILD], before)

    def test_da3_two_adjacent_closed_set_removes_exactly_the_requested(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([CHILD, LEAF])
        assert_deleted_exactly(self, ped, [CHILD, LEAF], ok=ok)

    def test_da3_open_adjacent_parent_set_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [ORDINARY, ORDINARY2], before)

    def test_da4_three_closed_deletions_remove_exactly_the_requested_set(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([ORDINARY, CHILD, LEAF])
        assert_deleted_exactly(self, ped, [ORDINARY, CHILD, LEAF], ok=ok)

    def test_da4_open_three_set_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(
            self, ped, [HALF_SIRE, ORDINARY, CHILD], before)

    def test_da5_request_order_does_not_change_the_deleted_set(self):
        first = isolate(rows_to_ped())
        second = isolate(rows_to_ped())
        ok_ab = first.delete_animals([ORDINARY2, LEAF])
        ok_ba = second.delete_animals([LEAF, ORDINARY2])
        self.assertTrue(ok_ab)
        self.assertTrue(ok_ba)
        self.assertEqual(set(originals(first)), set(originals(second)))
        self.assertEqual(set(ALL) - {ORDINARY2, LEAF}, set(originals(first)))
        self.assertEqual(stable_graph(first), stable_graph(second))

    def test_da6_unrelated_survivors_are_preserved(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.delete_animals([ORDINARY2, LEAF]))
        self.assertNotIn(ORDINARY2, originals(ped))
        self.assertNotIn(LEAF, originals(ped))
        self.assertIn(FOUNDER, originals(ped))
        self.assertIn(HALF_SIRE, originals(ped))
        self.assertIn(HALF_DAM, originals(ped))
        self.assertIn(ORDINARY, originals(ped))
        self.assertIn(CHILD, originals(ped))

    def test_da7_retained_known_parent_edges_are_preserved(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.delete_animals([ORDINARY2, LEAF]))
        graph = stable_graph(ped)
        self.assertEqual(("MISSING", "MISSING"), graph[FOUNDER])
        self.assertEqual(("MISSING", FOUNDER), graph[HALF_DAM])
        self.assertEqual((ORDINARY, HALF_DAM), graph[CHILD])

    def test_da13_no_duplicate_or_stale_id_map_entries_after_batch(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.delete_animals([ORDINARY2, LEAF]))
        remaining = set(originals(ped))
        self.assertEqual(remaining, {int(k) for k in ped.idmap})
        self.assertEqual(remaining, {int(v) for v in ped.backmap.values()})
        self.assertNotIn(ORDINARY2, ped.idmap)
        self.assertNotIn(LEAF, ped.idmap)
        self.assertEqual(len(ped.idmap), len(ped.pedigree))
        self.assertEqual(len(ped.backmap), len(ped.pedigree))
        self.assertEqual(len(set(ped.idmap.values())), len(ped.idmap))

    def test_two_trailing_animals_are_both_removed(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([CHILD, LEAF])
        assert_deleted_exactly(self, ped, [CHILD, LEAF], ok=ok)


class TestAlreadyCorrectAdjacentDescending(unittest.TestCase):
    def test_da3_descending_open_adjacent_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [ORDINARY2, ORDINARY], before)

    def test_da2_descending_open_nonadjacent_is_refused(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [CHILD, HALF_SIRE], before)

    def test_da2_descending_closed_nonadjacent_removes_the_requested_set(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF, ORDINARY2])
        assert_deleted_exactly(self, ped, [LEAF, ORDINARY2], ok=ok)


# ===========================================================================
# DA-8 --  supersedes dangling success.
# ===========================================================================
class TestOrphanContract(unittest.TestCase):
    def test_da8_referenced_parent_is_refused_not_left_dangling(self):
        """
        delete_animals still does not rewrite sireID/damID. A surviving
        child that names the requested parent now blocks the delete.
        """
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [ORDINARY], before)
        child = next(a for a in ped.pedigree if int(a.originalID) == CHILD)
        self.assertEqual(PRE_AID[ORDINARY], child.sireID)
        self.assertNotEqual(ped.kw["missing_parent"], child.sireID)

    def test_da8_descendants_are_not_cascade_deleted(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [ORDINARY], before)
        self.assertIn(CHILD, originals(ped))
        self.assertIn(LEAF, originals(ped))

    def test_parent_delete_with_default_renumber_refuses_before_mutation(self):
        """
        : referenced-parent deletion is refused before the
        loop.  never sees a dangling intermediate.
        """
        ped = rows_to_ped()
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [HALF_SIRE], before)
        self.assertIn(HALF_SIRE, originals(ped))
        self.assertEqual((HALF_SIRE, HALF_DAM), stable_graph(ped)[ORDINARY2])


# ===========================================================================
# DA-11 / DA-12 missing_parent. Isolated loop does not rewrite slots.
#  is the fast_reorder wrapper, tested separately.
# ===========================================================================
class TestMissingParentBoundary(unittest.TestCase):
    def test_da12_default_zero_is_the_shipped_sentinel(self):
        ped = rows_to_ped()
        self.assertEqual(0, ped.kw["missing_parent"])

    def test_da12_isolated_single_delete_does_not_hardcode_a_new_sentinel(self):
        ped = isolate(rows_to_ped())
        ped.delete_animals([LEAF])
        self.assertEqual(0, ped.kw["missing_parent"])
        self.assertTrue(all(
            a.sireID != CUSTOM and a.damID != CUSTOM for a in ped.pedigree))

    def test_da11_custom_sentinel_open_batch_is_refused(self):
        ped = isolate(rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM))
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        before = snapshot_ids(ped)
        assert_refused_unchanged(self, ped, [HALF_SIRE, CHILD], before)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        founder = next(a for a in ped.pedigree if int(a.originalID) == FOUNDER)
        self.assertEqual(CUSTOM, founder.sireID)
        self.assertEqual(CUSTOM, founder.damID)

    def test_da11_custom_sentinel_closed_batch_removes_the_requested_set(self):
        ped = isolate(rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM))
        ok = ped.delete_animals([ORDINARY2, LEAF])
        remaining, edges = expected_after({ORDINARY2, LEAF})
        self.assertTrue(ok)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        self.assertEqual(remaining, set(originals(ped)))
        self.assertEqual(edges, stable_graph(ped))

    def test_da11_custom_single_leaf_preserves_minus_999_slots(self):
        ped = isolate(rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM))
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        founder = next(a for a in ped.pedigree if int(a.originalID) == FOUNDER)
        self.assertEqual(CUSTOM, founder.sireID)
        half_sire = next(
            a for a in ped.pedigree if int(a.originalID) == HALF_SIRE)
        self.assertEqual(CUSTOM, half_sire.damID)

    def test_medium3_fast_wrapper_still_binds_the_configured_sentinel(self):
        """FR-5 control. One leaf, so  arithmetic is not in play."""
        from unittest import mock
        from PyPedal import pyp_utils

        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        ped.kw["reorder"] = 1
        ped.kw["renumber"] = 0
        ped.kw["slow_reorder"] = False
        with mock.patch(
            "PyPedal.pyp_newclasses.pyp_utils.fast_reorder",
            wraps=pyp_utils.fast_reorder,
        ) as spy:
            ped.delete_animals([LEAF])
        self.assertTrue(spy.called)
        args, kwargs = spy.call_args
        bound = inspect.signature(pyp_utils.fast_reorder).bind(*args, **kwargs)
        bound.apply_defaults()
        self.assertEqual(CUSTOM, bound.arguments["missingparent"])


# ===========================================================================
# DA-14 / DA-15 metadata and offspring -- only what the method promises.
# ===========================================================================
class TestMetadataAndOffspring(unittest.TestCase):
    def test_da14_successful_delete_refreshes_metadata_record_count(self):
        ped = isolate(rows_to_ped())
        self.assertEqual(7, ped.metadata.num_records)
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertEqual(6, ped.metadata.num_records)
        self.assertEqual(6, len(ped.pedigree))

    def test_da14_duplicate_request_refreshes_metadata_once(self):
        """``[A, A]`` is set semantics: delete A once and return True."""
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF, LEAF])
        self.assertTrue(ok)
        self.assertEqual(6, ped.metadata.num_records)
        self.assertEqual(6, len(ped.pedigree))

    def test_da15_isolated_path_rebuilds_offspring_lists(self):
        ped = isolate(rows_to_ped())
        before = next(a for a in ped.pedigree if int(a.originalID) == ORDINARY2)
        self.assertIn(PRE_AID[LEAF], before.unks)
        ped.delete_animals([LEAF])
        after = next(a for a in ped.pedigree if int(a.originalID) == ORDINARY2)
        self.assertNotIn(PRE_AID[LEAF], after.unks)

    def test_da15_successful_default_renumber_rebuilds_offspring(self):
        ped = rows_to_ped()
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        leftover = [
            key
            for animal in ped.pedigree
            for key in list(animal.sons) + list(animal.daus) + list(animal.unks)
        ]
        self.assertNotIn(PRE_AID[LEAF], leftover)
        self.assertNotIn(7, leftover)


# ===========================================================================
# DA-16 / DA-17 / empty / mixed -- pin current behaviour.
# ===========================================================================
class TestUnresolvedAndDuplicateRequests(unittest.TestCase):
    def test_da16_duplicate_request_deletes_once_and_returns_true(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF, LEAF])
        self.assertTrue(ok)
        self.assertNotIn(LEAF, originals(ped))
        self.assertEqual(6, len(ped.pedigree))
        self.assertEqual(6, ped.metadata.num_records)

    def test_da17_nonexistent_id_raises_and_mutates_nothing(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([777])
        self.assertEqual(before, snapshot_ids(ped))

    def test_empty_list_returns_true_and_keeps_every_animal(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([])
        self.assertTrue(ok)
        self.assertEqual(list(ALL), originals(ped))
        self.assertEqual(7, ped.metadata.num_records)

    def test_valid_then_invalid_refuses_without_mutation(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF, 777])
        self.assertEqual(before, snapshot_ids(ped))

    def test_invalid_then_valid_mutates_nothing(self):
        ped = isolate(rows_to_ped())
        before = snapshot_ids(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([777, LEAF])
        self.assertEqual(before, snapshot_ids(ped))


class TestDa18NoUnexpectedConfigurationMutation(unittest.TestCase):
    def test_isolated_success_does_not_rewrite_kw_missing_parent(self):
        ped = isolate(rows_to_ped())
        before = dict(ped.kw)
        ped.delete_animals([LEAF])
        self.assertEqual(before["missing_parent"], ped.kw["missing_parent"])
        self.assertEqual(before["renumber"], ped.kw["renumber"])
        self.assertEqual(before["reorder"], ped.kw["reorder"])
        self.assertEqual(before["pedformat"], ped.kw["pedformat"])


class TestUnrenumberedLoadIsCollateral(unittest.TestCase):
    """
    If the pedigree was never renumbered, ``renumberedID`` is still the
    initialisation sentinel and ``backmap[renumberedID]`` KeyError's
    after the animal has been located. ``delanimal()`` already warns
    that unrenumbered deletion is unsafe. Characterise only; this is
    not a  product claim.
    """

    def test_unrenumbered_single_leaf_raises_without_mutation(self):
        ped = rows_to_ped(
            renumber=False,
            pedigree_is_renumbered=False,
        )
        self.assertEqual(LEAF, ped.idmap[LEAF])
        self.assertEqual(-999, ped.pedigree[-1].renumberedID)
        before = snapshot_ids(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF])
        self.assertEqual(before, snapshot_ids(ped))
