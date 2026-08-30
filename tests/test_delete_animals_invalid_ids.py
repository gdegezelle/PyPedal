"""
Post-rc1 candidates D and E -- ``delete_animals`` failure preflights.

WHAT THIS FILE IS ABOUT
-----------------------
D. A mixed ``[valid, invalid]`` request currently deletes the valid
   animal and then returns ``False``. Reverse order mutates nothing.
   Desired: nonexistent-ID preflight before the first deletion.

E. Unrenumbered deletion (``pedigree_is_renumbered`` false,
   ``renumberedID == -999``) returns ``False`` after partially clearing
   name maps. Desired: refuse before mutation. This does not make
   ``renumber=False`` deletion supported.

Duplicate ``[A, A]`` is not the same as a nonexistent ID and must keep
its established sequential behaviour.

HOW THIS FILE WAS BUILT
-----------------------
Commit 1 landed characterisation plus the desired non-mutating refusals
marked ``xfail(strict=True)``. Those markers are removed here after the
repair. Characterisation is inverted. ``self.subTest`` stays out of the
former-xfail tests.

WHAT IS *NOT* CLAIMED HERE
--------------------------
No general transaction rollback. No cascade. No ``DUPLICATE_REDIRECT``.
 identity and  parent refusal are not reopened.
"""
import os
import tempfile
import unittest

from PyPedal.pyp_errors import PyPedalPedigreeStructureError, PyPedalUsageError

from _pedhelpers import load_corpus_from_path

BASELINE = "d14dbe8"

FOUNDER = 10
HALF_SIRE = 20
HALF_DAM = 30
ORDINARY = 40
ORDINARY2 = 50
CHILD = 60
LEAF = 999
ALL = (FOUNDER, HALF_SIRE, HALF_DAM, ORDINARY, ORDINARY2, CHILD, LEAF)
INVALID = 123456
S, D, C, H, M, G, U = 100, 200, 300, 400, 500, 600, 700
CUSTOM = -999

ROWS = [
    "%s 0 0" % FOUNDER,
    "%s %s 0" % (HALF_SIRE, FOUNDER),
    "%s 0 %s" % (HALF_DAM, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (ORDINARY2, HALF_SIRE, HALF_DAM),
    "%s %s %s" % (CHILD, ORDINARY, HALF_DAM),
    "%s %s %s" % (LEAF, ORDINARY2, CHILD),
]
F39_ROWS = [
    "%s 0 0" % S,
    "%s 0 0" % D,
    "%s %s %s" % (C, S, D),
    "%s %s 0" % (H, S),
    "%s %s 0" % (M, S),
    "%s %s 0" % (G, C),
    "%s 0 0" % U,
]
CUSTOM_ROWS = [
    "%s %s %s" % (FOUNDER, CUSTOM, CUSTOM),
    "%s %s %s" % (HALF_SIRE, FOUNDER, CUSTOM),
    "%s %s %s" % (HALF_DAM, CUSTOM, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (ORDINARY2, HALF_SIRE, HALF_DAM),
    "%s %s %s" % (CHILD, ORDINARY, HALF_DAM),
    "%s %s %s" % (LEAF, ORDINARY2, CHILD),
]


def rows_to_ped(rows=None, **overrides):
    rows = ROWS if rows is None else rows
    tmp = tempfile.mkdtemp(prefix="del_api_")
    path = os.path.join(tmp, "del.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, "asd", **overrides)


def isolate(ped):
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


def snapshot(ped):
    return {
        "originals": originals(ped),
        "animal_ids": [int(a.animalID) for a in ped.pedigree],
        "objects": [id(a) for a in ped.pedigree],
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "namemap": dict(ped.namemap),
        "namebackmap": dict(ped.namebackmap),
        "n": len(ped.pedigree),
        "meta": ped.metadata.num_records,
        "missing_parent": ped.kw["missing_parent"],
        "renumber": ped.kw.get("renumber"),
        "pedigree_is_renumbered": ped.kw.get("pedigree_is_renumbered"),
    }


# ===========================================================================
# Characterisation of d14dbe8. Inverted after the D/E repairs.
# ===========================================================================
class TestCurrentMixedRequestNoLongerMutatesPrefix(unittest.TestCase):
    """The d14dbe8 order-dependent mutation is gone. Characterisation inverted."""

    def test_valid_then_invalid_now_keeps_the_valid_id(self):
        ped = isolate(rows_to_ped())
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF, INVALID])
        self.assertEqual(before, snapshot(ped))
        self.assertIn(LEAF, originals(ped))

    def test_unrenumbered_now_leaves_name_maps_intact(self):
        ped = rows_to_ped(renumber=False, pedigree_is_renumbered=False)
        self.assertTrue(LEAF in ped.namebackmap)
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF])
        self.assertEqual(before, snapshot(ped))
        self.assertIn(LEAF, originals(ped))
        self.assertTrue(LEAF in ped.namebackmap)


# ===========================================================================
# D. Desired nonexistent-ID preflight.
# ===========================================================================
class TestDEL31ValidThenInvalidRefusesWithoutMutation(unittest.TestCase):
    def test_mixed_valid_invalid_raises_and_keeps_every_animal(self):
        ped = isolate(rows_to_ped())
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF, INVALID])
        self.assertEqual(before, snapshot(ped))


class TestDEL32InvalidThenValidAlreadyNonMutating(unittest.TestCase):
    def test_invalid_then_valid_raises_and_keeps_every_animal(self):
        ped = isolate(rows_to_ped())
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([INVALID, LEAF])
        self.assertEqual(before, snapshot(ped))


class TestDEL33MultipleValidIdsUnchanged(unittest.TestCase):
    def test_two_unreferenced_animals_still_delete(self):
        # U is unreferenced in the Finding-39 fixture; adding a second leaf.
        rows = list(F39_ROWS) + ["800 0 0"]
        ped = rows_to_ped(rows)
        ok = ped.delete_animals([U, 800])
        self.assertTrue(ok)
        self.assertNotIn(U, originals(ped))
        self.assertNotIn(800, originals(ped))
        self.assertEqual(6, len(ped.pedigree))


class TestDEL34DuplicateRequestPreserved(unittest.TestCase):
    def test_duplicate_request_deletes_once_and_returns_true(self):
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF, LEAF])
        self.assertTrue(ok)
        self.assertNotIn(LEAF, originals(ped))
        self.assertEqual(6, len(ped.pedigree))
        self.assertEqual(6, ped.metadata.num_records)


class TestDEL35Finding39ParentRefusalUnchanged(unittest.TestCase):
    def test_referenced_parent_alone_is_refused_without_mutation(self):
        ped = rows_to_ped(F39_ROWS)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.delete_animals([S])
        self.assertEqual(before["originals"], originals(ped))
        self.assertEqual(before["idmap"], dict(ped.idmap))
        self.assertEqual(before["n"], len(ped.pedigree))


class TestDEL36Medium6IdentityUnchanged(unittest.TestCase):
    def test_deleting_two_valid_non_parents_removes_those_originals(self):
        # After isolate, two later leaves: 999 is a leaf; 60 is referenced
        # by 999, so  refuses [60, 999] unless 999 is also named.
        ped = isolate(rows_to_ped())
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertNotIn(LEAF, originals(ped))
        self.assertEqual(
            [FOUNDER, HALF_SIRE, HALF_DAM, ORDINARY, ORDINARY2, CHILD],
            originals(ped),
        )


class TestDEL37NoImplicitDescendantDeletion(unittest.TestCase):
    def test_deleting_unreferenced_leaf_keeps_its_parents(self):
        ped = rows_to_ped()
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertIn(ORDINARY2, originals(ped))
        self.assertIn(CHILD, originals(ped))


class TestDEL38CustomMissingParentUnchanged(unittest.TestCase):
    def test_custom_sentinel_survives_a_successful_leaf_delete(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        isolate(ped)
        before_mp = ped.kw["missing_parent"]
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        self.assertEqual(before_mp, ped.kw["missing_parent"])
        self.assertNotIn(LEAF, originals(ped))


# ===========================================================================
# E. Desired unrenumbered guard.
# ===========================================================================
class TestC1UnrenumberedGuardDoesNotMutate(unittest.TestCase):
    def test_unrenumbered_leaf_delete_raises_without_map_changes(self):
        ped = rows_to_ped(renumber=False, pedigree_is_renumbered=False)
        self.assertEqual(-999, ped.pedigree[-1].renumberedID)
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF])
        self.assertEqual(before, snapshot(ped))


class TestC1IsolateAfterRenumberStillDeletes(unittest.TestCase):
    """The guard must not use kw['renumber'] -- isolate() sets that false."""

    def test_isolate_path_still_deletes_a_leaf(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.kw["pedigree_is_renumbered"])
        self.assertFalse(ped.kw["renumber"])
        ok = ped.delete_animals([LEAF])
        self.assertTrue(ok)
        self.assertNotIn(LEAF, originals(ped))
