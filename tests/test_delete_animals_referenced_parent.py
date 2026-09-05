"""
 -- 4.0 supported ``delete_animals`` parent-deletion boundary.

WHAT THIS FILE IS ABOUT
-----------------------
A successful deletion must not leave a dangling known-parent ID.
DUPLICATE_REDIRECT cannot be implemented after a Finding-37 load (no
live same-originalID pair; request identity is originalID). That future
policy is deferred.

The 4.0 supported contract is:

    A. non-parent deletion proceeds;
    B. a parent may be deleted only when every animal that references
       it is also explicitly in animal_list;
    C. a parent still referenced by a survivor is refused BEFORE any
       mutation, returning False;
    D. descendants are never implicitly added to the request.

Request identities are idmap keys (originalID after load+renumber).
Parent slots are current animalIDs. Tests use originals 100..700 so
those domains are not equal.

HOW THIS FILE WAS BUILT
-----------------------
Commit 1 landed the supported refusal/success contract with
``xfail(strict=True)``. Those markers are removed here after the
preflight. ``self.subTest`` stays out of the former-xfail tests.

WHAT IS *NOT* CLAIMED HERE
--------------------------
No orphaning to missing_parent. No cascade. No duplicate redirect.
 identity,  sentinel binding,  refusal of
genuinely dangling input, C1, and general C3 are not reopened.
"""
import copy
import os
import unittest

from PyPedal.pyp_errors import PyPedalPedigreeStructureError

from _pedhelpers import owned_temp_dir, load_corpus_from_path

BASELINE = "a8b62e4"

S, D, C, H, M, G, U = 100, 200, 300, 400, 500, 600, 700
ALL = (S, D, C, H, M, G, U)
CUSTOM = -999

ROWS = [
    "%s 0 0" % S,
    "%s 0 0" % D,
    "%s %s %s" % (C, S, D),
    "%s %s 0" % (H, S),
    "%s %s 0" % (M, S),
    "%s %s 0" % (G, C),
    "%s 0 0" % U,
]
CUSTOM_ROWS = [
    "%s %s %s" % (S, CUSTOM, CUSTOM),
    "%s %s %s" % (D, CUSTOM, CUSTOM),
    "%s %s %s" % (C, S, D),
    "%s %s %s" % (H, S, CUSTOM),
    "%s %s %s" % (M, S, CUSTOM),
    "%s %s %s" % (G, C, CUSTOM),
    "%s %s %s" % (U, CUSTOM, CUSTOM),
]


def rows_to_ped(rows=None, **overrides):
    rows = ROWS if rows is None else rows
    tmp = owned_temp_dir(prefix="f39_")
    path = os.path.join(tmp, "f39.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, "asd", **overrides)


def originals(ped):
    return [int(a.originalID) for a in ped.pedigree]


def by_original(ped, oid):
    return next(a for a in ped.pedigree if int(a.originalID) == oid)


def snapshot(ped):
    return {
        "originals": originals(ped),
        "animal_ids": [int(a.animalID) for a in ped.pedigree],
        "objects": [id(a) for a in ped.pedigree],
        "slots": {
            int(a.originalID): (a.sireID, a.damID, a.animalID, a.renumberedID)
            for a in ped.pedigree
        },
        "offspring": {
            int(a.originalID): (
                dict(a.sons), dict(a.daus), dict(a.unks),
            )
            for a in ped.pedigree
        },
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "namemap": dict(ped.namemap),
        "namebackmap": dict(ped.namebackmap),
        "n": len(ped.pedigree),
        "meta": ped.metadata.num_records,
        "kw": copy.deepcopy(dict(ped.kw)),
    }


def no_dangling(ped):
    missing = ped.kw["missing_parent"]
    live = {int(a.animalID) for a in ped.pedigree}
    bad = []
    for animal in ped.pedigree:
        for side, pid in (("sire", animal.sireID), ("dam", animal.damID)):
            if pid == missing or str(pid) == str(missing):
                continue
            if int(pid) not in live:
                bad.append((int(animal.originalID), side, pid))
    return bad


def assert_unchanged(test, before, ped):
    after = snapshot(ped)
    test.assertEqual(before, after)


def assert_parent_refused(test, ped, requested, before):
    with test.assertRaises(PyPedalPedigreeStructureError):
        ped.delete_animals(requested)
    assert_unchanged(test, before, ped)


# ===========================================================================
# Anti-vacuity. Green on the broken tree and the repaired tree.
# ===========================================================================
class TestFixtureIsNotVacuous(unittest.TestCase):
    def test_load_renumbers_so_original_and_animal_ids_differ(self):
        ped = rows_to_ped()
        self.assertEqual(sorted(ALL), sorted(originals(ped)))
        oids = [int(a.originalID) for a in ped.pedigree]
        aids = [int(a.animalID) for a in ped.pedigree]
        self.assertNotEqual(oids, aids)
        self.assertEqual(list(range(1, 8)), sorted(aids))
        sire = by_original(ped, S)
        child = by_original(ped, C)
        self.assertNotEqual(int(sire.originalID), int(sire.animalID))
        self.assertEqual(sire.animalID, child.sireID)
        self.assertNotEqual(sire.originalID, child.sireID)

    def test_s_is_referenced_by_c_h_m_and_not_by_g_or_u(self):
        ped = rows_to_ped()
        sire_aid = by_original(ped, S).animalID
        dam_aid = by_original(ped, D).animalID
        child_aid = by_original(ped, C).animalID
        self.assertEqual(sire_aid, by_original(ped, C).sireID)
        self.assertEqual(dam_aid, by_original(ped, C).damID)
        self.assertEqual(sire_aid, by_original(ped, H).sireID)
        self.assertEqual(sire_aid, by_original(ped, M).sireID)
        self.assertEqual(child_aid, by_original(ped, G).sireID)
        self.assertNotEqual(sire_aid, by_original(ped, G).sireID)
        self.assertNotEqual(sire_aid, by_original(ped, U).sireID)
        self.assertEqual(0, ped.kw["missing_parent"])


# ===========================================================================
# Desired 4.0 contract. Was xfail(strict=True) on a8b62e4.
# ===========================================================================
class TestDp42SireRefusal(unittest.TestCase):
    def test_dp42_delete_sire_with_surviving_child_refuses_before_mutation(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        child = by_original(ped, C)
        sire = by_original(ped, S)
        self.assertEqual(sire.animalID, child.sireID)
        self.assertNotIn(C, [S])
        assert_parent_refused(self, ped, [S], before)


class TestDp43DamRefusal(unittest.TestCase):
    def test_dp43_delete_dam_with_surviving_child_refuses_before_mutation(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        self.assertEqual(by_original(ped, D).animalID, by_original(ped, C).damID)
        assert_parent_refused(self, ped, [D], before)


class TestDp44MultipleOffspringRefusal(unittest.TestCase):
    def test_dp44_delete_sire_of_c_h_m_refuses_before_mutation(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        self.assertEqual(set(ALL), set(originals(ped)))


class TestDp45BothParentsRefusal(unittest.TestCase):
    def test_dp45_delete_sire_and_dam_while_c_survives_refuses(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S, D], before)


class TestDp46PartialDependentsRefusal(unittest.TestCase):
    def test_dp46_delete_s_and_c_while_h_survives_refuses(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        self.assertIn(H, originals(ped))
        self.assertEqual(by_original(ped, S).animalID, by_original(ped, H).sireID)
        assert_parent_refused(self, ped, [S, C], before)


class TestDp47ExplicitClosureAllowed(unittest.TestCase):
    """Already green: no survivor names a requested parent."""

    def test_dp47_delete_c_and_g_succeeds(self):
        ped = rows_to_ped()
        ok = ped.delete_animals([C, G])
        self.assertTrue(ok)
        self.assertEqual({S, D, H, M, U}, set(originals(ped)))
        self.assertEqual([], no_dangling(ped))

    def test_dp47_delete_s_and_everyone_who_would_dangle_succeeds(self):
        ped = rows_to_ped()
        ok = ped.delete_animals([S, C, H, M, G])
        self.assertTrue(ok)
        self.assertEqual({D, U}, set(originals(ped)))
        self.assertEqual([], no_dangling(ped))


class TestDp49RefusalSnapshot(unittest.TestCase):
    def test_dp49_referenced_parent_refusal_is_a_full_state_noop(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)


class TestDp410CustomSentinelRefusal(unittest.TestCase):
    def test_dp410_custom_missing_parent_sire_delete_refuses_unchanged(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        founder = by_original(ped, S)
        self.assertEqual(CUSTOM, founder.sireID)
        self.assertEqual(CUSTOM, founder.damID)


class TestDp411DefaultSentinelRefusal(unittest.TestCase):
    def test_dp411_default_zero_sire_delete_refuses_unchanged(self):
        ped = rows_to_ped()
        self.assertEqual(0, ped.kw["missing_parent"])
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)


class TestDp412IdDomain(unittest.TestCase):
    def test_dp412_blocked_parent_id_is_the_current_animal_id(self):
        ped = rows_to_ped()
        sire = by_original(ped, S)
        child = by_original(ped, C)
        self.assertNotEqual(int(sire.originalID), int(sire.animalID))
        self.assertEqual(sire.animalID, child.sireID)
        self.assertNotEqual(sire.originalID, child.sireID)
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        self.assertEqual(sire.animalID, by_original(ped, C).sireID)


class TestDp415NoDanglingOnSuccess(unittest.TestCase):
    def test_dp415_non_parent_and_closed_parent_sets_have_no_dangling(self):
        leaf = rows_to_ped()
        self.assertTrue(leaf.delete_animals([U]))
        self.assertEqual([], no_dangling(leaf))
        closed = rows_to_ped()
        self.assertTrue(closed.delete_animals([C, G]))
        self.assertEqual([], no_dangling(closed))


class TestDp416Finding37AcceptsSuccess(unittest.TestCase):
    def test_dp416_default_renumber_succeeds_after_supported_deletes(self):
        ped = rows_to_ped()
        self.assertTrue(ped.delete_animals([U]))
        self.assertEqual([], no_dangling(ped))
        ped2 = rows_to_ped()
        self.assertTrue(ped2.delete_animals([C, G]))
        self.assertEqual([], no_dangling(ped2))
        self.assertEqual({S, D, H, M, U}, set(originals(ped2)))


class TestDp48RefusalDoesNotCascade(unittest.TestCase):
    def test_dp48_dp19_refused_sire_delete_does_not_remove_c_h_m_g(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        self.assertEqual(set(ALL), set(originals(ped)))


class TestDp48ExplicitChildIsNotCascade(unittest.TestCase):
    def test_dp47_explicit_g_is_removed_only_because_requested(self):
        ped = rows_to_ped()
        self.assertTrue(ped.delete_animals([C, G]))
        self.assertNotIn(G, originals(ped))
        self.assertIn(H, originals(ped))
        self.assertIn(M, originals(ped))
        self.assertIn(U, originals(ped))


class TestDp417Dp418Edges(unittest.TestCase):
    def test_dp418_successful_u_delete_leaves_c_edges(self):
        ped = rows_to_ped()
        before_c = (by_original(ped, C).sireID, by_original(ped, C).damID)
        self.assertTrue(ped.delete_animals([U]))
        self.assertEqual(
            before_c,
            (by_original(ped, C).sireID, by_original(ped, C).damID),
        )


# ===========================================================================
# Already-correct neighbours. Not xfailed.
# ===========================================================================
class TestDp41NonParentAlreadyWorks(unittest.TestCase):
    def test_dp41_delete_unrelated_u_succeeds(self):
        ped = rows_to_ped()
        sire_aid = by_original(ped, S).animalID
        ok = ped.delete_animals([U])
        self.assertTrue(ok)
        self.assertEqual({S, D, C, H, M, G}, set(originals(ped)))
        self.assertEqual([], no_dangling(ped))
        self.assertEqual(sire_aid, by_original(ped, C).sireID)

    def test_dp41_default_path_also_succeeds_for_u(self):
        ped = rows_to_ped()
        self.assertTrue(ped.kw["renumber"])
        self.assertTrue(ped.delete_animals([U]))
        self.assertEqual([], no_dangling(ped))


class TestDp413Medium6NonParentIdentity(unittest.TestCase):
    def test_dp413_two_non_parent_leaves_on_the_m6_fixture(self):
        from test_delete_animals import (
            CHILD, LEAF, ALL as M6_ALL, isolate, rows_to_ped as m6_load,
            originals as m6_orig, assert_deleted_exactly,
        )
        ped = isolate(m6_load())
        ok = ped.delete_animals([CHILD, LEAF])
        assert_deleted_exactly(self, ped, [CHILD, LEAF], ok=ok)
        self.assertEqual(set(M6_ALL) - {CHILD, LEAF}, set(m6_orig(ped)))


class TestDp414ClosedParentSetHasNoDangling(unittest.TestCase):
    def test_dp414_parent_plus_all_dependents_leaves_no_dangling_survivor(self):
        ped = rows_to_ped()
        self.assertTrue(ped.delete_animals([S, C, H, M, G]))
        self.assertEqual({D, U}, set(originals(ped)))
        self.assertEqual([], no_dangling(ped))


class TestDp417NoSireDamSwap(unittest.TestCase):
    def test_dp417_closed_delete_does_not_swap_remaining_parent_slots(self):
        ped = rows_to_ped()
        before = {
            int(a.originalID): (a.sireID, a.damID)
            for a in ped.pedigree
        }
        self.assertTrue(ped.delete_animals([C, G]))
        for oid in (S, D, H, M, U):
            animal = by_original(ped, oid)
            self.assertEqual(before[oid], (animal.sireID, animal.damID))


class TestCurrentMasterCharacterisation(unittest.TestCase):
    """
    DP4-20. Isolated and default parent delete now refuse before
    mutation. Python 2's dangling / false-success path is not
    preserved. Deliberate 4.0 correctness divergence.
    """

    def test_isolated_sire_delete_refuses_without_mutation(self):
        ped = rows_to_ped()
        ped.kw["renumber"] = False
        ped.kw["reorder"] = 0
        ped.kw["set_offspring"] = False
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        self.assertIn(S, originals(ped))

    def test_default_sire_delete_refuses_without_mutation(self):
        ped = rows_to_ped()
        before = snapshot(ped)
        assert_parent_refused(self, ped, [S], before)
        self.assertIn(S, originals(ped))
