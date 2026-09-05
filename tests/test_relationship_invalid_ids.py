"""
Post-rc1 candidate A -- ``relationship()`` unresolved IDs must raise.

WHAT THIS FILE IS ABOUT
-----------------------
``pyp_metrics.relationship()`` documents ``anim_a`` / ``anim_b`` as current
renumbered animal IDs. Unresolved IDs currently return ``0.0`` after a
logged warning. ``0.0`` is also the coefficient for genuinely unrelated
animals, so the failure is ambiguous.

Desired 4.0 contract: IDs outside the current pedigree raise
``PyPedalUsageError``. Valid related and genuinely-zero pairs are
unchanged. Original IDs are not translated through ``idmap``.

HOW THIS FILE WAS BUILT
-----------------------
Commit 1 landed characterisation of the ``d14dbe8`` ``0.0`` failure plus
the desired refusal marked ``xfail(strict=True)``. Those markers are
removed here after the repair. Characterisation is inverted: unresolved
IDs now raise. ``self.subTest`` stays out of the former-xfail tests.

WHAT IS *NOT* CLAIMED HERE
--------------------------
The coefficient formula is not changed. Lacy, ``a_coefficients``,
``delete_animals``, and mating COI are out of scope.
"""
import os
import unittest

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

from _pedhelpers import owned_temp_dir, chdir_tmp, load_corpus, load_corpus_from_path

BASELINE = "d14dbe8"


STUD_ROWS = [
    "100 0 0",
    "200 0 0",
    "300 100 200",
    "400 0 0",
]


def studbook(**overrides):
    tmp = owned_temp_dir(prefix="rel_api_")
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
    }


# ===========================================================================
# Anti-vacuity and current-master characterisation (green on d14dbe8)
# ===========================================================================
class TestFixtureIsNotVacuous(unittest.TestCase):
    def test_studbook_original_and_animal_ids_differ(self):
        ped = studbook()
        oids = [int(a.originalID) for a in ped.pedigree]
        aids = [int(a.animalID) for a in ped.pedigree]
        self.assertEqual([100, 200, 400, 300], oids)
        self.assertEqual([1, 2, 3, 4], aids)
        self.assertNotEqual(oids, aids)
        self.assertEqual({100: 1, 200: 2, 400: 3, 300: 4}, ped.idmap)

    def test_mrode_pins_the_established_related_pair(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertAlmostEqual(pyp_metrics.relationship(4, 3, ped), 0.25, places=3)


class TestCurrentUnresolvedNoLongerReturnsZero(unittest.TestCase):
    """The d14dbe8 0.0 failure is gone. Characterisation inverted."""

    def test_both_unresolved_now_raise(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(99999, 88888, ped)

    def test_original_ids_on_a_studbook_now_raise(self):
        ped = studbook()
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(300, 100, ped)


# ===========================================================================
# Desired contract. REL-1 / REL-2 / REL-6 are already true (not xfailed).
# ===========================================================================
class TestREL1ValidRelatedPairUnchanged(unittest.TestCase):
    def test_mrode_4_vs_3_is_one_quarter(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.relationship(4, 3, ped)
        self.assertAlmostEqual(0.25, got, places=12)
        self.assertEqual(before, snapshot(ped))


class TestREL2ValidZeroRelatedPairRemainsZero(unittest.TestCase):
    def test_mrode_founders_are_unrelated(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertEqual(0.0, pyp_metrics.relationship(1, 2, ped))

    def test_studbook_unrelated_founder_is_zero_on_current_ids(self):
        ped = studbook()
        with chdir_tmp():
            self.assertEqual(0.0, pyp_metrics.relationship(3, 1, ped))


class TestREL3FirstIdUnresolvedRaises(unittest.TestCase):
    def test_first_id_unresolved_raises_usage_error(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(99999, 3, ped)
        self.assertEqual(before, snapshot(ped))


class TestREL4SecondIdUnresolvedRaises(unittest.TestCase):
    def test_second_id_unresolved_raises_usage_error(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(3, 88888, ped)
        self.assertEqual(before, snapshot(ped))


class TestREL5BothUnresolvedRaise(unittest.TestCase):
    def test_both_unresolved_raise_usage_error(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(99999, 88888, ped)
        self.assertEqual(before, snapshot(ped))


class TestREL6CallerPedigreeUnchangedOnSuccess(unittest.TestCase):
    def test_self_relationship_does_not_mutate(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with chdir_tmp():
            self.assertEqual(1.0, pyp_metrics.relationship(5, 5, ped))
        self.assertEqual(before, snapshot(ped))


class TestREL7OriginalIdsAreNotTheDomain(unittest.TestCase):
    def test_studbook_original_ids_raise_usage_error(self):
        ped = studbook()
        before = snapshot(ped)
        self.assertNotEqual(300, ped.idmap[300])
        with chdir_tmp():
            with self.assertRaises(PyPedalUsageError):
                pyp_metrics.relationship(300, 100, ped)
        self.assertEqual(before, snapshot(ped))


class TestREL7CurrentIdsOnTheStudbookStillCompute(unittest.TestCase):
    def test_current_ids_on_the_same_studbook_still_compute(self):
        ped = studbook()
        with chdir_tmp():
            self.assertAlmostEqual(0.5, pyp_metrics.relationship(4, 1, ped), places=12)
