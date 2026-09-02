"""
Post-rc1 candidate B -- Lacy unexpected failure must not return -999.9.

WHAT THIS FILE IS ABOUT
-----------------------
The supported Lacy success path is already scientifically validated
(``f_e ≈ 2.909090909`` on ``new_lacy.ped``). Unexpected exceptions in
``a_effective_founders_lacy`` / ``effective_founders_lacy`` currently
return ``{"fa_effective_founders": -999.9}``, which is not a legal
``f_e``. ``effective_founders_lacy`` already re-raises
``PyPedalValidationError``; ``a_effective_founders_lacy`` swallows it.

Desired 4.0 contract: preserve existing typed ``PyPedalError`` subclasses;
wrap only an unexpected ``Exception`` as ``PyPedalError``. Success schema
and founder mathematics are unchanged.

HOW THIS FILE WAS BUILT
-----------------------
Commit 1 landed characterisation of the token plus the desired refusal
marked ``xfail(strict=True)``. Those markers are removed here after the
repair. Characterisation is inverted: unexpected failure now raises.
``self.subTest`` stays out of the former-xfail tests. Candidate C is
not tested here.

WHAT IS *NOT* CLAIMED HERE
--------------------------
The Lacy calculation is not changed. ``relationship()`` and
``a_coefficients`` are out of scope.
"""
import os
import tempfile
import unittest

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalError, PyPedalUsageError, PyPedalValidationError

from _pedhelpers import chdir_tmp, load_corpus, load_corpus_from_path

BASELINE = "d14dbe8"
LACY_PUBLISHED = 2.909090909090909


def snapshot(ped):
    return {
        "n": len(ped.pedigree),
        "objects": [id(animal) for animal in ped.pedigree],
        "animal_ids": [int(animal.animalID) for animal in ped.pedigree],
        "original_ids": [animal.originalID for animal in ped.pedigree],
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
    }


def founders_only():
    tmp = tempfile.mkdtemp(prefix="lacy_fail_")
    path = os.path.join(tmp, "founders.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("1 0 0\n2 0 0\n")
    return load_corpus_from_path(path, "asd")


# ===========================================================================
# Current-master characterisation (green on d14dbe8)
# ===========================================================================
class TestCurrentTokenNoLongerReturned(unittest.TestCase):
    """The d14dbe8 -999.9 token is gone. Characterisation inverted."""

    def test_a_effective_founders_lacy_raises_on_none_pedigree(self):
        ped = load_corpus("new_lacy.ped")
        ped.pedigree = None
        with chdir_tmp():
            with self.assertRaises(PyPedalError):
                pyp_metrics.a_effective_founders_lacy(ped)

    def test_effective_founders_lacy_raises_on_none_pedigree(self):
        ped = load_corpus("new_lacy.ped")
        ped.pedigree = None
        with chdir_tmp():
            with self.assertRaises(PyPedalError):
                pyp_metrics.effective_founders_lacy(ped)


# ===========================================================================
# Desired contract. LFS-1 / LFS-2 / LFS-4 already hold (not xfailed).
# ===========================================================================
class TestLFS1EstablishedSuccessUnchanged(unittest.TestCase):
    def test_a_effective_founders_lacy_matches_lacy_appendix_a(self):
        ped = load_corpus("new_lacy.ped")
        with chdir_tmp():
            got = pyp_metrics.a_effective_founders_lacy(ped)
        self.assertAlmostEqual(LACY_PUBLISHED, got["fa_effective_founders"], places=9)
        self.assertEqual(
            {"fa_animal_count", "fa_founder_count",
             "fa_descendant_count", "fa_effective_founders"},
            set(got),
        )

    def test_effective_founders_lacy_matches_lacy_appendix_a(self):
        ped = load_corpus("new_lacy.ped")
        with chdir_tmp():
            got = pyp_metrics.effective_founders_lacy(ped)
        self.assertAlmostEqual(LACY_PUBLISHED, got["fa_effective_founders"], places=9)


class TestLFS2KnownTypedErrorsPreserved(unittest.TestCase):
    def test_bad_mode_raises_usage_error_on_both_entry_points(self):
        ped = load_corpus("new_lacy.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.a_effective_founders_lacy(ped, mode="not-a-mode")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.effective_founders_lacy(ped, mode="not-a-mode")

    def test_effective_founders_lacy_founders_only_raises_validation_error(self):
        ped = founders_only()
        with chdir_tmp():
            with self.assertRaises(PyPedalValidationError):
                pyp_metrics.effective_founders_lacy(ped)


class TestLFS4CallerStateUnchangedOnKnownFailure(unittest.TestCase):
    def test_bad_mode_does_not_mutate(self):
        ped = load_corpus("new_lacy.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.a_effective_founders_lacy(ped, mode="not-a-mode")
        self.assertEqual(before, snapshot(ped))


class TestLFS3UnexpectedFailureIsTyped(unittest.TestCase):
    def test_none_pedigree_raises_pypedal_error_not_the_token(self):
        ped = load_corpus("new_lacy.ped")
        before = snapshot(ped)
        ped.pedigree = None
        with chdir_tmp():
            with self.assertRaises(PyPedalError):
                pyp_metrics.a_effective_founders_lacy(ped)
            with self.assertRaises(PyPedalError):
                pyp_metrics.effective_founders_lacy(ped)
        # The animal list is gone because the test replaced it; maps must
        # still be the pre-call maps (the routine must not rebuild them).
        self.assertEqual(before["idmap"], dict(ped.idmap))
        self.assertEqual(before["backmap"], dict(ped.backmap))


class TestLFS3FoundersOnlyAEffectiveRaises(unittest.TestCase):
    def test_a_effective_founders_lacy_founders_only_raises(self):
        ped = founders_only()
        before = snapshot(ped)
        with chdir_tmp():
            with self.assertRaises(PyPedalError):
                pyp_metrics.a_effective_founders_lacy(ped)
        self.assertEqual(before, snapshot(ped))

