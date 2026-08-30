"""Known numerical differences from PyPedal 2.0.4, pinned on the Python 3 side.

PyPedal 2.0.4 return values for these routines are recorded here so the
suite does not need the historical interpreter. A value drifting away from
its registered difference fails.

The registry's rule is that anything in it is intentional and cited, and
anything not in it that differs is a regression. This file is the half of that
rule a test can enforce.
"""
import unittest
import warnings

from _pedhelpers import load_corpus

from PyPedal import pyp_errors, pyp_metrics

# pedigree -> pedformat
CORPUS = {"mrode.ped": "asd", "new_lacy.ped": "asd",
          "generations.ped": "asdbx", "hartlandclark.ped": "asdb",
          "boichard2a.ped": "asdg", "userfield.ped": "asdu"}

# Measured under pypedal-py27:audit. Recorded, not recomputed.
PYPEDAL_2_0_4 = {
    ("mrode.ped", "a_effective_founders_lacy"): 3.230283911671924,
    ("mrode.ped", "effective_founders_lacy"): 2.0,
    ("new_lacy.ped", "a_effective_founders_lacy"): 2.909090909090909,
    ("generations.ped", "a_effective_founders_lacy"): 4.612612612612613,
    ("hartlandclark.ped", "a_effective_founders_lacy"): 7.204573214944647,
    ("hartlandclark.ped", "effective_founders_lacy"): 3.0,
    ("boichard2a.ped", "a_effective_founders_boichard"): 4.0,
    ("boichard2a.ped", "a_effective_ancestors_definite"): 2.0,
    ("boichard2a.ped", "a_effective_founders_lacy"): 4.0,
    ("userfield.ped", "a_effective_founders_lacy"): 2.0,
}


def _lacy(routine, name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        return getattr(pyp_metrics, routine)(
            load_corpus(name, CORPUS[name]))["fa_effective_founders"]


class TestClassABoichardRefusesWhereTwoPointZeroFourReturnedInfinity(unittest.TestCase):
    """
    Class (A). PyPedal 2.0.4 returns **inf** from
    ``a_effective_ancestors_definite`` on five of six corpus pedigrees, because
    it has no generations precondition and computes from missing-value
    sentinels. Refusing is the repair, and it predates this phase.
    """

    WITHOUT_GENERATIONS = ("mrode.ped", "new_lacy.ped", "generations.ped",
                           "hartlandclark.ped", "userfield.ped")

    def test_the_boichard_routines_refuse_without_generations(self):
        for name in self.WITHOUT_GENERATIONS:
            for routine in ("a_effective_founders_boichard",
                            "a_effective_ancestors_definite",
                            "a_effective_ancestors_indefinite"):
                with self.subTest(pedigree=name, routine=routine):
                    with self.assertRaises(pyp_errors.PyPedalError):
                        getattr(pyp_metrics, routine)(
                            load_corpus(name, CORPUS[name]))

    def test_they_do_not_refuse_where_generations_exist(self):
        """Guard on the guard: the precondition must not refuse everything."""
        ped = load_corpus("boichard2a.ped", "asdg")
        self.assertEqual(4.0, pyp_metrics.a_effective_founders_boichard(ped))


class TestClassBEffectiveFoundersLacyDivergenceIsPreExisting(unittest.TestCase):
    """
    Class (B). ``effective_founders_lacy`` already differed from 2.0.4 on
    pedigrees with NO half-founder, and this phase did not change those values.

    2.0.4's figures there look like founder counts rather than effective
    numbers (3.0 on a three-founder pedigree). That is recorded as an open
    question, not repaired -- and the current value on ``new_lacy.ped`` is
    confirmed against Lacy's own published 2.91, so the current tree is right
    whatever 2.0.4 is doing.
    """

    def test_new_lacy_is_the_published_value_not_two_point_zero_fours(self):
        self.assertAlmostEqual(2.909090909090909,
                               _lacy("effective_founders_lacy", "new_lacy.ped"),
                               places=9)
        self.assertNotAlmostEqual(3.0,
                                  _lacy("effective_founders_lacy", "new_lacy.ped"),
                                  places=6)

    def test_generations_is_unchanged_by_this_phase(self):
        self.assertAlmostEqual(4.612612612612613,
                               _lacy("effective_founders_lacy", "generations.ped"),
                               places=9)


class TestClassCPhantomFounders(unittest.TestCase):
    """
    Class (C), INTRODUCED. Lacy (1989) p.113. Confined to pedigrees with a
    half-founder.
    """

    CHANGED = {"mrode.ped": 2.797814207650273,
               "hartlandclark.ped": 5.831988609}

    UNCHANGED = {"new_lacy.ped": 2.909090909090909,
                 "generations.ped": 4.612612612612613,
                 "boichard2a.ped": 4.0,
                 "userfield.ped": 2.0}

    def test_the_changed_values(self):
        for name, want in self.CHANGED.items():
            for routine in ("a_effective_founders_lacy", "effective_founders_lacy"):
                with self.subTest(pedigree=name, routine=routine):
                    self.assertAlmostEqual(want, _lacy(routine, name), places=6)

    def test_the_changed_values_are_no_longer_two_point_zero_fours(self):
        for (name, routine), old in PYPEDAL_2_0_4.items():
            if name not in self.CHANGED or "lacy" not in routine:
                continue
            with self.subTest(pedigree=name, routine=routine):
                self.assertNotAlmostEqual(old, _lacy(routine, name), places=6)

    def test_half_founder_free_pedigrees_stay_bit_exact_with_2_0_4(self):
        """
        The evidence that the change did not disturb the ordinary case. These
        are exact equality, not approximate.
        """
        for name, want in self.UNCHANGED.items():
            with self.subTest(pedigree=name):
                self.assertEqual(
                    want, _lacy("a_effective_founders_lacy", name))
                self.assertEqual(
                    PYPEDAL_2_0_4[(name, "a_effective_founders_lacy")],
                    _lacy("a_effective_founders_lacy", name))


class TestClassDTheBounds(unittest.TestCase):
    """
    Class (D), INTRODUCED. Boichard pp.9-10. 2.0.4 returns (NaN, 0.6667) on
    boichard2a.ped; both of those are wrong, the upper bound impossibly so.
    """

    def test_the_bounds_are_real_numbers_at_or_above_one(self):
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            load_corpus("boichard2a.ped", "asdg"))
        self.assertEqual((2.0, 2.0), (f_l, f_u))
        self.assertGreaterEqual(f_l, 1.0)


class TestWhatDidNotDiverge(unittest.TestCase):
    """
    The  and  repairs moved Python 3 INTO agreement with 2.0.4, against
    this phase's own prediction. Pinned as exact equality, because "close" would
    not be the claim being made.
    """

    def test_boichard_routines_are_bit_exact_with_2_0_4(self):
        ped_name = "boichard2a.ped"
        self.assertEqual(
            PYPEDAL_2_0_4[(ped_name, "a_effective_founders_boichard")],
            pyp_metrics.a_effective_founders_boichard(
                load_corpus(ped_name, "asdg")))
        self.assertEqual(
            PYPEDAL_2_0_4[(ped_name, "a_effective_ancestors_definite")],
            pyp_metrics.a_effective_ancestors_definite(
                load_corpus(ped_name, "asdg")))


if __name__ == "__main__":
    unittest.main()
