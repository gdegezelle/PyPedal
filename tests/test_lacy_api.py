"""
Public Lacy effective-founder domain: phantom only.

Lacy (1989) p.113 specifies phantom founders. PyPedal 4 supports that
treatment on both entry points and refuses historical strict/absorb and
half= variants rather than returning a non-probability f_e.
"""
import unittest
import warnings

from _pedhelpers import chdir_tmp, corpus, load_corpus
from oracles import lacy_f_e

from PyPedal import pyp_errors, pyp_metrics

ENTRY_POINTS = (pyp_metrics.a_effective_founders_lacy,
                pyp_metrics.effective_founders_lacy)

CORPUS = (("new_lacy.ped", "asd"), ("generations.ped", "asdbx"),
          ("mrode.ped", "asd"), ("hartlandclark.ped", "asdb"))

HALF_FOUNDER_PEDIGREES = ("mrode.ped", "hartlandclark.ped")


def _f_e(fn, name, pedformat, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        with chdir_tmp():
            return fn(load_corpus(name, pedformat), **kw)["fa_effective_founders"]


class TestTheDefaultIsThePapersRule(unittest.TestCase):

    def test_default_mode_is_phantom(self):
        self.assertEqual("phantom", pyp_metrics.LACY_DEFAULT_MODE)
        self.assertEqual(("phantom",), pyp_metrics.LACY_MODES)

    def test_both_entry_points_match_the_independent_oracle(self):
        for name, pedformat in CORPUS:
            want, gate = lacy_f_e(corpus(name))
            self.assertTrue(gate, f"{name}: oracle q must sum to one unforced")
            for fn in ENTRY_POINTS:
                with self.subTest(pedigree=name, routine=fn.__name__):
                    self.assertAlmostEqual(want, _f_e(fn, name, pedformat),
                                           places=9)

    def test_the_two_entry_points_agree(self):
        for name, pedformat in CORPUS:
            with self.subTest(pedigree=name):
                self.assertAlmostEqual(
                    _f_e(pyp_metrics.a_effective_founders_lacy, name, pedformat),
                    _f_e(pyp_metrics.effective_founders_lacy, name, pedformat),
                    places=9)

    def test_half_founder_free_pedigrees_are_unchanged(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                self.assertAlmostEqual(
                    2.909090909090909, _f_e(fn, "new_lacy.ped", "asd"), places=9)
                self.assertAlmostEqual(
                    4.612612613, _f_e(fn, "generations.ped", "asdbx"), places=6)


class TestPhantomValues(unittest.TestCase):

    PHANTOM = {
        "mrode.ped": 2.797814208,
        "hartlandclark.ped": 5.831988609,
    }

    def test_new_values_are_the_papers(self):
        for name, want in self.PHANTOM.items():
            pedformat = dict(CORPUS)[name]
            for fn in ENTRY_POINTS:
                with self.subTest(pedigree=name, routine=fn.__name__):
                    self.assertAlmostEqual(want, _f_e(fn, name, pedformat),
                                           places=6)

    def test_f_e_never_exceeds_the_number_of_founder_sources(self):
        for name, pedformat in CORPUS:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                sources = (sum(1 for a in ped.pedigree if a.founder == "y")
                           + len(pyp_metrics.lacy_phantom_slots(ped)))
                self.assertLessEqual(_f_e(pyp_metrics.a_effective_founders_lacy,
                                          name, pedformat), sources + 1e-9)


class TestPublicLacyDomain(unittest.TestCase):
    """Both entry points accept phantom and refuse historical variants."""

    def _ped(self, name="mrode.ped"):
        return load_corpus(name, dict(CORPUS)[name])

    def test_default_equals_explicit_phantom_on_both_entry_points(self):
        for fn in ENTRY_POINTS:
            for name, pedformat in CORPUS:
                with self.subTest(routine=fn.__name__, pedigree=name):
                    default = _f_e(fn, name, pedformat)
                    explicit = _f_e(fn, name, pedformat, mode="phantom")
                    self.assertAlmostEqual(default, explicit, places=9)

    def test_mode_phantom_succeeds(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                got = fn(self._ped(), mode="phantom")
                self.assertAlmostEqual(
                    2.797814208, got["fa_effective_founders"], places=6)

    def test_mode_strict_raises_usage_error(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError) as raised:
                    fn(self._ped(), mode="strict")
                self.assertIn("phantom", str(raised.exception))

    def test_mode_absorb_raises_usage_error(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError) as raised:
                    fn(self._ped(), mode="absorb")
                self.assertIn("phantom", str(raised.exception))

    def test_half_true_raises_usage_error(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    fn(self._ped(), half=True)

    def test_half_false_raises_usage_error(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    fn(self._ped(), half=False)

    def test_invalid_mode_raises_usage_error(self):
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    fn(self._ped(), mode="not-a-mode")

    def test_conflicting_mode_and_half_raise_usage_error(self):
        for fn in ENTRY_POINTS:
            for mode in ("phantom", "strict", "absorb"):
                for half in (True, False):
                    with self.subTest(routine=fn.__name__, mode=mode, half=half):
                        with self.assertRaises(pyp_errors.PyPedalUsageError):
                            fn(self._ped(), mode=mode, half=half)

    def test_the_usage_error_is_also_a_value_error(self):
        with self.assertRaises(ValueError):
            pyp_metrics.a_effective_founders_lacy(self._ped(), mode="nonsense")

    def test_default_call_does_not_warn(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pyp_metrics.a_effective_founders_lacy(self._ped())
        self.assertEqual([], [w for w in caught
                              if issubclass(w.category, DeprecationWarning)])

    def test_half_founder_pedigree_is_refused_the_same_way(self):
        ped = self._ped("hartlandclark.ped")
        for fn in ENTRY_POINTS:
            with self.subTest(routine=fn.__name__):
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    fn(ped, mode="absorb")
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    fn(ped, half=False)


class TestHistoricalPartitionsAreNotProbabilityVectors(unittest.TestCase):
    """Independent oracle evidence; not a supported public calculation."""

    def test_oracle_strict_and_absorb_fail_the_sum_to_one_gate(self):
        for name in HALF_FOUNDER_PEDIGREES:
            path = corpus(name)
            with self.subTest(pedigree=name):
                _fe, phantom_ok = lacy_f_e(path, mode="phantom")
                self.assertTrue(phantom_ok)
                _fe, strict_ok = lacy_f_e(path, mode="lacy")
                self.assertFalse(strict_ok)
                _fe, absorb_ok = lacy_f_e(path, mode="half")
                self.assertFalse(absorb_ok)


if __name__ == "__main__":
    unittest.main()
