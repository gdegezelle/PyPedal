"""Opt-in override of the Boichard reference-population antichain refusal."""
import inspect
import os
import unittest

from PyPedal import pyp_errors, pyp_metrics

from _pedhelpers import load_corpus, load_corpus_from_path, owned_temp_dir

NON_ANTICHAIN = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 5 4 2\n"
ANTICHAIN = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 3 4 2\n"

ROUTINES = (
    pyp_metrics.a_effective_founders_boichard,
    pyp_metrics.a_effective_ancestors_definite,
    pyp_metrics.a_effective_ancestors_indefinite,
)


def _write(text):
    tmp = owned_temp_dir(prefix="pypedal_related_ref_")
    path = os.path.join(tmp, "fixture.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return load_corpus_from_path(path, "asdg")


class TestAllowRelatedReferenceSignature(unittest.TestCase):
    def test_default_is_false_and_keyword_only(self):
        for routine in ROUTINES:
            with self.subTest(routine=routine.__name__):
                parameters = inspect.signature(routine).parameters
                self.assertIn("allow_related_reference", parameters)
                param = parameters["allow_related_reference"]
                self.assertEqual(inspect.Parameter.KEYWORD_ONLY, param.kind)
                self.assertIs(False, param.default)


class TestDefaultStillRefuses(unittest.TestCase):
    def test_definite_refuses_a_related_reference(self):
        ped = _write(NON_ANTICHAIN)
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertIn("antichain", str(caught.exception))
        self.assertIn("allow_related_reference=True", str(caught.exception))

    def test_indefinite_refuses_a_related_reference(self):
        ped = _write(NON_ANTICHAIN)
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            pyp_metrics.a_effective_ancestors_indefinite(ped, n=1)
        self.assertIn("antichain", str(caught.exception))
        self.assertIn("allow_related_reference=True", str(caught.exception))


class TestOverrideComputes(unittest.TestCase):
    def test_definite_override_returns_finite_f_a(self):
        ped = _write(NON_ANTICHAIN)
        got = pyp_metrics.a_effective_ancestors_definite(
            ped, output=False, allow_related_reference=True)
        self.assertGreaterEqual(got, 1.0)
        self.assertTrue(got < float("inf"))

    def test_indefinite_override_returns_finite_bounds(self):
        ped = _write(NON_ANTICHAIN)
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            ped, n=1, output=False, allow_related_reference=True)
        self.assertGreaterEqual(f_l, 1.0)
        self.assertGreaterEqual(f_u, f_l)


class TestOverrideDoesNotChangeMathematics(unittest.TestCase):
    def test_valid_antichain_is_identical_with_either_flag(self):
        ped = _write(ANTICHAIN)
        off = pyp_metrics.a_effective_ancestors_definite(
            ped, output=False, allow_related_reference=False)
        on = pyp_metrics.a_effective_ancestors_definite(
            ped, output=False, allow_related_reference=True)
        self.assertEqual(off, on)
        off_b = pyp_metrics.a_effective_ancestors_indefinite(
            ped, n=1, output=False, allow_related_reference=False)
        on_b = pyp_metrics.a_effective_ancestors_indefinite(
            ped, n=1, output=False, allow_related_reference=True)
        self.assertEqual(off_b, on_b)

    def test_founders_flag_does_not_change_a_related_reference(self):
        ped = _write(NON_ANTICHAIN)
        off = pyp_metrics.a_effective_founders_boichard(
            ped, output=False, allow_related_reference=False)
        on = pyp_metrics.a_effective_founders_boichard(
            ped, output=False, allow_related_reference=True)
        self.assertEqual(off, on)
        self.assertGreaterEqual(off, 1.0)

    def test_corpus_antichain_is_identical_with_either_flag(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        off = pyp_metrics.a_effective_ancestors_definite(
            ped, output=False, allow_related_reference=False)
        on = pyp_metrics.a_effective_ancestors_definite(
            ped, output=False, allow_related_reference=True)
        self.assertEqual(off, on)


if __name__ == "__main__":
    unittest.main()
