"""
``a_effective_ancestors_definite`` against the paper, the oracle, and 2.0.4.

Three independent checks on the same numbers, because agreeing with any one of
them alone would be weak evidence:

* the values Boichard, Maignel & Verrier (1997) PRINT in Tables I and II;
* the PyPedal-free oracle built from Appendices A and B;
* PyPedal 2.0.4, whose figures were measured by running it in the audit's
  Python 2.7 container and are recorded here as constants.

's repair also introduces two refusals -- R2 (half-founders) and R3
(reference population not an antichain) -- and a refusal that never fires is
indistinguishable from one that is broken, so both are exercised.
"""
import os
import unittest

from _pedhelpers import owned_temp_dir, corpus, load_corpus, load_corpus_from_path
from oracles import boichard_f_a, boichard_read

from PyPedal import pyp_errors, pyp_metrics

# Two-decimal published values; tolerance is the rounding half-width.
TWO_DP = 5e-2 + 1e-9

# Measured by running probe_boichard.py under pypedal-py27:audit against the
# read-only 2.0.4 checkout. Recorded rather than recomputed so the suite does
# not require Docker.
PYPEDAL_2_0_4 = {
    "boichard2a.ped": 2.0,
    "boichard_fig1.ped": 4.444444444444445,
    "boichard_fig2.ped": 2.9411764705882346,
}


class TestAgainstThePublishedTables(unittest.TestCase):

    def test_table_i_f_a(self):
        """Boichard Table I, p.9: marginal contributions imply f_a = 1/0.225."""
        ped = load_corpus("boichard_fig1.ped", "asdg")
        got = pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertAlmostEqual(1.0 / 0.225, got, places=9)

    def test_table_ii_f_a(self):
        """Boichard Table II, p.12, "Total" row: published f_a = 2.94."""
        ped = load_corpus("boichard_fig2.ped", "asdg")
        got = pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertAlmostEqual(2.94, got, delta=TWO_DP)

    def test_f_a_is_never_below_one(self):
        """
        The invariant that exposed the defect. f_a is the reciprocal of a sum of
        squared probabilities and cannot be below 1; the routine used to return
        0.0816.
        """
        for name in PYPEDAL_2_0_4:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, "asdg")
                self.assertGreaterEqual(
                    pyp_metrics.a_effective_ancestors_definite(ped), 1.0)


class TestAgainstTheIndependentOracle(unittest.TestCase):
    """
    Production and the PyPedal-free oracle implement the same appendix by
    different routes -- the oracle parses the pedigree itself and shares no code
    with PyPedal. Agreement to machine precision is the check.
    """

    def test_production_matches_the_oracle(self):
        for name in ("boichard_fig1.ped", "boichard_fig2.ped", "boichard2a.ped"):
            with self.subTest(pedigree=name):
                rows, gens = boichard_read(corpus(name), gen_col=3)
                newest = max(gens.values(), key=float)
                reference = [a for a, _, _ in rows if gens[a] == newest]
                want, _, _ = boichard_f_a(rows, reference)

                ped = load_corpus(name, "asdg")
                got = pyp_metrics.a_effective_ancestors_definite(ped)
                self.assertAlmostEqual(want, got, places=12)


class TestAgainstPyPedal204(unittest.TestCase):
    """
    MEASURED, and the result is the opposite of what was predicted.

    The adjudication expected the repair to introduce a permanent divergence
    from 2.0.4, on the grounds that 2.0.4's Appendix B carries index-domain
    defects. It does carry them -- but they live in its half-founder branch,
    and  now REFUSES half-founder pedigrees outright. On everything the
    repaired routine will actually accept, 2.0.4 was already right.

    So commit 4 introduced no divergence at all: it moved Python 3 INTO
    agreement with both 2.0.4 and the paper, from a value (0.0816) that was
    impossible.
    """

    def test_bit_exact_with_2_0_4(self):
        for name, want in PYPEDAL_2_0_4.items():
            with self.subTest(pedigree=name):
                ped = load_corpus(name, "asdg")
                got = pyp_metrics.a_effective_ancestors_definite(ped)
                self.assertEqual(
                    want, got,
                    f"{name}: expected bit-exact agreement with PyPedal 2.0.4")


class TestRefusals(unittest.TestCase):
    """
    R3 is a refusal, not a computation, and an untested refusal is
    indistinguishable from a broken one.

    R2 USED TO BE ONE TOO, and the two tests that pinned it are now the two
    value tests below -- rewritten rather than deleted, and deliberately kept
    on the SAME fixture, so that what changed about this pedigree is visible in
    one place instead of the old expectation simply vanishing.
    """

    #: The pedigree 's refusal was demonstrated on. Animal 4 has a known
    #: sire and an unknown dam.
    HALF_FOUNDER = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n"

    @staticmethod
    def _write(text, pedformat="asdg"):
        tmp = owned_temp_dir(prefix="pypedal_boichard_")
        path = os.path.join(tmp, "fixture.ped")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return load_corpus_from_path(path, pedformat)

    def test_r2_computes_a_half_founder_pedigree(self):
        """
        WAS test_r2_refuses_a_half_founder_pedigree.

        Appendix A step 4 halves a half-founder's contribution and Appendix B
        never restates it, so PyPedal refused rather than picking a reading.
        The reading was then adjudicated -- Reading C, phantom completion ahead
        of an unchanged Appendix B, mathematically implied and independently
        supported, NOT source-explicit Appendix-B text -- and implemented.

        f_a = 1.6 here; the derivation is written out in
        tests/test_boichard_half_founder.py.
        """
        ped = self._write(self.HALF_FOUNDER)
        self.assertAlmostEqual(
            1.6, pyp_metrics.a_effective_ancestors_definite(ped), places=12)

    def test_r2_no_longer_raises_a_typed_refusal(self):
        """
        WAS test_r2_names_the_offending_animal.

        The blocked item is closed, so nothing may still raise it. Asserting
        the absence explicitly matters because a refusal reintroduced by a
        merge would otherwise show up only as a distant test failing.
        """
        ped = self._write(self.HALF_FOUNDER)
        try:
            pyp_metrics.a_effective_ancestors_definite(ped)
        except pyp_errors.PyPedalNotImplementedError as exc:  # pragma: no cover
            self.fail("half-founder pedigrees must not be refused here "
                      f"(blocked_item={exc.blocked_item!r})")

    def test_r3_refuses_a_reference_population_that_is_not_an_antichain(self):
        """
        Animal 5 and its own offspring 6 are both in generation 2, so excluding
        reference-population members from selection would silently discard a
        real contribution rather than an uninteresting one.
        """
        ped = self._write("1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 5 4 2\n")
        with self.assertRaises(pyp_errors.PyPedalError) as ctx:
            pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertIn("antichain", str(ctx.exception))

    def test_an_ordinary_pedigree_is_not_refused(self):
        """
        Guard on the guards: refusals that fire on everything would be a
        regression dressed as caution.
        """
        ped = self._write("1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 3 4 2\n")
        self.assertGreaterEqual(
            pyp_metrics.a_effective_ancestors_definite(ped), 1.0)


class TestEffectiveFoundersAppendixA(unittest.TestCase):
    """
    . ``a_effective_founders_boichard`` against Appendix A and eq. 1.

    Two defects were repaired, both in the half-founder path and both shared
    with PyPedal 2.0.4: the halving was applied before propagation rather than
    after it, and half-founders were excluded from the founder set even though
    article p.7 says "when an animal has only one known parent, the animal is
    considered as a founder".
    """

    # Measured under pypedal-py27:audit; identical in the repaired Python 3,
    # because none of these pedigrees contains a half-founder.
    PYPEDAL_2_0_4 = {
        "boichard2a.ped": 4.0,
        "boichard_fig1.ped": 5.031446540880503,
        "boichard_fig2.ped": 5.5555555555555545,
    }

    def test_table_ii_f_e(self):
        """Boichard Table II, p.12, "Total" row: published f_e = 5.6."""
        ped = load_corpus("boichard_fig2.ped", "asdg")
        self.assertAlmostEqual(
            5.6, pyp_metrics.a_effective_founders_boichard(ped), delta=TWO_DP)

    def test_half_founder_free_pedigrees_are_unchanged(self):
        """
        The repair touches only the half-founder path, so every pedigree
        without one must be bit-exact with 2.0.4 -- which is also the evidence
        that the repair did not disturb the ordinary case.
        """
        for name, want in self.PYPEDAL_2_0_4.items():
            with self.subTest(pedigree=name):
                ped = load_corpus(name, "asdg")
                self.assertEqual(
                    want, pyp_metrics.a_effective_founders_boichard(ped))

    def test_founder_contributions_sum_to_one_on_a_half_founder_pedigree(self):
        """
        The defining check, and the one that failed before.

        Appendix A step 4 says dividing by N is what makes founder
        contributions sum to one. On this pedigree -- 1 x 2 -> 3, 3 x unknown
        -> 4, 3 x 4 -> {5, 6} -- the old code summed to 0.625, because animal 4
        was not counted as a founder and its retained contribution went
        nowhere. Nothing normalises the vector, so 1.0 here is a real result.
        """
        ped = TestRefusals._write(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        reference = [int(a.animalID) for a in ped.pedigree if str(a.gen) == "2"]
        q, founders = pyp_metrics.boichard_probabilities_of_gene_origin(
            ped, reference)
        self.assertEqual([1, 2, 4], sorted(founders),
                         "the half-founder must be in the founder set (p.7)")
        self.assertAlmostEqual(1.0, sum(q[f] for f in founders), places=12)

    def test_half_founder_contributions_are_the_paper_s(self):
        """
        Hand-derivable, so it is asserted exactly. The known parent of the
        half-founder receives half of its FULL contribution -- that is what the
        step 3 / step 4 ordering decides -- and the half-founder keeps the
        other half as its own founder contribution.
        """
        ped = TestRefusals._write(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        reference = [int(a.animalID) for a in ped.pedigree if str(a.gen) == "2"]
        q, _ = pyp_metrics.boichard_probabilities_of_gene_origin(ped, reference)
        for animal, want in {1: 0.375, 2: 0.375, 3: 0.75, 4: 0.25}.items():
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, q[animal], places=12)

    def test_f_e_on_a_half_founder_pedigree_is_possible(self):
        """
        f_e <= k. The old code returned 5.12 while counting only two founders,
        which is impossible; the repaired routine returns 2.909091 over three
        founder sources.
        """
        ped = TestRefusals._write(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        got = pyp_metrics.a_effective_founders_boichard(ped)
        self.assertAlmostEqual(2.909090909090909, got, places=9)
        self.assertLessEqual(got, 3.0)

    def test_half_founders_are_not_refused_here(self):
        """
         is settled by Appendix A, so the founders routine ACCEPTS
        half-founder pedigrees. Only the ancestors routines refuse them, and
        only because Appendix B is silent where Appendix A is explicit. The two
        must not be conflated.
        """
        ped = TestRefusals._write(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        self.assertGreater(pyp_metrics.a_effective_founders_boichard(ped), 0.0)


class TestCallerPedigreeIsNotConsumed(unittest.TestCase):
    """
    Appendix B step 3 deletes the parents of each selected ancestor. That is a
    legitimate step of the algorithm and it must not reach the caller (audit
    ). The engine applies it to local arrays, so this holds by
    construction -- asserted anyway, since it held "by construction" before too.
    """

    def test_repeated_calls_agree(self):
        ped = load_corpus("boichard_fig1.ped", "asdg")
        first = pyp_metrics.a_effective_ancestors_definite(ped)
        second = pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertEqual(first, second)

    def test_parentage_survives(self):
        ped = load_corpus("boichard_fig1.ped", "asdg")
        before = [(int(a.animalID), int(a.sireID), int(a.damID)) for a in ped.pedigree]
        pyp_metrics.a_effective_ancestors_definite(ped)
        after = [(int(a.animalID), int(a.sireID), int(a.damID)) for a in ped.pedigree]
        self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()
