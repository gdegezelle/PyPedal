"""
The paper-derived Boichard oracle, checked against the paper's published tables.

This validates the ORACLE, not PyPedal. Nothing here imports ``pyp_metrics``.
The point of the sequencing is that the instrument is calibrated against
published numbers before it is used to judge production code, so that a later
disagreement is evidence about production rather than about the oracle.

What the published examples DO establish: Appendix A, Appendix B, eq. 1, the
sum-to-one properties, and the exact f_e / f_a values of Tables I and II.

What they do NOT establish, and what is measured here instead of assumed: the
three residual semantics R1 (tie-breaking), R2 (half-founders inside Appendix
B) and R3 (reference-population members as candidates). Neither published
example contains a half-founder, and in both the population under study is
disjoint from the numbered ancestors, so R2 and R3 are simply not exercised.
R1 *is* exercised, and the result is more interesting than expected -- see
``TestTieBreakingIsNotDeterminedByTheSource``.
"""
import unittest

from _pedhelpers import corpus
from oracles import (
    BOICHARD_TIE_BREAKS,
    boichard_bounds,
    boichard_check_topological,
    boichard_f_a,
    boichard_f_e,
    boichard_founders,
    boichard_half_founders,
    boichard_marginal_contributions,
    boichard_read,
)

# Published values are printed to 2-3 significant decimals; the comparison
# tolerance is the corresponding rounding half-width, inclusive.
TWO_DP = 5e-2 + 1e-9
THREE_DP = 5e-4 + 1e-9


def _fig1():
    rows, gens = boichard_read(corpus("boichard_fig1.ped"), gen_col=3)
    reference = [a for a, _, _ in rows if gens[a] == "2"]
    return rows, reference


def _fig2(reference):
    rows, _ = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
    return rows, reference


class TestFixturesAreWellFormed(unittest.TestCase):

    def test_topologically_ordered(self):
        for name in ("boichard_fig1.ped", "boichard_fig2.ped"):
            with self.subTest(pedigree=name):
                rows, _ = boichard_read(corpus(name), gen_col=3)
                self.assertTrue(boichard_check_topological(rows))

    def test_neither_published_example_has_a_half_founder(self):
        """
        This is why the published examples cannot settle  or 's
        half-founder handling: they contain no animal with exactly one known
        parent. Reproducing Table I and Table II is therefore no evidence at
        all about the halving rule, and must never be cited as such.
        """
        for name in ("boichard_fig1.ped", "boichard_fig2.ped"):
            with self.subTest(pedigree=name):
                rows, _ = boichard_read(corpus(name), gen_col=3)
                self.assertEqual([], boichard_half_founders(rows))


class TestAppendixATableI(unittest.TestCase):
    """
    Boichard Table I, article p.9. The "Step 1" column is the raw probability
    of gene origin, which is exactly Appendix A.
    """

    # Table I, Step 1 -- every animal, not just the founders.
    PUBLISHED_STEP_1 = {1: 0.100, 2: 0.225, 3: 0.125, 4: 0.100,
                        5: 0.250, 6: 0.150, 7: 0.300, 8: 0.300}

    def test_founder_contributions_match_step_1(self):
        rows, reference = _fig1()
        _, q, _ = boichard_f_e(rows, reference)
        for animal in boichard_founders(rows):
            with self.subTest(animal=animal):
                self.assertAlmostEqual(self.PUBLISHED_STEP_1[animal], q[animal],
                                       delta=THREE_DP)

    def test_founder_contributions_sum_to_one_without_being_normalised(self):
        """
        Appendix A step 4 says dividing by N is what makes founder
        contributions sum to one. An UNFORCED 1.0 is therefore a real check on
        the half-founder handling and the propagation, not a tautology.
        """
        rows, reference = _fig1()
        _, _, q_sum = boichard_f_e(rows, reference)
        self.assertAlmostEqual(1.0, q_sum, places=12)

    def test_founder_set_includes_every_parentless_animal(self):
        rows, _ = _fig1()
        self.assertEqual([1, 2, 3, 4, 6, 8], sorted(boichard_founders(rows)))


class TestAppendixBTableI(unittest.TestCase):
    """Boichard Table I again -- the marginal contributions and f_a."""

    # The "Marginal contribution" column, article p.9.
    PUBLISHED_MARGINAL = {1: 0.100, 2: 0.150, 3: 0.000, 4: 0.100,
                          5: 0.050, 6: 0.000, 7: 0.300, 8: 0.300}
    PUBLISHED_F_A = 1.0 / 0.225      # = 4.444..., implied by the column

    def test_f_a_matches_the_published_column(self):
        rows, reference = _fig1()
        got, _, _ = boichard_f_a(rows, reference)
        self.assertAlmostEqual(self.PUBLISHED_F_A, got, places=9)

    def test_marginal_contributions_sum_to_one_without_being_normalised(self):
        """p.8: "the marginal contributions over all ancestors sum to one"."""
        rows, reference = _fig1()
        _, _, p_sum = boichard_f_a(rows, reference)
        self.assertAlmostEqual(1.0, p_sum, places=12)

    def test_number_of_contributors_does_not_exceed_the_founder_count(self):
        """p.8: "The number of ancestors with a positive contribution is less
        than or equal to the total number of founders"."""
        rows, reference = _fig1()
        _, order, _ = boichard_f_a(rows, reference)
        positive = [a for a, v in order if v > 0]
        self.assertLessEqual(len(positive), len(boichard_founders(rows)))

    def test_the_multiset_of_contributions_is_the_published_one(self):
        """
        Which ANIMAL receives the 0.050 depends on a tie the paper does not
        adjudicate (see the next class), so the per-animal column is not
        tie-break-invariant. The multiset of values is.
        """
        rows, reference = _fig1()
        for tie_break in BOICHARD_TIE_BREAKS:
            with self.subTest(tie_break=tie_break):
                _, order, _ = boichard_f_a(rows, reference, tie_break=tie_break)
                got = sorted(round(v, 6) for _, v in order)
                want = sorted(round(v, 6)
                              for v in self.PUBLISHED_MARGINAL.values() if v > 0)
                self.assertEqual(want, got)


class TestTieBreakingIsNotDeterminedByTheSource(unittest.TestCase):
    """
    , measured on the paper's own example.

    Table I contains three ties. The paper resolves them INCONSISTENTLY with
    respect to animal ID:

        round 1:  7 vs 8 at 0.300  -> picks 7   (lower id)
        round 4:  1 vs 4 at 0.100  -> picks 1   (lower id)
        round 6:  3 vs 5 at 0.050  -> picks 5   (HIGHER id)

    So no ID-based rule reproduces the paper's selection sequence, and the
    earlier supposition that "Table I breaks its tie the same way as our
    convention" was drawn from round 1 alone and does not survive the other
    two. The convention is a convention.

    What the ties do NOT affect is the statistic: f_a is identical under both
    conventions here, because the tie is between two candidates holding equal
    contributions. The choice changes only which animal is credited.
    """

    def test_f_a_is_invariant_under_the_tie_break_convention(self):
        rows, reference = _fig1()
        values = set()
        for tie_break in BOICHARD_TIE_BREAKS:
            got, _, _ = boichard_f_a(rows, reference, tie_break=tie_break)
            values.add(round(got, 12))
        self.assertEqual(1, len(values),
                         f"f_a should not depend on the tie-break here, got {values}")

    def test_the_conventions_credit_different_animals(self):
        """
        If this ever stopped being true the convention would have become
        irrelevant, and the R1 caveat could be retired. It is true today.
        """
        rows, reference = _fig1()
        credited = {}
        for tie_break in BOICHARD_TIE_BREAKS:
            _, order, _ = boichard_f_a(rows, reference, tie_break=tie_break)
            credited[tie_break] = {a for a, v in order if v > 0}
        self.assertNotEqual(credited["lowest_id"], credited["highest_id"])

    def test_neither_convention_reproduces_the_papers_selection_order(self):
        """
        Published order: 7, 8, 2, 1, 4, 5, 3.

        ``lowest_id`` matches the first five picks and then diverges (it takes
        3, the paper takes 5). ``highest_id`` takes 5 but reverses rounds 1 and
        4. Recorded so that neither is mistaken for the paper's rule.
        """
        rows, reference = _fig1()
        published_order = [7, 8, 2, 1, 4, 5, 3]
        for tie_break in BOICHARD_TIE_BREAKS:
            with self.subTest(tie_break=tie_break):
                _, order, _ = boichard_f_a(rows, reference, tie_break=tie_break)
                self.assertNotEqual(published_order, [a for a, _ in order])


class TestAppendixTableII(unittest.TestCase):
    """
    Boichard Table II, article p.12 -- three populations over one pedigree.

    f_e and f_a are deterministic and are checked. N_g (2.5, 1.8, 1.1) is a
    Monte-Carlo gene-drop figure and is deliberately NOT checked here.
    """

    # population -> (reference animals, published f_e, published f_a)
    ROWS = {
        "total": (list(range(7, 15)) + [19, 20], 5.6, 2.94),
        "family 1": (list(range(7, 15)), 4.0, 2.0),
        "family 2": ([19, 20], 2.0, 2.0),
    }

    def test_published_f_e(self):
        for label, (reference, want, _) in self.ROWS.items():
            with self.subTest(population=label):
                rows, ref = _fig2(reference)
                got, _, _ = boichard_f_e(rows, ref)
                self.assertAlmostEqual(want, got, delta=TWO_DP)

    def test_published_f_a(self):
        for label, (reference, _, want) in self.ROWS.items():
            with self.subTest(population=label):
                rows, ref = _fig2(reference)
                got, _, _ = boichard_f_a(rows, ref)
                self.assertAlmostEqual(want, got, delta=TWO_DP)

    def test_both_vectors_sum_to_one_on_every_population(self):
        for label, (reference, _, _) in self.ROWS.items():
            with self.subTest(population=label):
                rows, ref = _fig2(reference)
                _, _, q_sum = boichard_f_e(rows, ref)
                _, _, p_sum = boichard_f_a(rows, ref)
                self.assertAlmostEqual(1.0, q_sum, places=12)
                self.assertAlmostEqual(1.0, p_sum, places=12)

    def test_f_a_never_exceeds_f_e(self):
        """
        Stated at p.11: the effective number of ancestors is an intermediate
        estimate, at most the effective number of founders, because it accounts
        for bottlenecks the founder count ignores.
        """
        for label, (reference, _, _) in self.ROWS.items():
            with self.subTest(population=label):
                rows, ref = _fig2(reference)
                fe, _, _ = boichard_f_e(rows, ref)
                fa, _, _ = boichard_f_a(rows, ref)
                self.assertLessEqual(fa, fe + 1e-9)


class TestAVectorInvariant(unittest.TestCase):
    """
    ``0 <= a(i) <= 1`` is structural, and the oracle raises rather than clamps.

    Exercised across every fixture and both tie-break conventions; the engine
    asserts internally on each round, so simply driving it to exhaustion is the
    test.
    """

    def test_engine_completes_without_violating_the_invariant(self):
        cases = [(*_fig1(),)] + [_fig2(r) for r, _, _ in
                                 (TestAppendixTableII.ROWS["total"],
                                  TestAppendixTableII.ROWS["family 1"],
                                  TestAppendixTableII.ROWS["family 2"])]
        for rows, reference in cases:
            for tie_break in BOICHARD_TIE_BREAKS:
                with self.subTest(n=len(rows), tie_break=tie_break):
                    order = list(boichard_marginal_contributions(
                        rows, reference, tie_break=tie_break))
                    self.assertTrue(order)


class TestBoundsAgainstTheExactValue(unittest.TestCase):
    """
    . The bounds on pp.9-10 are a truncation of the same sequence, so
    ``f_l <= f_a <= f_u`` must hold at every n, including both endpoints.
    """

    def test_bracket_the_exact_value_at_every_n(self):
        rows, reference = _fig1()
        exact, order, _ = boichard_f_a(rows, reference)
        for n in range(1, len(order) + 1):
            with self.subTest(n=n):
                f_l, f_u, meta = boichard_bounds(rows, reference, n)
                self.assertLessEqual(f_l, exact + 1e-9,
                                     f"lower bound above the exact value at n={n}")
                self.assertLessEqual(exact, f_u + 1e-9,
                                     f"exact value above the upper bound at n={n}")
                self.assertLessEqual(f_l, f_u + 1e-9)

    def test_the_final_endpoint_is_exact_not_a_division_by_zero(self):
        """
        Once the residual mass is gone the truncation has reached the exact
        answer, and neither residual term may be evaluated -- (f - n) is zero
        there. This is the endpoint safeguard, not a scientific choice.
        """
        rows, reference = _fig1()
        exact, order, _ = boichard_f_a(rows, reference)
        f_l, f_u, meta = boichard_bounds(rows, reference, len(order))
        self.assertTrue(meta["exact"])
        self.assertAlmostEqual(exact, f_l, places=12)
        self.assertAlmostEqual(exact, f_u, places=12)

    def test_bounds_converge_as_more_ancestors_are_taken(self):
        rows, reference = _fig1()
        _, order, _ = boichard_f_a(rows, reference)
        widths = []
        for n in range(1, len(order) + 1):
            f_l, f_u, _ = boichard_bounds(rows, reference, n)
            widths.append(f_u - f_l)
        self.assertAlmostEqual(0.0, widths[-1], places=12)
        self.assertLess(widths[-1], widths[0])


if __name__ == "__main__":
    unittest.main()
