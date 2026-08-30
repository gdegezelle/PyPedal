"""
The independent Boichard N_g oracle, checked against published tables.

This validates the ORACLE, not PyPedal. Nothing here imports ``pyp_metrics``.
The sequencing is the point: the instrument is calibrated against published
numbers before it is used to judge production, so that a later disagreement is
evidence about production rather than about the oracle.

What the published material DOES establish here:

* Lacy (1989) Table 1 column H -- 10 000 gene drops, article p.116 -- validates
  this module's SEGREGATION AND GENE-COUNTING MECHANICS.
* Boichard (1997) Table II, article p.12 -- validates eq. 2 itself, and its
  Family 2 row discriminates between the two candidate estimators.

What it does NOT establish, and what is measured rather than assumed:

* that Lacy's f_g equals Boichard's N_g. It does not, and
  ``TestLacyFgIsNotBoichardNg`` pins the disagreement.
* which averaging Boichard used. Only the Family 2 row separates the two
  readings; Total and Family 1 are consistent with both, and the tests say so.
"""
import unittest
from fractions import Fraction

import pytest
from _pedhelpers import corpus
from oracles import (
    NgExactEnumerationTooLarge,
    boichard_read,
    exact_ng,
    lacy_f_g,
    mc_ng,
    ng_check_bounds,
    ng_distribution,
    ng_distribution_bruteforce,
    ng_founder_genes,
    ng_n_binary_choices,
    renumber,
)
from oracles import (
    read_pedigree as ml_read_pedigree,
)

# Published values are printed to 2-3 decimals; tolerance is the corresponding
# rounding half-width, inclusive.
TWO_DP = 5e-2 + 1e-9
THREE_DP = 5e-4 + 1e-9


def _lacy_table1_pedigree(n1, n2):
    """
    Lacy Table 1: two founder PAIRS with n1 and n2 first-generation offspring
    respectively; the descendant population is those offspring.
    """
    rows = [(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]
    nid, ref = 5, []
    for _ in range(n1):
        rows.append((nid, 1, 2))
        ref.append(nid)
        nid += 1
    for _ in range(n2):
        rows.append((nid, 3, 4))
        ref.append(nid)
        nid += 1
    return rows, ref


#: (n1, n2, published f_e, published f_g, published H) -- Lacy Table 1, p.116.
LACY_TABLE_1 = [
    (2, 2, 4.00, 3.00, .812),
    (3, 3, 4.00, 3.50, .833),
    (4, 4, 4.00, 3.75, .844),
    (2, 4, 3.60, 3.21, .819),
    (4, 8, 3.60, 3.54, .840),
]


class TestSegregationMechanicsAgainstLacysPublishedH(unittest.TestCase):
    """
    Lacy's H column came from 10 000 gene drops. Reproducing it by exhaustive
    enumeration validates that this module drops genes and counts them the way
    the literature does.

    It says NOTHING about whether f_g equals N_g -- see
    ``TestLacyFgIsNotBoichardNg`` immediately below, which uses these very same
    pedigrees to show that it does not.
    """

    def test_published_h(self):
        for n1, n2, _fe, _fg, h in LACY_TABLE_1:
            with self.subTest(offspring=(n1, n2)):
                rows, ref = _lacy_table1_pedigree(n1, n2)
                got = exact_ng(rows, ref)["h"]
                self.assertAlmostEqual(h, got, delta=THREE_DP)


class TestLacyFgIsNotBoichardNg(unittest.TestCase):
    """
    Boichard p.11 calls N_g the quantity "called 'founder genome equivalent' by
    Lacy, 1989". That is a naming remark; the paper derives no equivalence.

    Lacy p.115 devalues each founder by an allele RETENTION probability r_i;
    Boichard eq. 2 uses realised gene FREQUENCIES. On Lacy's own Table 1
    pedigrees the two disagree well outside the published rounding, so no
    implementation may treat one as an oracle for the other.
    """

    #: Measured relative gaps (f_g - N_g)/f_g on Lacy's own Table 1 pedigrees:
    #: 8.5%, 12.6%, 13.5%, 12.1%, 10.9%. The published f_g values are printed to
    #: 2 dp, i.e. to within 0.5-1.6% -- so the disagreement is an order of
    #: magnitude larger than the rounding, on every row.
    MIN_RELATIVE_GAP = 0.05

    def test_published_f_g_and_exact_n_g_disagree(self):
        for n1, n2, _fe, f_g, _h in LACY_TABLE_1:
            with self.subTest(offspring=(n1, n2)):
                rows, ref = _lacy_table1_pedigree(n1, n2)
                n_g = exact_ng(rows, ref)["n_g"]
                self.assertLess(n_g, f_g)
                self.assertGreater((f_g - n_g) / f_g, self.MIN_RELATIVE_GAP)

    def test_new_lacy_ped_published_f_g_is_not_its_n_g(self):
        """
        ``new_lacy.ped`` IS Lacy's Appendix A worked example, published f_g 2.18.
        Its exact Boichard N_g over the same descendant set is materially lower.
        """
        rows, _gens = boichard_read(corpus("new_lacy.ped"))
        descendants = [a for a, s, d in rows if not (s == 0 and d == 0)]
        n_g = exact_ng(rows, descendants)["n_g"]
        self.assertAlmostEqual(1.8418, n_g, delta=1e-3)

        ped, _back = renumber(ml_read_pedigree(corpus("new_lacy.ped")))
        f_g, _r = lacy_f_g(ped)
        self.assertAlmostEqual(2.18, f_g, delta=TWO_DP)

        self.assertGreater(f_g - n_g, 0.3)


class TestBoichardTableII(unittest.TestCase):
    """
    Boichard Figure 2 / Table II, article pp.11-12. The fixture is the already
    committed ``corpus/boichard_fig2.ped``, whose ``g`` column encodes the
    paper's own population under study.
    """

    def _rows_and_reference(self):
        rows, gens = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
        return rows, [a for a, _s, _d in rows if gens[a] == "2"]

    def test_family_2_exactly_reproduces_published_1_1(self):
        rows, reference = self._rows_and_reference()
        f2 = [r for r in rows if r[0] >= 15]
        r2 = [a for a in reference if a >= 15]
        got = exact_ng(f2, r2)
        self.assertAlmostEqual(1.1, got["n_g"], delta=TWO_DP)

    def test_family_1_reproduces_published_1_8(self):
        rows, reference = self._rows_and_reference()
        f1 = [r for r in rows if r[0] <= 14]
        r1 = [a for a in reference if a <= 14]
        self.assertAlmostEqual(1.8, exact_ng(f1, r1)["n_g"], delta=TWO_DP)

    def test_total_reproduces_published_2_5(self):
        """
        The whole pedigree's definitional state space is 2**28. It is computed
        exactly anyway, because the factorisation collapses the eight full sibs
        in each family -- no sampling is involved.
        """
        rows, reference = self._rows_and_reference()
        got = exact_ng(rows, reference)
        self.assertAlmostEqual(2.5, got["n_g"], delta=TWO_DP)
        self.assertEqual(1 << 28, got["state_space"])

    def test_total_agrees_with_an_independent_convolution(self):
        """
        The two families use DISJOINT founder genes, so the pooled statistic can
        also be obtained by convolving the two family distributions. That is a
        different computation from the factorised enumeration above, and it must
        give the same exact answer.
        """
        rows, reference = self._rows_and_reference()
        f1 = [r for r in rows if r[0] <= 14]
        f2 = [r for r in rows if r[0] >= 15]
        r1 = [a for a in reference if a <= 14]
        r2 = [a for a in reference if a >= 15]
        d1 = ng_distribution(f1, r1)
        d2 = ng_distribution(f2, r2)
        n1, n2 = 2 * len(r1), 2 * len(r2)
        n = n1 + n2
        e_ng = Fraction(0)
        for s1, p1 in d1.items():
            for s2, p2 in d2.items():
                s = (s1 * n1 * n1 + s2 * n2 * n2) / Fraction(n * n)
                e_ng += p1 * p2 / (2 * s)
        self.assertEqual(exact_ng(rows, reference)["n_g_exact"], e_ng)


class TestWhichAveragingBoichardUsed(unittest.TestCase):
    """
    Boichard says only that the procedure "is replicated to obtain an accurate
    estimate of the parameter of interest", which does not say whether the
    replicates average N_g or average SUM f_k^2.

    Only the Family 2 row separates the two readings. This is stated as a test
    rather than as prose so that the narrowness of the evidence is visible.
    """

    def _family_2(self):
        rows, gens = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
        f2 = [r for r in rows if r[0] >= 15]
        r2 = [a for a, _s, _d in rows if a >= 15 and gens[a] == "2"]
        return exact_ng(f2, r2)

    def test_mean_of_per_replicate_ng_matches_the_published_value(self):
        self.assertAlmostEqual(1.1, self._family_2()["n_g"], delta=TWO_DP)

    def test_inverting_the_mean_of_sum_f2_does_not(self):
        got = self._family_2()["inv_mean_s"]
        self.assertAlmostEqual(1.0, got, delta=1e-12)
        self.assertGreater(abs(1.1 - got), TWO_DP)

    def test_the_other_two_rows_do_not_discriminate(self):
        """
        Both readings round to the published Total 2.5 and Family 1 1.8, so the
        adjudication rests on Family 2 alone. Recorded, not glossed over.
        """
        rows, gens = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
        reference = [a for a, _s, _d in rows if gens[a] == "2"]
        f1 = [r for r in rows if r[0] <= 14]
        r1 = [a for a in reference if a <= 14]
        for subject, want in ((exact_ng(f1, r1), 1.8),
                              (exact_ng(rows, reference), 2.5)):
            self.assertAlmostEqual(want, subject["n_g"], delta=TWO_DP)
            self.assertAlmostEqual(want, subject["inv_mean_s"], delta=TWO_DP)


class TestDerivedBounds(unittest.TestCase):
    """
    The two bounds this repository is entitled to assert, and the one it is not.

    ``N_g >= 1/2``  because SUM f_k = 1 forces SUM f_k^2 <= 1. Equality is
                    FIXATION of a single founder gene and is a legitimate
                    result, so ``N_g < 1`` must never be treated as invalid.

    ``N_g <= f_e``  per replicate, against f_e computed from the SAME realised
                    frequencies, by (a+b)^2 <= 2(a^2+b^2) with p_i = f_i1 + f_i2.

    Nothing here asserts a relation between E[N_g] and an NRM-based expected
    f_e; that does not follow from the derivation.
    """

    FIXTURES = {
        "lone founder, R = itself": ([(1, 0, 0)], [1]),
        "two founders, one offspring": ([(1, 0, 0), (2, 0, 0), (3, 1, 2)], [3]),
        "full sibs": ([(1, 0, 0), (2, 0, 0), (4, 1, 2), (5, 1, 2)], [4, 5]),
        "selfing": ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 3)], [4]),
        "half-founder": ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0)], [4]),
        "founder inside R": ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 0, 0)],
                             [3, 4, 5]),
    }

    def test_every_enumerated_outcome_satisfies_both_bounds(self):
        for name, (rows, reference) in self.FIXTURES.items():
            with self.subTest(fixture=name):
                ok, worst, slack = ng_check_bounds(rows, reference)
                self.assertTrue(ok)
                self.assertGreaterEqual(worst, 0.5)
                self.assertGreaterEqual(slack, 0.0)

    def test_fixation_gives_exactly_one_half_and_that_is_legal(self):
        """
        Selfing lets animal 4 draw the same gene twice, fixing it in R = {4}.
        N_g = 0.5 exactly. A postcondition demanding N_g >= 1 would reject a
        correct answer, which is why the plan's original ``N_g >= 1`` was wrong.
        """
        rows, reference = self.FIXTURES["selfing"]
        dist = ng_distribution(rows, reference)
        attainable = sorted(float(1 / (2 * s)) for s in dist)
        self.assertEqual(0.5, attainable[0])
        self.assertIn(0.5, attainable)

    def test_one_is_not_a_lower_bound(self):
        rows, reference = self.FIXTURES["selfing"]
        self.assertLess(exact_ng(rows, reference)["min_ng"], 1.0)

    def test_lone_founder_reference_is_exactly_one(self):
        """
        FG-1, stated so its expected value is unambiguous: R is the founder
        itself, its two genes sit at frequency 1/2 each, SUM f_k^2 = 1/2 and
        N_g = 1 with no drift involved at all.
        """
        got = exact_ng([(1, 0, 0)], [1])
        self.assertEqual(1.0, got["n_g"])
        self.assertEqual(1, got["state_space"])


class TestHalfFounderRepresentationsAreExactlyEquivalent(unittest.TestCase):
    """
    FG-10. The two research representations of an unknown parental slot are:

    ``phantom``  a distinct dummy founder carrying two unique genes, one of
                 which the half-founder samples (Baumung et al. 2015 p.102).
    ``slot``     one unique gene placed directly in the slot.

    They are EXACTLY distributionally equivalent for SUM f_k^2, because a
    phantom founder has exactly one offspring and therefore exactly one
    transmission opportunity: its untransmitted gene has frequency zero in R and
    contributes zero to the sum, and which of the two symmetric genes was
    transmitted is a relabelling.

    The comparison is of the full Fraction-valued distributions, not of means.
    """

    CASES = {
        "FG-6 known sire, unknown dam":
            ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0)], [4]),
        "FG-7 unknown sire, known dam":
            ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 0, 3)], [4]),
        "half-founder with descendants":
            ([(1, 0, 0), (2, 0, 0), (3, 1, 0), (4, 2, 3), (5, 2, 3)], [4, 5]),
        "half-founder itself in R":
            ([(1, 0, 0), (2, 1, 0)], [2]),
        "FG-9 several distinct unknown slots":
            ([(1, 0, 0), (2, 1, 0), (3, 0, 1), (4, 1, 0), (5, 2, 3), (6, 4, 5)],
             [5, 6]),
        "FG-8 both parents unknown (control, no slots)":
            ([(1, 0, 0), (2, 0, 0), (3, 1, 2)], [3]),
        "mrode.ped topology, R = all non-founders":
            ([(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 1, 0), (5, 3, 2), (6, 1, 5)],
             [4, 5, 6]),
    }

    def test_distributions_are_identical_not_merely_close(self):
        for name, (rows, reference) in self.CASES.items():
            with self.subTest(fixture=name):
                self.assertEqual(ng_distribution(rows, reference, "phantom"),
                                 ng_distribution(rows, reference, "slot"))

    def test_expectations_are_identical_as_exact_rationals(self):
        for name, (rows, reference) in self.CASES.items():
            with self.subTest(fixture=name):
                self.assertEqual(exact_ng(rows, reference, "phantom")["n_g_exact"],
                                 exact_ng(rows, reference, "slot")["n_g_exact"])

    def test_the_two_representations_really_do_differ_structurally(self):
        """
        Anti-vacuity: if the modes produced identical state spaces the equality
        above would be trivial. Wherever there IS an unknown slot, ``phantom``
        names two extra genes and one extra binary choice per slot.
        """
        rows, reference = self.CASES["FG-9 several distinct unknown slots"]
        self.assertEqual(8, len(ng_founder_genes(rows, "phantom")))
        self.assertEqual(5, len(ng_founder_genes(rows, "slot")))
        self.assertGreater(ng_n_binary_choices(rows, "phantom"),
                           ng_n_binary_choices(rows, "slot"))

    def test_control_case_has_no_slots_at_all(self):
        rows, _ = self.CASES["FG-8 both parents unknown (control, no slots)"]
        self.assertEqual(ng_founder_genes(rows, "phantom"),
                         ng_founder_genes(rows, "slot"))


class TestTheFastPathIsTheDefinition(unittest.TestCase):
    """
    ``distribution()`` factorises exchangeable full sibs into binomials;
    ``distribution_bruteforce()`` enumerates all 2**b outcomes and IS the
    definition. Anywhere brute force can run, the two must agree exactly --
    Fraction for Fraction, key for key -- or the fast path is not the thing it
    claims to compute.
    """

    CASES = [
        ([(1, 0, 0)], [1]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2)], [3]),
        ([(1, 0, 0), (2, 0, 0), (4, 1, 2), (5, 1, 2)], [4, 5]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 3)], [4]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0)], [4]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 0, 3)], [4]),
        ([(1, 0, 0), (2, 1, 0), (3, 0, 1), (4, 1, 0), (5, 2, 3), (6, 4, 5)], [5, 6]),
        ([(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 1, 0), (5, 3, 2), (6, 1, 5)],
         [4, 5, 6]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 3, 4), (6, 3, 4)], [5, 6]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 0, 0)], [3, 4, 5]),
        ([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 1, 2), (6, 3, 4)],
         [4, 5, 6]),
    ]

    def test_factorised_equals_brute_force_exactly(self):
        for rows, reference in self.CASES:
            for mode in ("phantom", "slot"):
                with self.subTest(rows=rows, mode=mode):
                    self.assertEqual(
                        ng_distribution_bruteforce(rows, reference, mode),
                        ng_distribution(rows, reference, mode))

    def test_lacy_table_1_row_agrees_with_brute_force(self):
        """The (2,2) row is 2**8 outcomes, small enough to check both ways."""
        rows, reference = _lacy_table1_pedigree(2, 2)
        self.assertEqual(ng_distribution_bruteforce(rows, reference),
                         ng_distribution(rows, reference))


class TestOracleRefusesRatherThanSampling(unittest.TestCase):

    def test_exact_enumeration_refuses_past_its_cap(self):
        """
        A chain of 40 non-exchangeable animals: nothing to factorise, 2**80
        outcomes. The oracle must refuse rather than quietly sample.
        """
        rows = [(1, 0, 0), (2, 0, 0)]
        for a in range(3, 43):
            rows.append((a, a - 2, a - 1))
        with self.assertRaises(NgExactEnumerationTooLarge):
            exact_ng(rows, [42])

    def test_brute_force_refuses_past_its_cap(self):
        rows, gens = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
        reference = [a for a, _s, _d in rows if gens[a] == "2"]
        with self.assertRaises(NgExactEnumerationTooLarge):
            ng_distribution_bruteforce(rows, reference)

    def test_monte_carlo_is_locally_seeded_and_repeatable(self):
        rows = [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 3, 4), (6, 3, 4)]
        a = mc_ng(rows, [5, 6], 2000, seed=11)
        b = mc_ng(rows, [5, 6], 2000, seed=11)
        c = mc_ng(rows, [5, 6], 2000, seed=12)
        self.assertEqual(a["n_g"], b["n_g"])
        self.assertNotEqual(a["n_g"], c["n_g"])

    def test_monte_carlo_agrees_with_exact_enumeration(self):
        rows = [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 3, 4), (6, 3, 4)]
        want = exact_ng(rows, [5, 6])["n_g"]
        got = mc_ng(rows, [5, 6], 40000, seed=2026)
        self.assertLess(abs(got["n_g"] - want), 4 * got["sem"])

    def test_empty_reference_population_is_refused(self):
        with self.assertRaises(ValueError):
            exact_ng([(1, 0, 0), (2, 0, 0), (3, 1, 2)], [])


class TestOracleDoesNotImportProduction(unittest.TestCase):
    """
    Checked over the parsed IMPORT statements, not over the raw text: the module
    docstring names both ``pyp_metrics`` and ``oracle_lacy`` on purpose, to say
    what it is not. A substring check would forbid the explanation along with
    the dependency.
    """

    def _imported_modules(self):
        import ast

        import oracle_founder_genomes as mod
        with open(mod.__file__) as fh:
            tree = ast.parse(fh.read())
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.add(node.module or "")
        return names

    def test_module_never_imports_pypedal(self):
        for name in self._imported_modules():
            self.assertFalse(name.split(".")[0] in ("PyPedal", "pyp_metrics"),
                             "unexpected production import: %s" % name)

    def test_module_does_not_use_lacy_f_g_as_its_oracle(self):
        self.assertNotIn("oracle_lacy", self._imported_modules())

    def test_the_import_check_can_actually_fail(self):
        """Anti-vacuity: the parser really does see this module's imports."""
        self.assertIn("oracle_boichard", self._imported_modules())


if __name__ == "__main__":
    pytest.main([__file__])
