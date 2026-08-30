"""Half-founders in gene dropping: the gene source of a missing parental side.

THE RULE, and where it comes from
---------------------------------
Suwanlee et al. (2007) p.490 drops alleles from TWO parents and never says what a
missing parental side contributes. That silence is a historical fact and nothing
here revises it: it is why ``dropped_ancestral_inbreeding`` refused half-founders
under blocked item **b/HF**.

The rule is settled instead by **Baumung R., Farkas J., Boichard D., Meszaros G.,
Solkner J. & Curik I. (2015)**, "GRain: a computer program to calculate ancestral
and partial inbreeding coefficients using a gene dropping approach", *J Anim Breed
Genet* 132:100-108, **p.102, "Computational strategy"**, source-explicit:

    "For each half-founder, that is an animal with just one parent known, a dummy
    founder is created and the unknown second parent is assigned an artificially
    created new identification number and also provided with two unique alleles."

Corroborated p.104. Evidence status: **SOURCE-EXPLICIT**, but a *method statement,
not a derivation* -- GRAIN offers no justification.

Operationally: each unknown parental slot becomes its own ordinary founder, with
two unique alleles, never flagged, transmitting one gamete at p = 0.5, unrelated
to every recorded founder and to every other dummy. Two half-founders that both
record the sentinel 0 are **not** the same individual -- the pedigree carries no
information that they are.

WHAT ANCHORS THESE TESTS
------------------------
The primary anchors are the ANALYTIC FIXTURES below, whose expectations are
derived by hand from the rule and stated in each docstring before any code runs.
Agreement between production and ``the independent oracle`` is
CORROBORATION, not proof: they are two implementations of the same published
rule, and two structurally similar implementations can be wrong together. Where a
value is exact it is asserted exactly, with no tolerance at all.

GRAIN Figure 3 / Table 1 is a **candidate oracle whose validation is DEFERRED**:
the 45-animal pedigree exists only as a rendered figure and no reproducible
topology is available in this repository. It is deliberately not reconstructed
from the image or back-solved from the published values, and nothing here claims
to reproduce it.

ID COORDINATES -- read this before adding a fixture
---------------------------------------------------
``load_corpus_from_path`` renumbers, and renumbering SORTS FOUNDERS FIRST, so the
animal called 3 in a fixture file is generally not ``animalID`` 3. Every
expectation below is therefore keyed on the **originalID** -- the number as
written in the fixture text -- and results are translated with ``by_original``.
Asserting on raw ``animalID`` literals would silently test a different animal.
"""
import math
import os
import sys
import tempfile
import unittest

from PyPedal import pyp_metrics

from _pedhelpers import corpus, load_corpus, load_corpus_from_path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "oracles"))
from oracle_suwanlee import gene_drop, read_pedigree      # noqa: E402

# ---------------------------------------------------------------------------
# Predeclared experiment constants. Written down BEFORE the first run and not
# adjusted afterwards. Taken from the pattern already established by
# tests/test_suwanlee_production.py rather than invented here.
# ---------------------------------------------------------------------------

#: For the stochastic expectations (fixtures E and G).
ROUNDS, LOCI, SEED = 40, 5000, 20260820
N = ROUNDS * LOCI

#: For the EXACT expectations. Any (rounds, loci) gives the same answer -- the
#: value is exactly zero at every locus of every replicate -- so these are small
#: on purpose. A tolerance would be meaningless here and none is used.
EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED = 20, 500, 424242

#: Independent seeds for the oracle, which computes one replicate block per call.
ORACLE_SEEDS = (1, 2, 3, 4, 5, 6, 7, 8)
ORACLE_LOCI = 25_000


def half_width(expected, p_flagged):
    """
    Five standard errors of the mean over N loci, for a per-locus contribution
    that is 0 or 1/2 with probability ``p_flagged`` of being 1/2.

    Five sigma, not three: predeclared, and wide enough that this suite cannot
    flake on a legitimate implementation while still being far narrower than the
    gap to any rival reading of the rule.
    """
    variance = p_flagged * 0.25 - expected ** 2
    return 5.0 * math.sqrt(variance) / math.sqrt(N)


# ---------------------------------------------------------------------------
# Fixtures. Each is stated in its own coordinates, with the derivation.
# ---------------------------------------------------------------------------

#: A -- known sire, unknown dam. The minimal half-founder.
FIXTURE_A = "1 0 0\n2 0 0\n3 1 2\n4 3 0\n"

#: B -- unknown sire, known dam. The mirror of A.
FIXTURE_B = "1 0 0\n2 0 0\n3 1 2\n4 0 3\n"

#: C -- two INDEPENDENT half-founders, 3 and 4, whose unknown dams must not be
#: confused with one another.
FIXTURE_C = "1 0 0\n2 0 0\n3 1 0\n4 2 0\n5 3 4\n6 0 0\n7 5 6\n"

#: E -- a flag that ARISES after the half-founder and reaches a descendant.
#: 4 and 5 are full sibs of (3, 2); 6 is their offspring; 8 observes 6's flags.
FIXTURE_E = "1 0 0\n2 0 0\n3 1 0\n4 3 2\n5 3 2\n6 4 5\n7 0 0\n8 6 7\n"

#: G -- a mixed pedigree: founders, full-parent descendants and half-founders
#: together, so the three cases are exercised in one drop rather than separately.
FIXTURE_MIXED = ("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 0 2\n"
                 "6 3 4\n7 3 5\n8 6 7\n")


def write_pedigree(text, name="fixture.ped"):
    """Materialise a fixture. Output is confined to a tmpdir by the loader."""
    tmp = tempfile.mkdtemp(prefix="pypedal_hf_")
    path = os.path.join(tmp, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path


def load_fixture(text):
    return load_corpus_from_path(write_pedigree(text), "asd")


def by_original(ped, result):
    """
    Re-key a production result from ``animalID`` onto the originalID written in
    the fixture file, which is the coordinate system every derivation uses.
    """
    return {int(animal.originalID): float(result[animal.animalID])
            for animal in ped.pedigree}


def production(ped, rounds=ROUNDS, loci=LOCI, seed=SEED):
    return by_original(ped, pyp_metrics.dropped_ancestral_inbreeding(
        ped, rounds=rounds, loci=loci, seed=seed))


def oracle(rows, loci=ORACLE_LOCI, seeds=ORACLE_SEEDS):
    """
    Mean over independent seeds of the INDEPENDENT oracle. Already keyed on the
    IDs as written in the file, because the oracle does not renumber.
    """
    runs = [gene_drop(rows, loci=loci, seed=seed) for seed in seeds]
    return {animal: sum(run[animal] for run in runs) / len(runs)
            for animal in runs[0]}


def rows_of(text):
    return read_pedigree(write_pedigree(text))


# ---------------------------------------------------------------------------
# The independent oracle, against the analytic fixtures. Green immediately:
# the oracle implements the GRAIN rule and does not refuse.
# ---------------------------------------------------------------------------

class TestOracleAgainstTheAnalyticFixtures(unittest.TestCase):
    """
    The oracle earns its place here by satisfying the hand derivations on its
    own, with no production involvement at all. It is checked BEFORE production
    exists so that it cannot later be tuned into agreement with production.
    """

    def test_a_known_sire_unknown_dam_is_exactly_zero(self):
        got = oracle(rows_of(FIXTURE_A))
        self.assertEqual(0.0, got[4])

    def test_b_unknown_sire_known_dam_is_exactly_zero(self):
        got = oracle(rows_of(FIXTURE_B))
        self.assertEqual(0.0, got[4])

    def test_c_independent_dummies_never_share_identity(self):
        got = oracle(rows_of(FIXTURE_C))
        self.assertEqual(0.0, got[7])

    def test_e_the_derived_descendant_expectation(self):
        got = oracle(rows_of(FIXTURE_E))
        self.assertLessEqual(abs(got[8] - 0.125), 0.01,
                             "oracle f_a(8) = %.6f, derived 0.125" % got[8])

    def test_mrode_animal_six(self):
        rows = read_pedigree(corpus("mrode.ped"))
        got = oracle(rows)
        for animal in (1, 2, 3, 4, 5):
            self.assertEqual(0.0, got[animal])
        self.assertLessEqual(abs(got[6] - 0.0625), 0.01,
                             "oracle f_a(6) = %.6f, derived 0.0625" % got[6])

    def test_the_oracle_no_longer_refuses_a_half_founder(self):
        """
        The behaviour this commit changes. The oracle used to raise ValueError
        here, correctly, because Suwanlee did not define the case; GRAIN p.102
        does, so it now computes.
        """
        got = oracle(rows_of(FIXTURE_A))
        self.assertEqual({1, 2, 3, 4}, set(got))


class TestOracleDummyIdentityDiscrimination(unittest.TestCase):
    """
    Fixture C is only meaningful if the wrong rule would give a different answer.
    Measured here rather than asserted from theory, so the discriminating power
    of the exact-zero assertion is on the record.

    Merging the two unknown dams into one shared dummy makes animal 5 autozygous
    with probability 1/8 per locus, and animal 7 then inherits a flagged allele:

        P(5 autozygous) = P(3 passes its D allele) * P(4 passes its D allele)
                          * P(both are the same copy)
                        = 1/2 * 1/2 * 1/2 = 1/8
        E[f_a(7)]       = (1/8) / 2 = 0.0625

    So the correct rule and the merge bug are separated by 0.0625, against an
    assertion with zero tolerance.
    """

    def test_a_merged_dummy_would_be_visible(self):
        rows = rows_of(FIXTURE_C)
        correct = oracle(rows)
        self.assertEqual(0.0, correct[7])

        # The same pedigree with the two unknown dams named as ONE real animal:
        # the pedigree-level equivalent of merging the dummies.
        merged = oracle(read_pedigree(write_pedigree(
            "1 0 0\n2 0 0\n8 0 0\n3 1 8\n4 2 8\n5 3 4\n6 0 0\n7 5 6\n")))
        self.assertGreater(merged[7], 0.03,
                           "the merge bug must be visible, else fixture C "
                           "proves nothing")
        self.assertLessEqual(abs(merged[7] - 0.0625), 0.01)


class TestOraclePhantomEquivalence(unittest.TestCase):
    """
    Fixture D, oracle side. GRAIN's rule says the unknown side is an ORDINARY
    founder, so a pedigree with an implicit dummy and the same pedigree with that
    founder written out explicitly must be the same experiment.

    Exact rather than statistical: minting a dummy consumes no RNG, and each
    non-founder consumes exactly two draws in both pedigrees, so the draw
    sequence is identical. Label VALUES differ, but results depend only on the
    equality pattern of labels, which is isomorphic.
    """

    PAIRS = (
        ("A", FIXTURE_A, "1 0 0\n2 0 0\n9 0 0\n3 1 2\n4 3 9\n", (1, 2, 3, 4)),
        ("E", FIXTURE_E,
         "1 0 0\n2 0 0\n9 0 0\n3 1 9\n4 3 2\n5 3 2\n6 4 5\n7 0 0\n8 6 7\n",
         (1, 2, 3, 4, 5, 6, 7, 8)),
    )

    def test_implicit_and_explicit_completion_agree_exactly(self):
        for name, implicit, explicit, originals in self.PAIRS:
            with self.subTest(fixture=name):
                for seed in ORACLE_SEEDS:
                    a = gene_drop(rows_of(implicit), loci=2000, seed=seed)
                    b = gene_drop(rows_of(explicit), loci=2000, seed=seed)
                    for animal in originals:
                        self.assertEqual(
                            a[animal], b[animal],
                            "animal %s diverged at seed %s: implicit %r vs "
                            "explicit %r" % (animal, seed, a[animal], b[animal]))

    def test_mrode_implicit_and_explicit_completion_agree_exactly(self):
        implicit = read_pedigree(corpus("mrode.ped"))
        explicit = read_pedigree(write_pedigree(
            "1 0 0\n2 0 0\n7 0 0\n3 1 2\n4 1 7\n5 4 3\n6 5 2\n"))
        for seed in ORACLE_SEEDS:
            a = gene_drop(implicit, loci=2000, seed=seed)
            b = gene_drop(explicit, loci=2000, seed=seed)
            for animal in range(1, 7):
                self.assertEqual(a[animal], b[animal],
                                 "mrode animal %s diverged at seed %s"
                                 % (animal, seed))


class TestOracleIndependence(unittest.TestCase):
    """
    The oracle is only evidence if it is genuinely independent. Asserted on the
    source text, not on convention, because the failure mode -- someone factoring
    out a "shared helper" to remove duplication -- looks like an improvement.
    """

    ORACLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "oracles", "oracle_suwanlee.py")
    PRODUCTION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "PyPedal", "pyp_metrics.py")

    @staticmethod
    def _import_lines(path):
        """
        The IMPORT statements only. Prose in a docstring may name the other side
        -- production's comments cite the oracle by path, and should -- so a
        substring search over the whole file would fail on documentation rather
        than on a dependency.
        """
        with open(path, encoding="utf-8") as handle:
            return [line.strip() for line in handle
                    if line.strip().startswith(("import ", "from "))]

    def test_the_oracle_does_not_import_pypedal(self):
        for line in self._import_lines(self.ORACLE):
            self.assertNotIn("PyPedal", line, "oracle imports production: %r" % line)

    def test_production_does_not_import_the_oracle(self):
        for line in self._import_lines(self.PRODUCTION):
            self.assertNotIn("oracle", line, "production imports an oracle: %r" % line)
            self.assertNotIn("difftest", line,
                             "production imports from difftest: %r" % line)

    def test_neither_side_shares_a_dummy_founder_helper(self):
        """
        Both implement the same published rule; neither may call the same code to
        do it. If a helper were shared, agreement between them would measure
        nothing.
        """
        with open(self.ORACLE, encoding="utf-8") as handle:
            oracle_source = handle.read()
        with open(self.PRODUCTION, encoding="utf-8") as handle:
            production_source = handle.read()
        for shared in ("def _dummy_founder", "def _mint_dummy",
                       "def _half_founder_state"):
            self.assertNotIn(shared, oracle_source)
            self.assertNotIn(shared, production_source)


# ---------------------------------------------------------------------------
# Production. Every expectation below is the adjudicated behaviour, asserted
# BEFORE it exists. The markers come off in the commit that implements it.
# ---------------------------------------------------------------------------

class TestProductionExactZeroCases(unittest.TestCase):
    """
    Fixtures A, B, C and F -- the cases where the GRAIN rule gives a value that
    is exactly zero, at every locus of every replicate, for every seed.

    Derivation, common to all three: a dummy founder is a founder, so it carries
    no flag; and a fresh dummy label can never equal any other label, so no
    autozygosity arises through it. An animal whose parents are a founder (or an
    unflagged descendant) and a dummy therefore has nothing arriving flagged.

    These are asserted with NO tolerance. Zero is not a rounding of anything.
    """

    def test_a_known_sire_unknown_dam(self):
        ped = load_fixture(FIXTURE_A)
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(0.0, got[4])

    def test_b_unknown_sire_known_dam(self):
        ped = load_fixture(FIXTURE_B)
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(0.0, got[4])

    def test_f_the_two_orientations_agree_everywhere(self):
        """
        F, the symmetry fixture. A and B differ only in which parental slot holds
        the sentinel. Sire and dam are symmetric under Mendelian segregation, so
        swapping them cannot change the science -- and because both pedigrees
        have the same shape, the same seed must give the same numbers animal for
        animal, not merely the same distribution.
        """
        a = production(load_fixture(FIXTURE_A), EXACT_ROUNDS, EXACT_LOCI,
                       EXACT_SEED)
        b = production(load_fixture(FIXTURE_B), EXACT_ROUNDS, EXACT_LOCI,
                       EXACT_SEED)
        self.assertEqual(a, b)

    def test_c_independent_half_founders_share_no_identity(self):
        """
        The identity test. Animals 3 and 4 are half-founders with DIFFERENT
        unknown dams, so animal 5 draws from disjoint origin sets and can never
        be autozygous; nothing is ever flagged and animal 7 is exactly zero.

        Merging the two unknown sides -- the obvious wrong implementation, since
        both record the sentinel 0 -- would make animal 7 approximately 0.0625.
        See TestOracleDummyIdentityDiscrimination, where that is measured.
        """
        ped = load_fixture(FIXTURE_C)
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(0.0, got[7])
        for animal in (1, 2, 3, 4, 5, 6):
            self.assertEqual(0.0, got[animal])


class TestProductionDerivedExpectations(unittest.TestCase):
    """
    Fixture E -- the case where a flag ARISES downstream of a half-founder and
    must reach a descendant.

    Animals 4 and 5 are full sibs of the non-inbred, mutually unrelated pair
    (3, 2), so animal 6 is the offspring of a full-sib mating:

        P(6 autozygous at a locus) = 1/4     exactly
        f_a(6)                     = 0.0     exactly -- nothing ARRIVED flagged
        E[f_a(8)]                  = (1/4)/2 = 0.125

    The middle line is the convention that defines the estimator: a flag must not
    count for the animal that creates it, and must count for its descendants.

    The dummy is load-bearing here. Of that 1/4, the dummy allele accounts for
    1/16. Under a rival reading in which the unknown side contributes a
    never-matching null, P(6 autozygous) = 3/16 and E[f_a(8)] = 0.09375 -- about
    13 half-widths away, so this fixture discriminates the rules rather than
    merely observing one.
    """

    def test_the_creating_animal_is_exactly_zero(self):
        got = production(load_fixture(FIXTURE_E))
        self.assertEqual(0.0, got[6])

    def test_the_descendant_matches_the_derived_expectation(self):
        got = production(load_fixture(FIXTURE_E))
        width = half_width(0.125, 0.25)
        self.assertAlmostEqual(0.0024206, width, places=7)
        self.assertLessEqual(
            abs(got[8] - 0.125), width,
            "f_a(8) = %.6f is outside 0.125 +/- %.7f" % (got[8], width))

    def test_the_null_dummy_reading_is_excluded(self):
        """
        Stated as its own assertion so the discrimination is executable, not just
        described: 0.09375 must be outside the acceptance interval.
        """
        got = production(load_fixture(FIXTURE_E))
        width = half_width(0.125, 0.25)
        self.assertGreater(abs(0.09375 - 0.125), width)
        self.assertGreater(
            abs(got[8] - 0.09375), width,
            "f_a(8) = %.6f is compatible with the null-dummy reading" % got[8])


class TestProductionMixedPedigree(unittest.TestCase):
    """
    G, first half -- founders, full-parent descendants and half-founders in ONE
    pedigree, so the three domain cases are exercised together rather than in
    isolation. Animal 4 is a half-founder with a known sire, animal 5 with a
    known dam, and animal 3 has both parents.
    """

    def test_the_mixed_pedigree_computes_for_every_animal(self):
        ped = load_fixture(FIXTURE_MIXED)
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(8, len(got))
        for animal, value in got.items():
            self.assertGreaterEqual(value, 0.0, "animal %s" % animal)
            self.assertLessEqual(value, 1.0, "animal %s" % animal)


class TestProductionOnTheFormerRefusalFixtures(unittest.TestCase):
    """
    G, second half -- the two corpus pedigrees that carry half-founders, and were
    the EXPECTED_REFUSAL fixtures.

    mrode.ped, renumbered, is:

        1 0 0 | 2 0 0 | 3 2 1 | 4 2 0  <- half-founder | 5 4 3 | 6 5 1

    (originalID 1 and 2 swap places under renumbering; the derivation below is in
    ORIGINAL coordinates, where the half-founder is animal 4 with known sire 1.)

        f_a(1..5) = 0.0 exactly. Animal 4's sire is a founder and its dam is a
        dummy, so nothing arrives flagged, and a founder label can never equal a
        dummy label so animal 4 is never autozygous.

        P(5 autozygous) = 1/2 * 1/2 * 1/2 = 1/8 -- animals 4 and 3 must each pass
        their allele descended from the shared founder, and it must be the same
        copy. This equals mrode.ped's published inbreeding coefficient for animal
        5, F = 0.125, which is an independent check on the derivation.

        E[f_a(6)] = (1/8) / 2 = 0.0625.

    ON THE VALUE 0.0 FOR ANIMAL 4. This is derived HERE, afresh, from the
    source-explicit GRAIN rule. An early revision of
    the algorithm notes asserted 0.0 on different
    reasoning and that claim was WITHDRAWN as presuming phantom semantics no
    approved source then established. The withdrawn hypothesis is not evidence
    and is not used as one. That the two coincide is a result, not a premise.
    """

    def test_mrode_is_no_longer_refused(self):
        ped = load_corpus("mrode.ped", "asd")
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(6, len(got))

    def test_mrode_animals_one_to_five_are_exactly_zero(self):
        ped = load_corpus("mrode.ped", "asd")
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        for animal in (1, 2, 3, 4, 5):
            self.assertEqual(0.0, got[animal], "original animal %s" % animal)

    def test_mrode_animal_six_matches_the_derived_expectation(self):
        ped = load_corpus("mrode.ped", "asd")
        got = production(ped)
        width = half_width(0.0625, 0.125)
        self.assertAlmostEqual(0.001848775, width, places=9)
        self.assertLessEqual(
            abs(got[6] - 0.0625), width,
            "f_a(6) = %.6f is outside 0.0625 +/- %.7f" % (got[6], width))

    def test_hartlandclark_is_no_longer_refused(self):
        ped = load_corpus("hartlandclark.ped", "asdb")
        got = production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(15, len(got))
        for animal, value in got.items():
            self.assertGreaterEqual(value, 0.0, "animal %s" % animal)
            self.assertLessEqual(value, 1.0, "animal %s" % animal)


class TestProductionPhantomEquivalence(unittest.TestCase):
    """
    Fixture D, production side -- the sharpest statement of the rule. A pedigree
    with an implicit dummy must be the same experiment as the same pedigree with
    that founder written out explicitly.

    Attempted as EXACT same-seed equality first. Minting a dummy consumes no RNG
    and each non-founder consumes exactly two draws in both pedigrees, so the
    draw sequences coincide provided reorder/renumber preserves the relative
    order of the non-founders.

    If it does not, the plan's predeclared fallback applies and this becomes a
    distributional comparison. The fixture is NOT to be reshaped to force the
    exact form: that would be fitting the evidence to the conclusion.
    """

    def test_mrode_implicit_and_explicit_completion_agree(self):
        implicit = production(load_corpus("mrode.ped", "asd"),
                              EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        explicit = production(
            load_fixture("1 0 0\n2 0 0\n7 0 0\n3 1 2\n4 1 7\n5 4 3\n6 5 2\n"),
            EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        for animal in range(1, 7):
            self.assertEqual(
                implicit[animal], explicit[animal],
                "animal %s: implicit %r vs explicit %r"
                % (animal, implicit[animal], explicit[animal]))


class TestProductionStateAndReproducibility(unittest.TestCase):
    """
    H and I, on half-founder input specifically. The existing suite already pins
    both on half-founder-free pedigrees; the dummy-founder path is new code and
    must not be the place where a caller's pedigree starts being mutated or the
    global RNG starts moving.
    """

    def test_the_same_seed_reproduces_exactly(self):
        first = production(load_fixture(FIXTURE_E), EXACT_ROUNDS, EXACT_LOCI,
                           EXACT_SEED)
        second = production(load_fixture(FIXTURE_E), EXACT_ROUNDS, EXACT_LOCI,
                            EXACT_SEED)
        self.assertEqual(first, second)

    def test_a_different_seed_gives_a_different_answer(self):
        """
        Only asserted where variation is scientifically possible. Fixtures A, B
        and C are exactly zero at every seed by construction, so asserting
        variation there would be asserting a defect.
        """
        first = production(load_fixture(FIXTURE_E), EXACT_ROUNDS, EXACT_LOCI, 11)
        second = production(load_fixture(FIXTURE_E), EXACT_ROUNDS, EXACT_LOCI, 12)
        self.assertNotEqual(first[8], second[8])

    def test_the_process_global_numpy_stream_is_not_disturbed(self):
        import numpy as np
        np.random.seed(1234)
        before = np.random.get_state()
        production(load_fixture(FIXTURE_E), EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        after = np.random.get_state()
        self.assertEqual(before[0], after[0])
        self.assertTrue((before[1] == after[1]).all())
        self.assertEqual(before[2:], after[2:])

    def test_the_callers_animals_and_options_are_untouched(self):
        import copy
        ped = load_fixture(FIXTURE_E)
        animals = copy.deepcopy([vars(a).copy() for a in ped.pedigree])
        options = copy.deepcopy(ped.kw)
        production(ped, EXACT_ROUNDS, EXACT_LOCI, EXACT_SEED)
        self.assertEqual(animals, [vars(a).copy() for a in ped.pedigree])
        self.assertEqual(options, ped.kw)


class TestProductionAgreesWithTheIndependentOracle(unittest.TestCase):
    """
    Corroboration, deliberately last and deliberately loose. Two implementations
    of one published rule agreeing is weaker evidence than either satisfying an
    analytic fixture, and this must never become the primary anchor.
    """

    def test_production_and_oracle_agree_on_mrode(self):
        got = production(load_corpus("mrode.ped", "asd"))
        want = oracle(read_pedigree(corpus("mrode.ped")))
        for animal in range(1, 7):
            self.assertLessEqual(
                abs(got[animal] - want[animal]), 0.01,
                "animal %s: production %.6f vs oracle %.6f"
                % (animal, got[animal], want[animal]))

    def test_production_and_oracle_agree_on_hartlandclark(self):
        """
        The heavier half-founder fixture: 15 animals, FOUR half-founders (8, 9,
        10 and 11), so four independent dummy founders are minted per replicate
        and several descend into the same lines. mrode exercises one.
        """
        got = production(load_corpus("hartlandclark.ped", "asdb"))
        want = oracle(read_pedigree(corpus("hartlandclark.ped")))
        worst, where = 0.0, None
        for animal in sorted(want):
            deviation = abs(got[animal] - want[animal])
            if deviation > worst:
                worst, where = deviation, animal
        self.assertLessEqual(
            worst, 0.01,
            "worst production-vs-oracle deviation %.5f at animal %s"
            % (worst, where))

    def test_the_half_founders_of_hartlandclark_are_not_all_zero(self):
        """
        Guards the test above against being vacuous. If every value were zero,
        agreement would prove only that both sides return zeros.
        """
        got = production(load_corpus("hartlandclark.ped", "asdb"))
        self.assertGreater(sum(1 for v in got.values() if v > 0.0), 0)


if __name__ == "__main__":
    unittest.main()
