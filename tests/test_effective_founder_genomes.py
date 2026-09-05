"""
 -- the contract of ``pyp_metrics.effective_founder_genomes``.

TWO KINDS OF TEST LIVE HERE, AND THE DIFFERENCE MATTERS
-------------------------------------------------------
**Permanent tests** state facts that are true now and must stay true after the
repair: the routine estimates Boichard eq. 2, it reproduces the published
Table II value on a clean pedigree, and its results obey the derived bounds.
They are ordinary tests and are never removed.

**Strict-xfail tests** state the CORRECT contract for behaviour the adjudicated
repair is expected to change. They fail today, deliberately, and they are marked
``xfail(strict=True)`` so that the day production starts satisfying them the
suite says so loudly instead of quietly passing.

    When the repair lands: REMOVE THE MARKERS, KEEP THE TESTS.

Deleting them would delete the only executable statement of what the repair was
for. The evidence of *current* behaviour lives in
``the independent oracle`` and in the Finding-31
document, not here -- so nothing is lost by these tests turning green.

This research branch therefore carries strict xfails on purpose. The release
implementation must return the suite to zero xfails.

Tier 3 (a ``reference=`` API) is DEFERRED by operator decision and is
deliberately not tested here, not even as an xfail.
"""
import os
import random
import unittest

import pytest
from _pedhelpers import owned_temp_dir, chdir_tmp, corpus, load_corpus, load_corpus_from_path
from oracles import boichard_read, exact_ng, ng_distribution

from PyPedal import pyp_metrics
from PyPedal.pyp_errors import (
    PyPedalError, PyPedalPedigreeStructureError, PyPedalUsageError)



def _xfail(defect, what):
    return pytest.mark.xfail(
        strict=True,
        reason="effective founder genomes: %s -- %s" % (defect, what),
    )


# ---------------------------------------------------------------------------
# PREDECLARED STOCHASTIC PROTOCOL
#
# Fixed here, from the ORACLE's exact variance, before any acceptance run --
# never tuned to an observed production result and never seed-shopped.
#
#   fixture     corpus/boichard_fig2.ped, R = the g == '2' cohort
#   target      E[N_g] = 2.522425 exactly (factorised enumeration, 2**28 space)
#   Var[N_g]    0.051623 per replicate  ->  sd 0.227207
#   estimator   the EWMA has effective sample size 3 regardless of `rounds`,
#               so a single call has sd = sqrt(Var/3) = 0.131178
#   protocol    ROUNDS rounds, REPLICATES replicates, accept if the replicate
#               mean lies within 5 sigma of the target:
#                   5 * sqrt(Var / 3 / REPLICATES) = 0.1037
#   The band is deliberately computed for the CURRENT, worse estimator. The
#   repaired arithmetic mean has smaller variance, so the same band still holds.
# ---------------------------------------------------------------------------
PROTOCOL = {
    "fixture": "boichard_fig2.ped",
    "rounds": 300,
    "replicates": 40,
    "oracle_target": 2.522425,
    "oracle_variance": 0.051623,
    "ewma_effective_sample_size": 3.0,
    "band": 0.1037,
}


def _fig2_reference():
    rows, gens = boichard_read(corpus("boichard_fig2.ped"), gen_col=3)
    return rows, [a for a, _s, _d in rows if gens[a] == "2"]


def _replicates(ped_factory, rounds, n):
    out = []
    for _ in range(n):
        out.append(float(pyp_metrics.effective_founder_genomes(
            ped_factory(), rounds=rounds, quiet=True)))
    return out


def _mean_sd(xs):
    m = sum(xs) / float(len(xs))
    return m, (sum((x - m) ** 2 for x in xs) / float(len(xs))) ** 0.5


# ===========================================================================
# PERMANENT -- true now, and required to stay true after the repair
# ===========================================================================

class TestItEstimatesBoichardEquation2(unittest.TestCase):
    """
    The routine's *mathematics* is not the defect. On a clean pedigree with no
    half-founder, no founder inside R and single-digit generation labels, it
    reproduces Boichard's published Table II value.

    This is what makes the finding a repair rather than a rewrite, and it is the
    control that must not regress while the defects around it are fixed.
    """

    def test_reproduces_published_table_ii_total(self):
        rows, reference = _fig2_reference()
        target = exact_ng(rows, reference)["n_g"]
        self.assertAlmostEqual(PROTOCOL["oracle_target"], target, delta=1e-6)

        with chdir_tmp():
            values = _replicates(lambda: load_corpus("boichard_fig2.ped"),
                                 PROTOCOL["rounds"], PROTOCOL["replicates"])
        mean, _sd = _mean_sd(values)
        self.assertAlmostEqual(target, mean, delta=PROTOCOL["band"])
        # ... and the paper's own printed value, to its printed precision
        self.assertAlmostEqual(2.5, mean, delta=5e-2 + PROTOCOL["band"])

    def test_the_reference_population_is_the_papers_own(self):
        """Anti-vacuity: the g column really does select Boichard's ellipse."""
        _rows, reference = _fig2_reference()
        self.assertEqual(list(range(7, 15)) + [19, 20], sorted(reference))


class TestDerivedBoundsHoldOnProductionOutput(unittest.TestCase):
    """
    ``N_g >= 0.5``, not ``>= 1``: SUM f_k = 1 forces SUM f_k^2 <= 1. Equality is
    fixation of a single founder gene, a legitimate gene-drop outcome.

    A production value below 1 is therefore NOT evidence of a defect, and this
    test exists partly to stop a future postcondition from being written too
    strong.
    """

    def test_clean_pedigree_results_are_at_least_one_half(self):
        with chdir_tmp():
            for name in ("boichard_fig2.ped", "new_lacy.ped"):
                for rounds in (1, 25):
                    with self.subTest(pedigree=name, rounds=rounds):
                        v = pyp_metrics.effective_founder_genomes(
                            load_corpus(name), rounds=rounds, quiet=True)
                        self.assertGreaterEqual(float(v), 0.5)

    def test_values_below_one_are_attainable_and_legal(self):
        """
        A selfing pedigree can fix a founder gene in R, giving exactly 0.5. The
        oracle enumerates it; nothing may reject it as impossible.
        """
        attainable = [float(1 / (2 * s)) for s in
                      ng_distribution([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 3)],
                                      [4])]
        self.assertEqual(0.5, min(attainable))


class TestRoundsValidation(unittest.TestCase):
    """
    Each round is one Monte Carlo replicate, so `rounds < 1` is not a
    calculation with a defensible answer. It used to be coerced to 1 and logged;
    it is now a typed caller-contract refusal.
    """

    def test_non_positive_rounds_is_refused(self):
        with chdir_tmp():
            for rounds in (0, -5):
                with self.subTest(rounds=rounds):
                    with self.assertRaises(PyPedalUsageError):
                        pyp_metrics.effective_founder_genomes(
                            load_corpus("new_lacy.ped"), rounds=rounds, quiet=True)


# ===========================================================================
# STRICT XFAIL -- the adjudicated contract the repair must deliver.
# Remove the markers when it lands. Do not delete the tests.
# ===========================================================================

#: Number of fresh founders declared in generation '12'. Their N_g is EXACT --
#: see ``TestD1ReferencePopulationSelection``.
D1_FRESH_FOUNDERS = 6


def _deep_gen_text():
    """
    Twelve declared generations, so that '9' sorts above '12' as a string.

    Generations 1-11 are a narrow chain of full-sib pairs descending from two
    founders; drift there drives N_g down to about 0.63 by generation 9.
    Generation 12 is instead six FRESH unrelated founders, whose N_g is exactly
    6 with no Monte Carlo variance at all -- nothing segregates.

    That asymmetry is deliberate: it makes "which generation was selected" a
    deterministic question with two answers that cannot be confused, instead of
    a comparison of two noisy stochastic means.
    """
    rows, prev, nid = ["1 0 0 1", "2 0 0 1"], (1, 2), 3
    for g in range(2, 12):
        rows.append("%d %d %d %d" % (nid, prev[0], prev[1], g))
        rows.append("%d %d %d %d" % (nid + 1, prev[0], prev[1], g))
        prev, nid = (nid, nid + 1), nid + 2
    for _ in range(D1_FRESH_FOUNDERS):
        rows.append("%d 0 0 12" % nid)
        nid += 1
    return "\n".join(rows) + "\n"


def _write(tmpdir, name, text):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


class TestD1ReferencePopulationSelection(unittest.TestCase):
    """
    D1 -- the unrepaired residue of audit .

    ``NewAnimal.gen`` is a string, and the routine picks R with
    ``sorted(set(gens))[-1]``, which orders labels lexicographically. On a
    pedigree with ten or more declared generations that selects generation '9'.

    ``pyp_metrics._most_recent_generation`` already exists and fixes exactly this
    for the three Boichard routines; this one was missed.
    """

    def _fixture(self):
        tmp = owned_temp_dir(prefix="fg31_d1_")
        return _write(tmp, "deep.ped", _deep_gen_text())

    def test_the_fixture_really_does_sort_wrongly_as_strings(self):
        """Anti-vacuity: without ten or more generations there is no defect."""
        ped = load_corpus_from_path(self._fixture(), "asdg")
        gens = sorted({a.gen for a in ped.pedigree})
        self.assertEqual("9", gens[-1])
        self.assertEqual("12", max(gens, key=float))

    def test_selects_the_numerically_most_recent_generation(self):
        """
        Generation '12' is six fresh founders, so selecting it gives exactly
        6.0 with no variance. Generation '9' -- the lexicographic answer -- is
        two heavily drifted full sibs, worth about 0.63. Only one of those two
        numbers can come out, so this needs no stochastic tolerance.
        """
        got = pyp_metrics.effective_founder_genomes(
            load_corpus_from_path(self._fixture(), "asdg"), rounds=3, quiet=True)
        self.assertAlmostEqual(float(D1_FRESH_FOUNDERS), float(got), places=9)

    def test_non_numeric_generation_labels_are_refused_not_guessed(self):
        """
        ``_most_recent_generation`` refuses labels that are not ordered
        quantities rather than falling back to string order. Same contract as
        the three Boichard routines.
        """
        tmp = owned_temp_dir(prefix="fg31_d1b_")
        path = _write(tmp, "alpha.ped",
                      "1 0 0 early\n2 0 0 early\n3 1 2 late\n4 1 2 late\n")
        with self.assertRaises(PyPedalError):
            pyp_metrics.effective_founder_genomes(
                load_corpus_from_path(path, "asdg"), rounds=2, quiet=True)


class TestD2FoundersInsideTheReferencePopulation(unittest.TestCase):
    """
    D2 -- the ``continue`` that skips a founder's transmission also skips its
    frequency tally, so a founder belonging to R contributes no gene copies.

    R is the analyst's population under study (Boichard App. A/B step 1). Lacy's
    okapi and Goeldi's-monkey analyses are exactly populations containing
    wild-caught founders, so "founders are never in R" is not available as a
    defence.
    """

    #: founders 1, 2 in generation 1; 3, 4 their offspring in generation 2;
    #: founder 5 also declared in generation 2 -- a wild-caught import.
    TEXT = "1 0 0 1\n2 0 0 1\n3 1 2 2\n4 1 2 2\n5 0 0 2\n"

    def test_founder_gene_copies_are_counted(self):
        tmp = owned_temp_dir(prefix="fg31_d2_")
        path = _write(tmp, "founder_in_r.ped", self.TEXT)
        target = exact_ng([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 0, 0)],
                          [3, 4, 5])["n_g"]
        values = _replicates(lambda: load_corpus_from_path(path, "asdg"), 300, 40)
        mean, _sd = _mean_sd(values)
        self.assertAlmostEqual(target, mean, delta=0.15)


class TestD3HalfFounders(unittest.TestCase):
    """
    D3 -- 's original symptom, re-measured on this branch.

    A missing parent is the sentinel 0, and ``pedigree[0 - 1]`` is
    ``pedigree[-1]``. What that resolves to decides the outcome: an unrelated
    animal with no alleles yet trips the guard and the routine returns the
    sentinel ``0``; anything else is silently adopted as a parent.

    The correct contract is the phantom-founder rule -- Lacy p.113, Boichard
    App. A step 4, Baumung et al. 2015 p.102 -- which FG-10 shows is exactly
    equivalent to putting one unique gene in the unknown slot.
    """

    #: Two half-founders sharing the known parent, with R = both of them.
    #:
    #: A SINGLE half-founder with R = {itself} does not discriminate: whatever
    #: the wrapped lookup returns, the animal ends up with two distinct genes at
    #: frequency 1/2 and N_g is 1 either way. Two of them do discriminate,
    #: because the phantom rule gives each its OWN unique unknown-side gene
    #: whereas the wrapped lookup lets them share one. Oracle 1.667 against a
    #: measured 1.166.
    #: The ``g`` column is present so that the oracle's R and production's R
    #: are provably the same set of animals. Without it every animal shares the
    #: missing-generation default and R silently becomes the whole pedigree,
    #: which would compare two different quantities.
    CASES = {
        "known sire, unknown dam": (
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 2\n5 3 0 2\n",
            [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0), (5, 3, 0)], [4, 5]),
        "unknown sire, known dam": (
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 0 3 2\n5 0 3 2\n",
            [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 0, 3), (5, 0, 3)], [4, 5]),
    }

    def test_half_founder_agrees_with_the_source_backed_oracle(self):
        tmp = owned_temp_dir(prefix="fg31_d3_")
        # No subTest: a strict xfail must fail as ONE outcome, and a subTest
        # that passes while a sibling fails muddies which contract is unmet.
        measured = {}
        for name, (text, rows, reference) in sorted(self.CASES.items()):
            path = _write(tmp, name.replace(", ", "_").replace(" ", "_") + ".ped", text)
            target = exact_ng(rows, reference)["n_g"]
            values = _replicates(lambda: load_corpus_from_path(path, "asdg"), 200, 30)
            mean, _sd = _mean_sd(values)
            measured[name] = (target, mean)
        for name, (target, mean) in sorted(measured.items()):
            self.assertAlmostEqual(target, mean, delta=0.15,
                                   msg="%s: oracle %.4f, production %.4f"
                                       % (name, target, mean))

    def test_mrode_agrees_with_the_oracle_instead_of_returning_the_sentinel(self):
        """
        ``mrode.ped`` returned the sentinel 0 because its half-founder's dam
        lookup wrapped to animal 6, whose alleles were still ``['','']``. It has
        no generation column, so R is the whole pedigree.
        """
        target = exact_ng([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 0),
                           (5, 4, 3), (6, 5, 2)], [1, 2, 3, 4, 5, 6])["n_g"]
        self.assertAlmostEqual(1.9664, target, delta=1e-3)
        with chdir_tmp():
            values = _replicates(lambda: load_corpus("mrode.ped"), 200, 40)
        mean, _sd = _mean_sd(values)
        self.assertAlmostEqual(target, mean, delta=0.15)

    def test_the_empty_string_is_no_longer_a_founder_gene(self):
        """
        The wrapped lookup could resolve to the animal being processed, after
        which ``random.choice(['', 'X__2'])`` put the EMPTY STRING into eq. 2's
        denominator as though it were a founder gene. Nothing may produce a
        falsy gene label now.
        """
        tmp = owned_temp_dir(prefix="fg31_d3b_")
        path = _write(tmp, "hf_last.ped", "1 0 0\n2 0 0\n3 1 2\n4 3 0\n")
        ped = load_corpus_from_path(path, "asd")
        most_recent = pyp_metrics._most_recent_generation(
            [a.gen for a in ped.pedigree], "t")
        plan = pyp_metrics._build_gene_drop_plan(ped, most_recent, "t")
        labels = [g for pair in plan.founder_genes if pair for g in pair]
        labels += [g for pair in plan.slot_genes for g in pair if g]
        self.assertTrue(all(labels), labels)
        self.assertEqual(len(labels), len(set(labels)))


class TestD11FounderGeneDenominator(unittest.TestCase):
    """
    D11 -- the report said ``2 * count(founder == 'y')``, which counts only
    animals with BOTH parents unknown and so understates 2f wherever a
    half-founder exists.

    The slot optimisation materialises one gene per unknown slot, but the
    conceptual founder count is unchanged: each unknown parental slot is a
    diploid phantom founder contributing TWO founder genes, one of which simply
    has frequency zero in R. The optimisation must not silently redefine f.
    """

    def _plan(self, text, pedformat="asd"):
        tmp = owned_temp_dir(prefix="fg31_d11_")
        ped = load_corpus_from_path(_write(tmp, "d11.ped", text), pedformat)
        most_recent = pyp_metrics._most_recent_generation(
            [a.gen for a in ped.pedigree], "t")
        return pyp_metrics._build_gene_drop_plan(ped, most_recent, "t")

    def test_two_conceptual_genes_per_phantom_founder(self):
        plan = self._plan("1 0 0\n2 0 0\n3 1 2\n4 3 0\n5 3 0\n")
        self.assertEqual(2, plan.n_true_founders)
        self.assertEqual(2, plan.n_slots)
        # the old count would have been 2 * 2 = 4
        self.assertEqual(8, plan.n_founder_genes)

    def test_matches_the_oracle_2f_for_the_phantom_model(self):
        from oracles import ng_founder_genes
        rows = [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 0), (5, 4, 3), (6, 5, 2)]
        plan = self._plan("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
        self.assertEqual(len(ng_founder_genes(rows, "phantom")),
                         plan.n_founder_genes)

    def test_no_half_founders_leaves_the_count_unchanged(self):
        plan = self._plan("1 0 0\n2 0 0\n3 1 2\n4 1 2\n")
        self.assertEqual(0, plan.n_slots)
        self.assertEqual(4, plan.n_founder_genes)


class TestD12StructuralSafety(unittest.TestCase):
    """
    D12 -- ``pedigree[int(parentID) - 1]`` was the parent lookup. It is an
    accident that works only while position == animalID - 1, and it turns the
    missing-parent sentinel 0 into ``pedigree[-1]`` rather than an error.

    Parents now resolve through an explicit ID map, validated once before any
    replicate and before the report file is opened.  makes the normal
    loader produce a canonically ordered pedigree; a direct caller can still
    build something invalid, and these are the refusals it gets. 
    itself is not touched.
    """

    def _load_unchecked(self, text, **overrides):
        tmp = owned_temp_dir(prefix="fg31_d12_")
        options = dict(renumber=False, reorder=False, pedigree_is_renumbered=True)
        options.update(overrides)
        return load_corpus_from_path(_write(tmp, "d12.ped", text), "asd", **options)

    def test_offspring_before_parent_is_refused(self):
        ped = self._load_unchecked("1 0 0\n2 3 4\n3 0 0\n4 0 0\n")
        with self.assertRaises(PyPedalPedigreeStructureError):
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)

    def test_unresolvable_known_parent_is_refused_not_treated_as_missing(self):
        ped = self._load_unchecked("1 0 0\n2 0 0\n3 1 99\n")
        with self.assertRaises(PyPedalPedigreeStructureError):
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)

    def test_self_parent_is_refused(self):
        ped = self._load_unchecked("1 0 0\n2 0 0\n3 3 1\n")
        with self.assertRaises(PyPedalPedigreeStructureError):
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)

    def test_no_untyped_exception_escapes_these_paths(self):
        """None of the refusals may surface as IndexError/KeyError/ZeroDivision."""
        for text in ("1 0 0\n2 3 4\n3 0 0\n4 0 0\n",
                     "1 0 0\n2 0 0\n3 1 99\n",
                     "1 0 0\n2 0 0\n3 3 1\n"):
            ped = self._load_unchecked(text)
            try:
                pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)
            except PyPedalError:
                pass
            except Exception as exc:               # noqa: BLE001 - that is the point
                self.fail("untyped %s: %s" % (type(exc).__name__, exc))


class TestD4RoundsMustBuyPrecision(unittest.TestCase):
    """
    D4 -- ``summary_freqs['n_g'] = (_ng + n_g)/2`` is an exponentially weighted
    moving average, not a mean. Its weights are 2**-(R-r), so
    ``SUM w^2 -> 1/3`` and the effective sample size is 3 no matter how large
    ``rounds`` is. The documented parameter buys no precision past ~3 rounds.

    Per-replicate inversion is NOT the defect and must be kept: Boichard's
    published Family 2 value of 1.1 is E[1/(2S)] = 1.10625, whereas
    1/(2 E[S]) = 1.0 exactly.

    Replicate counts are predeclared in ``PROTOCOL``; the ratio threshold is
    0.5, which the arithmetic mean clears by a factor of ten at these rungs
    while the EWMA sits at 1.0.
    """

    LOW_ROUNDS, HIGH_ROUNDS, REPLICATES = 3, 300, 120

    def test_the_ewma_weights_are_what_the_algebra_says(self):
        """
        Permanent, and the reason the ratio test below is not merely empirical.

        The recurrence is ``m_1 = x_1`` and ``m_r = (x_r + m_(r-1)) / 2``. Push
        unit vectors through it to read off the weights: they are
        ``1/2, 1/4, ..., 1/2**(R-1), 1/2**(R-1)``, they sum to 1, and their
        squares sum to 1/3 in the limit -- an effective sample size of 3.
        """
        for rounds in (1, 2, 5, 30):
            with self.subTest(rounds=rounds):
                weights = []
                for r in range(rounds):
                    weights = ([w / 2.0 for w in weights] + [0.5]) if r else [1.0]
                self.assertAlmostEqual(1.0, sum(weights), places=9)
                self.assertEqual(rounds, len(weights))
        ess = 1.0 / sum(w * w for w in weights)
        self.assertAlmostEqual(3.0, ess, places=6)

    def test_more_rounds_reduces_spread(self):
        with chdir_tmp():
            low = _replicates(lambda: load_corpus("boichard_fig2.ped"),
                              self.LOW_ROUNDS, self.REPLICATES)
            high = _replicates(lambda: load_corpus("boichard_fig2.ped"),
                               self.HIGH_ROUNDS, self.REPLICATES)
        _m1, sd_low = _mean_sd(low)
        _m2, sd_high = _mean_sd(high)
        self.assertLess(sd_high, 0.5 * sd_low,
                        "sd(%d rounds)=%.4f vs sd(%d rounds)=%.4f"
                        % (self.HIGH_ROUNDS, sd_high, self.LOW_ROUNDS, sd_low))


class TestD5AllFounderReferencePopulationIsValid(unittest.TestCase):
    """
    D5, corrected during the repair phase.

    An all-founder R used to raise a bare ``ZeroDivisionError``, and the
    research phase initially filed that under "undefined calculation". It is
    not. The exception existed only because production skipped founders before
    tallying them (D2), leaving the frequency denominator at zero.

    With D2 repaired the calculation is perfectly well defined, and its value is
    exact: *f* unrelated diploid founders each contributing both gene copies
    give 2f founder genes at frequency 1/(2f), so

        SUM f_k^2 = 2f * (1/2f)^2 = 1/(2f)   and   N_g = f

    with **no** Monte Carlo variance -- every replicate gives the same answer,
    because no transmission happens at all.
    """

    def _all_founder(self, f):
        tmp = owned_temp_dir(prefix="fg31_d5_")
        text = "".join("%d 0 0\n" % i for i in range(1, f + 1))
        return _write(tmp, "allfounder%d.ped" % f, text)

    def test_all_founder_reference_population_gives_exactly_f(self):
        for f in (1, 3, 7):
            with self.subTest(founders=f):
                path = self._all_founder(f)
                rows = [(i, 0, 0) for i in range(1, f + 1)]
                self.assertEqual(float(f), exact_ng(rows, list(range(1, f + 1)))["n_g"])
                got = pyp_metrics.effective_founder_genomes(
                    load_corpus_from_path(path, "asd"), rounds=3, quiet=True)
                self.assertAlmostEqual(float(f), float(got), places=9)

    def test_it_is_deterministic_because_nothing_segregates(self):
        path = self._all_founder(5)
        for rounds in (1, 2, 10):
            with self.subTest(rounds=rounds):
                got = pyp_metrics.effective_founder_genomes(
                    load_corpus_from_path(path, "asd"), rounds=rounds, quiet=True)
                self.assertAlmostEqual(5.0, float(got), places=9)


class TestD6ReproducibilityAndGlobalRng(unittest.TestCase):
    """
    D6 -- the routine calls the argument-less ``random.seed()`` once per round,
    reseeding a PROCESS-GLOBAL generator from OS entropy. There is no ``seed``
    parameter, results are not reproducible, and any caller-side seeding is
    destroyed.

    Note for the differential record: PyPedal 2.0.4 has the same defect against
    a DIFFERENT generator -- see the Python-2 section of the Finding-31 document.
    """

    def test_a_seed_makes_the_result_reproducible(self):
        with chdir_tmp():
            a = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=20, seed=7, quiet=True)
            b = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=20, seed=7, quiet=True)
        self.assertEqual(a, b)

    def test_the_global_random_state_is_left_alone(self):
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            random.seed(20260822)
            before = random.getstate()
            pyp_metrics.effective_founder_genomes(ped, rounds=3, quiet=True)
            self.assertEqual(before, random.getstate())

    def test_the_numpy_global_random_state_is_left_alone(self):
        """
        D6a. PyPedal 2.0.4's ``from numpy import random`` shadowed its standard
        library import, so the legacy routine consumed and reset NUMPY's global
        generator while the Python 3 port used the standard library's. The two
        implementations polluted different global streams -- the one measured
        Python 2 / Python 3 divergence in this finding. The repair supersedes
        both: neither is touched.
        """
        import numpy
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            numpy.random.seed(20260822)
            before = numpy.random.get_state()
            pyp_metrics.effective_founder_genomes(ped, rounds=3, quiet=True)
            after = numpy.random.get_state()
        self.assertEqual(before[0], after[0])
        self.assertTrue((before[1] == after[1]).all())
        self.assertEqual(before[2:], after[2:])

    def test_an_unseeded_call_still_varies(self):
        """Anti-vacuity: reproducibility must come from the seed, not from the
        routine having quietly become deterministic."""
        with chdir_tmp():
            values = {float(pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=5, quiet=True))
                for _ in range(12)}
        self.assertGreater(len(values), 1)

    def test_different_seeds_give_different_draws(self):
        with chdir_tmp():
            a = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=5, seed=1, quiet=True)
            b = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=5, seed=2, quiet=True)
        self.assertNotEqual(a, b)

    def test_the_same_pedobj_can_be_reused_without_carrying_state(self):
        """
        The old routine wrote each replicate onto the caller's animals, so a
        second call started from the first call's leftovers. A seeded call must
        now give the same answer on a REUSED pedigree as on a fresh one.
        """
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            first = pyp_metrics.effective_founder_genomes(ped, rounds=8, seed=99,
                                                         quiet=True)
            second = pyp_metrics.effective_founder_genomes(ped, rounds=8, seed=99,
                                                          quiet=True)
            fresh = pyp_metrics.effective_founder_genomes(
                load_corpus("mrode.ped"), rounds=8, seed=99, quiet=True)
        self.assertEqual(first, second)
        self.assertEqual(first, fresh)


class TestD7CallerStateIsNotMutated(unittest.TestCase):
    """
    D7 -- the routine writes simulation scratch into ``animal.alleles`` on the
    caller's own ``NewAnimal`` records, so a second call starts from the first
    call's leftovers and any downstream reader of ``alleles`` sees them.
    """

    #: Every NewAnimal field the routine could plausibly disturb.
    FIELDS = ("alleles", "gen", "igen", "founder", "sireID", "damID", "animalID",
              "originalID", "renumberedID", "sex", "ancestor", "sons",
              "daus", "unks", "pedcomp", "fa", "paddedID")

    def _snapshot(self, ped):
        return ([[repr(getattr(a, f, None)) for f in self.FIELDS] for a in ped.pedigree],
                [str(a.originalID) for a in ped.pedigree],
                dict(ped.kw))

    def test_pedigree_alleles_are_unchanged_by_a_call(self):
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            before = [list(a.alleles) for a in ped.pedigree]
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)
            self.assertEqual(before, [list(a.alleles) for a in ped.pedigree])

    def test_nothing_else_on_the_pedigree_moves_either(self):
        """Fields, pedigree ORDER and kw, on a pedigree that has half-founders."""
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = self._snapshot(ped)
            pyp_metrics.effective_founder_genomes(ped, rounds=4, seed=5, quiet=True)
            self.assertEqual(before, self._snapshot(ped))

    def test_half_founder_placeholder_alleles_survive_the_call(self):
        """
        ``NewAnimal.__init__`` seeds a half-founder with ``['', 'X__2']``. The
        old routine overwrote both slots on the caller's record; the repaired
        one reads neither and writes nothing.
        """
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            half = [a for a in ped.pedigree if "" in a.alleles]
            self.assertTrue(half, "fixture no longer contains a half-founder")
            before = {a.animalID: list(a.alleles) for a in half}
            pyp_metrics.effective_founder_genomes(ped, rounds=3, quiet=True)
            self.assertEqual(before, {a.animalID: list(a.alleles) for a in half})


class TestD8ReportSideEffect(unittest.TestCase):
    """
    D8 -- ``<filetag>_gene_drop.out`` is written into the working directory on
    every call, with no way to turn it off. ``pyp_nrm.inbreeding`` already
    establishes the ``output=`` convention this routine should follow.
    """

    @staticmethod
    def _report(ped):
        """
        The report path is derived from ``kw['filetag']``, which is the pedfile
        path -- NOT the working directory. Asserting on ``os.listdir('.')``
        looks right and is vacuous, because the file never lands there.
        """
        return "%s_gene_drop.out" % ped.kw["filetag"]

    def test_output_can_be_suppressed(self):
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True,
                                                  output=False)
            self.assertFalse(os.path.exists(self._report(ped)))

    def test_output_defaults_to_on(self):
        """The historical side effect stays on unless a caller opts out."""
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)
            self.assertTrue(os.path.exists(self._report(ped)))

    def test_the_report_names_the_corrected_founder_gene_total(self):
        """
        D11 reaches the file, not just the plan. ``mrode.ped`` has two ordinary
        founders and one unknown parental slot, so 2 * (2 + 1) = 6 conceptual
        founder genes -- where the old count of 2 * 2 = 4 ignored the phantom.
        The phantom's transmitted gene is labelled and legible in the report;
        the empty string that used to appear as ``Allele : `` is gone.
        """
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            pyp_metrics.effective_founder_genomes(ped, rounds=2, seed=3, quiet=True)
            with open(self._report(ped)) as fh:
                text = fh.read()
        self.assertIn("Number of distinct founder alleles: 6", text)
        self.assertIn("__d:", text)
        self.assertNotIn("Allele : ", text)

    def test_suppressing_output_does_not_change_the_science(self):
        with chdir_tmp():
            on = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=6, seed=77, quiet=True)
            off = pyp_metrics.effective_founder_genomes(
                load_corpus("boichard_fig2.ped"), rounds=6, seed=77, quiet=True,
                output=False)
        self.assertEqual(on, off)

    def test_a_refused_call_writes_no_report_at_all(self):
        """
        Preconditions are validated before the file is opened, so a refusal
        cannot leave a partial report that looks like a successful run.
        """
        with chdir_tmp() as tmp:
            ped = load_corpus("new_lacy.ped")
            with self.assertRaises(PyPedalError):
                pyp_metrics.effective_founder_genomes(ped, rounds=2,
                                                      chrometype="sex", quiet=True)
            self.assertFalse(os.path.exists(self._report(ped)))
            self.assertEqual([], [f for f in os.listdir(tmp)
                                  if f.endswith("_gene_drop.out")])


class TestD9SexChromosomeIsNotImplemented(unittest.TestCase):
    """
    D9 -- ``chrometype`` and ``heterogametic`` are validated and then never
    read. ``chrometype='sex'`` silently returns an autosomal N_g, which is a
    wrong answer to the question the caller asked.

    The repair is to refuse, not to invent sex-chromosome gene-drop
    mathematics.
    """

    def test_the_parameters_are_never_read_in_the_autosomal_path(self):
        """
        Permanent, mechanical. ``chrometype='sex'`` now refuses outright, and
        for the only implemented mode -- autosomal -- neither parameter has any
        read site at all. ``heterogametic`` is genuinely irrelevant there, which
        is what the docstring now says instead of implying it is honoured.
        """
        source = open(pyp_metrics.__file__.replace(".pyc", ".py")).read()
        body = source[source.index("def effective_founder_genomes"):]
        body = body[:body.index("\ndef ", 1)]
        anchor = "heterogametic = 'm'"
        tail = body[body.rindex(anchor) + len(anchor):]
        self.assertIn("rng.choice", tail, "anti-vacuity: the tail is the "
                      "calculation proper, not an empty slice")
        self.assertEqual(0, tail.count("chrometype"))
        self.assertEqual(0, tail.count("heterogametic"))

    def test_sex_chromosome_refuses_loudly(self):
        with chdir_tmp():
            with self.assertRaises(PyPedalError):
                pyp_metrics.effective_founder_genomes(
                    load_corpus("new_lacy.ped"), rounds=2, chrometype="sex",
                    quiet=True)


if __name__ == "__main__":
    pytest.main([__file__])
