"""
The Suwanlee gene-drop oracle, the denominator adjudication, and the size of the
gap to production.

**No production code is changed by anything here.** b's repair is deferred
to its own implementation boundary; this commit's job is to build the
instrument, run the predeclared experiment, and MEASURE the gap so the deferral
is a quantified decision rather than an assertion.
"""
import os
import sys
import unittest
import warnings

from _pedhelpers import corpus, load_corpus

from PyPedal import pyp_metrics

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "oracles"))
from oracle_suwanlee import (                 # noqa: E402
    CANDIDATES,
    LOCI,
    SEEDS,
    adjudicate,
    gene_drop,
    measure,
    read_pedigree,
    tolerance,
)

# Suwanlee Figure 1, article p.491, the f_a_g column.
PUBLISHED_GENE_DROP = {11: 0.000, 12: 0.000, 13: 0.125, 14: 0.125,
                       15: 0.375, 16: 0.375,
                       17: 0.585, 18: 0.585, 19: 0.585, 20: 0.585,
                       21: 0.728, 22: 0.728, 23: 0.822}


def _rows():
    return read_pedigree(corpus("suwanlee_fig1.ped"))


class TestTheAdjudication(unittest.TestCase):
    """
    The predeclared D1-vs-D2 procedure, run in full. The criterion lives in the
    oracle module and is not adjusted here -- if these constants ever move, the
    experiment has changed and the verdict has to be re-earned.
    """

    @classmethod
    def setUpClass(cls):
        cls.report = adjudicate(_rows(), PUBLISHED_GENE_DROP)

    def test_the_criterion_is_the_predeclared_one(self):
        self.assertEqual(200_000, LOCI)
        self.assertEqual((1, 2, 3, 4, 5), SEEDS)
        self.assertEqual(1.96, self.report["z"])
        self.assertEqual(5e-4, self.report["rounding_halfwidth"])
        self.assertEqual(1_000_000, self.report["published_replicates"])

    def test_the_six_published_values_expand_to_thirteen_checks(self):
        """
        Figure 1 prints ONE row of values per generation, labelling one to four
        animals each. Six distinct published values, thirteen animal-level
        checks. The seven extra checks are the same six values at co-labelled
        siblings -- independent Monte-Carlo draws of one expectation, so
        corroborating rather than additional source evidence.

        An earlier report said "13 published rows". It was a counting error;
        this pins the distinction so it cannot recur.
        """
        self.assertEqual(13, len(PUBLISHED_GENE_DROP))
        self.assertEqual(6, len(set(PUBLISHED_GENE_DROP.values())))

    def test_d1_satisfies_every_distinct_published_value(self):
        """
        The criterion is stated over the six values, so it is checked over the
        six -- not merely over the thirteen expanded targets, where a value
        could in principle pass on one sibling and fail on another.
        """
        groups = {}
        for animal, value in PUBLISHED_GENE_DROP.items():
            groups.setdefault(value, []).append(animal)
        rows = {r["animal"]: r for r in self.report["candidates"]["D1"]["rows"]}
        for value, animals in groups.items():
            with self.subTest(published=value):
                self.assertTrue(all(rows[a]["inside"] for a in animals))

    def test_d2_fails_five_of_the_six_distinct_values(self):
        groups = {}
        for animal, value in PUBLISHED_GENE_DROP.items():
            groups.setdefault(value, []).append(animal)
        rows = {r["animal"]: r for r in self.report["candidates"]["D2"]["rows"]}
        failed = [v for v, animals in groups.items()
                  if not all(rows[a]["inside"] for a in animals)]
        self.assertEqual(5, len(failed))
        # Only the trivial zero row fits, which is no evidence at all.
        self.assertEqual([0.0], sorted(set(groups) - set(failed)))

    def test_both_candidates_were_evaluated_on_every_row(self):
        """A verdict from a favourable subset would not be a verdict."""
        for candidate in CANDIDATES:
            with self.subTest(candidate=candidate):
                rows = self.report["candidates"][candidate]["rows"]
                self.assertEqual(sorted(PUBLISHED_GENE_DROP),
                                 sorted(r["animal"] for r in rows))

    def test_the_verdict_is_d1(self):
        self.assertEqual("D1", self.report["verdict"])

    def test_d1_fits_every_published_row(self):
        self.assertTrue(self.report["candidates"]["D1"]["fits_every_row"])

    def test_d2_is_rejected_decisively_not_marginally(self):
        """
        D2 does not merely lose on points. It misses eleven of thirteen rows by
        0.086 to 0.171 -- two orders of magnitude outside tolerance -- so the
        verdict does not rest on the choice of tolerance.
        """
        rows = self.report["candidates"]["D2"]["rows"]
        self.assertFalse(self.report["candidates"]["D2"]["fits_every_row"])
        missed = [r for r in rows if not r["inside"]]
        self.assertGreaterEqual(len(missed), 10)
        for row in missed:
            with self.subTest(animal=row["animal"]):
                self.assertGreater(abs(row["deviation"]), 20 * row["tolerance"])

    def test_the_residual_bias_is_recorded(self):
        """
        Every non-zero D1 deviation is positive. That is inside the predeclared
        tolerance, so the criterion is met -- but independent seeds agreeing on
        the sign is a systematic difference, not noise, and it is asserted here
        so that it cannot quietly disappear from the record.
        """
        rows = [r for r in self.report["candidates"]["D1"]["rows"]
                if r["published"] > 0]
        self.assertTrue(all(r["deviation"] > 0 for r in rows))
        self.assertTrue(all(abs(r["deviation"]) < 0.002 for r in rows))


class TestTheMechanism(unittest.TestCase):
    """Properties of the estimator that hold independently of the denominator."""

    def test_founders_have_no_ancestral_inbreeding(self):
        got = gene_drop(_rows(), loci=2000, seed=11)
        for animal in (1, 2, 3, 4):
            with self.subTest(animal=animal):
                self.assertEqual(0.0, got[animal])

    def test_the_first_inbred_generation_still_has_none(self):
        """
        Animals 11 and 12 are the first inbred individuals (f = 0.125), and
        their own f_a is 0: an allele identical by descent for the first time
        "contributes to ancestral inbreeding coefficient OF OFFSPRING".
        """
        got = gene_drop(_rows(), loci=5000, seed=12)
        self.assertEqual(0.0, got[11])
        self.assertEqual(0.0, got[12])

    def test_it_is_a_probability(self):
        got = gene_drop(_rows(), loci=2000, seed=13)
        for animal, value in got.items():
            with self.subTest(animal=animal):
                self.assertGreaterEqual(value, 0.0)
                self.assertLessEqual(value, 1.0)

    def test_d2_is_never_below_d1(self):
        """
        A locus with one flagged allele counts a half under D1 and a whole one
        under D2, so the ordering is structural. It is what makes the two
        distinguishable at all.
        """
        d1 = gene_drop(_rows(), loci=5000, seed=14, candidate="D1")
        d2 = gene_drop(_rows(), loci=5000, seed=14, candidate="D2")
        for animal in d1:
            with self.subTest(animal=animal):
                self.assertLessEqual(d1[animal], d2[animal] + 1e-12)

    def test_a_half_founder_pedigree_uses_a_dummy_founder(self):
        """
        This test previously asserted a REFUSAL. Suwanlee's drop assumes every
        non-founder has two parents, and p.490 does not say what a missing
        parental side contributes, so raising was the correct behaviour on the
        evidence then available.

        Baumung et al. (2015) p.102 supplies the rule -- one dummy founder per
        unknown parental slot, given two unique alleles and dropped through by
        ordinary Mendelian segregation -- so the oracle now computes instead of
        refusing. **Suwanlee remains silent**; the authority here is GRAIN.

        Note what is still true: the drop never indexes into "whatever happens to
        be there". The dummy is constructed explicitly, which is what the old
        refusal was protecting against.
        """
        got = gene_drop([(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0)], loci=1000)
        self.assertEqual({1, 2, 3, 4}, set(got))
        self.assertEqual(0.0, got[4])


class TestTheGapToProduction(unittest.TestCase):
    """
    b, now REPAIRED. This class recorded the gap while it existed and now
    guards its closure.

    Before the repair ``dropped_ancestral_inbreeding`` implemented neither
    candidate: it had no allele flag at all, and set a per-locus indicator when
    an IMMEDIATE PARENT was autozygous, so ancestral inbreeding more than one
    generation back was invisible to it. The MEASURED table below is kept as the
    historical record of that state -- it is data about a past implementation,
    not an expectation of the current one, and the sign test that reads it is
    still a true statement about those numbers.
    """

    # Measured: published, adjudicated oracle (D1), PRE-REPAIR production.
    # The third element is history. Nothing here re-runs the old mechanism, so
    # these values can never be reproduced from this tree; they are the reason
    # the repair had to replace the mechanism rather than adjust it.
    MEASURED = {13: (0.125, 0.1247, 0.2120),
                15: (0.375, 0.3746, 0.4905),
                17: (0.585, 0.5861, 0.6188),
                21: (0.728, 0.7284, 0.7112),
                23: (0.822, 0.8211, 0.8192)}

    def test_production_is_now_the_suwanlee_estimator(self):
        """
        The inverse of the assertion this class used to make. It demanded a gap
        of more than 0.05 against the adjudicated oracle; the repaired routine
        agrees with it instead.

        The threshold is loose on purpose. Both sides are Monte Carlo, and this
        runs at a locus count chosen to keep the default suite fast, so it is a
        coarse "same estimator" check. The predeclared statistical comparison
        against the published values lives in
        tests/test_suwanlee_production.py.
        """
        ped = load_corpus("suwanlee_fig1.ped", "asd")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            production = pyp_metrics.dropped_ancestral_inbreeding(
                ped, rounds=20, loci=200, seed=1)
        idmap = ped.idmap
        worst = 0.0
        for animal, (_published, oracle_value, _pre_repair) in self.MEASURED.items():
            renumbered = idmap.get(animal, idmap.get(str(animal)))
            worst = max(worst, abs(float(production[renumbered]) - oracle_value))
        self.assertLess(
            worst, 0.05,
            "production should now agree with the adjudicated D1 estimator; "
            "worst deviation was %.4f" % worst)

    def test_the_pre_repair_difference_had_no_consistent_sign(self):
        """
        Why the repair replaced the mechanism instead of correcting its output.

        The old routine ran HIGH on early generations (0.212 against 0.125) and
        LOW on later ones (0.711 against 0.728). It was a different quantity,
        not a biased version of the same one, so no scaling factor, constant
        offset, denominator patch or clipping could have reconciled them.

        This reads only the recorded table, so it remains a true statement about
        the pre-repair implementation.
        """
        highs = [a for a, (_p, oracle, seen) in self.MEASURED.items() if seen > oracle]
        lows = [a for a, (_p, oracle, seen) in self.MEASURED.items() if seen < oracle]
        self.assertTrue(highs, "expected the old routine to overshoot somewhere")
        self.assertTrue(lows, "expected the old routine to undershoot somewhere")

    def test_production_returns_a_probability_for_every_animal(self):
        """
        The shape of the result, unchanged by the repair: one entry per animal,
        every value a probability. Not clamped -- a violation here is a defect.
        """
        ped = load_corpus("suwanlee_fig1.ped", "asd")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=50)
        self.assertEqual(len(ped.pedigree), len(got))
        for value in got.values():
            self.assertGreaterEqual(float(value), 0.0)
            self.assertLessEqual(float(value), 1.0)


class TestTolerance(unittest.TestCase):

    def test_rounding_dominates_where_the_estimate_is_precise(self):
        """
        With no sampling error the acceptance half-width is exactly the
        publication rounding interval -- the two are combined additively, not
        in quadrature, because they are different kinds of uncertainty.
        """
        self.assertAlmostEqual(5e-4, tolerance(0.0, 0.0), places=12)

    def test_sampling_error_widens_it(self):
        self.assertGreater(tolerance(0.5, 1e-3), tolerance(0.5, 0.0))


if __name__ == "__main__":
    unittest.main()
