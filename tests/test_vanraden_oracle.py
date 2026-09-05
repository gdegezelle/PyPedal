"""
VanRaden (2008) Method 1: the oracle, its hand fixture, and the adjudicated
properties.

**The paper publishes no worked G.** Unlike Boichard, Lacy and Ballou there is
no table of numbers to reproduce, so this oracle is *definition-derived*: it
encodes the equations from article p.4416 and is validated against a fixture
small enough to verify mentally plus the algebra those equations imply. That is
a weaker footing than the pedigree oracles have, and saying so is part of the
point -- C3's status is honestly lower than 's.

Production is checked against the oracle here too, since the SNP path already
exists; only the guards and the missing APIs are added in later commits.
"""
import os
import sys
import unittest

import numpy as np

from PyPedal import pyp_newclasses, pyp_snp
from _pedhelpers import owned_temp_dir

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "oracles"))
from oracle_vanraden import (                 # noqa: E402
    HAND_COUNTS,
    HAND_F_G,
    HAND_G,
    HAND_P,
    RARE_COUNTS,
    RARE_F_G,
    RARE_P,
    genomic_inbreeding,
    grm_method_1,
    scaling_denominator,
    wright_relationships,
    z_matrix,
)


def _load_with_genotypes(pedigree_rows, genotype_rows):
    tmp = owned_temp_dir(prefix="pypedal_vr_")
    pedfile = os.path.join(tmp, "ped.ped")
    snpfile = os.path.join(tmp, "geno.txt")
    with open(pedfile, "w", encoding="utf-8") as fh:
        fh.write("\n".join(pedigree_rows) + "\n")
    with open(snpfile, "w", encoding="utf-8") as fh:
        fh.write("\n".join(genotype_rows) + "\n")
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        return pyp_newclasses.loadPedigree(options={
            "pedfile": pedfile, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0, "renumber": True,
            "snpfile": snpfile})
    finally:
        os.chdir(cwd)


class TestTheHandFixture(unittest.TestCase):
    """
    Three individuals, four loci, every frequency 0.5. P is all ones, so
    Z = counts - 1 and the denominator is 2. Every entry is checkable mentally.
    """

    def test_z_is_counts_minus_twice_p(self):
        want = np.array([[-1.0, 0.0, 1.0, 0.0],
                         [1.0, 0.0, -1.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0]])
        np.testing.assert_allclose(want, z_matrix(HAND_COUNTS, HAND_P))

    def test_the_denominator(self):
        self.assertAlmostEqual(2.0, scaling_denominator(HAND_P), places=12)

    def test_g(self):
        np.testing.assert_allclose(np.array(HAND_G),
                                   grm_method_1(HAND_COUNTS, HAND_P))

    def test_genomic_inbreeding_is_the_diagonal_minus_one(self):
        np.testing.assert_allclose(np.array(HAND_F_G),
                                   genomic_inbreeding(HAND_COUNTS, HAND_P))

    def test_an_all_heterozygote_sits_exactly_on_the_lower_bound(self):
        """
        Individual 3 is heterozygous at every locus, so it carries no deviation
        from the base population: G_33 = 0 and F_g = -1 exactly. That the lower
        bound is attainable is why it must be `>= -1`, not `> -1`.
        """
        self.assertEqual(-1.0, genomic_inbreeding(HAND_COUNTS, HAND_P)[2])


class TestTheAdjudicatedProperties(unittest.TestCase):
    """
    The six claims C3 was asked to confirm or reject, each checked on random
    inputs as well as the fixture so that a property is not mistaken for a
    coincidence of one example.
    """

    @staticmethod
    def _random_case(rng, n=6, m=25):
        counts = rng.integers(0, 3, size=(n, m))
        p = rng.uniform(0.05, 0.95, size=m)
        return counts, p

    def test_g_is_symmetric(self):
        rng = np.random.default_rng(7)
        for _ in range(20):
            counts, p = self._random_case(rng)
            g = grm_method_1(counts, p)
            np.testing.assert_allclose(g, g.T, atol=1e-12)

    def test_the_denominator_is_positive_when_a_locus_is_polymorphic(self):
        self.assertGreater(scaling_denominator([0.0, 0.5, 1.0]), 0.0)
        self.assertEqual(0.0, scaling_denominator([0.0, 1.0, 0.0]))

    def test_the_diagonal_is_never_negative(self):
        """G_ii = ||z_i||^2 / c with c > 0."""
        rng = np.random.default_rng(8)
        for _ in range(20):
            counts, p = self._random_case(rng)
            self.assertTrue(np.all(np.diag(grm_method_1(counts, p)) >= -1e-12))

    def test_genomic_inbreeding_is_never_below_minus_one(self):
        rng = np.random.default_rng(9)
        for _ in range(20):
            counts, p = self._random_case(rng)
            self.assertTrue(np.all(genomic_inbreeding(counts, p) >= -1.0 - 1e-12))

    def test_there_is_no_generic_finite_upper_bound(self):
        """
        As p -> 0 with the individual homozygous for the rare allele, the
        numerator tends to 4 while its own denominator term tends to 0. The
        paper concurs qualitatively: "the genomic inbreeding coefficient is
        greater if the individual is homozygous for rare alleles than if
        homozygous for common alleles."

        So the pedigree-inbreeding validators must never be applied to F_g.
        """
        self.assertAlmostEqual(RARE_F_G,
                               genomic_inbreeding(RARE_COUNTS, RARE_P)[0],
                               places=6)
        previous = 0.0
        for p in (0.1, 0.01, 0.001, 0.0001):
            got = genomic_inbreeding([[2]], [p])[0]
            self.assertGreater(got, previous)
            previous = got
        self.assertGreater(previous, 1000.0)

    def test_wright_scaled_relationships_are_the_correlation_form(self):
        g = grm_method_1(HAND_COUNTS, HAND_P)
        w = wright_relationships(HAND_COUNTS, HAND_P)
        self.assertAlmostEqual(g[0][1] / np.sqrt(g[0][0] * g[1][1]), w[0][1],
                               places=12)
        self.assertTrue(np.isnan(w[2][2]),
                        "a zero diagonal leaves the Wright form undefined, and "
                        "that must be visible rather than reported as zero")


class TestMonomorphicInput(unittest.TestCase):

    def test_a_monomorphic_locus_is_harmless(self):
        """
        It contributes 0 to ZZ' and 0 to the denominator, so Method 1 tolerates
        it. Adding one must not change G at all.
        """
        counts = [[0, 1, 2, 1, 0], [2, 1, 0, 1, 0], [1, 1, 1, 1, 0]]
        p = HAND_P + [0.0]
        np.testing.assert_allclose(np.array(HAND_G), grm_method_1(counts, p),
                                   atol=1e-12)

    def test_all_monomorphic_input_is_undefined(self):
        with self.assertRaises(ZeroDivisionError):
            grm_method_1([[0, 2], [0, 2]], [0.0, 1.0])


class TestFrequenciesAreValidated(unittest.TestCase):
    """
    ``p`` is an argument, not an estimate, and it is not blindly trusted --
    source-faithful and validated are different things.
    """

    def test_one_frequency_per_locus_is_required(self):
        with self.assertRaises(ValueError):
            z_matrix(HAND_COUNTS, [0.5, 0.5])

    def test_frequencies_must_be_finite(self):
        with self.assertRaises(ValueError):
            z_matrix(HAND_COUNTS, [0.5, 0.5, float("nan"), 0.5])

    def test_frequencies_must_be_probabilities(self):
        for bad in ([0.5, 0.5, 1.5, 0.5], [0.5, 0.5, -0.1, 0.5]):
            with self.subTest(p=bad):
                with self.assertRaises(ValueError):
                    z_matrix(HAND_COUNTS, bad)


class TestProductionMatchesTheOracle(unittest.TestCase):
    """
    ``form_grm_from_snp`` at its default is VanRaden Method 1 with the sample
    frequency estimate. Given the same frequencies, it must equal the oracle.
    """

    PEDIGREE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
    GENOTYPES = ["1 chip1 10 0120120120",
                 "2 chip1 10 1201201201",
                 "3 chip1 10 2012012012",
                 "4 chip1 10 0000000000"]

    def _counts_and_sample_p(self):
        counts = np.array([[int(c) for c in row.split()[3]]
                           for row in self.GENOTYPES], dtype=float)
        return counts, counts.sum(axis=0) / (2.0 * counts.shape[0])

    def test_default_grm_equals_method_1_with_the_sample_estimate(self):
        ped = _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)
        got = pyp_snp.form_grm_from_snp(ped)
        counts, p = self._counts_and_sample_p()
        np.testing.assert_allclose(grm_method_1(counts, p), np.asarray(got),
                                   atol=1e-12)

    def test_the_default_is_symmetric_and_its_diagonal_is_non_negative(self):
        ped = _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)
        g = np.asarray(pyp_snp.form_grm_from_snp(ped))
        np.testing.assert_allclose(g, g.T, atol=1e-12)
        self.assertTrue(np.all(np.diag(g) >= -1e-12))


if __name__ == "__main__":
    unittest.main()


class TestTheC3Guards(unittest.TestCase):
    """
    The two refusals. Both are derived algebraically from Method 1, not quoted
    from VanRaden -- the paper is silent on each because there is nothing for it
    to say, which is a different thing from the paper leaving a choice open.
    """

    PEDIGREE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
    GENOTYPES = ["1 chip1 4 0120", "2 chip1 4 1201",
                 "3 chip1 4 2012", "4 chip1 4 0000"]
    MONOMORPHIC = ["1 chip1 4 0000", "2 chip1 4 0000",
                   "3 chip1 4 0000", "4 chip1 4 0000"]

    def test_scale_m_false_is_refused(self):
        """
        It produced Z = counts - p rather than counts - 2p, so G_jj - 1 was not
        a genomic inbreeding coefficient. An option in name only.
        """
        from PyPedal import pyp_errors as errors
        ped = _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)
        with self.assertRaises(errors.PyPedalUsageError):
            pyp_snp.form_grm_from_snp(ped, scale_m=False)
        with self.assertRaises(errors.PyPedalUsageError):
            pyp_snp.form_m_matrix_from_snp(ped, scale_m=False)

    def test_an_unimplemented_method_is_refused_not_silently_substituted(self):
        from PyPedal import pyp_errors as errors
        ped = _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)
        for method in (2, 3):
            with self.subTest(method=method):
                with self.assertRaises(errors.PyPedalUsageError):
                    pyp_snp.form_grm_from_snp(ped, method=method)

    def test_all_monomorphic_input_raises_rather_than_dividing_by_zero(self):
        from PyPedal import pyp_validate
        ped = _load_with_genotypes(self.PEDIGREE, self.MONOMORPHIC)
        with self.assertRaises(pyp_validate.PyPedalValidationError):
            pyp_snp.form_grm_from_snp(ped)

    def test_a_single_monomorphic_locus_is_still_fine(self):
        """
        Guard on the guard. Method 1 tolerates monomorphic loci individually --
        only the degenerate all-monomorphic case is undefined -- so a guard that
        fired on any of them would be a regression dressed as caution.
        """
        genotypes = ["1 chip1 4 0120", "2 chip1 4 1200",
                     "3 chip1 4 2010", "4 chip1 4 0000"]
        ped = _load_with_genotypes(self.PEDIGREE, genotypes)
        g = np.asarray(pyp_snp.form_grm_from_snp(ped))
        self.assertEqual((4, 4), g.shape)
        self.assertTrue(np.all(np.isfinite(g)))


class TestBasePopulationFrequencies(unittest.TestCase):
    """
    VanRaden p.4416 requires P from the unselected base population; PyPedal
    estimated it from the genotyped sample and had no way to be told otherwise.
    The sample estimate stays as a documented fallback -- it is the paper's own
    "simple estimate" -- but it is an estimate OF the named quantity, not that
    quantity, and p.4419 says estimation biases genomic inbreeding coefficients.
    """

    PEDIGREE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
    GENOTYPES = ["1 chip1 4 0120", "2 chip1 4 1201",
                 "3 chip1 4 2012", "4 chip1 4 0000"]

    def _ped(self):
        return _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)

    def _counts(self):
        return np.array([[int(c) for c in row.split()[3]]
                         for row in self.GENOTYPES], dtype=float)

    def test_supplied_frequencies_are_used_verbatim(self):
        supplied = [0.5, 0.5, 0.5, 0.5]
        got = pyp_snp.form_grm_from_snp(self._ped(), base_frequencies=supplied)
        np.testing.assert_allclose(grm_method_1(self._counts(), supplied),
                                   np.asarray(got), atol=1e-12)

    def test_supplying_frequencies_changes_g(self):
        """
        The point of the argument. If base frequencies made no difference there
        would be nothing for the paper to warn about.
        """
        sample = np.asarray(pyp_snp.form_grm_from_snp(self._ped()))
        base = np.asarray(pyp_snp.form_grm_from_snp(
            self._ped(), base_frequencies=[0.3, 0.4, 0.6, 0.2]))
        self.assertFalse(np.allclose(sample, base))

    def test_a_mapping_is_aligned_by_locus_identity(self):
        """
        Not by iteration order. A dict written in a different order must give
        the same answer -- silently trusting caller order is exactly the failure
        this validation exists to prevent.
        """
        forward = pyp_snp.form_grm_from_snp(
            self._ped(), base_frequencies={0: 0.3, 1: 0.4, 2: 0.6, 3: 0.2})
        shuffled = pyp_snp.form_grm_from_snp(
            self._ped(), base_frequencies={3: 0.2, 1: 0.4, 0: 0.3, 2: 0.6})
        np.testing.assert_allclose(np.asarray(forward), np.asarray(shuffled),
                                   atol=1e-12)

    def test_the_wrong_number_of_frequencies_is_refused(self):
        from PyPedal import pyp_errors as errors
        for bad in ([0.5, 0.5], [0.5] * 5):
            with self.subTest(n=len(bad)):
                with self.assertRaises(errors.PyPedalUsageError):
                    pyp_snp.form_grm_from_snp(self._ped(), base_frequencies=bad)

    def test_a_mapping_missing_a_locus_is_refused(self):
        from PyPedal import pyp_errors as errors
        with self.assertRaises(errors.PyPedalUsageError):
            pyp_snp.form_grm_from_snp(self._ped(),
                                      base_frequencies={0: 0.3, 1: 0.4, 2: 0.6})

    def test_a_mapping_for_a_different_panel_is_refused(self):
        from PyPedal import pyp_errors as errors
        with self.assertRaises(errors.PyPedalUsageError):
            pyp_snp.form_grm_from_snp(
                self._ped(),
                base_frequencies={0: 0.3, 1: 0.4, 2: 0.6, 3: 0.2, 9: 0.5})

    def test_non_finite_and_out_of_range_frequencies_are_refused(self):
        from PyPedal import pyp_errors as errors
        for bad in ([0.5, 0.5, float("nan"), 0.5],
                    [0.5, 0.5, 1.5, 0.5],
                    [0.5, 0.5, -0.1, 0.5]):
            with self.subTest(p=bad):
                with self.assertRaises(errors.PyPedalUsageError):
                    pyp_snp.form_grm_from_snp(self._ped(), base_frequencies=bad)

    def test_the_sample_fallback_is_unchanged(self):
        """Omitting the argument must behave exactly as before."""
        counts = self._counts()
        sample_p = counts.sum(axis=0) / (2.0 * counts.shape[0])
        np.testing.assert_allclose(
            grm_method_1(counts, sample_p),
            np.asarray(pyp_snp.form_grm_from_snp(self._ped())), atol=1e-12)


class TestTheRestoredGenomicApis(unittest.TestCase):
    """
    The two functions missing since the port (audit section 3.1), which left
    PyPedal able to PARSE the ``G`` and ``Y`` pedformat columns while unable to
    compute either.

    Their bounds contracts are deliberately different, and that is the point of
    testing them together: homozygosity is a proportion and really is bounded in
    [0, 1]; F_g is bounded below by -1 and not bounded above at all.
    """

    PEDIGREE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
    GENOTYPES = ["1 chip1 4 0120", "2 chip1 4 1201",
                 "3 chip1 4 2012", "4 chip1 4 0000"]

    def _ped(self):
        return _load_with_genotypes(self.PEDIGREE, self.GENOTYPES)

    def test_genomic_inbreeding_is_the_diagonal_minus_one(self):
        ped = self._ped()
        g = np.asarray(pyp_snp.form_grm_from_snp(ped))
        got = pyp_snp.compute_genomic_inbreeding_from_grm(ped, grm=g)
        for row in range(len(ped.snp)):
            animal_id = ped.snp.iloc[row, 0]
            with self.subTest(animal=animal_id):
                self.assertAlmostEqual(g[row][row] - 1.0, got[animal_id],
                                       places=12)

    def test_it_matches_the_independent_oracle(self):
        ped = self._ped()
        counts = np.array([[int(c) for c in row.split()[3]]
                           for row in self.GENOTYPES], dtype=float)
        p = counts.sum(axis=0) / (2.0 * counts.shape[0])
        want = genomic_inbreeding(counts, p)
        got = pyp_snp.compute_genomic_inbreeding_from_grm(ped)
        for row in range(len(ped.snp)):
            with self.subTest(row=row):
                self.assertAlmostEqual(want[row], got[ped.snp.iloc[row, 0]],
                                       places=12)

    def test_it_populates_the_animals_and_sets_g_computed(self):
        """
        Closes the gap where the G pedformat code parsed but nothing computed
        it.
        """
        ped = self._ped()
        pyp_snp.compute_genomic_inbreeding_from_grm(ped)
        self.assertTrue(ped.kw["g_computed"])
        self.assertTrue(any(a.genomic_inbreeding != 0.0 for a in ped.pedigree))

    def test_a_matrix_of_the_wrong_size_is_refused(self):
        from PyPedal import pyp_errors as errors
        ped = self._ped()
        with self.assertRaises(errors.PyPedalUsageError):
            pyp_snp.compute_genomic_inbreeding_from_grm(ped, grm=np.eye(3))

    def test_homozygosity_is_a_proportion_of_typed_loci(self):
        ped = self._ped()
        got = pyp_snp.compute_genomic_homozygosity_from_snp(ped)
        # "0120": 0 and 2 and 0 are homozygous, 1 is not -> 3/4
        # "0000": all homozygous -> 1.0
        expected = {"0120": 0.75, "1201": 0.5, "2012": 0.75, "0000": 1.0}
        for row in range(len(ped.snp)):
            genotype = ped.snp.iloc[row, 3]
            with self.subTest(genotype=genotype):
                self.assertAlmostEqual(expected[genotype],
                                       got[ped.snp.iloc[row, 0]], places=12)

    def test_homozygosity_is_bounded_in_zero_one(self):
        got = pyp_snp.compute_genomic_homozygosity_from_snp(self._ped())
        for animal_id, value in got.items():
            with self.subTest(animal=animal_id):
                self.assertGreaterEqual(value, 0.0)
                self.assertLessEqual(value, 1.0)

    def _with_genotypes(self, genotypes):
        """
        Substitute genotype strings directly.

        The loader rejects anything but 0, 1 and 2 -- correctly, since a dosage
        cannot be anything else -- so a missing code cannot arrive through it.
        ``missing_code`` exists for genotype frames assembled by other means,
        and testing it means bypassing the reader rather than weakening it.
        """
        ped = self._ped()
        for row, genotype in enumerate(genotypes):
            ped.snp.iloc[row, 3] = genotype
        return ped

    def test_an_untyped_individual_is_missing_not_zero(self):
        """
        Zero typed loci is a measurement that was never made. Reporting 0.0
        from a 0/0 would be a missing result wearing the clothes of a complete
        one -- and 0.0 is a perfectly plausible homozygosity, so nothing
        downstream could tell them apart.
        """
        ped = self._with_genotypes(["0120", "1201", "2012", "NNNN"])
        got = pyp_snp.compute_genomic_homozygosity_from_snp(ped, missing_code="N")
        self.assertEqual(ped.kw["missing_homozygosity"], got[ped.snp.iloc[3, 0]])

    def test_partially_typed_individuals_use_only_typed_loci(self):
        ped = self._with_genotypes(["01N0", "1201", "2012", "0000"])
        got = pyp_snp.compute_genomic_homozygosity_from_snp(ped, missing_code="N")
        # "01N0" -> typed 0, 1, 0; homozygous 0 and 0 -> 2/3
        self.assertAlmostEqual(2.0 / 3.0, got[ped.snp.iloc[0, 0]], places=12)

    def test_the_loader_still_refuses_a_missing_code(self):
        """
        missing_code is a computation-time accommodation, not a licence to
        loosen the reader. A dosage really cannot be N.
        """
        from PyPedal import pyp_errors as errors
        with self.assertRaises(errors.PyPedalInputError):
            _load_with_genotypes(self.PEDIGREE,
                                 ["1 chip1 4 0120", "2 chip1 4 1201",
                                  "3 chip1 4 2012", "4 chip1 4 NNNN"])

    def test_the_two_contracts_are_not_interchangeable(self):
        """
        F_g may exceed 1 and may be negative; homozygosity may do neither. A
        single "genomic" bounds check would be wrong for one of them.
        """
        from PyPedal import pyp_validate
        ped = self._ped()
        pyp_validate.check_genomic_inbreeding(ped, {1: 50.0, 2: -0.5},
                                              "test")   # must not raise
        with self.assertRaises(pyp_validate.PyPedalValidationError):
            pyp_validate.check_genomic_inbreeding(ped, {1: -1.5}, "test")
