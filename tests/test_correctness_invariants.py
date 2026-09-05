"""
Class A -- CORRECTNESS tests.

Every assertion here is backed by an independent
mathematical derivation, a published textbook value, or a definitional invariant
of the quantity being computed. None of them is justified by "Python 2 does it
too" and none is justified by an existing PyPedal expected value.

The distinction matters. ``tests/test_legacy_numeric_compatibility.py`` pins values that
the two implementations agree on; agreement is evidence of a faithful migration,
not of mathematical correctness. This file is the one that would have caught the
defects the audit found.

Scope discipline: an invariant is asserted only for the class of algorithm it is
proven for. ``F in [0, 1]`` holds for a *pedigree* inbreeding coefficient, which
is a probability; it does not hold for a genomic estimator, which may be
negative. ``A[i,i] == 1 + F_i`` holds for the *pedigree* NRM, not for a
SNP-derived genomic relationship matrix. Preconditions are asserted in the test
body rather than assumed.
"""
import os
import sys
import unittest

from PyPedal import pyp_errors, pyp_nrm, pyp_metrics, pyp_utils, pyp_validate

from _pedhelpers import owned_temp_dir, corpus, load_corpus, load_corpus_from_path, nrm_value
from oracles import ballou_f_a, lacy_f_e, lacy_n_half_founders, oracle_inbreeding

COI_METHODS = ("tabular", "vanraden", "meu_luo", "mod_meu_luo")


def coi(ped, method):
    result = pyp_nrm.inbreeding(ped, method=method, output=False)
    if isinstance(result, tuple):
        result = result[0]
    return {int(k): float(v) for k, v in result["fx"].items()}


class TestInbreedingAgainstIndependentOracle(unittest.TestCase):
    """
    The oracle in ``the independent oracle`` does not import
    PyPedal. It implements A = LDL' from first principles and cross-checks three
    traversals of the same mathematics against each other, so a PyPedal defect
    cannot propagate into the expected value.
    """

    def test_mrode_animal_5_is_one_eighth(self):
        """Mrode (2005) Table 2.1: F(5) = 0.125, a published textbook value."""
        ped = load_corpus("mrode.ped")
        self.assertEqual(6, len(ped.pedigree))
        expected = oracle_inbreeding(corpus("mrode.ped"))
        self.assertAlmostEqual(0.125, expected[5], places=12)
        for method in COI_METHODS:
            with self.subTest(method=method):
                self.assertAlmostEqual(0.125, coi(ped, method)[5], places=12)

    def test_mrode_every_animal_matches_oracle(self):
        ped = load_corpus("mrode.ped")
        expected = oracle_inbreeding(corpus("mrode.ped"))
        got = coi(ped, "tabular")
        for animal_id, want in expected.items():
            with self.subTest(animal=animal_id):
                self.assertAlmostEqual(want, got[animal_id], places=12)

    def test_hartlandclark_matches_oracle_on_all_methods(self):
        ped = load_corpus("hartlandclark.ped")
        self.assertEqual(15, len(ped.pedigree))
        expected = oracle_inbreeding(corpus("hartlandclark.ped"))
        self.assertAlmostEqual(0.14453125, expected[15], places=12)
        for method in COI_METHODS:
            with self.subTest(method=method):
                got = coi(ped, method)
                for animal_id, want in expected.items():
                    self.assertAlmostEqual(want, got[animal_id], places=12,
                                           msg="animal %d" % animal_id)


class TestFinding3DuplicateAncestorGuard(unittest.TestCase):
    """
    , repaired in Commit 4.

    ``inbreeding_modified_meuwissen_luo`` lost the ``not in ancs`` / ``not in
    ancd`` membership guard in the migration, so an ancestor reachable by more
    than one path was pushed onto the ancestor list repeatedly and its
    contribution counted once per path.

    Three regression cases, deliberately chosen to fail in different ways:

      hartlandclark.ped  a PLAUSIBLE wrong answer -- 0.20703125 against a true
                         0.14453125, 43% high, inside [0, 1] and therefore
                         invisible to any range check;
      doug.ped           an EXTREME violation -- F(41) = 4.72, impossible for a
                         probability;
      new_ids.ped        a wrong answer that disagrees with PyPedal 2.0.4 --
                         F = 1.0 against 0.5 on 7 of 17 values.

    Together they show why an output-range postcondition alone is insufficient
    (audit Findings 21 vs 22): it catches doug and new_ids but not
    hartlandclark. What catches all three is cross-method agreement, since the
    four algorithms compute one quantity.
    """

    def test_hartlandclark_animal_15_matches_the_oracle(self):
        ped = load_corpus("hartlandclark.ped")
        self.assertAlmostEqual(0.14453125, coi(ped, "mod_meu_luo")[15], places=12)

    def test_all_four_methods_agree_across_the_corpus(self):
        """One quantity, four algorithms: they must all produce the same F."""
        for name, pedformat in (("mrode.ped", "asd"),
                                ("new_lacy.ped", "asd"),
                                ("hartlandclark.ped", "asdb"),
                                ("generations.ped", "asdbx"),
                                ("boichard2a.ped", "asdg")):
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                expected = oracle_inbreeding(corpus(name))
                for method in COI_METHODS:
                    got = coi(ped, method)
                    for animal_id, want in expected.items():
                        self.assertAlmostEqual(
                            want, got[animal_id], places=12,
                            msg="%s / %s / animal %d" % (name, method, animal_id))

    def test_string_pedigrees_agree_across_methods(self):
        """
        doug.ped and new_ids.ped have string IDs, so the integer-parsing oracle
        cannot read them. Cross-method agreement is still available and is the
        evidence that caught the extreme violations here.
        """
        for name, pedformat, sepchar in (("doug.ped", "ASDx", " "),
                                         ("new_ids.ped", "ASD", " "),
                                         ("horse.ped", "ASD", ",")):
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat, sepchar=sepchar)
                reference = coi(ped, "tabular")
                for method in ("vanraden", "meu_luo", "mod_meu_luo"):
                    got = coi(ped, method)
                    for animal_id, want in reference.items():
                        self.assertAlmostEqual(
                            want, got[animal_id], places=12,
                            msg="%s / %s / animal %d" % (name, method, animal_id))

    def test_repeated_common_ancestor_pedigree(self):
        """
        A purpose-built worst case: a full-sib mating whose parents share both
        grandparents, so the common ancestors are reachable by four paths each.
        This is precisely the structure the dropped guard mishandled.

            1, 2        unrelated founders
            3, 4        full sibs, both out of 1 x 2
            5, 6        full sibs, both out of 3 x 4
            7           out of 5 x 6
        """
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4", "6 3 4", "7 5 6"]
        tmp = owned_temp_dir(prefix="pypedal_sib_")
        path = os.path.join(tmp, "fullsib.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(rows) + "\n")

        expected = oracle_inbreeding(path)
        # sanity: this pedigree must actually be inbred, or it proves nothing
        self.assertGreater(expected[5], 0.0)
        self.assertGreater(expected[7], expected[5])

        ped = load_corpus_from_path(path, "asd")
        for method in COI_METHODS:
            with self.subTest(method=method):
                got = coi(ped, method)
                for animal_id, want in expected.items():
                    self.assertAlmostEqual(want, got[animal_id], places=12,
                                           msg="%s / animal %d" % (method, animal_id))


class TestLacyEffectiveFounders(unittest.TestCase):
    """
    Lacy (1989) f_e = 1 / sum(q_k^2). The oracle's validity gate is that the
    founder contribution vector q sums to 1 **without being normalised** --
    PyPedal divides by the total, which makes sum-to-one true by construction
    and therefore useless as a check.

    Only pedigrees that pass the gate are asserted here. ``mrode.ped`` and
    ``hartlandclark.ped`` contain half-founders the oracle does not model and
    remain unadjudicated (audit ).
    """

    def _assert_fe(self, pedfile, pedformat):
        want, gate = lacy_f_e(corpus(pedfile))
        self.assertTrue(gate, "oracle precondition failed: q does not sum to 1 "
                              "unforced, so this pedigree is not adjudicable")
        ped = load_corpus(pedfile, pedformat)
        got = pyp_metrics.effective_founders_lacy(ped)["fa_effective_founders"]
        self.assertAlmostEqual(want, got, places=9)

    def test_new_lacy(self):
        self._assert_fe("new_lacy.ped", "asd")

    def test_generations(self):
        self._assert_fe("generations.ped", "asdbx")


class TestBallouAncestralInbreeding(unittest.TestCase):
    """
    Stage A5a. ``ballou_ancestral_inbreeding`` had **no caller and no test
    anywhere in the repository** before this class existed.

    Ballou's recursion, **equation 1 of Ballou (1997) article p.170**::

        f_a(x) = [ f_a(s) + (1 - f_a(s)) * F_s
                 + f_a(d) + (1 - f_a(d)) * F_d ] / 2

    **Provenance ceiling lifted.** This was previously grounded only in the
    transcription at ``the historical manual`` line 88, a secondary
    source. The paper has since been read and the transcription was correct, so
    these assertions are now paper-validated rather than
    consistency-with-a-transcription. The published values themselves are
    checked in ``tests/test_published_fixtures.py``.

    **What is deliberately absent.** No test in this class, or anywhere else,
    compares Ballou against ``dropped_ancestral_inbreeding``. They are different
    estimators, and this is no longer a caution but a published result: Suwanlee
    et al. (2007) Figure 1 prints both on one pedigree and they differ from the
    third generation onward. A systematic gap is the EXPECTED result, so
    agreement is not a correctness criterion and disagreement is not a defect.
    Nor may the reverse ordering become an invariant -- Suwanlee proves no
    theorem. See audit item .
    """

    # Hand-derived from the recurrence. Full-sib mating of unrelated non-inbred
    # parents gives F = 0.25, which is standard; everything else follows.
    HAND_DERIVED = {
        # 5 is the first inbred animal (F=0.25) and so has f_a = 0: its own
        # inbreeding is not "in the past" from its own point of view.
        "ballou_simple.ped": {5: 0.0, 7: 0.125, 9: 0.0625},
        # 7 = 5 x 6, both parents inbred at 0.25:
        #   f_a(7) = [0 + (1-0)(0.25) + 0 + (1-0)(0.25)] / 2 = 0.25
        "ballou_two_gen.ped": {5: 0.0, 6: 0.0, 7: 0.25},
        # 6 is inbred (F=0.125, half sibs via founder 1) but both parents are
        # non-inbred, so every f_a in the pedigree is 0.
        "ballou_noninbred.ped": {4: 0.0, 5: 0.0, 6: 0.0},
    }

    CORPUS = (("mrode.ped", "asd"), ("hartlandclark.ped", "asdb"),
              ("new_lacy.ped", "asd"), ("generations.ped", "asdbx"),
              ("boichard2a.ped", "asdg"), ("ballou_simple.ped", "asd"),
              ("ballou_two_gen.ped", "asd"), ("ballou_noninbred.ped", "asd"))

    def test_oracle_reproduces_hand_derived_values(self):
        """The oracle must be right before it is used to judge anything."""
        for name, expected in self.HAND_DERIVED.items():
            got = ballou_f_a(corpus(name))
            for animal_id, want in expected.items():
                with self.subTest(pedigree=name, animal=animal_id):
                    self.assertAlmostEqual(want, got[animal_id], places=12)

    def test_pypedal_agrees_with_the_oracle(self):
        """
        Both are deterministic, so exact agreement is the right bar -- there is
        no sampling error to allow for.
        """
        for name, pedformat in self.CORPUS:
            with self.subTest(pedigree=name):
                want = ballou_f_a(corpus(name))
                ped = load_corpus_from_path(corpus(name), pedformat)
                got = pyp_metrics.ballou_ancestral_inbreeding(ped)
                by_original = {
                    int(ped.pedigree[int(k) - 1].originalID): float(v)
                    for k, v in got.items()
                }
                self.assertEqual(set(want), set(by_original))
                for animal_id, value in want.items():
                    self.assertAlmostEqual(value, by_original[animal_id], places=12)

    def test_f_a_is_a_probability(self):
        """Definitional: f_a is the probability of having inherited an allele
        that has undergone inbreeding at least once."""
        for name, pedformat in self.CORPUS:
            with self.subTest(pedigree=name):
                ped = load_corpus_from_path(corpus(name), pedformat)
                for value in pyp_metrics.ballou_ancestral_inbreeding(ped).values():
                    self.assertGreaterEqual(float(value), 0.0)
                    self.assertLessEqual(float(value), 1.0)

    def test_founders_have_zero_ancestral_inbreeding(self):
        """An animal with no known parents has no ancestors to have been inbred."""
        for name, pedformat in self.CORPUS:
            with self.subTest(pedigree=name):
                ped = load_corpus_from_path(corpus(name), pedformat)
                result = pyp_metrics.ballou_ancestral_inbreeding(ped)
                missing = ped.kw["missing_parent"]
                for animal in ped.pedigree:
                    if (str(animal.sireID) == str(missing)
                            and str(animal.damID) == str(missing)):
                        self.assertEqual(0.0, float(result[animal.animalID]))


class TestLacyOracleHalfFounderHandling(unittest.TestCase):
    """
    Audit item . These assert properties of the **oracle**, not of PyPedal:
    they establish that the oracle is fit to adjudicate a half-founder pedigree
    before anything is adjudicated with it.

    The validity gate is that the contribution vector q sums to 1 without being
    normalised. An animal with exactly one known parent carries founder genome
    through its unknown side, so under the strict definition ('lacy') that genome
    belongs to no contributor and q falls short of 1. Naming each unknown parent
    as its own founder ('phantom') is faithful to the definition -- a founder is
    an animal whose parents are unknown -- and restores the gate.

    Satisfying the gate is a **necessary** condition, not a proof of agreement
    with Lacy (1989).
    """

    HALF_FOUNDER_FREE = (("new_lacy.ped", 2.909090909090909),
                         ("generations.ped", 4.612612612612613),
                         ("boichard2a.ped", 4.0))
    WITH_HALF_FOUNDERS = (("mrode.ped", 1), ("hartlandclark.ped", 4))

    def test_half_founder_counts_are_as_documented(self):
        for name, _fe in self.HALF_FOUNDER_FREE:
            with self.subTest(pedigree=name):
                self.assertEqual(0, lacy_n_half_founders(corpus(name)))
        for name, expected in self.WITH_HALF_FOUNDERS:
            with self.subTest(pedigree=name):
                self.assertEqual(expected, lacy_n_half_founders(corpus(name)))

    def test_phantom_mode_is_a_noop_without_half_founders(self):
        """With nothing to complete, all three modes must coincide."""
        for name, expected in self.HALF_FOUNDER_FREE:
            with self.subTest(pedigree=name):
                values = {}
                for mode in ("lacy", "half", "phantom"):
                    f_e, gate = lacy_f_e(corpus(name), mode=mode)
                    self.assertTrue(gate, f"{mode} gate failed on {name}")
                    values[mode] = f_e
                for mode, f_e in values.items():
                    self.assertAlmostEqual(expected, f_e, places=9, msg=mode)

    def test_strict_mode_gate_fails_exactly_on_half_founder_pedigrees(self):
        """The gate must discriminate, or it is not doing any work."""
        for name, _fe in self.HALF_FOUNDER_FREE:
            with self.subTest(pedigree=name, mode="lacy"):
                self.assertTrue(lacy_f_e(corpus(name), mode="lacy")[1])
        for name, _n in self.WITH_HALF_FOUNDERS:
            with self.subTest(pedigree=name, mode="lacy"):
                self.assertFalse(lacy_f_e(corpus(name), mode="lacy")[1])

    def test_phantom_mode_restores_the_gate_on_half_founder_pedigrees(self):
        for name, _n in self.WITH_HALF_FOUNDERS:
            with self.subTest(pedigree=name):
                f_e, gate = lacy_f_e(corpus(name), mode="phantom")
                self.assertTrue(
                    gate,
                    "phantom completion must make q a probability vector; "
                    "without that the oracle cannot check Lacy's probability vector",
                )
                self.assertGreaterEqual(f_e, 1.0)

    def test_pypedal_half_flag_semantics_do_not_yield_a_probability_vector(self):
        """
        Characterisation, not a defect claim about the oracle: mirroring
        PyPedal's ``half=True`` -- promoting the whole half-founder into the
        contributor set -- overshoots, because the true founders upstream of it
        are credited as well. Recorded so the measurement is not repeated.
        """
        for name, _n in self.WITH_HALF_FOUNDERS:
            with self.subTest(pedigree=name):
                self.assertFalse(lacy_f_e(corpus(name), mode="half")[1])


class TestEffectiveFounderNumberIsBounded(unittest.TestCase):
    """
    f_e = 1 / sum(q_k^2) over a contribution vector q of length k that is
    non-negative and sums to 1. Both bounds follow from that:

      f_e >= 1   because sum(q_k^2) <= 1 for any probability vector;
      f_e <= k   because equal contributions minimise sum(q_k^2) at 1/k.

    In words: the effective number of founders is at least one and never
    exceeds the actual number of founders, reaching it only when every founder
    contributes equally. The upper bound is the sharper of the two -- it is
    what distinguishes a merely implausible result from an impossible one.

    **k counts founder SOURCES, and this test once got that wrong.** The caveat
    that used to stand here warned that ``count(founder == 'y')`` understates k
    on a half-founder pedigree -- an animal's unknown side is a further,
    independent source -- and predicted that "a future failure would need
    checking against the k question before being treated as a defect".

    That failure duly arrived when  was settled. Lacy (1989) p.113 makes the
    unknown parent a founder in its own right, so on ``mrode.ped`` there are
    three sources, not two, and f_e = 2.798 is comfortably inside 3 while
    violating the old bound of 2. The test now counts sources the way the paper
    does; the mathematics never changed, only the k being fed to it.

    This also reinstates the bound withdrawn in ``d1129c4``, which was withdrawn
    because it could not be justified with the k then available.

    Each routine is called on a FRESHLY LOADED pedigree. That matters: the
    effective-ancestor routines mutate the pedigree they are given (see
    test_known_defects.py), so a value measured after one of them has run says
    nothing about the routine that produced it.
    """

    def test_lacy_f_e_is_within_bounds(self):
        for name, pedformat in (("new_lacy.ped", "asd"),
                                ("generations.ped", "asdbx"),
                                ("mrode.ped", "asd"),
                                ("hartlandclark.ped", "asdb")):
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                # Founder SOURCES: animals with no known parent, plus one per
                # unknown parent of an animal that has one (Lacy p.113).
                sources = (sum(1 for a in ped.pedigree
                               if getattr(a, "founder", "n") == "y")
                           + len(pyp_metrics.lacy_phantom_slots(ped)))
                self.assertGreater(sources, 0)
                f_e = pyp_metrics.effective_founders_lacy(ped)["fa_effective_founders"]
                self.assertGreaterEqual(f_e, 1.0)
                self.assertLessEqual(f_e, sources)

    def test_boichard_f_e_is_within_bounds(self):
        """
        ``boichard2a.ped`` is the only corpus pedigree carrying a generation
        column, which this routine's documented precondition requires.
        """
        ped = load_corpus("boichard2a.ped", "asdg")
        founders = sum(1 for a in ped.pedigree
                       if getattr(a, "founder", "n") == "y")
        self.assertEqual(4, founders)
        f_e = float(pyp_metrics.a_effective_founders_boichard(ped))
        self.assertGreaterEqual(f_e, 1.0)
        self.assertLessEqual(f_e, founders)


class TestNrmDefinitionalProperties(unittest.TestCase):
    """
    Properties true of *any* numerator relationship matrix, by construction. The
    NRM is a covariance-structure matrix, so it is symmetric; and its diagonal
    is A_ii = 1 + F_i by the definition of the inbreeding coefficient.

    These are asserted on small pedigrees only. A full O(n^2) scan is a test
    activity, not something a production call should ever do -- see the layering
    rule in ``PyPedal/pyp_validate.py``.
    """

    SMALL = (("mrode.ped", "asd"), ("new_lacy.ped", "asd"),
             ("hartlandclark.ped", "asdb"), ("generations.ped", "asdbx"))

    def test_nrm_is_symmetric(self):
        for name, fmt in self.SMALL:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, fmt)
                a = pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="dense")
                n = len(ped.pedigree)
                for i in range(n):
                    for j in range(i + 1, n):
                        self.assertAlmostEqual(
                            nrm_value(a, i, j), nrm_value(a, j, i), places=12,
                            msg="A[%d,%d] != A[%d,%d] in %s" % (i, j, j, i, name))

    def test_nrm_diagonal_is_one_plus_inbreeding(self):
        """Pedigree NRM only. A SNP-derived GRM does not satisfy this."""
        for name, fmt in self.SMALL:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, fmt)
                a = pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="dense")
                fx = coi(ped, "tabular")
                for i in range(len(ped.pedigree)):
                    # renumbered pedigree: index i corresponds to animalID i+1
                    self.assertAlmostEqual(
                        1.0 + fx[i + 1], nrm_value(a, i, i), places=12,
                        msg="A[%d,%d] != 1 + F(%d) in %s" % (i, i, i + 1, name))


class TestOrderingPrecondition(unittest.TestCase):
    """
    ``fast_a_matrix`` computes the NRM in a single forward pass, which is only
    valid when every parent precedes its offspring in the list. That is a
    precondition of the algorithm, not an incidental property, so it is asserted
    rather than assumed.
    """

    def test_parent_precedes_offspring_after_reorder_and_renumber(self):
        for name in ("mrode.ped", "new_lacy.ped", "hartlandclark.ped",
                     "generations.ped", "boichard2a.ped"):
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                position = {int(a.animalID): idx
                            for idx, a in enumerate(ped.pedigree)}
                violations = []
                for idx, animal in enumerate(ped.pedigree):
                    for parent in (int(animal.sireID), int(animal.damID)):
                        if parent and parent in position and position[parent] > idx:
                            violations.append((int(animal.animalID), parent))
                self.assertEqual([], violations,
                                 "offspring precedes parent in %s" % name)

    def test_renumbered_index_matches_animal_id(self):
        """The ``index == animalID - 1`` contract the analysis layer relies on."""
        for name in ("mrode.ped", "hartlandclark.ped", "boichard2a.ped"):
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                for idx, animal in enumerate(ped.pedigree):
                    self.assertEqual(idx + 1, int(animal.animalID))


class TestStringIdLoadRoundTrip(unittest.TestCase):
    """
    , repaired in Commit 3.

    The uppercase pedformat codes ``A``/``S``/``D`` mark a column as holding a
    string ID that must be hashed to an integer. ``NewAnimal.__init__`` lost
    that branch and applied ``int()`` unconditionally, and ``preprocess()``
    swallowed the resulting ``ValueError``, so every string pedigree loaded as
    an **empty** pedigree with no error reported.

    The invariant asserted here is a load round-trip, independent of any other
    implementation: a file with N data rows cannot produce fewer than N
    animals, and every ID in the file must be present afterwards. It is a
    stronger and more useful statement than "matches Python 2", and in
    particular it makes the silent-empty-pedigree failure impossible to
    reintroduce. The exact record counts are pinned separately, as behavioural
    compatibility, in ``test_legacy_numeric_compatibility.py``.
    """

    STRING_PEDIGREES = (
        ("doug.ped", "ASDx", " "),
        ("new_ids.ped", "ASD", " "),
        ("horse.ped", "ASD", ","),
    )

    @staticmethod
    def _data_rows(path, sepchar):
        rows = []
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#") and not line.startswith("%"):
                    rows.append(line.split(sepchar))
        return rows

    def test_string_pedigree_never_loads_empty(self):
        """A load failure must never present itself as a successful empty load."""
        for name, pedformat, sepchar in self.STRING_PEDIGREES:
            with self.subTest(pedigree=name):
                rows = self._data_rows(corpus(name), sepchar)
                self.assertGreater(len(rows), 0, "fixture has no data rows")
                ped = load_corpus(name, pedformat, sepchar=sepchar)
                self.assertGreater(
                    len(ped.pedigree), 0,
                    "%s has %d data rows but loaded 0 animals" % (name, len(rows)))

    def test_every_animal_in_the_file_is_present_after_load(self):
        for name, pedformat, sepchar in self.STRING_PEDIGREES:
            with self.subTest(pedigree=name):
                rows = self._data_rows(corpus(name), sepchar)
                ped = load_corpus(name, pedformat, sepchar=sepchar)
                # implicit parents may be added, so the count can only grow
                self.assertGreaterEqual(len(ped.pedigree), len(rows))
                loaded = {str(a.name).strip() for a in ped.pedigree}
                for row in rows:
                    self.assertIn(row[0].strip(), loaded,
                                  "animal %r from %s is missing after load"
                                  % (row[0], name))

    def test_string_ids_hash_to_distinct_integers(self):
        """Distinct input IDs must remain distinct animals after hashing."""
        for name, pedformat, sepchar in self.STRING_PEDIGREES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat, sepchar=sepchar)
                ids = [int(a.animalID) for a in ped.pedigree]
                self.assertEqual(len(ids), len(set(ids)))


class TestStringIdCollisionDetection(unittest.TestCase):
    """
    String IDs are mapped to integers by md5 modulo 2**63 - 1. A collision is
    vanishingly unlikely but would silently merge two animals into one and
    corrupt every relationship computed afterwards. Detection costs one dict
    lookup per ID, so it is done rather than argued about.

    The registry is scoped to a single pedigree and covers animal, sire and dam
    IDs together: a collision between an animal ID and a sire ID merges records
    just as effectively as one between two animal IDs.
    """

    def test_distinct_strings_colliding_are_rejected(self):
        from PyPedal import pyp_newclasses
        from PyPedal.pyp_errors import PyPedalStringIDCollisionError

        kw = {}
        first = pyp_newclasses.hashed_string_id("Doug", kw, "animal")
        # Feed the registry a forged entry so a genuine collision is simulated
        # without needing to find a real md5 collision.
        kw[pyp_newclasses._STRING_ID_REGISTRY][first] = "SomeOtherAnimal"
        with self.assertRaises(PyPedalStringIDCollisionError):
            pyp_newclasses.hashed_string_id("Doug", kw, "sire")

    def test_repeating_the_same_string_is_not_a_collision(self):
        from PyPedal import pyp_newclasses

        kw = {}
        first = pyp_newclasses.hashed_string_id("Doug", kw, "animal")
        second = pyp_newclasses.hashed_string_id("Doug", kw, "sire")
        self.assertEqual(first, second)

    def test_registry_is_scoped_to_one_pedigree(self):
        from PyPedal import pyp_newclasses

        first_kw, second_kw = {}, {}
        pyp_newclasses.hashed_string_id("Doug", first_kw, "animal")
        self.assertNotIn(pyp_newclasses._STRING_ID_REGISTRY, second_kw)


class TestAnalysisDoesNotMutateTheCallersPedigree(unittest.TestCase):
    """
    , repaired as a state-integrity fix.

    ``a_effective_ancestors_definite`` and ``_indefinite`` used to erase
    parentage from the pedigree they were given. Both built
    ``list(reversed(pedobj.pedigree))`` -- a shallow copy holding the *same*
    NewAnimal objects -- and then assigned ``missing_parent`` to ``sireID`` and
    ``damID``. Clearing an ancestor's parents once its marginal contribution is
    counted is a legitimate step of the Boichard algorithm, but it was being
    applied to the live pedigree, so the destruction outlived the call. On
    ``boichard2a.ped``, ``_definite`` destroyed 1 of 14 animals' parentage and
    ``_indefinite`` 7 of 14.

    Analysis reads a pedigree; it does not consume it. That is the invariant
    here, and it is independent of the Boichard mathematics.

    Both exit paths are covered, and WHICH fixture supplies each has changed
    TWICE. These routines used to raise ``PyPedalValidationError`` on
    ``boichard2a.ped`` because the value they computed was impossible; since
     was repaired against Appendix B, ``_definite`` returns 2.0 there. The
    raise path then moved to a half-founder pedigree, which  refused.
     is now implemented, so that pedigree returns a value too, and the
    raise path has moved again -- to 's non-antichain guard, which is
    still a refusal and is expected to stay one.

    An early exit through a raise is the path most likely to leave state
    half-modified, so it has to keep being exercised by SOMETHING. That it has
    had to be re-homed twice is exactly why the guard-on-the-guard below
    exists.
    """

    ROUTINES = ("a_effective_ancestors_definite",
                "a_effective_ancestors_indefinite",
                "a_effective_founders_boichard")

    # Two known founders, one animal with a single known parent, and a
    # generation column. Under  this used to reach a refusal; it now
    # computes, and is kept because a newly supported path is exactly where a
    # state leak would be new.
    HALF_FOUNDER_PED = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n"

    # Animal 5 and its own offspring 6 are both in generation 2, so the
    # reference population is not an antichain and  refuses. This is the
    # current supplier of the raise path.
    NON_ANTICHAIN_PED = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 5 4 2\n"

    @staticmethod
    def _write(text, prefix):
        tmp = owned_temp_dir(prefix=prefix)
        path = os.path.join(tmp, "fixture.ped")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return load_corpus_from_path(path, "asdg")

    def _half_founder_pedigree(self):
        return self._write(self.HALF_FOUNDER_PED, "pypedal_r2_")

    def _refusing_pedigree(self):
        return self._write(self.NON_ANTICHAIN_PED, "pypedal_r3_")

    @staticmethod
    def _parentage(ped):
        return [(int(a.animalID), int(a.sireID), int(a.damID))
                for a in ped.pedigree]

    def test_pedigree_survives_the_call(self):
        for name in self.ROUTINES:
            with self.subTest(routine=name):
                ped = load_corpus("boichard2a.ped", "asdg")
                before = self._parentage(ped)
                raised = False
                try:
                    getattr(pyp_metrics, name)(ped)
                except pyp_errors.PyPedalError:
                    raised = True
                self.assertEqual(
                    before, self._parentage(ped),
                    "%s mutated the caller's pedigree (raised=%s)" % (name, raised))

    def test_validation_raise_path_is_covered(self):
        """
        Guard on the guard: if these routines ever stop raising here, the test
        below silently stops exercising the early-exit path and must be given a
        new fixture that does.

        This has fired TWICE, which is the point of having it. 's repair
        made ``boichard2a.ped`` return 2.0 instead of raising, and the raise
        path moved to the half-founder refusal (). 's
        implementation then made the half-founder pedigree return 1.6, and the
        raise path moved again, to 's non-antichain guard.
        """
        ped = self._refusing_pedigree()
        with self.assertRaises(pyp_errors.PyPedalError) as ctx:
            pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertIn("antichain", str(ctx.exception))

    def test_pedigree_survives_a_refusal(self):
        """
        The refusal path must leave the pedigree untouched too -- a routine that
        validated late could already have modified state before raising.
        """
        ped = self._refusing_pedigree()
        before = self._parentage(ped)
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertEqual(before, self._parentage(ped))

    def test_pedigree_survives_the_newly_supported_half_founder_path(self):
        """
         completes the pedigree before running Appendix B. That
        completion is analysis-local, so it must be invisible here: no phantom
        joins the pedigree, and no parental slot is filled in on the caller's
        animals.
        """
        for name in self.ROUTINES:
            with self.subTest(routine=name):
                ped = self._half_founder_pedigree()
                before = self._parentage(ped)
                getattr(pyp_metrics, name)(ped)
                self.assertEqual(before, self._parentage(ped))
                self.assertEqual(6, len(ped.pedigree))

    def test_repeated_calls_are_reproducible(self):
        """
        The sharpest consequence of the defect: calling twice gave a different
        answer the second time, because the first call had eaten the pedigree.
        """
        for name in ("a_effective_ancestors_definite",
                     "a_effective_ancestors_indefinite"):
            with self.subTest(routine=name):
                ped = load_corpus("boichard2a.ped", "asdg")
                routine = getattr(pyp_metrics, name)

                def call():
                    try:
                        return ("value", routine(ped))
                    except pyp_errors.PyPedalError as exc:
                        return (type(exc).__name__, str(exc))

                self.assertEqual(call(), call())

    def test_subsequent_analysis_is_unaffected(self):
        """A COI computed after these routines must match one computed before."""
        ped = load_corpus("boichard2a.ped", "asdg")
        before = coi(ped, "tabular")
        for name in self.ROUTINES:
            try:
                getattr(pyp_metrics, name)(ped)
            except pyp_errors.PyPedalError:
                pass
        after = coi(ped, "tabular")
        self.assertEqual(before, after)


class TestBoichardFoundersGenerationArgument(unittest.TestCase):
    """
    Stage A4, step 1. ``a_effective_founders_boichard(gen=...)`` selects the
    reference population, and the docstring has always documented it that way.
    It was assigned and then ignored -- the loop hardcoded the most recent
    generation -- so a caller asking about generation 2 silently received an
    analysis of generation 3.

    This covers the ARGUMENT only. The one-parent-known branch of the same
    routine is  and stays untouched; ``generation_split.ped`` deliberately
    gives every animal two known parents so that branch is not exercised here.
    """

    # gen 2 = {5,6}, both out of 1 x 2: founders 1,2 contribute 0.5 each and
    #   3,4 contribute nothing.  sum q^2 = 0.5  ->  f_e = 2
    # gen 3 = {7,8}: 7 out of 5 x 6, 8 out of 3 x 4, so all four contribute
    #   0.25.        sum q^2 = 0.25 ->  f_e = 4
    HAND_DERIVED = {"2": 2.0, "3": 4.0}

    def _load(self):
        return load_corpus_from_path(corpus("generation_split.ped"), "asdg")

    def test_gen_selects_the_reference_population(self):
        for generation, expected in self.HAND_DERIVED.items():
            with self.subTest(gen=generation):
                got = pyp_metrics.a_effective_founders_boichard(
                    self._load(), gen=int(generation))
                self.assertAlmostEqual(expected, got, places=12)

    def test_default_is_the_most_recent_generation_and_is_unchanged(self):
        """The default path must not move: gen=None still means gens[0]."""
        self.assertAlmostEqual(
            self.HAND_DERIVED["3"],
            pyp_metrics.a_effective_founders_boichard(self._load()),
            places=12)

    def test_int_and_str_generations_are_equivalent(self):
        """
        ``NewAnimal.gen`` holds a string while the signature annotates gen as
        Optional[int]. Comparing 3 to '3' would select no animals at all, so the
        comparison is string-normalised and both spellings must work.
        """
        self.assertEqual(
            pyp_metrics.a_effective_founders_boichard(self._load(), gen=2),
            pyp_metrics.a_effective_founders_boichard(self._load(), gen="2"))

    def test_an_absent_generation_raises_rather_than_analysing_something_else(self):
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_metrics.a_effective_founders_boichard(self._load(), gen=99)


class TestGenerationIDsAreOrderedNumerically(unittest.TestCase):
    """
    Audit . ``NewAnimal.gen`` holds a string, and the Boichard
    routines chose their default reference population with a lexicographic
    sort, in which ``'9' > '12'``. On any pedigree with ten or more generations
    that silently analysed **generation 9**.

    The documented contract is unambiguous. ``pyp-metrics.tex:14-38`` says, for
    all three Boichard routines: *"By default the most recent generation -- the
    generation with the largest generation ID -- will be used as the reference
    population."*

    ``deep_generations.ped`` has 12 generations, which is the smallest count
    that exposes the defect. Every animal in it has two known parents or none,
    so the disputed  branch is deliberately not exercised.
    """

    def _deep(self):
        return load_corpus_from_path(corpus("deep_generations.ped"), "asdg")

    def test_the_fixture_actually_exposes_the_defect(self):
        """
        Guard on the guard: if the fixture ever drops below ten generations, or
        its generation IDs stop being ambiguous under a lexicographic sort, the
        tests below stop testing anything.
        """
        gens = {str(a.gen) for a in self._deep().pedigree}
        self.assertGreaterEqual(len(gens), 10)
        self.assertNotEqual(
            max(gens),                       # lexicographic
            max(gens, key=float),            # numeric
            "fixture no longer distinguishes lexicographic from numeric order")

    def test_default_reference_population_is_the_numerically_largest(self):
        default = pyp_metrics.a_effective_founders_boichard(self._deep())
        numeric_max = pyp_metrics.a_effective_founders_boichard(self._deep(), gen=12)
        lexicographic_max = pyp_metrics.a_effective_founders_boichard(self._deep(), gen=9)

        self.assertAlmostEqual(numeric_max, default, places=12)
        self.assertNotAlmostEqual(
            lexicographic_max, default, places=9,
            msg="generation 9 and generation 12 give the same answer here, so "
                "this fixture cannot detect the defect")

    def test_helper_returns_the_original_value_not_the_parsed_number(self):
        """
        Downstream code compares against ``individual.gen`` directly, so the
        helper must hand back the value as it was stored.
        """
        got = pyp_metrics._most_recent_generation(["1", "9", "12", "10"], "t")
        self.assertEqual("12", got)

    def test_non_numeric_generations_raise_rather_than_guessing_an_order(self):
        """
        "Largest generation ID" has no meaning for labels that are not ordered
        quantities. Falling back to lexicographic or insertion order would be
        inventing a contract PyPedal has never documented.
        """
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            pyp_metrics._most_recent_generation(["1", "2", "spring"], "t")
        self.assertIn("not numeric", str(caught.exception))

    def test_numeric_strings_and_numbers_are_both_accepted(self):
        self.assertEqual(3, pyp_metrics._most_recent_generation([1, 2, 3], "t"))
        self.assertEqual("3.0", pyp_metrics._most_recent_generation(["3.0", "1"], "t"))


class TestIndefiniteBoundsAreReal(unittest.TestCase):
    """
    , audit  -- **the refusal has been retired and replaced**,
    which is what the class it supersedes said had to happen.

    ``a_effective_ancestors_indefinite`` claimed to return approximate lower and
    upper bounds for f_a. It did not: ``f_l`` and ``f_u`` were one expression
    evaluated twice and assigned to two names, so no input could make them
    differ, and both fell below 1 on valid input (0.6514 on ``boichard2a.ped``),
    impossible for the reciprocal of a sum of squared probabilities. PyPedal
    2.0.4 returns ``NaN`` for the lower bound on the same pedigree.

    an earlier revision stage A3 replaced that with an explicit refusal, and recorded that
    the assertions "must be replaced -- not deleted -- by bound-validity
    assertions when the reference arrives". Boichard pp.9-10 have now been read,
    so that is what these are.

    The defining property is that the two bounds BRACKET the exact Appendix B
    value, which is only meaningful because both are computed from the same
    marginal-contribution sequence the exact routine consumes.
    """

    PEDIGREES = (("boichard2a.ped", "asdg"), ("boichard_fig1.ped", "asdg"),
                 ("boichard_fig2.ped", "asdg"))

    def test_the_bounds_bracket_the_exact_value(self):
        for name, pedformat in self.PEDIGREES:
            for n in (1, 2, 3, 25):
                with self.subTest(pedigree=name, n=n):
                    exact = pyp_metrics.a_effective_ancestors_definite(
                        load_corpus(name, pedformat))
                    f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(name, pedformat), n=n)
                    self.assertLessEqual(f_l, exact + 1e-9)
                    self.assertLessEqual(exact, f_u + 1e-9)

    def test_the_two_bounds_are_not_the_same_expression_twice(self):
        """
        The original defect in one assertion: with only the first ancestor
        taken there is real unexplained mass, so the bounds MUST differ. If they
        ever coincide there again, the two formulas have collapsed back into
        one.
        """
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            load_corpus("boichard_fig1.ped", "asdg"), n=1)
        self.assertNotAlmostEqual(f_l, f_u, places=6)
        self.assertLess(f_l, f_u)

    def test_bounds_are_never_below_one(self):
        for name, pedformat in self.PEDIGREES:
            with self.subTest(pedigree=name):
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_corpus(name, pedformat))
                self.assertGreaterEqual(f_l, 1.0)
                self.assertGreaterEqual(f_u, 1.0)

    def test_taking_every_ancestor_gives_the_exact_value(self):
        """
        The endpoint. Once no contribution is unexplained the truncation IS the
        exact answer, and neither residual term may be evaluated -- (f - n) is
        zero there. Reaching it without a ZeroDivisionError is the safeguard
        working.
        """
        for name, pedformat in self.PEDIGREES:
            with self.subTest(pedigree=name):
                exact = pyp_metrics.a_effective_ancestors_definite(
                    load_corpus(name, pedformat))
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_corpus(name, pedformat), n=10 ** 6)
                self.assertAlmostEqual(exact, f_l, places=9)
                self.assertAlmostEqual(exact, f_u, places=9)

    @staticmethod
    def _fixture(text):
        tmp = owned_temp_dir(prefix="pypedal_bl2_")
        path = os.path.join(tmp, "fixture.ped")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path

    def test_it_agrees_with_the_definite_routine_about_what_it_answers(self):
        """
        WAS test_it_refuses_exactly_when_the_definite_routine_refuses.

        The bounds inherit 's R1/R2/R3 because they consume the same
        engine, so the two routines must agree about what they will and will
        not answer. When this was written the only agreement available was a
        shared refusal on a half-founder pedigree;  is now implemented,
        so the same pedigree exercises the agreement from the ACCEPTING side --
        which is the stronger direction, because a bound computed on semantics
        the exact routine declined is the invisibly-wrong answer the refusal
        existed to prevent, and so is a refusal on input it accepts.

         still refuses, so both directions are covered.
        """
        half_founder = self._fixture(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        non_antichain = self._fixture(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 5 4 2\n")

        # Accepted by both, now that R2 is resolved.
        exact = pyp_metrics.a_effective_ancestors_definite(
            load_corpus_from_path(half_founder, "asdg"))
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            load_corpus_from_path(half_founder, "asdg"), n=10 ** 6)
        self.assertAlmostEqual(exact, f_l, places=9)
        self.assertAlmostEqual(exact, f_u, places=9)

        # Refused by both, for the same reason, through the same guard.
        with self.assertRaises(pyp_errors.PyPedalError) as definite:
            pyp_metrics.a_effective_ancestors_definite(
                load_corpus_from_path(non_antichain, "asdg"))
        with self.assertRaises(pyp_errors.PyPedalError) as bounded:
            pyp_metrics.a_effective_ancestors_indefinite(
                load_corpus_from_path(non_antichain, "asdg"))
        self.assertIn("antichain", str(definite.exception))
        self.assertIn("antichain", str(bounded.exception))

    def test_the_bounds_bracket_the_exact_value_on_a_half_founder_pedigree(self):
        """
         must inherit R2 through the shared generator rather than acquiring
        a half-founder interpretation of its own. If the bounds ever stopped
        bracketing here, they would be truncating a different sequence from the
        one the exact routine runs to exhaustion.
        """
        path = self._fixture(
            "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n")
        exact = pyp_metrics.a_effective_ancestors_definite(
            load_corpus_from_path(path, "asdg"))
        for n in (1, 2):
            with self.subTest(n=n):
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_corpus_from_path(path, "asdg"), n=n)
                self.assertLessEqual(f_l, exact + 1e-9)
                self.assertLessEqual(exact, f_u + 1e-9)

    def test_n_must_be_positive(self):
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            pyp_metrics.a_effective_ancestors_indefinite(
                load_corpus("boichard2a.ped", "asdg"), n=0)

    def test_the_callers_pedigree_survives(self):
        """
        The routine used to clear the parentage of 7 of 14 animals on this
        pedigree (). Still asserted now that it returns rather than
        raises -- the success path is the one that does the deleting.
        """
        ped = load_corpus("boichard2a.ped", "asdg")
        before = [(int(a.animalID), int(a.sireID), int(a.damID))
                  for a in ped.pedigree]
        pyp_metrics.a_effective_ancestors_indefinite(ped)
        after = [(int(a.animalID), int(a.sireID), int(a.damID))
                 for a in ped.pedigree]
        self.assertEqual(before, after)


class TestGeneDropAncestralInbreeding(unittest.TestCase):
    """
    , repaired in Commit 5.

    ``dropped_ancestral_inbreeding`` accumulated one estimate per replicate and
    never divided by ``rounds``, so it returned ``rounds`` times the
    coefficient -- up to 100 at the default ``rounds=100``, for a quantity
    bounded in [0, 1].

    Note carefully what is NOT asserted: that the result is identical at
    different ``rounds``. This is a Monte Carlo estimator, so that would be an
    invalid invariant, and the audit retracted it as such (R6). What must hold
    is that the value is a probability and does not scale with the replicate
    count. Both assertions are exact at a fixed seed, so this suite stays
    deterministic; the variance-and-convergence experiment, which is
    inherently statistical, lives in ``the independent oracle``.

    **Fixture: suwanlee_fig1.ped, not hartlandclark.ped.** suwanlee_fig1.ped is
    the published 23-animal pedigree and carries non-trivial ancestral
    inbreeding (to 0.822), so it exercises more of the estimator than the old
    fixture did.

    It was originally chosen because hartlandclark.ped could not be evaluated at
    all: four of its fifteen animals are half-founders, and the estimator refused
    that input under b/HF. **That refusal is gone** -- Baumung et al. (2015)
    p.102 supplies the rule, and hartlandclark.ped now computes. The fixture
    choice stands on its own merits regardless, and the half-founder numerics
    live in tests/test_half_founder_gene_drop.py against hand-derived analytic
    expectations rather than here.
    """

    def test_ancestral_inbreeding_is_a_probability(self):
        for rounds in (1, 5, 20):
            with self.subTest(rounds=rounds):
                ped = load_corpus("suwanlee_fig1.ped", "asd")
                result = pyp_metrics.dropped_ancestral_inbreeding(
                    ped, rounds=rounds, loci=20, seed=42)
                self.assertTrue(result)
                for animal_id, value in result.items():
                    self.assertGreaterEqual(value, 0.0)
                    self.assertLessEqual(
                        value, 1.0,
                        "F_a(%s) = %r at rounds=%d" % (animal_id, value, rounds))

    def test_result_does_not_scale_with_rounds(self):
        """
        The defect multiplied the estimate by `rounds`. A correct estimator
        differs between these two only by sampling noise.
        """
        def mean_at(rounds):
            ped = load_corpus("suwanlee_fig1.ped", "asd")
            result = pyp_metrics.dropped_ancestral_inbreeding(
                ped, rounds=rounds, loci=20, seed=42)
            return sum(result.values()) / len(result)

        one, twenty = mean_at(1), mean_at(20)
        self.assertGreater(one, 0.0)
        self.assertLess(twenty, 4.0 * one)
        self.assertGreater(twenty, 0.25 * one)

    def test_it_leaves_no_scratch_state_on_the_callers_animals(self):
        """
        Audit . The simulation used to write an undeclared
        ``ancestor_alleles`` attribute onto every NewAnimal in the caller's
        pedigree and leave the last replicate's state there after returning.

        ``ancestor_alleles`` is not a NewAnimal field -- it appears nowhere in
        ``pyp_newclasses`` -- so this was scratch storage borrowed from the
        caller's objects, not a result being published.
        """
        ped = load_corpus("suwanlee_fig1.ped", "asd")
        pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5, seed=42)
        leaked = [a.animalID for a in ped.pedigree
                  if hasattr(a, "ancestor_alleles")]
        self.assertEqual([], leaked)

    def test_repeated_calls_give_the_same_answer(self):
        """
        The consequence that matters: with scratch state left behind, a second
        call started from a mutated pedigree. Deterministic at a fixed seed, so
        any difference here is state leaking between calls.
        """
        ped = load_corpus("suwanlee_fig1.ped", "asd")
        first = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=3, loci=10, seed=7)
        second = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=3, loci=10, seed=7)
        self.assertEqual(first, second)

    def test_variance_falls_across_independent_seeds(self):
        """
        Stage A5b. The Monte Carlo property this estimator must have, and the
        one that distinguishes a genuine sampling estimator from an arbitrary
        function of ``rounds``: spread across independent seeds shrinks as the
        replicate count grows.

        Deliberately NOT asserted: strict monotonicity at every intermediate
        replicate count. Sample standard deviation over eight seeds is itself a
        random quantity, and requiring it to fall at each step would be flaky by
        construction (audit R6 makes the same point about a different invalid
        invariant). Instead: a FIXED seed set, so the test is deterministic and
        cannot flake at all, and only two WELL-SEPARATED counts.

        Measured behaviour on suwanlee_fig1.ped, loci=20, these eight seeds --
        sd 0.014432 (rounds=2), 0.008992 (10), 0.004174 (50), 0.002000 (200),
        which tracks the 1/sqrt(rounds) a mean-of-replicates estimator should
        show. The assertion below demands only a factor of 2 between rounds=2
        and rounds=50, where the measured factor is 3.5.

        The threshold is a REGRESSION AND CHARACTERISATION criterion, not a
        mathematical theorem, and it is deliberately not re-tuned to whatever
        the estimator happens to produce after a mechanism change. If a repaired
        estimator fails it, that is a finding to report rather than a number to
        adjust.

        This is a property of the gene-drop estimator alone. It says nothing
        about Ballou's recursion, which is deterministic and has no sampling
        variance -- see  on why the two are not compared.
        """
        import statistics

        seeds = (11, 23, 37, 51, 67, 83, 97, 113)

        def spread(rounds):
            means = []
            for seed in seeds:
                ped = load_corpus("suwanlee_fig1.ped", "asd")
                result = pyp_metrics.dropped_ancestral_inbreeding(
                    ped, rounds=rounds, loci=20, seed=seed)
                means.append(sum(result.values()) / len(result))
            return statistics.stdev(means)

        few, many = spread(2), spread(50)
        self.assertGreater(few, 0.0, "no spread at all suggests the seed is ignored")
        self.assertLess(
            many, few / 2.0,
            "spread across seeds did not fall between rounds=2 (sd=%.6f) and "
            "rounds=50 (sd=%.6f)" % (few, many))


class TestMendelianTransmission(unittest.TestCase):
    """
    , repaired in Commit 6.

    The gene-dropping simulation chose which parental allele to transmit with
    ``np.random.rand() < frequency``, where ``frequency`` was documented as a
    minor allele frequency and defaulted to 0.05. Mendelian segregation is 0.5
    by definition, and in this routine an allele frequency has no role at all:
    founder alleles are unique origin labels, so this is IBD gene dropping,
    which tracks descent rather than allelic state.

    This is a DOMAIN invariant, not an output-range postcondition -- audit
     as distinct from . Both 0.05 and 0.5 leave F_a inside
    [0, 1], so no range check can tell them apart. Only the transmission rate
    itself can, which is why the statistical test below exists and why it is
    kept separate from the F_a range test above.
    """

    def test_observed_transmission_rate_is_one_half(self):
        """
        The draw is exercised directly at the constant the routine now uses.
        Deterministic at a fixed seed: with 200,000 draws the tolerance below
        is roughly 10 binomial standard deviations, so this cannot flake.
        """
        import numpy as np

        self.assertEqual(0.5, pyp_metrics.MENDELIAN_TRANSMISSION_P)
        np.random.seed(20260817)
        draws = 200000
        hits = int(np.sum(np.random.rand(draws) < pyp_metrics.MENDELIAN_TRANSMISSION_P))
        self.assertAlmostEqual(0.5, hits / float(draws), places=2)

    def test_result_is_independent_of_frequency(self):
        import warnings

        def run(**kwargs):
            ped = load_corpus("suwanlee_fig1.ped", "asd")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                return pyp_metrics.dropped_ancestral_inbreeding(
                    ped, rounds=5, loci=20, seed=42, **kwargs)

        baseline = run()
        for frequency in (0.05, 0.5, 1.0):
            with self.subTest(frequency=frequency):
                self.assertEqual(baseline, run(frequency=frequency))

    def test_passing_frequency_warns(self):
        import warnings

        ped = load_corpus("suwanlee_fig1.ped", "asd")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pyp_metrics.dropped_ancestral_inbreeding(
                ped, rounds=1, loci=5, frequency=0.05, seed=42)
        self.assertTrue(any(issubclass(w.category, DeprecationWarning)
                            for w in caught))

    def test_not_passing_frequency_is_silent(self):
        """The deprecation must not become background noise for normal calls."""
        import warnings

        ped = load_corpus("suwanlee_fig1.ped", "asd")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pyp_metrics.dropped_ancestral_inbreeding(
                ped, rounds=1, loci=5, seed=42)
        self.assertEqual([], [w for w in caught
                              if issubclass(w.category, DeprecationWarning)])


class TestReorderKeywordBinding(unittest.TestCase):
    """
    , repaired in Commit 7.

    ``load()`` called ``reorder(self.pedigree, self.kw['missing_parent'],
    self.kw['reorder_max_rounds'])`` positionally, but the signature is
    ``reorder(myped, filetag, io, missingparent, debug, max_rounds)``. So
    ``missing_parent`` was bound to ``filetag`` and ``reorder_max_rounds`` to
    ``io``, while the parameters actually intended kept their defaults.

    At the defaults the two happen to coincide, which is why this was latent.
    It becomes a silent corruption the moment a pedigree uses a non-zero
    missing-parent sentinel: reorder would look for parent ID 0 while the data
    says -999, and would not recognise founders.
    """

    def test_binding_is_by_keyword(self):
        import inspect
        signature = inspect.signature(pyp_utils.reorder)
        parameters = list(signature.parameters)
        # The defect is only interesting because these are NOT adjacent to
        # myped; assert the shape that made positional binding wrong.
        self.assertEqual("filetag", parameters[1])
        self.assertEqual("io", parameters[2])
        self.assertEqual("missingparent", parameters[3])
        self.assertEqual("max_rounds", parameters[5])

    def test_non_default_missing_parent_is_honoured(self):
        """
        The regression case: a non-zero sentinel on the reorder path, with a
        non-default round limit. Founders must be identified against the
        configured sentinel, so parents still precede offspring afterwards.
        """
        rows = ["1 -999 -999", "2 -999 -999", "3 1 2", "4 1 2", "5 3 4"]
        tmp = owned_temp_dir(prefix="pypedal_sentinel_")
        path = os.path.join(tmp, "sentinel.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(reversed(rows)) + "\n")   # offspring first

        ped = load_corpus_from_path(
            path, "asd",
            missing_parent=-999, reorder=True, renumber=False,
            reorder_max_rounds=7, pedigree_is_renumbered=False)

        self.assertEqual(5, len(ped.pedigree))
        position = {int(a.animalID): i for i, a in enumerate(ped.pedigree)}
        for index, animal in enumerate(ped.pedigree):
            for parent in (int(animal.sireID), int(animal.damID)):
                if parent != -999 and parent in position:
                    self.assertLess(position[parent], index,
                                    "offspring precedes parent after reorder")


class TestFailuresAreLoud(unittest.TestCase):
    """
     (with Findings 1 and 14), repaired in Commit 8.

    The analysis layer used blanket ``except Exception`` handlers that returned
    a sentinel or a partial result. Three of them turned a failure into a
    confident-looking wrong answer:

      * ``preprocess()`` returned an EMPTY pedigree on any parse error -- and
        every caller ignores its return value, so the flag was never consulted;
      * ``recurse_pedigree()`` caught ``RecursionError`` and returned the
        ancestors found so far, silently truncating the ancestor list on a
        pedigree deeper than the recursion limit;
      * ``effective_founders_lacy()`` skipped a founder whose contribution
        failed, computing f_e from an incomplete contribution vector.

    Scope: exactly these handlers plus the two sentinel returns. This is not a
    sweep of every ``except Exception`` in the package.
    """

    def test_malformed_pedigree_raises_instead_of_loading_empty(self):
        """A pedigree that cannot be parsed must not present as an empty one."""
        from PyPedal.pyp_errors import PyPedalError

        tmp = owned_temp_dir(prefix="pypedal_bad_")
        path = os.path.join(tmp, "malformed.ped")
        with open(path, "w", encoding="utf-8") as handle:
            # 'asd' promises integer IDs; these are not integers, and the
            # pedformat has no uppercase codes to license strings.
            handle.write("alpha beta gamma\ndelta beta gamma\n")

        with self.assertRaises(PyPedalError):
            load_corpus_from_path(path, "asd")

    def test_recursion_limit_is_not_swallowed(self):
        """
        A pedigree deeper than the recursion limit must fail, not return a
        truncated ancestor list. Uses a chain longer than the limit rather than
        lowering the limit globally, which would destabilise other tests.
        """

        depth = sys.getrecursionlimit() + 200
        rows = ["1 0 0"] + ["%d %d 0" % (i, i - 1) for i in range(2, depth + 1)]
        tmp = owned_temp_dir(prefix="pypedal_deep_")
        path = os.path.join(tmp, "deep.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(rows) + "\n")

        ped = load_corpus_from_path(path, "asd")
        self.assertEqual(depth, len(ped.pedigree))
        with self.assertRaises(RecursionError):
            pyp_nrm.recurse_pedigree(ped, depth, [])

    def test_shallow_pedigree_still_recurses_normally(self):
        """The guard must not make ordinary traversal fail."""
        ped = load_corpus("mrode.ped", "asd")
        ancestors = pyp_nrm.recurse_pedigree(ped, 6, [])
        self.assertGreater(len(ancestors), 1)


class TestPedigreeInbreedingIsAProbability(unittest.TestCase):
    """
    F is the probability that an animal's two alleles at a locus are identical
    by descent, so 0 <= F <= 1. Asserted for pedigree-based COI only -- genomic
    inbreeding estimators are differently defined and may be negative.
    """

    def test_coi_lies_in_unit_interval(self):
        for name in ("mrode.ped", "new_lacy.ped", "hartlandclark.ped",
                     "generations.ped", "boichard2a.ped"):
            for method in ("tabular", "vanraden", "meu_luo"):
                with self.subTest(pedigree=name, method=method):
                    ped = load_corpus(name)
                    for animal_id, value in coi(ped, method).items():
                        self.assertGreaterEqual(value, 0.0,
                                                "F(%d) < 0 in %s" % (animal_id, name))
                        self.assertLessEqual(value, 1.0,
                                             "F(%d) > 1 in %s" % (animal_id, name))


if __name__ == "__main__":
    unittest.main()
