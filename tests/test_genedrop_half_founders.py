"""
 -- the half-founder path in ``dropped_ancestral_inbreeding``.

 had two outcomes, and conflating them would lose one of them. Both
are now closed, but they closed SEPARATELY and for different reasons, which is
the thing this file exists to keep straight.

**The corruption mechanism is CLOSED (unchanged).** A missing parent could reach
``pedigree[int(sireID) - 1]``, which with the sentinel ``missing_parent == 0``
is ``pedigree[-1]`` -- the last animal in the pedigree, adopted silently as a
parent with no IndexError. Separately, allele slots were strings initialised to
``""``, and ``"" == ""`` is True, so any animal whose pair had not been assigned
yet read as homozygous by descent at every locus. Both are repaired, and both
repairs are general rather than specific to any fixture. **Nothing here weakens
them**, and they do not depend on which inputs are accepted.

**The scientific treatment is now SUPPORTED (this is what changed).** Suwanlee
et al. (2007) p.490 is still silent on the gene source of a missing parental
side -- that has not changed and is not revised. The rule comes from **Baumung
et al. (2015) p.102**, source-explicit: one dummy founder per unknown parental
slot, given two unique alleles and dropped through under ordinary Mendelian
segregation. So ``b/HF`` is closed and the typed refusal is gone.

The accepted domain is now all three parentage cases:

    0 known parents  -> founder, allowed
    1 known parent   -> half-founder, allowed (dummy founder, GRAIN p.102)
    2 known parents  -> descendant, allowed

Stating all three still matters, for the same reason it did when one was
refused: the obvious wrong implementation -- "require two known parents" --
would reject every founder in every pedigree, and a truthiness or ``in`` test on
the parent slots would not distinguish the three cases at all.

The scientific content of the half-founder case lives in
tests/test_half_founder_gene_drop.py, against hand-derived analytic fixtures.
What lives HERE is the structural half: that no sentinel can index a real
animal, that no placeholder can read as autozygous, and that the values are the
ones the rule gives rather than the ones the defect gave.

On the value 0.0 for mrode animal 4: it is derived in the fixtures file from the
GRAIN rule, afresh. An early revision of
the algorithm notes asserted 0.0 on reasoning
that presumed phantom semantics no approved source then established, and that
claim was withdrawn. The withdrawn hypothesis is not evidence and is not used as
one; that the two coincide is a result, not a premise.
"""
import copy
import unittest

import numpy as np

from PyPedal import pyp_errors, pyp_metrics

from _pedhelpers import load_corpus, load_corpus_from_path

import os
import tempfile


def _write_pedigree(text, pedformat="asd"):
    """A purpose-built pedigree, loaded with its output confined to a tmpdir."""
    tmp = tempfile.mkdtemp(prefix="pypedal_f30_")
    path = os.path.join(tmp, "fixture.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return load_corpus_from_path(path, pedformat)


class TestTheAcceptedDomain(unittest.TestCase):
    """
    The domain, pinned case by case. All three parentage cases are legitimate
    and must keep working.
    """

    def test_zero_known_parents_is_an_allowed_founder(self):
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n")
        founders = [a.animalID for a in ped.pedigree if a.founder == "y"]
        self.assertEqual([1, 2], sorted(int(a) for a in founders))
        got = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=2, loci=5, seed=11)
        self.assertEqual(0.0, float(got[1]))
        self.assertEqual(0.0, float(got[2]))

    def test_two_known_parents_is_an_allowed_descendant(self):
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=2, loci=5, seed=11)
        self.assertIn(3, got)
        self.assertGreaterEqual(float(got[3]), 0.0)
        self.assertLessEqual(float(got[3]), 1.0)

    def test_exactly_one_known_parent_is_supported(self):
        """
        The case that used to raise ``PyPedalNotImplementedError`` under
        ``b/HF``. It now computes, because Baumung et al. (2015) p.102
        defines the gene source the refusal was standing in for.
        """
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 3 0\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=2, loci=5, seed=11)
        self.assertEqual(4, len(got))
        for value in got.values():
            self.assertGreaterEqual(float(value), 0.0)
            self.assertLessEqual(float(value), 1.0)


class TestTheHalfFounderIsComputedNotRefused(unittest.TestCase):
    """
    Both parental orientations, and the properties that follow directly from
    the dummy-founder rule.

    The defect was symmetric -- either slot could hold the sentinel -- so the
    support has to be too.
    """

    def test_a_missing_dam_is_computed(self):
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 3 0\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5)
        self.assertIn(4, got)
        self.assertEqual(0.0, float(got[4]))

    def test_a_missing_sire_is_computed(self):
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 0 3\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5)
        self.assertIn(4, got)
        self.assertEqual(0.0, float(got[4]))

    def test_the_result_contains_no_dummy_founders(self):
        """
        The dummy founder is simulation-local. GRAIN emits dummy founders in its
        own output files (p.104); PyPedal does not, because the public return
        contract is a software-design question and not part of the rule. The
        mapping is keyed on the pedigree's own animals and nothing else.
        """
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 0\n4 2 0\n5 3 4\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5)
        self.assertEqual({int(a.animalID) for a in ped.pedigree},
                         {int(k) for k in got})

    def test_neither_half_founder_refusal_survives_but_they_closed_separately(self):
        """
        WAS test_boichards_half_founder_refusal_is_untouched, which asserted
        that b/HF closing left 's refusal in place. Both are now
        closed, and the guard inverts rather than being deleted -- its subject
        was never "one of them is still refusing", it was "these two are not
        the same question and must not be closed as though they were".

        They were not:

          b/HF  Baumung et al. (2015) p.102, SOURCE-EXPLICIT. The paper
                    says in so many words that a dummy founder is created for
                    each half-founder. Closed 2026-08-21.

             Boichard (1997) Appendix B, SILENT. Appendix A step 4 is
                    source-explicit, but carrying it in as pedigree completion
                    ahead of an unchanged Appendix B is MATHEMATICALLY IMPLIED
                    / INDEPENDENTLY SUPPORTED, selected against the paper's own
                    p.7, p.8 and p.17 invariants. NOT source-explicit
                    Appendix-B text.

        Different papers, different evidence classes, different phases. The
        mechanism they share -- give an unknown parental side its own founder
        -- is a coincidence of shape, not a shared justification, and the two
        implementations deliberately share no code: gene dropping mints an
        anonymous pair of allele arrays with no ID at all, while Appendix B
        completes the pedigree with numbered phantoms it must be able to
        select as ancestors.
        """
        self.assertFalse(hasattr(pyp_metrics, "_boichard_require_no_half_founders"))
        self.assertFalse(hasattr(pyp_metrics, "_suwanlee_require_no_half_founders"))

    def test_the_two_half_founder_rules_share_no_implementation(self):
        """
        The other half of the same guard. A later "let us not duplicate this"
        refactor merging the two would silently make one paper's rule stand in
        for the other's.
        """
        with open(os.path.abspath(pyp_metrics.__file__), encoding="utf-8") as fh:
            source = fh.read()
        for shared in ("def _half_founder_completion(",
                       "def _mint_dummy(",
                       "def _shared_phantom_rule("):
            self.assertNotIn(shared, source)
        # Gene dropping must not reach for the Appendix-B completion layer.
        drop = source[source.index("def dropped_ancestral_inbreeding"):]
        self.assertNotIn("_boichard_completed_arrays", drop)


class TestTheCorpusFixtures(unittest.TestCase):
    """
    The two corpus pedigrees that carry half-founders. These were the
    EXPECTED_REFUSAL fixtures and are now supported-domain fixtures.
    """

    HALF_FOUNDER_PEDIGREES = (("mrode.ped", "asd", [4]),
                              ("hartlandclark.ped", "asdb", [8, 9, 10, 11]))

    def test_the_half_founders_are_where_we_think_they_are(self):
        """
        The premise of these fixtures, checked rather than assumed. If a corpus
        pedigree changed, the tests below would still pass while testing
        nothing.
        """
        for name, pedformat, expected in self.HALF_FOUNDER_PEDIGREES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                missing = str(ped.kw["missing_parent"])
                half = [int(a.animalID) for a in ped.pedigree
                        if (str(a.sireID) == missing) != (str(a.damID) == missing)]
                self.assertEqual(expected, half)

    def test_both_are_computed(self):
        for name, pedformat, _expected in self.HALF_FOUNDER_PEDIGREES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                got = pyp_metrics.dropped_ancestral_inbreeding(
                    ped, rounds=5, loci=7, seed=11)
                self.assertEqual(len(ped.pedigree), len(got))
                for value in got.values():
                    self.assertGreaterEqual(float(value), 0.0)
                    self.assertLessEqual(float(value), 1.0)

    def test_mrode_animal_four_is_zero_and_not_one(self):
        """
        The symptom  was registered on: a deterministic ``f_a = 1.0``
        -- maximum possible ancestral inbreeding -- for an animal reached only
        through ``pedigree[-1]``.

        The half-founder is animal 4 in the FILE, which renumbering moves; the
        assertion is made in original coordinates so it cannot drift onto a
        different animal.

        Under the GRAIN rule its value is 0.0 exactly: its known parent is a
        founder and its unknown side is a dummy founder, so nothing arrives
        flagged, and a founder label can never equal a dummy label so it is
        never autozygous. 1.0 was the defect; 0.0 is the rule's answer.
        """
        ped = load_corpus("mrode.ped", "asd")
        got = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=50, loci=50, seed=11)
        by_original = {int(a.originalID): float(got[a.animalID])
                       for a in ped.pedigree}
        self.assertEqual(0.0, by_original[4])
        self.assertNotEqual(1.0, by_original[4])

    def test_ballou_still_computes_on_the_same_pedigrees(self):
        """
        Ballou's recursion defines a half-founder's unknown side explicitly --
        f_s = f_as = 0 -- so it was unaffected by the refusal and is unaffected
        by its removal. The two remain different estimators and must not be
        repaired, or compared, as one.
        """
        for name, pedformat, _expected in self.HALF_FOUNDER_PEDIGREES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                got = pyp_metrics.ballou_ancestral_inbreeding(ped)
                self.assertEqual(len(ped.pedigree), len(got))
                for value in got.values():
                    self.assertGreaterEqual(float(value), 0.0)
                    self.assertLessEqual(float(value), 1.0)


class TestNoSentinelCanIndex(unittest.TestCase):
    """
    The structural half of the repair, which holds regardless of the domain
    policy -- and which now matters MORE than it did, because the half-founder
    path is executed rather than refused.

    The refusal used to narrow the input domain so the defect was unreachable.
    It is reachable now, so these tests are what stands between the sentinel and
    ``pedigree[-1]``.
    """

    def test_the_sentinel_is_never_used_as_an_index(self):
        """
        ``pedigree[int(0) - 1]`` is ``pedigree[-1]``, which is a valid Python
        index and therefore fails silently. Asserted against the pedigree that
        used to trigger it: if any code path still resolved a missing parent
        positionally, animal 4 would adopt the last animal in the pedigree as
        its dam -- which is its own descendant, and unprocessed, so it read as
        autozygous at every locus and returned a deterministic 1.0.
        """
        ped = load_corpus("mrode.ped", "asd")
        missing = int(ped.kw["missing_parent"])
        self.assertEqual(0, missing)
        # The property under test, stated explicitly: the sentinel minus one is
        # a legal index into a non-empty list, which is why this needed an
        # explicit construction rather than an exception.
        self.assertEqual(ped.pedigree[missing - 1], ped.pedigree[-1])

        got = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=50, loci=50, seed=11)
        by_original = {int(a.originalID): float(got[a.animalID])
                       for a in ped.pedigree}
        # Deterministic 1.0 was the signature of the wrap. It is not merely
        # "different now" -- it is exactly zero, at every one of the 2500 loci.
        self.assertEqual(0.0, by_original[4])

    def test_no_last_animal_dependency(self):
        """
        The sharpest form of the wrap test, and the one that does not depend on
        knowing what the right answer is.

        ``pedigree[-1]`` is whatever happens to be last, so under the defect a
        half-founder's value was a function of an unrelated animal at the end of
        the file. Appending an animal that is nobody's ancestor must therefore
        leave every original coefficient EXACTLY unchanged. It does not merely
        change the answer under the defect -- it changes which animal is
        consulted -- so this catches a positional lookup even if the value it
        produced happened to look plausible.
        """
        base = "1 0 0\n2 0 0\n3 1 2\n4 3 0\n"
        appended = base + "5 1 2\n"

        first = pyp_metrics.dropped_ancestral_inbreeding(
            _write_pedigree(base), rounds=20, loci=50, seed=7)
        second = pyp_metrics.dropped_ancestral_inbreeding(
            _write_pedigree(appended), rounds=20, loci=50, seed=7)

        self.assertEqual(4, len(first))
        self.assertEqual(5, len(second))
        for animal in first:
            self.assertEqual(
                float(first[animal]), float(second[animal]),
                "animal %s changed when an unrelated animal was appended: "
                "%r -> %r" % (animal, first[animal], second[animal]))

    def test_no_placeholder_can_read_as_autozygous(self):
        """
        The second half of the original defect: allele slots were strings
        initialised to ``""``, and ``"" == ""`` is True, so any animal not yet
        processed read as homozygous by descent at every locus.

        Asserted structurally -- allele identity is an integer allocated from a
        counter, so there is no value in the domain that means "unset" and could
        compare equal to itself. The pedigree below has two half-founders whose
        dummies must not collide with each other or with anything unset; every
        coefficient is exactly zero, which no placeholder collision permits.
        """
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 0\n4 2 0\n5 3 4\n6 0 0\n7 5 6\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=20, loci=100,
                                                       seed=3)
        for animal, value in got.items():
            self.assertEqual(0.0, float(value), "animal %s" % animal)

    def test_an_unprocessed_known_parent_raises_an_internal_error(self):
        """
        The contract for the one case that is genuinely impossible: a KNOWN
        parent that has not been dropped through yet. That is a broken invariant
        -- the pedigree is not ordered parents-before-offspring -- and it must
        surface as PyPedalInternalError naming both animals, not as a bare
        KeyError leaking the state dict as the contract.

        Critically, this path must NOT be reachable via the sentinel. The
        missing-parent case is an explicit construction placed before this
        check, so a half-founder can never arrive here; if it could, the
        "impossible" branch would be doing routine work.
        """
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 3 0\n")
        # Reverse the pedigree in place, so offspring precede their parents.
        ped.pedigree = list(reversed(ped.pedigree))
        with self.assertRaises(pyp_errors.PyPedalInternalError) as ctx:
            pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5)
        message = str(ctx.exception)
        self.assertIn("has not been dropped through yet", message)

    def test_the_sentinel_never_raises_the_internal_error(self):
        """
        The complement of the test above, and the reason it is worth having
        both: a half-founder must be handled, not diagnosed. If the sentinel
        ever fell through to the state lookup it would raise here rather than
        computing, which is how a fall-through would most likely present.
        """
        ped = _write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 3 0\n5 0 3\n")
        got = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=2, loci=5)
        self.assertEqual(5, len(got))


class TestTheScientificInvariant(unittest.TestCase):
    """
    ``f_a`` is a probability, so it lies in [0, 1] -- asserted on the RAW
    returned values, with no clamping anywhere in the routine.

     is the standing proof that this bound is necessary but not
    sufficient: the corrupted value was exactly 1.0, inside the interval, and
    the postcondition passed while the output was wrong. It is kept because it
    is cheap and true, not because it is evidence of correctness.
    """

    PEDIGREES = (("mrode.ped", "asd"), ("hartlandclark.ped", "asdb"))

    def test_every_coefficient_is_a_probability(self):
        for name, pedformat in self.PEDIGREES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                got = pyp_metrics.dropped_ancestral_inbreeding(
                    ped, rounds=20, loci=50, seed=5)
                for animal, value in got.items():
                    self.assertGreaterEqual(float(value), 0.0, "animal %s" % animal)
                    self.assertLessEqual(float(value), 1.0, "animal %s" % animal)


class TestStateIsolationOnHalfFounders(unittest.TestCase):
    """
    The dummy-founder path is new code, and it is the obvious place for a
    caller's pedigree to start being mutated or for the global RNG to start
    moving. Both are pinned here on half-founder input specifically -- the
    existing suite pins them only where there is no half-founder.
    """

    def _ped(self):
        return load_corpus("hartlandclark.ped", "asdb")

    def test_the_callers_animals_are_untouched(self):
        """
        A DEEP snapshot. A shallow one would compare attribute containers by
        identity and pass even if their contents had been rewritten in place,
        which is exactly how the historical ancestor_alleles mutation hid.
        """
        ped = self._ped()
        before = copy.deepcopy([vars(a).copy() for a in ped.pedigree])
        pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=20, seed=9)
        self.assertEqual(before, [vars(a).copy() for a in ped.pedigree])

    def test_no_undeclared_scratch_attribute_is_attached(self):
        ped = self._ped()
        before = [set(vars(a)) for a in ped.pedigree]
        pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=20, seed=9)
        self.assertEqual(before, [set(vars(a)) for a in ped.pedigree])

    def test_the_options_dictionary_is_untouched(self):
        ped = self._ped()
        before = copy.deepcopy(ped.kw)
        pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=20, seed=9)
        self.assertEqual(before, ped.kw)

    def test_the_process_global_numpy_stream_is_not_disturbed(self):
        """
        The routine builds its own default_rng(seed). Nothing it does may move
        numpy's process-global state, because that would silently change results
        for every other caller in the process.
        """
        ped = self._ped()
        np.random.seed(1234)
        before = np.random.get_state()
        pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=20, seed=9)
        after = np.random.get_state()
        self.assertEqual(before[0], after[0])
        self.assertTrue((before[1] == after[1]).all())
        self.assertEqual(before[2:], after[2:])

    def test_the_same_seed_reproduces_exactly_on_a_half_founder_pedigree(self):
        first = pyp_metrics.dropped_ancestral_inbreeding(
            self._ped(), rounds=5, loci=20, seed=9)
        second = pyp_metrics.dropped_ancestral_inbreeding(
            self._ped(), rounds=5, loci=20, seed=9)
        self.assertEqual({k: float(v) for k, v in first.items()},
                         {k: float(v) for k, v in second.items()})


if __name__ == "__main__":
    unittest.main()
