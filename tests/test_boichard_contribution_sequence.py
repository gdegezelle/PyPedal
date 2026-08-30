"""
 --  must refuse where the exact calculation refuses.

WHAT THIS FILE IS ABOUT
-----------------------
Boichard, Maignel & Verrier (1997), *Genet Sel Evol* 29:5-23. The bounding
scheme on article pp.9-10 is a TRUNCATION of the Appendix-B sequence, not a
second method:

    f_u = 1 / [ sum(p_i^2) + (1-c)^2 / (f - n) ]
    f_l = 1 / [ sum(p_i^2) + m * p_n^2 ],    m = (1 - c) / p_n

Both are built on ``1 - c``, which is the unexplained mass ONLY because
article p.8 says the marginal contributions over all ancestors sum to one. If
the full sequence sums to ``S != 1`` then the unexplained mass is ``S - c`` and
``1 - c`` bounds nothing. Worse, the bracketing claim ``f_l <= f_a <= f_u``
names ``f_a`` -- exactly the quantity the exact routine refuses to produce on
such a sequence. The invariant has no middle term.

The paper's own bracketing claim settled this before either
routine was written:

    "a_effective_ancestors_indefinite must not return bounds for a pedigree on
    which a_effective_ancestors_definite refuses -- a bound built on unresolved
    semantics is exactly the invisibly-wrong answer the  refusal was
    written to prevent."

``a_effective_ancestors_definite`` called ``check_contribution_vector`` and
``a_effective_ancestors_indefinite`` did not. That was , and this
file is the contract the repair satisfies. See
``the algorithm notes``.

WHAT IS *NOT* CLAIMED HERE
--------------------------
Nothing in this file says a reference population is biologically invalid. An
analyst-defined R may be perfectly meaningful while the current Boichard/R3/
 implementation domain cannot produce a valid contribution sequence for
it. R3 is a declared, source-silent PyPedal convention and is neither weakened
nor redefined by this finding. The griffon reference populations used below
are TEST-SPECIFIED characterisation sets; neither is proposed as the correct
scientific griffon cohort.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with every test
asserting the repaired behaviour marked ``xfail(strict=True)``, so the contract
was written down before any code satisfied it and nothing could pass by
accident. Those markers are removed in the commit that implements the fix; a
strict xfail left in place after it starts passing is an error.

``self.subTest`` is deliberately kept out of the xfailed tests: pytest 9
handles the subTest/xfail interaction differently from pytest 8, a difference
this project has already been bitten by
(``docs/-R2-HALF-FOUNDER-VERIFICATION.md`` section 2.4). Plain loops are
used there instead.
"""
import os
import tempfile
import unittest

from _pedhelpers import load_corpus, load_corpus_from_path
from PyPedal import pyp_errors, pyp_metrics

# `griffon_cohort` is a TEST-ONLY helper -- see  section 12 for the
# audit confirming nothing of the sort ships in `PyPedal/`. It is imported
# rather than duplicated so there is exactly one definition of these sets.
from test_boichard_reference_population import griffon_cohort, load_griffon

# ---------------------------------------------------------------------------
# The synthetic reproducer.
#
# Sixteen animals, pedformat `asdg`. Generation 4 is the reference population
# and holds two animals: a descendant of eight founders, and animal 16, which
# is a FOUNDER SITTING INSIDE R.  zeroes reference-population members
# before ancestor selection, so the half of the gene pool animal 16 carries is
# credited to no ancestor at all and the sequence sums to 0.5.
#
# It matters that this fixture carries a `g` column and needs no `reference=`
# argument: it reaches the defect through the LEGACY generation path, which is
# what establishes that  predates  rather than being
# reachable only through the new API. Measured on pre-phase master `ed1fee4`
# as well as here -- see the finding document, section 5(b).
# ---------------------------------------------------------------------------
FOUNDER_IN_REFERENCE = """1 0 0 1
2 0 0 1
3 0 0 1
4 0 0 1
5 0 0 1
6 0 0 1
7 0 0 1
8 0 0 1
9 1 2 2
10 3 4 2
11 5 6 2
12 7 8 2
13 9 10 3
14 11 12 3
15 13 14 4
16 0 0 4
"""

#: Measured, not predicted. Renumbering puts founders first, so original 15 --
#: the only non-founder at generation 4 -- becomes animalID 16.
#:
#: Original 16 is a founder, and was animalID 1 until : the founder
#: pre-pass hoisted each founder to the front of the list in turn, which
#: reversed the founder block, and 16 is the last founder in the file. Founders
#: now keep their input order, so it is animalID 9, after founders 1..8. The
#: reference population is the same two animals either way -- original 15 and
#: 16 -- and every quantity computed from it is unchanged; only the internal
#: numbering moved. See the algorithm notes.
SYNTH_REFERENCE = [9, 16]
SYNTH_SEQUENCE = [(14, 0.25), (15, 0.25)]
SYNTH_SUM_P = 0.5
SYNTH_FOUNDERS = 9

#: Measured on this tree and on `ed1fee4`; identical on both.
GRIFFON_36_SUM_P = 0.3055555555555555
GRIFFON_36_SEQUENCE_LENGTH = 17
GRIFFON_36_FOUNDERS = 130

#: The valid contrast: the nine 1890 animals from known sire and dam.
GRIFFON_9_F_A = 12.461538461538462
GRIFFON_9_BOUNDS = {
    1: (9.0, 54.139896373057),
    2: (9.0, 33.99344262295082),
    3: (9.0, 24.66906474820144),
    5: (12.461538461538462, 18.359020852221214),
    25: (12.461538461538462, 12.461538461538462),
    10 ** 6: (12.461538461538462, 12.461538461538462),
}

#: Published controls. `f_a` and the bounds at every n, measured.
PUBLISHED = {
    "boichard2a.ped": {
        "f_a": 2.0,
        "bounds": {1: (2.0, 3.0), 2: (2.0, 2.0), 3: (2.0, 2.0),
                   25: (2.0, 2.0), 10 ** 6: (2.0, 2.0)},
    },
    "boichard_fig1.ped": {
        "f_a": 4.444444444444445,
        "bounds": {1: (3.333333333333333, 5.319148936170213),
                   2: (3.3333333333333335, 4.545454545454546),
                   3: (4.166666666666667, 4.477611940298508),
                   25: (4.444444444444445, 4.444444444444445),
                   10 ** 6: (4.444444444444445, 4.444444444444445)},
    },
    "boichard_fig2.ped": {
        "f_a": 2.9411764705882346,
        "bounds": {1: (2.5, 4.310344827586206),
                   2: (2.5, 3.03030303030303),
                   3: (2.9411764705882346, 2.999999999999999),
                   25: (2.9411764705882346, 2.9411764705882346),
                   10 ** 6: (2.9411764705882346, 2.9411764705882346)},
    },
}


def load_synthetic(**overrides):
    """Load the founder-in-R fixture from a temporary copy."""
    tmp = tempfile.mkdtemp(prefix="pypedal_test_")
    path = os.path.join(tmp, "f35_founder_in_reference.ped")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(FOUNDER_IN_REFERENCE)
    return load_corpus_from_path(path, "asdg", **overrides)


def griffon_36():
    """The 36-animal 1890 cohort. TEST-SPECIFIED, not an endorsed cohort."""
    return griffon_cohort(load_griffon(), {1890})


def griffon_9():
    """The nine 1890 animals from known sire and dam. Also TEST-SPECIFIED."""
    return griffon_cohort(load_griffon(), {1890}, complete_pedigree_only=True)


# ---------------------------------------------------------------------------
# 1. The measured state of the defect. These pass before AND after the repair;
#    they describe the input, not the behaviour under repair.
# ---------------------------------------------------------------------------

class TestTheInvalidSequenceIsWhatThisFindingSaysItIs(unittest.TestCase):
    """
    Pins the reproducer itself. If any of these drift, every conclusion drawn
    from the fixtures below is drawn from something else.
    """

    def test_the_griffon_cohort_sequence_does_not_sum_to_one(self):
        order = list(pyp_metrics.boichard_marginal_contributions(
            load_griffon(), griffon_36()))
        self.assertEqual(GRIFFON_36_SEQUENCE_LENGTH, len(order))
        self.assertEqual(GRIFFON_36_SUM_P, sum(value for _, value in order))

    def test_the_griffon_cohort_passes_the_r3_antichain_guard(self):
        """
        R3 is not what refuses this cohort, so a reader cannot dismiss the
        finding as an antichain problem already covered elsewhere.
        """
        ped = load_griffon()
        pyp_metrics._boichard_require_antichain(ped, griffon_36(), "finding-35")

    def test_the_mass_is_lost_to_founders_sitting_inside_the_reference(self):
        """
        The cause, measured. Boichard p.7: an animal with one unknown parent is
        a founder too, so the count that matters is 27, not the 25 with both
        parents unknown.
        """
        ped = load_griffon()
        cohort = set(griffon_36())
        missing = int(ped.kw["missing_parent"])
        by_id = {int(a.animalID): a for a in ped.pedigree}
        both = [i for i in cohort
                if int(by_id[i].sireID) == missing
                and int(by_id[i].damID) == missing]
        either = [i for i in cohort
                  if int(by_id[i].sireID) == missing
                  or int(by_id[i].damID) == missing]
        self.assertEqual(36, len(cohort))
        self.assertEqual(25, len(both))
        self.assertEqual(27, len(either), "Boichard p.7 founder definition")

    def test_the_founder_count_is_far_above_the_bounds_it_would_check(self):
        """
        Why `check_effective_number` cannot close the gap: it compares the
        bounds against the founder count, and on this pedigree that ceiling is
        nowhere near the values the bounded routine produces.
        """
        _i, _s, _d, _ph, n_founders = pyp_metrics._boichard_completed_arrays(
            load_griffon())
        self.assertEqual(GRIFFON_36_FOUNDERS, n_founders)

    def test_the_synthetic_fixture_reaches_the_defect_through_the_legacy_path(self):
        """
        No `reference=` anywhere. The reference population is generation 4,
        chosen by the legacy selector, which is what makes this fixture
        evidence that the defect predates .
        """
        ped = load_synthetic()
        reference = sorted(int(a.animalID) for a in ped.pedigree
                           if str(a.gen) == "4")
        self.assertEqual(SYNTH_REFERENCE, reference)
        order = list(pyp_metrics.boichard_marginal_contributions(
            load_synthetic(), reference))
        self.assertEqual(SYNTH_SEQUENCE,
                         [(int(i), float(v)) for i, v in order])
        self.assertEqual(SYNTH_SUM_P, sum(value for _, value in order))
        _i, _s, _d, phantoms, n_founders = \
            pyp_metrics._boichard_completed_arrays(load_synthetic())
        self.assertEqual(SYNTH_FOUNDERS, n_founders)
        self.assertEqual(0, len(phantoms), "no half-founders in this fixture")

    def test_the_exact_routine_already_refuses_both_fixtures(self):
        """The half of the contract that is already honoured."""
        with self.assertRaises(pyp_errors.PyPedalValidationError) as caught:
            pyp_metrics.a_effective_ancestors_definite(
                load_griffon(), reference=list(griffon_36()))
        self.assertIn("marginal contributions", str(caught.exception))
        self.assertIn("not 1", str(caught.exception))

        with self.assertRaises(pyp_errors.PyPedalValidationError) as caught:
            pyp_metrics.a_effective_ancestors_definite(load_synthetic())
        self.assertIn("sums to 0.5, not 1", str(caught.exception))


# ---------------------------------------------------------------------------
# 2. The contract under repair. Every test here is expected to FAIL until the
#    validation call lands, and the markers come off in that same commit.
# ---------------------------------------------------------------------------

class TestTheBoundedRoutineRefusesTheSameInvalidSequence(unittest.TestCase):

    def test_it_refuses_the_griffon_cohort_at_every_n(self):
        cohort = griffon_36()
        for n in (1, 5, 25, 10 ** 6):
            with self.assertRaises(pyp_errors.PyPedalValidationError):
                pyp_metrics.a_effective_ancestors_indefinite(
                    load_griffon(), n=n, reference=list(cohort))

    def test_it_refuses_the_synthetic_fixture_at_every_n(self):
        for n in (1, 2, 25, 10 ** 6):
            with self.assertRaises(pyp_errors.PyPedalValidationError):
                pyp_metrics.a_effective_ancestors_indefinite(
                    load_synthetic(), n=n)

    def test_both_routines_report_the_same_root_cause(self):
        """
        Not merely "both raise". The two diagnostics must name the same
        defective quantity and the same measured sum, differing only in the
        routine that reports it -- otherwise a reader cannot tell that the two
        refusals are the same refusal.
        """
        cohort = list(griffon_36())
        with self.assertRaises(pyp_errors.PyPedalValidationError) as exact:
            pyp_metrics.a_effective_ancestors_definite(
                load_griffon(), reference=cohort)
        with self.assertRaises(pyp_errors.PyPedalValidationError) as bounded:
            pyp_metrics.a_effective_ancestors_indefinite(
                load_griffon(), n=5, reference=cohort)

        exact_message = str(exact.exception)
        bounded_message = str(bounded.exception)
        self.assertEqual(
            exact_message.replace("a_effective_ancestors_definite", "R", 1),
            bounded_message.replace("a_effective_ancestors_indefinite", "R", 1),
            "the two routines must give the same diagnostic for the same "
            "invalid sequence, differing only in the routine name")
        self.assertIn("marginal contributions", bounded_message)
        self.assertIn(repr(GRIFFON_36_SUM_P), bounded_message)
        self.assertIn("not a valid probability vector", bounded_message)

    def test_the_diagnostic_does_not_blame_the_reference_population(self):
        """
        The message describes the contribution vector. An analyst-defined R may
        be perfectly meaningful; what has failed is this implementation's
        ability to produce a valid contribution sequence for it.
        """
        with self.assertRaises(pyp_errors.PyPedalValidationError) as caught:
            pyp_metrics.a_effective_ancestors_indefinite(
                load_griffon(), n=5, reference=list(griffon_36()))
        message = str(caught.exception).lower()
        for forbidden in ("biolog", "invalid reference", "invalid population",
                          "bad reference"):
            self.assertNotIn(forbidden, message)

    def test_it_writes_no_report_when_it_refuses(self):
        """
        The check precedes the `.dat` write, so a refused run leaves nothing
        behind claiming a bound it did not stand behind.

        The report path comes from `kw['filetag']`, which is derived from the
        pedfile (`pyp_newclasses.py:86`) and is NOT the working directory --
        asserting on the cwd instead would pass vacuously whatever the routine
        did.
        """
        ped = load_synthetic()
        report = "%s_fa_boichard_indefinite_.dat" % ped.kw["filetag"]
        self.assertFalse(os.path.exists(report), "stale file before the call")
        try:
            pyp_metrics.a_effective_ancestors_indefinite(ped, n=2)
        except pyp_errors.PyPedalValidationError:
            pass
        self.assertFalse(
            os.path.exists(report),
            "a refused run wrote %s" % os.path.basename(report))

    def test_the_two_routines_use_the_same_validator(self):
        """
        Structural, so the repair cannot be satisfied by a second, drifting
        implementation of the same criterion.
        """
        self.assertIn("check_contribution_vector",
                      pyp_metrics.a_effective_ancestors_definite.__code__.co_names)
        self.assertIn("check_contribution_vector",
                      pyp_metrics.a_effective_ancestors_indefinite.__code__.co_names)


# ---------------------------------------------------------------------------
# 3. Containment. Everything here passes before AND after the repair; a
#    failure means the fix reached past the invalid domain.
# ---------------------------------------------------------------------------

class TestValidReferencePopulationsAreUnaffected(unittest.TestCase):

    def test_the_published_controls_are_unchanged(self):
        for name, expected in sorted(PUBLISHED.items()):
            with self.subTest(pedigree=name):
                self.assertEqual(
                    expected["f_a"],
                    pyp_metrics.a_effective_ancestors_definite(
                        load_corpus(name, "asdg")))
                for n, bounds in sorted(expected["bounds"].items()):
                    self.assertEqual(
                        bounds,
                        pyp_metrics.a_effective_ancestors_indefinite(
                            load_corpus(name, "asdg"), n=n),
                        "%s bounds changed at n=%d" % (name, n))

    def test_the_published_sequences_sum_to_one(self):
        """
        The precondition the repair rests on, asserted on the valid domain so
        that domain is demonstrably unaffected by it.
        """
        for name in sorted(PUBLISHED):
            with self.subTest(pedigree=name):
                ped = load_corpus(name, "asdg")
                top = max((str(a.gen) for a in ped.pedigree), key=int)
                reference = sorted(int(a.animalID) for a in ped.pedigree
                                   if str(a.gen) == top)
                order = list(pyp_metrics.boichard_marginal_contributions(
                    load_corpus(name, "asdg"), reference))
                self.assertEqual(1.0, sum(value for _, value in order))

    def test_the_valid_griffon_reference_is_unchanged(self):
        reference = list(griffon_9())
        self.assertEqual(9, len(reference))
        self.assertEqual(
            GRIFFON_9_F_A,
            pyp_metrics.a_effective_ancestors_definite(
                load_griffon(), reference=reference))
        for n, bounds in sorted(GRIFFON_9_BOUNDS.items()):
            with self.subTest(n=n):
                self.assertEqual(
                    bounds,
                    pyp_metrics.a_effective_ancestors_indefinite(
                        load_griffon(), n=n, reference=reference))

    def test_the_bounds_still_bracket_the_exact_value(self):
        for name in sorted(PUBLISHED):
            exact = pyp_metrics.a_effective_ancestors_definite(
                load_corpus(name, "asdg"))
            for n in (1, 2, 3, 25):
                with self.subTest(pedigree=name, n=n):
                    f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(name, "asdg"), n=n)
                    self.assertLessEqual(f_l, exact + 1e-9)
                    self.assertLessEqual(exact, f_u + 1e-9)

    def test_a_truncated_prefix_is_not_what_gets_validated(self):
        """
        The distinction the repair turns on. At n=1 the prefix sums to p_1,
        which is far from 1 on every one of these pedigrees -- and that is
        correct, because bounding is the entire point. Only the FULL sequence
        carries the sum-to-one requirement.
        """
        for name in sorted(PUBLISHED):
            with self.subTest(pedigree=name):
                ped = load_corpus(name, "asdg")
                top = max((str(a.gen) for a in ped.pedigree), key=int)
                reference = sorted(int(a.animalID) for a in ped.pedigree
                                   if str(a.gen) == top)
                order = list(pyp_metrics.boichard_marginal_contributions(
                    load_corpus(name, "asdg"), reference))
                self.assertLess(order[0][1], 1.0)
                self.assertEqual(
                    PUBLISHED[name]["bounds"][1],
                    pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(name, "asdg"), n=1))

    def test_the_founders_routine_is_untouched_on_the_invalid_cohort(self):
        """
        Appendix A credits every founder directly, so `q` sums to one even
        here. This routine does not share the defect and must not acquire a
        refusal from the repair.
        """
        self.assertEqual(
            40.949888916316965,
            pyp_metrics.a_effective_founders_boichard(
                load_griffon(), reference=list(griffon_36())))


class TestValidationGatingIsUnchanged(unittest.TestCase):
    """
    `check_contribution_vector` respects `kw['validate']`. The repair must
    inherit that exactly, and must not introduce an unconditional refusal:
    with validation off, both routines compute as they always did.
    """

    def _griffon_without_validation(self):
        from test_boichard_reference_population import GRIFFON_OPTIONS
        from _pedhelpers import load_griffon_1871_1890
        options = dict(GRIFFON_OPTIONS)
        options["validate"] = False
        return load_griffon_1871_1890(options)

    def test_the_exact_routine_still_computes_with_validation_disabled(self):
        self.assertEqual(
            162.0,
            pyp_metrics.a_effective_ancestors_definite(
                self._griffon_without_validation(),
                reference=list(griffon_36())))

    def test_the_bounded_routine_still_computes_with_validation_disabled(self):
        self.assertEqual(
            (36.0, 102.14375788146279),
            pyp_metrics.a_effective_ancestors_indefinite(
                self._griffon_without_validation(), n=5,
                reference=list(griffon_36())))

    def test_the_synthetic_fixture_computes_with_validation_disabled(self):
        self.assertEqual(
            (4.0, 6.222222222222223),
            pyp_metrics.a_effective_ancestors_indefinite(
                load_synthetic(validate=False), n=2))


class TestTheEndpointGuardsAreUnchangedAndAccountedFor(unittest.TestCase):
    """
    Four documented endpoint behaviours of the contribution sequence
    must survive. Two of them are ALGEBRAICALLY
    UNREACHABLE once the contribution vector is known to be a probability
    vector, and that is recorded here rather than papered over with a
    manufactured fixture.

    ================================================================
    guard                                     status after 
    ================================================================
    residual <= tol -> f_l = f_u = f_a        REACHABLE, exercised below
    n == f, never divide by (f - n)           REACHABLE, exercised below,
                                              subsumed by the row above
    n_taken >= n_founders with residual > tol UNREACHABLE: k <= f (p.8) and
                                              c = 1 at k, so residual > tol
                                              implies n_taken < k <= f
    p_n <= tol with residual > tol            UNREACHABLE: the engine stops
                                              at `best <= tol`, so every
                                              yielded contribution exceeds it
    ================================================================

    Both unreachable guards were already unreachable on VALID input before
     -- the repair did not make them so. What changed is that the
    invalid sequences which could formerly reach them are now refused earlier,
    with a diagnostic naming the root cause instead of a symptom. They are
    retained as defensive assertions on an internally inconsistent state, and
    the two properties that make them unreachable are asserted below so the
    claim cannot rot.
    """

    #: The engine's own termination threshold and the routine's residual
    #: tolerance, both 1e-12 in `pyp_metrics`.
    TOL = 1e-12

    def _sequence(self, name):
        ped = load_corpus(name, "asdg")
        top = max((str(a.gen) for a in ped.pedigree), key=int)
        reference = sorted(int(a.animalID) for a in ped.pedigree
                           if str(a.gen) == top)
        order = list(pyp_metrics.boichard_marginal_contributions(
            load_corpus(name, "asdg"), reference))
        _i, _s, _d, _ph, n_founders = pyp_metrics._boichard_completed_arrays(
            load_corpus(name, "asdg"))
        return order, n_founders

    def test_the_zero_residual_endpoint_collapses_onto_the_exact_value(self):
        for name in sorted(PUBLISHED):
            with self.subTest(pedigree=name):
                exact = pyp_metrics.a_effective_ancestors_definite(
                    load_corpus(name, "asdg"))
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_corpus(name, "asdg"), n=10 ** 6)
                self.assertEqual(exact, f_l)
                self.assertEqual(exact, f_u)

    def test_taking_exactly_f_ancestors_does_not_divide_by_zero(self):
        """
        `boichard_fig1.ped` has k == f == 6, so n=6 drives `n_taken` to the
        founder count exactly. The residual is zero there, so the short-circuit
        fires and `(1-c)^2/(f-n)` is never evaluated -- which is the adjudicated
        behaviour, not an accident of this fixture.
        """
        order, n_founders = self._sequence("boichard_fig1.ped")
        self.assertEqual(6, len(order))
        self.assertEqual(6, n_founders)
        exact = pyp_metrics.a_effective_ancestors_definite(
            load_corpus("boichard_fig1.ped", "asdg"))
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            load_corpus("boichard_fig1.ped", "asdg"), n=n_founders)
        self.assertEqual(exact, f_l)
        self.assertEqual(exact, f_u)

    def test_no_yielded_contribution_can_trip_the_zero_probability_guard(self):
        """
        Why `p_n <= tol` is unreachable: the engine terminates at
        `best <= tol`, so nothing at or below the tolerance is ever yielded.
        """
        for name in sorted(PUBLISHED):
            with self.subTest(pedigree=name):
                order, _f = self._sequence(name)
                self.assertTrue(order)
                self.assertGreater(min(value for _, value in order), self.TOL)

    def test_the_contributor_count_never_exceeds_the_founder_count(self):
        """
        Why `n_taken >= n_founders` with mass outstanding is unreachable:
        Boichard p.8 gives k <= f, and the sequence is exhausted at k, where
        c = 1 and the residual short-circuit has already fired.
        """
        for name in sorted(PUBLISHED):
            with self.subTest(pedigree=name):
                order, n_founders = self._sequence(name)
                self.assertLessEqual(len(order), n_founders)

    def test_the_division_guards_are_still_present(self):
        """
        Structural. Unreachable is not the same as unnecessary: both guards
        stand between an internally inconsistent state and a division, and
        deleting them because no fixture reaches them would be removing the
        detector rather than the fault.
        """
        import inspect
        source = inspect.getsource(pyp_metrics.a_effective_ancestors_indefinite)
        self.assertIn("n_taken >= n_founders", source)
        self.assertIn("p_n <= tol", source)
        self.assertIn("residual <= tol", source)

    def test_the_bounds_interval_is_never_inverted(self):
        for name in sorted(PUBLISHED):
            for n in (1, 2, 3, 25, 10 ** 6):
                with self.subTest(pedigree=name, n=n):
                    f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(name, "asdg"), n=n)
                    self.assertLessEqual(f_l, f_u)

    def test_n_below_one_is_still_a_usage_error(self):
        for n in (0, -1):
            with self.subTest(n=n):
                with self.assertRaises(pyp_errors.PyPedalUsageError):
                    pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus("boichard2a.ped", "asdg"), n=n)


if __name__ == "__main__":
    unittest.main()
