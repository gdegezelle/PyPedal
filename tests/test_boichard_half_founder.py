"""
 -- half-founders in Boichard's Appendix B, under Reading C.

THE RULE, AND ITS EVIDENCE CLASS
--------------------------------
For every animal with exactly one known parent, the unknown parental SLOT
becomes a founder with a record of its own; Appendix B then runs on the
completed pedigree with no special case at all. Appendix A step 4 is NOT also
applied -- after completion there are no half-founders left for it to apply to.

The evidence class is the one the project adjudicated, and no higher:

    MATHEMATICALLY IMPLIED / INDEPENDENTLY SUPPORTED.
    NOT SOURCE-EXPLICIT APPENDIX-B TEXT.

What IS source-explicit is Appendix A step 4, article p.22:

    "(4) if an animal is a 'half founder' (ie, with one known parent and one
    unknown parent), multiply its contribution by 0.5. This is equivalent to
    considering the unknown parent as a founder. Divide the vector q by N, so
    that founder contributions sum to 1."

and the founder definition on p.7:

    "A founder is defined as an ancestor with unknown parents. Note that when
    an animal has only one known parent, the animal is considered as a
    founder."

**Appendix B never restates the halving, and never excludes it.** It is silent,
and silence is not an exclusion. Reading C was chosen by adjudication against
invariants the paper states about itself -- p.7 Sum q = 1, p.8 Sum p = 1 and
k <= f, p.17 f_a <= f_e -- which the two rival readings violate and Reading C
does not. Nothing here may be cited as Appendix B saying this.

WHAT ANCHORS THESE TESTS
------------------------
Analytic fixtures with the derivation written out, first. The independent
oracle corroborates; it is never the primary anchor, and it is checked BEFORE
production exists so that it cannot later be tuned into agreement with it.

The strongest single anchor is external to this phase: on FIXTURE_A the
Appendix-A halving that production has shipped since  gives
f_e = 2.909090909090909, and Reading C's completion gives the same number to
the last bit. Appendix A step 4 says in its own words that it should. That
value is already pinned independently in
tests/test_boichard_production.py::TestEffectiveFoundersAppendixA.

WHAT IS NOT SETTLED HERE
------------------------
R1 (tie-breaking) and R3 (reference-population candidacy) are untouched. Both
remain declared conventions, and lowest_id is still NOT Boichard's rule.
 is untouched and still open, so these fixtures carry an explicit
'g' column -- exactly as the adjudication probe does, and for the same reason.

ID COORDINATES -- read this before adding a fixture
---------------------------------------------------
These pedigrees are written so that load() renumbers them to themselves:
animals appear oldest-first and are already numbered 1..n, so
originalID == animalID == index + 1. Phantom IDs are assigned by production
ABOVE every real ID and are not pedigree members. If you add a fixture that
does not have this property, map through pedobj.idmap explicitly rather than
assuming it.
"""
import os
import sys
import tempfile
import unittest

from _pedhelpers import corpus, load_corpus, load_corpus_from_path
from oracles import (
    boichard_bounds,
    boichard_f_a,
    boichard_f_e,
    boichard_founders,
    boichard_marginal_contributions,
    boichard_phantom_complete,
    boichard_read,
    boichard_unknown_parent_slots,
    BOICHARD_PHANTOM_ENCODINGS,
)

from PyPedal import pyp_errors, pyp_metrics

# ---------------------------------------------------------------------------
# Fixtures. Every one carries an explicit 'g' column (pedformat 'asdg'),
# because pyp_utils.set_generation writes animal.igen while the Boichard
# routines read animal.gen -- , still open. Inferring generations
# here would make a scientific expectation depend on an unrelated defect.
#
# FIXTURE_A  known sire / unknown dam. Animal 4's dam is the only slot.
#            This is the pedigree the old R2 refusal was demonstrated on, kept
#            deliberately so the refusal tests become value tests on the SAME
#            input rather than on a friendlier one.
#
#            Completed: phantom P = dam of 4.  Founders {1, 2, P}, so f = 3.
#            Appendix A, reference {5, 6}, N = 2:
#              q(5)=q(6)=1 -> q(3)=1, q(4)=1 -> q(3)=1.5, q(P)=0.5
#                          -> q(1)=q(2)=0.75
#              /N:  q(1)=q(2)=0.375, q(P)=0.25.  Sum = 1.  f_e = 1/0.34375
#                                                              = 2.909090909...
#            Appendix B round 1: p = q, references zeroed, max 1.5 at animal 3
#              -> contribution 1.5/2 = 0.75.
#            Round 2: 3 becomes a pseudo-founder, a(3)=1, a(4)=0.5.
#              p(4) = 1*(1-0.5) = 0.5 and p(P) = 0.5*(1-0) = 0.5 -- A REAL
#              ANIMAL AND ITS OWN PHANTOM PARENT TIE. R1's lowest_id picks the
#              real animal 4, because phantoms are numbered above every real
#              ID. Contribution 0.5/2 = 0.25.
#            Round 3: everything is explained; p is zero throughout, stop.
#              Sum p = 1.0 exactly.  f_a = 1/(0.75^2 + 0.25^2) = 1/0.625 = 1.6
#
# FIXTURE_B  the mirror of A: unknown SIRE, known dam. Every published quantity
#            must be identical -- the paper's rule is about a missing side, not
#            about which side.
#
# FIXTURE_C  two independent half-founders, one slot each. Proves two sentinel
#            zeroes do not collapse into one shared parent: with two distinct
#            phantoms f_e = 4.0 and f_a = 2.0; sharing one would not.
#
# FIXTURE_D  a chain in which a PHANTOM is itself selected as an ancestor,
#            holding 0.5 of the total mass. Reading C makes phantoms ordinary
#            founders and R3 does not reach them, so this is forced, not
#            chosen: withholding that 0.5 would leave Sum p = 0.5 and break the
#            p.8 invariant Reading C was selected by.
#              Round 1: animal 3 and phantom(dam of 4) tie at 0.5; the real
#                       animal wins -> 0.5.
#              Round 2: only phantom(dam of 4) is left with mass -> 0.5.
#              Sum p = 1.0, f_a = 2.0, f_e = 2.909090909..., f_a <= f_e.
#
# FIXTURE_NONE  the control: no half-founders at all. Completion has nothing to
#            do, so every value must be bit-identical to today's.
# ---------------------------------------------------------------------------

FIXTURE_A = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 3 0 1\n5 3 4 2\n6 3 4 2\n"
FIXTURE_B = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 0 3 1\n5 3 4 2\n6 3 4 2\n"
FIXTURE_C = "1 0 0 1\n2 0 0 1\n3 1 0 1\n4 2 0 1\n5 3 4 2\n6 3 4 2\n"
FIXTURE_D = "1 0 0 1\n2 1 0 1\n3 2 0 1\n4 3 0 2\n"
FIXTURE_NONE = "1 0 0 1\n2 0 0 1\n3 1 2 1\n4 1 2 1\n5 3 4 2\n6 3 4 2\n"

#: (name, text, rows, reference) -- the rows/reference restated for the oracle,
#: which parses nothing from PyPedal. Getting one wrong yields a vacuously
#: agreeing comparison, so they are written out rather than derived.
FIXTURES = {
    "A": (FIXTURE_A,
          [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 0), (5, 3, 4), (6, 3, 4)],
          [5, 6]),
    "B": (FIXTURE_B,
          [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 0, 3), (5, 3, 4), (6, 3, 4)],
          [5, 6]),
    "C": (FIXTURE_C,
          [(1, 0, 0), (2, 0, 0), (3, 1, 0), (4, 2, 0), (5, 3, 4), (6, 3, 4)],
          [5, 6]),
    "D": (FIXTURE_D,
          [(1, 0, 0), (2, 1, 0), (3, 2, 0), (4, 3, 0)],
          [4]),
    "NONE": (FIXTURE_NONE,
             [(1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 1, 2), (5, 3, 4), (6, 3, 4)],
             [5, 6]),
}

#: Hand-derived above, then confirmed by the independent oracle. Exact values,
#: asserted with no tolerance at all.
EXPECTED = {
    "A": {"f_a": 1.6, "f_e": 2.909090909090909, "k": 2, "f": 3, "slots": 1},
    "B": {"f_a": 1.6, "f_e": 2.909090909090909, "k": 2, "f": 3, "slots": 1},
    "C": {"f_a": 2.0, "f_e": 4.0, "k": 2, "f": 4, "slots": 2},
    "D": {"f_a": 2.0, "f_e": 2.909090909090909, "k": 2, "f": 4, "slots": 3},
    "NONE": {"slots": 0},
}

HALF_FOUNDER_FIXTURES = ("A", "B", "C", "D")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def write_pedigree(text):
    """Materialise a fixture in a throwaway directory, never in the repo."""
    tmp = tempfile.mkdtemp(prefix="pypedal_bl1r2_")
    path = os.path.join(tmp, "fixture.ped")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def load_fixture(name):
    return load_corpus_from_path(write_pedigree(FIXTURES[name][0]), "asdg")


def oracle_rows(name):
    _text, rows, reference = FIXTURES[name]
    return rows, reference


# ---------------------------------------------------------------------------
# 1. The oracle, checked before production exists
# ---------------------------------------------------------------------------

class TestOracleAgainstTheAnalyticFixtures(unittest.TestCase):
    """
    The oracle shares no code with PyPedal and parses the pedigree itself. If
    it disagrees with the hand derivations in the header, the derivations are
    wrong and nothing below is worth running.
    """

    def test_f_a_matches_the_derivation(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                got, _order, _total = boichard_f_a(
                    rows, reference, half_founder="complete")
                self.assertAlmostEqual(EXPECTED[name]["f_a"], got, places=12)

    def test_f_e_matches_the_derivation(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                got, _q, _total = boichard_f_e(
                    rows, reference, half_founder="complete")
                self.assertAlmostEqual(EXPECTED[name]["f_e"], got, places=12)

    def test_marginal_contributions_sum_to_one_without_being_normalised(self):
        """
        Boichard p.8: "the marginal contributions over all ancestors sum to
        one". Appendix B step 8 divides by N and by nothing else, so a 1.0 here
        is unforced and is a real check -- this is the invariant that
        eliminated the in-place-halving reading.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                _f_a, _order, total = boichard_f_a(
                    rows, reference, half_founder="complete")
                self.assertAlmostEqual(1.0, total, places=12)

    def test_founder_contributions_sum_to_one_without_being_normalised(self):
        """Boichard p.7 -- the invariant that eliminated the no-halving reading."""
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                _f_e, _q, total = boichard_f_e(
                    rows, reference, half_founder="complete")
                self.assertAlmostEqual(1.0, total, places=12)

    def test_f_a_never_exceeds_f_e(self):
        """Boichard p.17: "Consequently f_a is always lower than or equal to f_e"."""
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                f_a, _o, _t = boichard_f_a(rows, reference, half_founder="complete")
                f_e, _q, _s = boichard_f_e(rows, reference, half_founder="complete")
                self.assertLessEqual(f_a, f_e + 1e-12)

    def test_ancestor_count_does_not_exceed_the_founder_count(self):
        """
        Boichard p.8: "The number of ancestors with a positive contribution is
        less than or equal to the total number of founders." ``f`` is counted
        on the COMPLETED pedigree, which is the one the engine ran on.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                completed, _phantoms = boichard_phantom_complete(rows)
                _f_a, order, _t = boichard_f_a(
                    rows, reference, half_founder="complete")
                self.assertLessEqual(len(order), len(boichard_founders(completed)))
                self.assertEqual(EXPECTED[name]["k"], len(order))
                self.assertEqual(EXPECTED[name]["f"],
                                 len(boichard_founders(completed)))


class TestOracleReadingCEqualsAppendixAOnFounders(unittest.TestCase):
    """
    Appendix A step 4's own words: halving a half-founder's contribution "is
    equivalent to considering the unknown parent as a founder". That is a
    source-explicit claim about Appendix A, and it is checkable: f_e computed
    by halving and f_e computed by completion must agree exactly.

    This is the bridge between the source-explicit part of the evidence and
    the adjudicated part. It does not make Reading C source-explicit for
    Appendix B, and must never be cited as though it did.
    """

    def test_halving_and_completion_agree_exactly(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                halved, _q1, _s1 = boichard_f_e(rows, reference, half_founder="halve")
                completed, _q2, _s2 = boichard_f_e(rows, reference,
                                                   half_founder="complete")
                self.assertAlmostEqual(halved, completed, places=12)


class TestOracleRivalReadingsAreExcluded(unittest.TestCase):
    """
    A test that only confirms the chosen reading observes rather than
    discriminates. The two rival readings are stated here as executable
    assertions, so that "Reading C passes the invariants" carries the weight it
    is supposed to carry.
    """

    def test_no_half_founder_handling_breaks_the_founder_sum(self):
        """Reading B: q stops being a probability of gene origin (p.7)."""
        broken = 0
        for name in HALF_FOUNDER_FIXTURES:
            rows, reference = oracle_rows(name)
            _f_e, _q, total = boichard_f_e(rows, reference, half_founder="none")
            if abs(total - 1.0) > 1e-9:
                broken += 1
        self.assertEqual(len(HALF_FOUNDER_FIXTURES), broken,
                         "the no-halving reading must fail Sum q = 1 on every "
                         "half-founder fixture; if it stopped doing so the "
                         "discrimination has been lost")

    def test_in_place_halving_breaks_the_ancestor_sum(self):
        """
        Reading A: halving truncates the half-founder's contribution and the
        unknown parent has no node to be credited to, so mass is simply lost.
        """
        broken = 0
        for name in HALF_FOUNDER_FIXTURES:
            rows, reference = oracle_rows(name)
            _f_a, _order, total = boichard_f_a(rows, reference,
                                               half_founder="halve")
            if abs(total - 1.0) > 1e-9:
                broken += 1
        self.assertGreater(broken, 0,
                           "the in-place-halving reading must fail Sum p = 1 "
                           "somewhere in this corpus, or it is not being "
                           "discriminated against")

    def test_the_control_cannot_tell_the_readings_apart(self):
        """
        On a pedigree with no half-founders all three readings have nothing to
        do. If they differ here, the scorer is biased and every result above is
        suspect.
        """
        rows, reference = oracle_rows("NONE")
        values = set()
        for rule in ("none", "halve", "complete"):
            f_a, _o, _t = boichard_f_a(rows, reference, half_founder=rule)
            f_e, _q, _s = boichard_f_e(rows, reference, half_founder=rule)
            values.add((round(f_a, 12), round(f_e, 12)))
        self.assertEqual(1, len(values))


class TestOraclePhantomEncodingDiscrimination(unittest.TestCase):
    """
    The phantom integers mean nothing, and that has to be MEASURED rather than
    argued from the fact that both encodings sort above the real IDs.

    Production uses contiguous max(real)+n; the adjudication probe used
    max(real)*100+n. If these two ever disagree on a scientific quantity, the
    numbering has become load-bearing through R1's tie-break and the correct
    response is to STOP and report -- not to pick whichever encoding gives the
    nicer answer.
    """

    def test_every_scientific_quantity_is_encoding_independent(self):
        for name in HALF_FOUNDER_FIXTURES:
            for encoding in BOICHARD_PHANTOM_ENCODINGS:
                with self.subTest(fixture=name, encoding=encoding):
                    rows, reference = oracle_rows(name)
                    f_a, order, total = boichard_f_a(
                        rows, reference, half_founder="complete",
                        phantom_encoding=encoding)
                    f_e, _q, q_total = boichard_f_e(
                        rows, reference, half_founder="complete",
                        phantom_encoding=encoding)
                    self.assertAlmostEqual(EXPECTED[name]["f_a"], f_a, places=12)
                    self.assertAlmostEqual(EXPECTED[name]["f_e"], f_e, places=12)
                    self.assertAlmostEqual(1.0, total, places=12)
                    self.assertAlmostEqual(1.0, q_total, places=12)
                    self.assertEqual(EXPECTED[name]["k"], len(order))

    def test_the_selected_sequence_agrees_after_mapping_phantoms_by_slot(self):
        """
        Stronger than agreeing on f_a: the two encodings must select the same
        ANIMALS in the same order with the same contributions, once each
        phantom is named by the parental slot it stands for rather than by its
        arbitrary integer.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                mapped = []
                for encoding in BOICHARD_PHANTOM_ENCODINGS:
                    _completed, phantoms = boichard_phantom_complete(rows, encoding)
                    order = list(boichard_marginal_contributions(
                        rows, reference, half_founder="complete",
                        phantom_encoding=encoding))
                    mapped.append([(phantoms.get(a, ("real", a)), round(v, 12))
                                   for a, v in order])
                self.assertEqual(mapped[0], mapped[1])

    def test_a_phantom_never_outranks_a_real_animal_on_a_tie(self):
        """
        R1 containment. Phantoms are numbered above every real ID under both
        encodings, so lowest_id still means "the oldest real animal" wherever a
        real animal is available. FIXTURE_A exercises exactly this: animal 4
        and its own phantom dam tie at 0.5 in round 2.
        """
        rows, reference = oracle_rows("A")
        _completed, phantoms = boichard_phantom_complete(rows)
        order = list(boichard_marginal_contributions(
            rows, reference, half_founder="complete"))
        self.assertEqual([(3, 0.75), (4, 0.25)],
                         [(a, round(v, 12)) for a, v in order])
        self.assertNotIn(list(phantoms)[0], [a for a, _ in order])

    def test_a_phantom_is_selected_when_it_holds_real_mass(self):
        """
        The other side of the same coin, and the reason phantoms cannot simply
        be excluded from candidacy: on FIXTURE_D a phantom holds 0.5 of the
        total, and dropping it would leave Sum p = 0.5.
        """
        rows, reference = oracle_rows("D")
        _completed, phantoms = boichard_phantom_complete(rows)
        order = list(boichard_marginal_contributions(
            rows, reference, half_founder="complete"))
        selected = [a for a, _ in order]
        chosen_phantoms = [a for a in selected if a in phantoms]
        self.assertEqual(1, len(chosen_phantoms))
        self.assertAlmostEqual(1.0, sum(v for _, v in order), places=12)


class TestOracleSlotIdentity(unittest.TestCase):
    """
    One phantom per missing parental SLOT. Two sentinel zeroes are two distinct
    unknown individuals, and nothing in the input licenses inferring otherwise.
    """

    def test_one_phantom_per_slot(self):
        for name, expected in EXPECTED.items():
            with self.subTest(fixture=name):
                rows, _reference = oracle_rows(name)
                slots = boichard_unknown_parent_slots(rows)
                _completed, phantoms = boichard_phantom_complete(rows)
                self.assertEqual(expected["slots"], len(slots))
                self.assertEqual(expected["slots"], len(phantoms))
                self.assertEqual(len(phantoms), len(set(phantoms)))

    def test_two_half_founders_do_not_share_a_parent(self):
        """
        FIXTURE_C's two half-founders both record 0. If they collapsed to one
        phantom the founder set would be {1, 2, shared} and f_e would not be
        4.0. The distinctness is structural -- the fill map is keyed on
        (animal, side) -- and this pins the consequence.
        """
        rows, reference = oracle_rows("C")
        completed, phantoms = boichard_phantom_complete(rows)
        self.assertEqual(2, len(phantoms))
        parents = {a: (s, d) for a, s, d in completed}
        self.assertNotEqual(parents[3][1], parents[4][1])
        f_e, _q, _t = boichard_f_e(rows, reference, half_founder="complete")
        self.assertAlmostEqual(4.0, f_e, places=12)

    def test_an_animal_with_both_parents_unknown_yields_no_slot(self):
        """A real founder is already a founder; completion must leave it alone."""
        rows, _reference = oracle_rows("A")
        completed, _phantoms = boichard_phantom_complete(rows)
        parents = {a: (s, d) for a, s, d in completed}
        self.assertEqual((0, 0), parents[1])
        self.assertEqual((0, 0), parents[2])

    def test_completion_is_a_no_op_without_half_founders(self):
        rows, _reference = oracle_rows("NONE")
        completed, phantoms = boichard_phantom_complete(rows)
        self.assertEqual([], list(phantoms))
        self.assertEqual(rows, completed)

    def test_phantom_ids_are_above_every_real_id(self):
        for name in HALF_FOUNDER_FIXTURES:
            for encoding in BOICHARD_PHANTOM_ENCODINGS:
                with self.subTest(fixture=name, encoding=encoding):
                    rows, _reference = oracle_rows(name)
                    _completed, phantoms = boichard_phantom_complete(rows, encoding)
                    top = max(a for a, _, _ in rows)
                    for pid in phantoms:
                        self.assertGreater(pid, top)
                        self.assertNotEqual(0, pid)

    def test_completion_preserves_parent_before_offspring(self):
        """Appendix B's single-pass steps 5 and 6 require it."""
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, _reference = oracle_rows(name)
                completed, _phantoms = boichard_phantom_complete(rows)
                seen = set()
                for a, s, d in completed:
                    for parent in (s, d):
                        if parent != 0:
                            self.assertIn(parent, seen)
                    seen.add(a)


class TestOracleExplicitCompletionEquivalence(unittest.TestCase):
    """
    If Reading C really is pedigree-level completion followed by an unchanged
    Appendix B, then analysing P under Reading C must equal analysing a
    hand-written P' -- a pedigree with REAL founder rows standing in for the
    unknown slots -- under no half-founder rule at all.

    P' is written out by hand rather than generated, so this compares Reading C
    against an independent statement of what it is supposed to mean, not
    against itself.
    """

    #: (fixture, explicit rows of P', {explicit stand-in id: (animal, side)})
    EXPLICIT = {
        "A": ([(90, 0, 0), (1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 3, 90),
               (5, 3, 4), (6, 3, 4)],
              {90: (4, "dam")}),
        "B": ([(90, 0, 0), (1, 0, 0), (2, 0, 0), (3, 1, 2), (4, 90, 3),
               (5, 3, 4), (6, 3, 4)],
              {90: (4, "sire")}),
        "C": ([(90, 0, 0), (91, 0, 0), (1, 0, 0), (2, 0, 0), (3, 1, 90),
               (4, 2, 91), (5, 3, 4), (6, 3, 4)],
              {90: (3, "dam"), 91: (4, "dam")}),
        "D": ([(90, 0, 0), (91, 0, 0), (92, 0, 0), (1, 0, 0), (2, 1, 90),
               (3, 2, 91), (4, 3, 92)],
              {90: (2, "dam"), 91: (3, "dam"), 92: (4, "dam")}),
    }

    def test_f_a_agrees_exactly(self):
        for name, (explicit, _names) in self.EXPLICIT.items():
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                want, _o, _t = boichard_f_a(explicit, reference)
                got, _o2, _t2 = boichard_f_a(rows, reference,
                                             half_founder="complete")
                self.assertAlmostEqual(want, got, places=12)

    def test_f_e_agrees_exactly(self):
        for name, (explicit, _names) in self.EXPLICIT.items():
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                want, _q, _t = boichard_f_e(explicit, reference,
                                            half_founder="none")
                got, _q2, _t2 = boichard_f_e(rows, reference,
                                             half_founder="complete")
                self.assertAlmostEqual(want, got, places=12)

    def test_the_ancestor_sequence_agrees_after_naming_the_stand_ins(self):
        for name, (explicit, names) in self.EXPLICIT.items():
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                _completed, phantoms = boichard_phantom_complete(rows)
                want = [(names.get(a, ("real", a)), round(v, 12))
                        for a, v in boichard_marginal_contributions(
                            explicit, reference)]
                got = [(phantoms.get(a, ("real", a)), round(v, 12))
                       for a, v in boichard_marginal_contributions(
                           rows, reference, half_founder="complete")]
                self.assertEqual(want, got)

    def test_bounds_agree_exactly(self):
        for name, (explicit, _names) in self.EXPLICIT.items():
            for n in (1, 2):
                with self.subTest(fixture=name, n=n):
                    rows, reference = oracle_rows(name)
                    want = boichard_bounds(explicit, reference, n)[:2]
                    got = boichard_bounds(rows, reference, n,
                                          half_founder="complete")[:2]
                    self.assertAlmostEqual(want[0], got[0], places=12)
                    self.assertAlmostEqual(want[1], got[1], places=12)


class TestOracleBoundsBracketTheExactValue(unittest.TestCase):
    """
    , pp.9-10. The bounds are a truncation of the same sequence, so they
    must bracket the exact value at every n and collapse to it at the end,
    without ever evaluating (1-c)^2/(f-n) at n == f.
    """

    def test_bounds_bracket_f_a_at_every_n(self):
        for name in HALF_FOUNDER_FIXTURES:
            rows, reference = oracle_rows(name)
            exact, order, _t = boichard_f_a(rows, reference,
                                            half_founder="complete")
            for n in range(1, len(order) + 1):
                with self.subTest(fixture=name, n=n):
                    f_l, f_u, _meta = boichard_bounds(
                        rows, reference, n, half_founder="complete")
                    self.assertLessEqual(f_l, exact + 1e-9)
                    self.assertLessEqual(exact, f_u + 1e-9)

    def test_the_final_endpoint_is_exact_not_a_division_by_zero(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                exact, order, _t = boichard_f_a(rows, reference,
                                                half_founder="complete")
                f_l, f_u, meta = boichard_bounds(
                    rows, reference, len(order), half_founder="complete")
                self.assertTrue(meta["exact"])
                self.assertAlmostEqual(exact, f_l, places=12)
                self.assertAlmostEqual(exact, f_u, places=12)


class TestOracleIndependence(unittest.TestCase):
    """
    Assertions about imports, not about numbers. The oracle is only evidence
    while it stays independent, and the cheapest way for that to be lost is a
    shared helper appearing later "to avoid duplication".
    """

    @staticmethod
    def _import_lines(path):
        """
        Only the import statements. Both files legitimately MENTION the other
        side in prose -- the oracle's own docstring opens "Does not import
        PyPedal", and pyp_metrics cites docs that live next to the harness --
        so a substring search over the whole file tests the comments rather
        than the dependencies.
        """
        with open(path, "r", encoding="utf-8") as fh:
            return [line.strip() for line in fh
                    if line.startswith(("import ", "from "))]

    def test_the_oracle_does_not_import_pypedal(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "oracles", "oracle_boichard.py")
        for line in self._import_lines(path):
            self.assertNotIn("PyPedal", line)

    def test_production_does_not_import_the_oracle(self):
        for line in self._import_lines(os.path.abspath(pyp_metrics.__file__)):
            self.assertNotIn("oracle", line)
            self.assertNotIn("difftest", line)


# ---------------------------------------------------------------------------
# 2. Production. Expected to fail until Reading C is implemented.
# ---------------------------------------------------------------------------

class TestProductionComputesTheAdjudicatedValues(unittest.TestCase):

    def test_f_a_matches_the_derivation(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                got = pyp_metrics.a_effective_ancestors_definite(
                    load_fixture(name))
                self.assertAlmostEqual(EXPECTED[name]["f_a"], got, places=12)

    def test_f_a_agrees_with_the_independent_oracle(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                want, _o, _t = boichard_f_a(rows, reference,
                                            half_founder="complete")
                got = pyp_metrics.a_effective_ancestors_definite(
                    load_fixture(name))
                self.assertAlmostEqual(want, got, places=12)

    def test_f_a_never_exceeds_f_e(self):
        """p.17, across the two production routines rather than inside one."""
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                f_a = pyp_metrics.a_effective_ancestors_definite(
                    load_fixture(name))
                f_e = pyp_metrics.a_effective_founders_boichard(
                    load_fixture(name))
                self.assertLessEqual(f_a, f_e + 1e-9)

    def test_appendix_a_halving_still_gives_the_same_f_e(self):
        """
        NOT an xfail: this passes TODAY. a_effective_founders_boichard already
        accepts half-founders --  is settled on Appendix A step 4, which is
        source-explicit -- and this phase does not change it.

        Pinned here, in the  file, because these are the numbers Reading
        C has to land on: step 4's own equivalence clause promises that
        completion and halving give the same f_e, so if this drifts, the
        anchor for every f_e expectation below has moved.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                got = pyp_metrics.a_effective_founders_boichard(
                    load_fixture(name))
                self.assertAlmostEqual(EXPECTED[name]["f_e"], got, places=12)

    def test_orientation_does_not_matter(self):
        """
        The rule is about a missing side, not about which side. A and B are
        mirrors and must agree bit for bit.
        """
        a = pyp_metrics.a_effective_ancestors_definite(
            load_fixture("A"))
        b = pyp_metrics.a_effective_ancestors_definite(
            load_fixture("B"))
        self.assertEqual(a, b)


class TestProductionMarginalContributions(unittest.TestCase):

    def test_contributions_sum_to_one(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                _text, _rows, reference = FIXTURES[name]
                order = list(pyp_metrics.boichard_marginal_contributions(
                    ped, reference))
                self.assertAlmostEqual(1.0, sum(v for _, v in order), places=12)

    def test_the_engine_agrees_with_the_oracle_after_mapping_by_slot(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                rows, reference = oracle_rows(name)
                _completed, oracle_phantoms = boichard_phantom_complete(rows)
                want = [(oracle_phantoms.get(a, ("real", a)), round(v, 12))
                        for a, v in boichard_marginal_contributions(
                            rows, reference, half_founder="complete")]

                ped = load_fixture(name)
                phantoms = pyp_metrics.boichard_phantom_ids(ped)
                got = [(phantoms.get(a, ("real", a)), round(v, 12))
                       for a, v in pyp_metrics.boichard_marginal_contributions(
                           ped, reference)]
                self.assertEqual(want, got)

    def test_phantom_ids_are_reported_one_per_slot_above_every_real_id(self):
        for name, expected in EXPECTED.items():
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                phantoms = pyp_metrics.boichard_phantom_ids(ped)
                self.assertEqual(expected["slots"], len(phantoms))
                top = max(int(a.animalID) for a in ped.pedigree)
                for pid, slot in phantoms.items():
                    self.assertGreater(pid, top)
                    self.assertIn(slot[1], ("sire", "dam"))

    def test_a_phantom_is_selected_when_it_holds_real_mass(self):
        ped = load_fixture("D")
        _text, _rows, reference = FIXTURES["D"]
        phantoms = pyp_metrics.boichard_phantom_ids(ped)
        order = list(pyp_metrics.boichard_marginal_contributions(ped, reference))
        self.assertEqual(1, len([a for a, _ in order if a in phantoms]))


class TestProductionBounds(unittest.TestCase):

    def test_bounds_bracket_the_exact_value(self):
        for name in HALF_FOUNDER_FIXTURES:
            exact = pyp_metrics.a_effective_ancestors_definite(
                load_fixture(name))
            for n in (1, 2):
                with self.subTest(fixture=name, n=n):
                    f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                        load_fixture(name), n=n)
                    self.assertLessEqual(f_l, exact + 1e-9)
                    self.assertLessEqual(exact, f_u + 1e-9)

    def test_the_endpoint_collapses_to_the_exact_value(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                exact = pyp_metrics.a_effective_ancestors_definite(
                    load_fixture(name))
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_fixture(name), n=EXPECTED[name]["k"])
                self.assertAlmostEqual(exact, f_l, places=12)
                self.assertAlmostEqual(exact, f_u, places=12)

    def test_the_two_routines_accept_the_same_pedigrees(self):
        """
         inherits R2 through the shared engine. A bound returned on a
        pedigree the exact routine declines would be a bound built on
        unresolved semantics.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                pyp_metrics.a_effective_ancestors_definite(
                    load_fixture(name))
                pyp_metrics.a_effective_ancestors_indefinite(
                    load_fixture(name), n=1)


class TestProductionStateIsolation(unittest.TestCase):
    """
    Completion is analysis-local. Nothing about it may reach the caller's
    pedigree. Scratch analysis state must not remain on caller objects.
    """

    def test_the_callers_pedigree_is_not_mutated(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                before = [dict(a.__dict__) for a in ped.pedigree]
                kw_before = dict(ped.kw)
                pyp_metrics.a_effective_ancestors_definite(ped)
                after = [dict(a.__dict__) for a in ped.pedigree]
                self.assertEqual(before, after)
                self.assertEqual(kw_before, dict(ped.kw))

    def test_no_phantom_joins_the_pedigree(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                n_before = len(ped.pedigree)
                phantoms = pyp_metrics.boichard_phantom_ids(ped)
                pyp_metrics.a_effective_ancestors_definite(ped)
                self.assertEqual(n_before, len(ped.pedigree))
                live = {int(a.animalID) for a in ped.pedigree}
                self.assertEqual(set(), live & set(phantoms))

    def test_repeated_calls_agree(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                first = pyp_metrics.a_effective_ancestors_definite(
                    ped)
                second = pyp_metrics.a_effective_ancestors_definite(
                    ped)
                self.assertEqual(first, second)

    def test_an_earlier_analysis_does_not_change_a_later_one(self):
        ped = load_fixture("A")
        alone = pyp_metrics.a_effective_ancestors_definite(
            load_fixture("A"))
        pyp_metrics.a_effective_founders_boichard(ped)
        pyp_metrics.a_effective_ancestors_indefinite(ped, n=1)
        after = pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertEqual(alone, after)


class TestProductionCrossEstimatorCorroboration(unittest.TestCase):
    """
    Lacy and Boichard founder-number estimators should agree on this
    half-founder fixture:

        "a_effective_founders_boichard and effective_founders_lacy(mode=
        'phantom') continue to agree on half-founder pedigrees (measured 4.0
        vs 4.0 on synth_hf.ped)."

    Two estimators, two papers, two mechanisms: Boichard halves the
    half-founder's own contribution (Appendix A step 4, p.22), Lacy gives the
    unknown parent a founder record (p.113). They agree exactly on a
    half-founder pedigree.

    CORROBORATION, NOT AUTHORITY. Lacy says nothing about Appendix B, and this
    must never be cited as evidence that Reading C is Boichard's text.
    """

    def test_the_two_estimators_agree_on_the_half_founder_fixture(self):
        """Passes today: neither routine has ever refused a half-founder."""
        boichard = pyp_metrics.a_effective_founders_boichard(
            load_corpus("synth_hf.ped", "asdg"))
        lacy = pyp_metrics.effective_founders_lacy(
            load_corpus("synth_hf.ped", "asdg"), mode="phantom")
        self.assertAlmostEqual(4.0, boichard, places=12)
        self.assertAlmostEqual(4.0, lacy["fa_effective_founders"], places=12)
        self.assertEqual(boichard, lacy["fa_effective_founders"])

    def test_the_control_has_no_half_founder_to_disagree_about(self):
        """
        Guard on the guard. If synth_nohf.ped ever grew a half-founder, or
        synth_hf.ped lost one, the pair would stop being a discriminating
        comparison and would keep passing.
        """
        rows, _gens = boichard_read(corpus("synth_hf.ped"), gen_col=3)
        control, _gens2 = boichard_read(corpus("synth_nohf.ped"), gen_col=3)
        self.assertEqual(1, len(boichard_unknown_parent_slots(rows)))
        self.assertEqual([], boichard_unknown_parent_slots(control))

    def test_f_a_does_not_exceed_f_e_on_the_pre_declared_fixture(self):
        f_a = pyp_metrics.a_effective_ancestors_definite(
            load_corpus("synth_hf.ped", "asdg"))
        self.assertLessEqual(f_a, 4.0 + 1e-9)


class TestProductionR1Containment(unittest.TestCase):
    """
     is NOT resolved by this phase and must behave exactly as it did.
    Reading C grows the candidate set on half-founder pedigrees; it does not
    change how a tie among candidates is broken, and lowest_id is still a
    PyPedal convention chosen for determinism rather than Boichard's rule --
    the paper resolves its own three Table-I ties inconsistently by ID.
    """

    def test_the_convention_is_unchanged(self):
        self.assertEqual("lowest_id", pyp_metrics.BOICHARD_TIE_BREAK)

    def test_f_a_is_invariant_under_the_tie_break_convention(self):
        """
        The R1 finding, re-measured on the newly supported domain: the choice
        changes WHICH animal is credited, never the statistic. If it ever
        changed f_a on a half-founder pedigree, R1 would have stopped being a
        free convention the moment R2 was implemented.
        """
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                _text, _rows, reference = FIXTURES[name]
                values = set()
                for convention in ("lowest_id", "highest_id"):
                    order = list(pyp_metrics.boichard_marginal_contributions(
                        load_fixture(name), reference, tie_break=convention))
                    ssq = sum(v * v for _, v in order)
                    values.add(round(1.0 / ssq, 12))
                    self.assertAlmostEqual(1.0, sum(v for _, v in order),
                                           places=12)
                self.assertEqual({round(EXPECTED[name]["f_a"], 12)}, values)

    def test_the_tie_break_is_still_validated(self):
        _text, _rows, reference = FIXTURES["A"]
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            list(pyp_metrics.boichard_marginal_contributions(
                load_fixture("A"), reference, tie_break="oldest"))

    def test_the_real_animal_wins_a_real_versus_phantom_tie(self):
        """
        The concrete R1 consequence of the phantom numbering. On fixture A,
        animal 4 and its own phantom dam both hold p = 0.5 in round 2. Under
        lowest_id the real animal is credited, because every phantom is
        numbered above every real ID.

        Under highest_id the phantom wins instead -- asserted here so that this
        is documented as a property of the CONVENTION rather than a hidden
        preference baked into the completion.
        """
        _text, _rows, reference = FIXTURES["A"]
        phantoms = pyp_metrics.boichard_phantom_ids(load_fixture("A"))
        phantom_id = next(iter(phantoms))

        lowest = [a for a, _ in pyp_metrics.boichard_marginal_contributions(
            load_fixture("A"), reference, tie_break="lowest_id")]
        highest = [a for a, _ in pyp_metrics.boichard_marginal_contributions(
            load_fixture("A"), reference, tie_break="highest_id")]

        self.assertNotIn(phantom_id, lowest)
        self.assertIn(phantom_id, highest)


class TestProductionR3Containment(unittest.TestCase):
    """
     is NOT resolved by this phase either. Reference-population members
    are still excluded from selection, the antichain guard still fires, and the
    reference population is still real animals only -- a phantom is a parent of
    an existing animal by construction, so it can never be a member.
    """

    def test_the_antichain_guard_still_fires_on_a_half_founder_pedigree(self):
        """
        The guard has to survive the newly supported domain. If completion had
        been placed before the R3 check, a half-founder pedigree with a
        non-antichain reference population would have started returning a
        number instead of refusing.
        """
        ped = load_corpus_from_path(write_pedigree(
            "1 0 0 1\n2 0 0 1\n3 1 0 1\n4 1 2 1\n5 3 4 2\n6 5 4 2\n"), "asdg")
        with self.assertRaises(pyp_errors.PyPedalError) as ctx:
            pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertIn("antichain", str(ctx.exception))

    def test_a_phantom_is_never_in_the_reference_population(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                _text, _rows, reference = FIXTURES[name]
                phantoms = pyp_metrics.boichard_phantom_ids(ped)
                self.assertEqual(set(), set(reference) & set(phantoms))

    def test_reference_members_are_still_excluded_from_selection(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                _text, _rows, reference = FIXTURES[name]
                order = list(pyp_metrics.boichard_marginal_contributions(
                    load_fixture(name), reference))
                self.assertEqual(set(), {a for a, _ in order} & set(reference))


class TestProductionReporting(unittest.TestCase):
    """
    A synthetic ID must never reach a human reader looking like a real animal.
    The report keeps the integer, because the structure is a list of IDs and
    dropping it would make the file untraceable, but never bare.
    """

    @staticmethod
    def _report(name, suffix):
        ped = load_fixture(name)
        if suffix == "definite":
            pyp_metrics.a_effective_ancestors_definite(ped)
        else:
            pyp_metrics.a_effective_ancestors_indefinite(ped, n=1)
        path = f"{ped.kw['filetag']}_fa_boichard_{suffix}_.dat"
        with open(path, encoding="utf-8") as fh:
            return fh.read()

    def test_a_selected_phantom_is_labelled_not_printed_bare(self):
        """
        Fixture D is the one where a phantom is actually selected, so this
        checks the labelling on a report that really contains one.
        """
        report = self._report("D", "definite")
        phantoms = pyp_metrics.boichard_phantom_ids(load_fixture("D"))
        selected = [a for a, _ in pyp_metrics.boichard_marginal_contributions(
            load_fixture("D"), FIXTURES["D"][2]) if a in phantoms]
        self.assertEqual(1, len(selected))
        self.assertIn("=phantom (", report)
        self.assertIn(f"{selected[0]}=phantom (dam of 4)", report)

    def test_both_reports_explain_the_phantoms(self):
        for suffix in ("definite", "indefinite"):
            with self.subTest(report=suffix):
                report = self._report("A", suffix)
                self.assertIn("phantom founder(s) were created", report)
                self.assertIn("NOT in the pedigree", report)
                self.assertIn("mathematically implied", report)

    def test_the_note_does_not_claim_appendix_b_says_this(self):
        """
        The evidence class is load-bearing and the report is where a reader
        without the docs will meet it. It must not read as though Appendix B
        states the rule.
        """
        report = self._report("A", "definite")
        self.assertIn("silent", report)
        self.assertNotIn("Appendix B states", report)
        self.assertNotIn("source-explicit", report)

    def test_a_half_founder_free_report_gains_no_note(self):
        """Non-vacuity in the other direction: the note is conditional."""
        report = self._report("NONE", "definite")
        self.assertNotIn("phantom", report)


class TestProductionDeepStateIsolation(unittest.TestCase):
    """
    Stronger than comparing parentage: the full attribute dictionary of every
    animal, the attribute key sets, kw, pedigree membership and order.
    """

    ROUTINES = ("a_effective_ancestors_definite",
                "a_effective_founders_boichard")

    @staticmethod
    def _snapshot(ped):
        return ([{k: repr(v) for k, v in sorted(a.__dict__.items())}
                 for a in ped.pedigree],
                {k: repr(v) for k, v in sorted(ped.kw.items())})

    def test_nothing_about_the_caller_changes(self):
        for name in HALF_FOUNDER_FIXTURES:
            for routine in self.ROUTINES:
                with self.subTest(fixture=name, routine=routine):
                    ped = load_fixture(name)
                    before = self._snapshot(ped)
                    getattr(pyp_metrics, routine)(ped)
                    self.assertEqual(before, self._snapshot(ped))

    def test_the_bounded_routine_leaves_the_caller_alone_too(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                before = self._snapshot(ped)
                pyp_metrics.a_effective_ancestors_indefinite(ped, n=1)
                self.assertEqual(before, self._snapshot(ped))

    def test_no_undeclared_attribute_is_added(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                keys = [set(a.__dict__) for a in ped.pedigree]
                pyp_metrics.a_effective_ancestors_definite(ped)
                self.assertEqual(keys, [set(a.__dict__) for a in ped.pedigree])

    def test_the_snapshot_would_notice_a_change(self):
        """
        Guard on the guard. A snapshot comparison that cannot fail proves
        nothing, so mutate one animal by hand and confirm it is detected.
        """
        ped = load_fixture("A")
        before = self._snapshot(ped)
        ped.pedigree[3].damID = 99
        self.assertNotEqual(before, self._snapshot(ped))

    def test_phantom_ids_are_stable_across_calls(self):
        for name in HALF_FOUNDER_FIXTURES:
            with self.subTest(fixture=name):
                ped = load_fixture(name)
                first = pyp_metrics.boichard_phantom_ids(ped)
                pyp_metrics.a_effective_ancestors_definite(ped)
                self.assertEqual(first, pyp_metrics.boichard_phantom_ids(ped))


class TestProductionContainmentOnTheOldDomain(unittest.TestCase):
    """
    Completion has nothing to do without half-founders, so the previously
    supported domain must be untouched. This is the in-suite half of
    containment; the cross-checkout half is
    the independent oracle
    """

    def test_the_control_pedigree_is_accepted_before_and_after(self):
        """Green today and must stay green -- not an xfail."""
        got = pyp_metrics.a_effective_ancestors_definite(
            load_fixture("NONE"))
        rows, reference = oracle_rows("NONE")
        want, _o, _t = boichard_f_a(rows, reference)
        self.assertAlmostEqual(want, got, places=12)

    def test_the_control_has_no_half_founders(self):
        """
        Guard on the guard: if this fixture ever acquires a half-founder the
        containment test above silently stops testing containment.
        """
        rows, _reference = oracle_rows("NONE")
        self.assertEqual([], boichard_unknown_parent_slots(rows))


if __name__ == "__main__":
    sys.exit(unittest.main())
