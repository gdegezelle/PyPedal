"""
Class C -- KNOWN-DEFECT characterisation tests.

Each test here asserts the **intended correct** behaviour of something that is
currently broken, and is marked expected-to-fail. That is deliberate: it means
the repository states, executably, what the right answer is before anyone tries
to produce it, so a fix cannot be declared complete by adjusting the
expectation.

Blocked entries use ``@pytest.mark.xfail(strict=True, reason=...)`` rather than
``@unittest.expectedFailure``. The reason is then visible in the test report
instead of only in a docstring, and **strict** means an unexpected PASS fails
the suite -- so a blocked item cannot turn green without someone noticing, which
is the whole point of the category.

Two kinds of entry live here, and they have different futures:

**Repairable.** The correct target value is established by independent evidence.
The commit that performs the repair removes the ``expectedFailure`` marker and
moves the assertion into ``tests/test_correctness_invariants.py`` in that same
commit. A repairable entry must not survive an earlier revision.

**Blocked.** The invariant proves the current output is impossible, but no
independent evidence establishes what the correct value *is*. These stay
expected failures until the reference that would adjudicate them is in hand. They are never deleted, never skipped, and
never "fixed" by making the algorithm satisfy the invariant -- doing that would
be guessing, and would replace a visibly wrong number with an invisibly wrong
one. See the blocked-items table in ``docs/python3-migration-audit.md``.
"""
import unittest

import pytest

from PyPedal import pyp_errors, pyp_nrm, pyp_metrics, pyp_validate

from _pedhelpers import corpus, load_corpus
from oracles import oracle_inbreeding


def coi(ped, method):
    result = pyp_nrm.inbreeding(ped, method=method, output=False)
    if isinstance(result, tuple):
        result = result[0]
    return {int(k): float(v) for k, v in result["fx"].items()}


# ---------------------------------------------------------------------------
# Repairable -- promoted to class A by the commit named on each test
# ---------------------------------------------------------------------------

#
#  (modified Meuwissen-Luo duplicate guard) was repaired in Commit 4.
# Its assertions were promoted to tests/test_correctness_invariants.py.
#


#
#  (gene-drop replicate normalisation) was repaired in Commit 5. Its
# assertions were promoted to tests/test_correctness_invariants.py.
#


#
#  (string-ID pedigrees) was repaired in Commit 3. Its assertions were
# promoted to tests/test_correctness_invariants.py (the load round-trip
# invariant and collision detection) and tests/test_legacy_numeric_compatibility.py (the
# record counts, which rest on agreement with PyPedal 2.0.4 rather than on an
# independent oracle).
#


#
#  (Mendelian transmission probability) was repaired in Commit 6. Its
# assertions were promoted to tests/test_correctness_invariants.py.
#


# ---------------------------------------------------------------------------
# Blocked -- expected failures for all of an earlier revision, by design
# ---------------------------------------------------------------------------

class TestBlockedEffectiveAncestors(unittest.TestCase):
    """
     / , audit Findings 13 and 17 -- **both now repaired against the
    source**, so this class no longer holds an expected failure.

    ``f_a`` is the reciprocal of a sum of squared probabilities, so ``f_a >= 1``
    unconditionally: a population cannot have fewer than one effective ancestor.
    Both routines used to violate it. ``a_effective_ancestors_definite``
    returned 0.0816, then 0.0625 after the index-domain repair;
    ``a_effective_ancestors_indefinite`` returned 0.6514, and PyPedal 2.0.4
    returns NaN for its lower bound.

    The cause was established in an earlier revision stage A0 and confirmed by the paper:
    the port had no ``a`` vector and therefore no ``q*(1-a)``
    marginal-contribution step, and its update pass traversed oldest-to-youngest
    where Appendix B step 5 requires youngest-to-oldest.

    Both are fixed, so the assertion below now PASSES and is kept as a
    regression test rather than deleted. The value it guards is no longer
    hypothetical: on ``boichard2a.ped`` the routine returns 2.0, which is what
    PyPedal 2.0.4 returns and what the definition gives by hand.
    """

    def test_definite_f_a_is_at_least_one(self):
        ped = load_corpus("boichard2a.ped")
        # precondition: the routine documents that generations must be defined
        generations = {getattr(a, "gen", None) for a in ped.pedigree}
        self.assertGreater(len(generations), 1,
                           "precondition not met; this would not be a valid test")
        f_a = pyp_metrics.a_effective_ancestors_definite(ped)
        self.assertGreaterEqual(float(f_a), 1.0)

    def test_definite_matches_the_hand_derived_value(self):
        """
        2.0 is not merely >= 1. It is the value hand-derived from the definition
        (ancestors 5 and 6 each explain half of the reference generation, so
        sum(p^2) = 0.5) and the value PyPedal 2.0.4 returns. Asserting the
        number, not just the bound, is what stops a future change from
        satisfying the invariant while getting the answer wrong.
        """
        ped = load_corpus("boichard2a.ped")
        self.assertAlmostEqual(
            2.0, float(pyp_metrics.a_effective_ancestors_definite(ped)), places=9)

    #
    # 's expected failure went first, in an earlier revision stage A3, when the routine
    # was made to REFUSE rather than return a number it could not compute. The
    # xfail count went 2 -> 1 then, and 1 -> 0 here. Neither was a blocked item
    # being quietly greened: A3 replaced a wrong value with an explicit refusal,
    # and this commit replaced the refusal's cause with the paper's algorithm.
    #
    # What guards the refusal itself is TestBlockedIndefiniteRefusesToRun in
    # tests/test_correctness_invariants.py.
    #


#
#  (Lacy on half-founder pedigrees) deliberately has NO expected failure.
#
# An earlier revision of this file added one, asserting
# f_e <= count(founder == 'y') for a_effective_founders_lacy, on the grounds
# that f_e = 1/sum(q_k^2) over a probability vector of length k is bounded above
# by k. That assertion was WRONG and has been withdrawn.
#
# The bound is valid only when k is the number of independent founder
# contribution sources the model represents. On a pedigree with half-founders
# the unknown parental side is an additional such source, so
# count(founder == 'y') is not the right k. The oracle demonstrates it: the
# phantom-completed mode is the only one whose q is a genuine probability
# vector, and it too exceeds count(founder == 'y') --
#
#     mrode.ped          count(founder=='y') = 2, k with phantoms = 3, f_e = 2.798
#     hartlandclark.ped  count(founder=='y') = 3, k with phantoms = 7, f_e = 5.832
#
# -- while satisfying f_e <= k against the phantom-completed k. So the values
# a_effective_founders_lacy returns (3.230 and 7.205) are not shown to be
# impossible. They are computed from a q that does not sum to 1, which makes
# them ill-defined rather than provably out of range, and no reference-free
# invariant refutes them.
#
# Establishing the correct founder-source count, and only then whether an upper
# bound applies, requires Lacy (1989). Until that is in hand there is nothing
# here that can honestly be asserted as intended behaviour, so nothing is
# asserted. The discriminating measurements live as passing oracle tests in
# tests/test_correctness_invariants.py::TestLacyOracleHalfFounderHandling.
#

#
#  (the effective-ancestor routines mutating the caller's pedigree)
# was repaired as a state-integrity fix. Its assertions were promoted to
# tests/test_correctness_invariants.py. The MATHEMATICS of those routines
# remains blocked -- see TestBlockedEffectiveAncestors above, whose expected
# failures are unaffected by that repair.
#


#
#  (half-founders corrupting the gene-drop parent lookup) was
# repaired by the refusal that closes its index defect. Its assertions were
# promoted to tests/test_genedrop_half_founders.py.
#
# It is worth recording WHY the expected failure could be removed while the
# science it guards stays blocked, because the two are usually the same thing
# here and in this case are not. The entry asserted a TYPED REFUSAL, not a
# value: the corruption mechanism (a sentinel reaching pedigree[-1], and ""
# == "" reading as autozygosity) is closed, while the scientific question --
# what gene source a missing parental side supplies under gene dropping -- is
# untouched and remains blocked as b/HF. Removing the marker therefore
# records a repair that happened, not a block that was quietly greened.
#
# The 0.0 once recorded as animal 4's correct value is NOT the target and was
# never asserted. See the note in tests/test_genedrop_half_founders.py.
#


if __name__ == "__main__":
    unittest.main()
