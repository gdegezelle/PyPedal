###############################################################################
# NAME: pyp_validate.py
# VERSION: see PyPedal.__version__
# LICENSE: LGPL
# Written for PyPedal 4.0; this module is not part of the original 2.0.4 sources.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
"""
Scientific validation for PyPedal results.

Why this module exists
----------------------
Several historical defects returned plausible-looking but wrong scientific
results, silently. Range checks on quantities with known mathematical bounds
catch those immediately: an ancestral inbreeding coefficient of 20 where the
quantity is a probability, or an effective ancestor count of 0.08 where the
quantity is the reciprocal of a sum of squared probabilities.

The controlling principle is that returning a plausible wrong number is worse
than raising. A caller who gets an exception knows the analysis failed. A caller
who gets 0.0816 for an effective ancestor count may publish it.

Two layers, deliberately separated
----------------------------------
**(a) Runtime**, enabled by ``kw['validate']`` (default True). Cheap checks on
values the calling routine has *already computed*: preconditions on its inputs,
finiteness, and scalar mathematical bounds. Everything here is O(1) or O(n) over
data already in hand.

This layer must never recompute an independent algorithm and must never change
a complexity class. PyPedal must stay practical on pedigrees of >100,000
animals, and a validation layer that turned an O(n) path into O(n^2), or that
ran a second COI implementation to cross-check the first, would destroy that for
no correctness gain that a test could not provide more cheaply.

**(b) Test and diagnostic only.** Cross-method agreement between COI algorithms,
full O(n^2) NRM symmetry and diagonal scans, Monte Carlo variance experiments,
and statistical Mendelian-transmission experiments all live in ``tests/``.
Full-matrix checking is available at runtime only through
``kw['diagnostics']``, which is off by default and is not implied by
``validate``.

What this module is NOT
-----------------------
A guard here is not a repair. ``check_effective_ancestors`` rejects an f_a
below 1. The Boichard Appendix B routines now implement the paper's
marginal-contribution step; the invariant still proves an output below 1 is
impossible. Refusing to return an impossible number is all this layer claims
to do.
"""

import logging
import math

from .pyp_errors import PyPedalError, PyPedalValidationError  # noqa: F401


def _enabled(pedobj, key="validate", default=True):
    """Read a validation flag off a pedigree, tolerating a bare options dict."""
    if pedobj is None:
        return default
    kw = getattr(pedobj, "kw", pedobj)
    try:
        return bool(kw.get(key, default))
    except AttributeError:
        return default


def _fail(pedobj, message):
    logging.error(message)
    raise PyPedalValidationError(message)


# ---------------------------------------------------------------------------
# Finiteness
# ---------------------------------------------------------------------------

def check_finite(pedobj, value, name):
    """
    Reject NaN and infinity. Both propagate silently through arithmetic and
    through most comparisons, so an unchecked NaN can travel a long way from
    where it was produced before anyone notices.
    """
    if not _enabled(pedobj):
        return value
    try:
        as_float = float(value)
    except (TypeError, ValueError):
        _fail(pedobj, "%s is not a number: %r" % (name, value))
    if math.isnan(as_float):
        _fail(pedobj, "%s is NaN" % name)
    if math.isinf(as_float):
        _fail(pedobj, "%s is infinite" % name)
    return value


# ---------------------------------------------------------------------------
# Scalar range postconditions
# ---------------------------------------------------------------------------

def check_probability(pedobj, value, name, tolerance=1e-9):
    """
    Assert 0 <= value <= 1 for a quantity that is a probability.

    Applies to pedigree-based coefficients: F is the probability that an
    animal's two alleles at a locus are identical by descent, and F_a is the
    probability of having inherited an allele that was once inbred.

    It does NOT apply to genomic estimators. A SNP-based inbreeding coefficient
    is a different quantity, defined relative to a reference allele frequency,
    and is legitimately negative when an animal is less homozygous than that
    reference. Do not call this on ``pyp_snp`` output.
    """
    if not _enabled(pedobj):
        return value
    check_finite(pedobj, value, name)
    as_float = float(value)
    if as_float < -tolerance or as_float > 1.0 + tolerance:
        _fail(pedobj, "%s = %r is outside [0, 1]; this quantity is a "
                      "probability, so the result is impossible" % (name, as_float))
    return value


def check_effective_number(pedobj, value, name, n_contributors=None):
    """
    Assert 1 <= value <= n_contributors for an *effective number* of founders
    or ancestors.

    Both bounds come from the same fact: the contributions q_k form a
    probability vector, non-negative and summing to 1, and the statistic is
    1 / sum(q_k^2).

      Lower bound. sum(q_k^2) <= 1 for any such vector, so f >= 1. A population
      cannot have fewer than one effective ancestor.

      Upper bound. For k contributors, sum(q_k^2) is minimised when the
      contributions are equal, giving sum(q_k^2) >= 1/k, so f <= k. The
      effective number of founders cannot exceed the actual number of
      founders, and equals it only when every founder contributes equally.

    The upper bound is passed explicitly rather than inferred, because it is
    only valid when `n_contributors` really is the size of the contribution
    vector the statistic was computed from. Callers that cannot establish that
    should omit it and get the lower bound alone.

    This is a fail-loud guard, not a repair. See the module docstring.
    """
    if not _enabled(pedobj):
        return value
    check_finite(pedobj, value, name)
    as_float = float(value)
    if as_float < 1.0:
        _fail(pedobj, "%s = %r is less than 1; an effective number of founders "
                      "or ancestors is the reciprocal of a sum of squared "
                      "probabilities and cannot be below 1. The value is "
                      "mathematically impossible and has been refused rather "
                      "than returned." % (name, as_float))
    if n_contributors is not None and as_float > float(n_contributors) + 1e-9:
        _fail(pedobj, "%s = %r exceeds the number of contributors (%d). The "
                      "contributions form a probability vector, so their "
                      "squared sum is at least 1/k and the effective number is "
                      "at most k, reached only when all contribute equally. "
                      "The value is mathematically impossible and has been "
                      "refused rather than returned."
                      % (name, as_float, n_contributors))
    return value


def check_contribution_vector(pedobj, contributions, name, tolerance=1e-6):
    """
    Assert that an already-computed contribution vector is a probability
    vector: non-negative and summing to 1.

    Only call this where the routine has computed the vector for its own
    purposes. Building one here in order to check it would be recomputing the
    algorithm, which this layer must not do.

    Note the sum-to-1 check is only meaningful on a vector that has *not* been
    normalised by its own total. Normalising makes it true by construction and
    the check vacuous.
    """
    if not _enabled(pedobj):
        return contributions
    values = list(contributions.values() if hasattr(contributions, "values")
                  else contributions)
    for index, value in enumerate(values):
        check_finite(pedobj, value, "%s[%d]" % (name, index))
        if float(value) < -tolerance:
            _fail(pedobj, "%s[%d] = %r is negative; contributions are "
                          "probabilities" % (name, index, value))
    total = sum(float(v) for v in values)
    if values and abs(total - 1.0) > tolerance:
        _fail(pedobj, "%s sums to %r, not 1; it is not a valid probability "
                      "vector" % (name, total))
    return contributions


# ---------------------------------------------------------------------------
# Preconditions on inputs
# ---------------------------------------------------------------------------

def require_generations(pedobj, routine):
    """
    Refuse to run a metric that documents "assumes that generations are
    defined" when the pedigree does not in fact define them.

    Previously such a call returned a number computed from missing-value
    sentinels. Refusing is a guard on the input, not a change to the algorithm.
    """
    if not _enabled(pedobj):
        return
    missing = pedobj.kw.get("missing_gen", -999.0)
    real = set()
    for animal in pedobj.pedigree:
        gen = getattr(animal, "gen", None)
        if gen is None:
            continue
        try:
            if float(gen) != float(missing):
                real.add(float(gen))
        except (TypeError, ValueError):
            continue
    if len(real) < 2:
        _fail(pedobj, "%s requires generations to be defined, but the pedigree "
                      "has %d distinct generation value(s). Set kw['set_generations'] "
                      "or supply a 'g' column in the pedigree format." % (routine, len(real)))


def require_renumbered(pedobj, routine):
    """
    Assert the ``index == animalID - 1`` contract that most of ``pyp_nrm`` and
    ``pyp_metrics`` silently assumes. A violation produces wrong answers by
    reading the wrong animal, with nothing to indicate it.
    """
    if not _enabled(pedobj):
        return
    for index, animal in enumerate(pedobj.pedigree):
        if int(animal.animalID) != index + 1:
            _fail(pedobj, "%s requires a renumbered pedigree (index == animalID - 1), "
                          "but pedigree[%d] has animalID %s. Load with "
                          "kw['renumber'] = True." % (routine, index, animal.animalID))


# ---------------------------------------------------------------------------
# Aggregate helpers used at the call sites
# ---------------------------------------------------------------------------

def check_inbreeding_coefficients(pedobj, fx, routine):
    """
    Validate a whole COI result dictionary. O(n) over values the caller already
    computed -- it does not recompute anything and does not touch the NRM.
    """
    if not _enabled(pedobj):
        return fx
    for animal_id, value in fx.items():
        check_probability(pedobj, value, "%s: F(%s)" % (routine, animal_id))
    return fx


def check_genomic_inbreeding(pedobj, values, routine):
    """
    Postcondition for VanRaden genomic inbreeding: finite, and ``>= -1``.

    **Deliberately NOT the pedigree-inbreeding contract.** ``F_g = G_jj - 1``
    where ``G_jj = ||z_j||^2 / (2 sum p(1-p)) >= 0``, so the lower bound is -1
    and is ATTAINABLE -- an individual heterozygous at every locus sits exactly
    on it. There is **no generic finite upper bound**: as an allele becomes rare
    the numerator tends to 4 while its own denominator term tends to 0, so an
    individual homozygous for rare alleles can have an arbitrarily large
    coefficient. VanRaden (2008) p.4416 says as much in words: "the genomic
    inbreeding coefficient is greater if the individual is homozygous for rare
    alleles than if homozygous for common alleles."

    Reusing ``check_inbreeding_coefficients`` here would therefore fire on
    perfectly legitimate values. The two quantities share a name and not a
    range.
    """
    if not _enabled(pedobj):
        return values
    for animal_id, value in values.items():
        name = '%s: F_g for animal %s' % (routine, animal_id)
        check_finite(pedobj, value, name)
        if float(value) < -1.0 - 1e-9:
            _fail(pedobj,
                  '%s = %r is below -1. F_g = G_jj - 1 and G_jj is a squared '
                  'norm divided by a positive scalar, so it cannot be negative '
                  'and F_g cannot be below -1.' % (name, float(value)))
    return values


def check_ancestral_inbreeding(pedobj, id2aic, routine):
    """
    Validate an ancestral-inbreeding result dictionary. O(n), no recomputation.

    Applies to BOTH ancestral-inbreeding routines, which are two different
    estimators of the same quantity: Ballou's deterministic recursion
    (``ballou_ancestral_inbreeding``) and the gene-dropping estimator
    attributable to Suwanlee et al. 2007 (``dropped_ancestral_inbreeding``).
    F_a is a probability under either, so the [0, 1] bound is legitimate for
    both. Nothing here compares the two against each other -- they have
    different mathematical properties and their agreement is not an invariant.
    """
    if not _enabled(pedobj):
        return id2aic
    for animal_id, value in id2aic.items():
        check_probability(pedobj, value, "%s: F_a(%s)" % (routine, animal_id))
    return id2aic
