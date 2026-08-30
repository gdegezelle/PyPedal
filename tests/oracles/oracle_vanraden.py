#!/usr/bin/env python3
"""
INDEPENDENT oracle for VanRaden (2008) Method 1.
Does not import PyPedal.

VanRaden PM (2008) "Efficient Methods to Compute Genomic Predictions",
*J Dairy Sci* 91:4414-4423, article p.4416:

    M   n individuals x m loci, elements -1, 0, 1 for the homozygote, the
        heterozygote and the other homozygote
    p_i frequency of the SECOND allele at locus i
    P   column i is 2(p_i - 0.5)
    Z   = M - P, "which sets mean values of the allele effects to 0"
    G   = ZZ' / (2 * sum(p_i (1 - p_i)))      <- Method 1
        "Division by 2*sum(p_i(1-p_i)) scales G to be analogous to the
         numerator relationship matrix A."
    F_g "The genomic inbreeding coefficient for individual j is simply G_jj - 1"

Genotypes are supplied here as ALLELE COUNTS 0/1/2, which is how PyPedal stores
them; ``M = counts - 1`` recovers the paper's coding, so ``Z = counts - 2p``.

WHAT THIS ORACLE CANNOT DO, STATED PLAINLY
------------------------------------------
**The paper publishes no worked G.** There is no table of numbers to reproduce,
as there is for Boichard, Lacy and Ballou. This oracle is therefore
*definition-derived*: it encodes the equations from p.4416 and is checked
against a fixture small enough to verify by hand, plus the algebraic properties
the equations imply. That is a weaker footing than the pedigree oracles have,
and it is recorded rather than glossed.

ALLELE FREQUENCIES ARE AN ARGUMENT, NOT AN ESTIMATE
---------------------------------------------------
p.4416: "Allele frequencies in P should be from the unselected base population
rather than from those that appear after selection or inbreeding. An earlier or
later base population can lead to greater or fewer relationships and to more or
less inbreeding." p.4419 adds that "frequency estimation can bias the genomic
inbreeding coefficients".

So ``p`` is required here. Estimating it from the sample is a modelling choice
the paper explicitly flags, and an oracle that made that choice silently could
not be used to test whether production makes it correctly.

THE MONOMORPHIC CASE IS ALGEBRA, NOT SEMANTICS
----------------------------------------------
A locus with p in {0, 1} contributes 0 to ZZ' and 0 to the denominator, so under
Method 1 it is harmless. If EVERY locus is monomorphic the denominator is 0 and
Method 1 is undefined -- not "undefined by convention", undefined. The paper is
silent because there is nothing to say. This raises.
"""
import argparse
import json

import numpy as np


def z_matrix(counts, p):
    """
    Z = M - P with M the paper's -1/0/1 coding and P column i = 2(p_i - 0.5).

    Equivalently ``counts - 2p`` for 0/1/2 allele counts, since
    ``(counts - 1) - (2p - 1) = counts - 2p``.
    """
    counts = np.asarray(counts, dtype=float)
    p = np.asarray(p, dtype=float)
    if counts.ndim != 2:
        raise ValueError("counts must be a 2-D array of individuals x loci")
    if p.shape != (counts.shape[1],):
        raise ValueError(
            f"expected exactly one allele frequency per locus: got {p.shape} "
            f"for {counts.shape[1]} loci")
    if not np.all(np.isfinite(p)):
        raise ValueError("allele frequencies must all be finite")
    if np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("allele frequencies must lie in [0, 1]")
    return counts - 2.0 * p


def scaling_denominator(p):
    """2 * sum(p_i (1 - p_i)). Zero exactly when every locus is monomorphic."""
    p = np.asarray(p, dtype=float)
    return float(2.0 * np.sum(p * (1.0 - p)))


def grm_method_1(counts, p):
    """G = ZZ' / (2 sum p(1-p)), VanRaden Method 1."""
    z = z_matrix(counts, p)
    denominator = scaling_denominator(p)
    if denominator <= 0.0:
        raise ZeroDivisionError(
            "2*sum(p(1-p)) is zero: every locus is monomorphic, so VanRaden "
            "Method 1 is mathematically undefined for this input. This is "
            "algebra, not a convention -- the paper is silent because there is "
            "nothing to state.")
    return (z @ z.T) / denominator


def genomic_inbreeding(counts, p):
    """F_g(j) = G_jj - 1, article p.4416."""
    return np.diag(grm_method_1(counts, p)) - 1.0


def wright_relationships(counts, p):
    """
    Off-diagonals on Wright's scale: G_jk / sqrt(G_jj * G_kk), p.4416.

    Undefined where a diagonal is zero -- an individual carrying no deviation
    from the base population at any locus -- so those entries are NaN rather
    than silently zero.
    """
    g = grm_method_1(counts, p)
    d = np.sqrt(np.diag(g))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.outer(d, d) > 0, g / np.outer(d, d), np.nan)


# --------------------------------------------------------------------------
# The hand fixture
# --------------------------------------------------------------------------

#: Three individuals, four loci, every frequency 0.5, so P is all ones and
#: Z = counts - 1. Chosen so that every entry of G can be checked mentally:
#:
#:     Z = [[-1, 0,  1, 0],
#:          [ 1, 0, -1, 0],
#:          [ 0, 0,  0, 0]]
#:     denominator = 2 * 4 * 0.25 = 2
#:     G = [[1, -1, 0], [-1, 1, 0], [0, 0, 0]]
#:     F_g = [0, 0, -1]
#:
#: Individual 3 is heterozygous everywhere: it carries no deviation from the
#: base population, so G_33 = 0 and F_g = -1, the lower bound exactly.
HAND_COUNTS = [[0, 1, 2, 1],
               [2, 1, 0, 1],
               [1, 1, 1, 1]]
HAND_P = [0.5, 0.5, 0.5, 0.5]
HAND_G = [[1.0, -1.0, 0.0],
          [-1.0, 1.0, 0.0],
          [0.0, 0.0, 0.0]]
HAND_F_G = [0.0, 0.0, -1.0]

#: One locus, a rare second allele, one individual homozygous for it. Shows
#: that F_g has no generic finite upper bound: as p -> 0 the numerator tends to
#: 4 while its own denominator term tends to 0.
#:
#:     Z = 2 - 2(0.01) = 1.98;  denominator = 2(0.01)(0.99) = 0.0198
#:     G = 1.98^2 / 0.0198 = 198;  F_g = 197
RARE_COUNTS = [[2]]
RARE_P = [0.01]
RARE_F_G = 197.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--hand", action="store_true",
                    help="print the hand fixture and its expected values")
    args = ap.parse_args()
    if args.hand:
        g = grm_method_1(HAND_COUNTS, HAND_P)
        print(json.dumps({
            "counts": HAND_COUNTS, "p": HAND_P,
            "denominator": scaling_denominator(HAND_P),
            "G": g.tolist(), "G_expected": HAND_G,
            "F_g": genomic_inbreeding(HAND_COUNTS, HAND_P).tolist(),
            "F_g_expected": HAND_F_G,
        }, indent=2, sort_keys=True))
        return
    ap.error("give --hand; this oracle has no pedigree input")


if __name__ == "__main__":
    main()
