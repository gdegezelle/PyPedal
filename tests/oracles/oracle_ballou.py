#!/usr/bin/env python3
"""
INDEPENDENT oracle for Ballou's ancestral inbreeding coefficient.
Does not import PyPedal. Reuses the independently cross-validated inbreeding
coefficients from oracle_meuwissen_luo.

    f_a(x) = [ f_a(s) + (1 - f_a(s)) * F_s
             + f_a(d) + (1 - f_a(d)) * F_d ] / 2

where F is the ordinary coefficient of inbreeding, s and d are the sire and dam,
and an unknown parent contributes f_a = 0 and F = 0. Animals are processed
oldest to youngest, so both parents' values are settled before a child's.

f_a is the probability that an individual has inherited an allele that has
undergone inbreeding at least once in the past. Note the consequence, which the
fixtures below exercise deliberately: the FIRST inbred animal in a line has
f_a = 0. Its own inbreeding is not "in the past" from its own point of view --
what matters is whether its PARENTS were inbred.

PROVENANCE -- PAPER-VALIDATED
----------------------------
The recurrence is **equation 1 of Ballou JD (1997), "Ancestral inbreeding only
minimally affects inbreeding depression in mammalian populations", J Hered
88:169-178, article p.170**, read directly:

    f_a = [ f_a(s) + (1 - f_a(s)) f_(s) + f_a(d) + (1 - f_a(d)) f_(d) ] / 2

with the paper's own gloss: "An individual's f_a is then the proportion of its
parent's genome that has been previously exposed to inbreeding (f_a of the
parent) plus the effect of the parent's inbreeding coefficient on the
proportion that has not been previously exposed (1 - f_a of parent), averaged
across both parents; it ranges from 0 to 1."

Independently corroborated twice more:

* Ballou Figure 1 (p.170) prints f and f_a beside each node of a small
  pedigree, and all of them reproduce -- see corpus/ballou_fig1.ped.
* Suwanlee et al. (2007) p.490 restates the same equation, from a different
  group, as their eq. (1).

The founder base case, f_a = 0, is read off Ballou Figure 1 rather than stated
in prose. That remains an inference, though a narrow one.

WHAT THIS ORACLE MUST NOT BE USED FOR
-------------------------------------
**Do not compare it against ``pyp_metrics.dropped_ancestral_inbreeding``, and do
not treat agreement or disagreement with that routine as a correctness signal.**
This was a caution before; it is now a published result. Suwanlee et al. (2007)
Figure 1 prints Ballou's coefficient and a gene-drop estimate SIDE BY SIDE on
one pedigree, and they differ from the third generation onward:

    generation      f_a (Ballou)    f_a (gene drop, 10^6 reps)
    15, 16              0.399              0.375
    17-20               0.662              0.585
    21, 22              0.847              0.728
    23                  0.944              0.822

Ballou's formula assumes f_a and f are independent (Suwanlee eq. 4); they show
the assumption does not hold. A systematic gap is the EXPECTED result.

Nor may the reverse become an invariant. Suwanlee proves no theorem: it reports
Ballou overestimating in the scenarios it simulated, and separately reports gene
dropping itself marginally overestimating under lethal-allele models. Use the
published rows as regression expectations, not "gene drop < Ballou" as a law.
"""
import argparse, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_meuwissen_luo import read_pedigree, renumber, inbreeding_tabular


def ancestral_inbreeding(ped, F=None):
    """
    Ballou's f_a for a topologically ordered pedigree of (id, sire, dam) triples
    with ids 1..n. ``F`` may be supplied to avoid recomputing; otherwise the
    independently cross-validated tabular coefficients are used.

    Returns a list indexed 0..n, with index 0 unused, matching ``F``'s shape.
    """
    if F is None:
        F = inbreeding_tabular(ped)
    n = len(ped)
    f_a = [0.0] * (n + 1)
    for i, s, d in ped:
        f_as = f_a[s] if s else 0.0
        f_s = F[s] if s else 0.0
        f_ad = f_a[d] if d else 0.0
        f_d = F[d] if d else 0.0
        f_a[i] = (f_as + (1.0 - f_as) * f_s + f_ad + (1.0 - f_ad) * f_d) / 2.0
    return f_a


def ancestral_inbreeding_from_file(path, sepchar=" ", animal=0, sire=1, dam=2):
    """Returns ``{original_id: f_a}``, keyed by the ids in the file."""
    rows = read_pedigree(path, sepchar=sepchar, animal=animal, sire=sire, dam=dam)
    ped, back = renumber(rows)
    f_a = ancestral_inbreeding(ped)
    return {back[i]: f_a[i] for i in range(1, len(ped) + 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pedfile")
    ap.add_argument("--animal", type=int, default=0)
    ap.add_argument("--sire", type=int, default=1)
    ap.add_argument("--dam", type=int, default=2)
    a = ap.parse_args()
    rows = read_pedigree(a.pedfile, animal=a.animal, sire=a.sire, dam=a.dam)
    ped, back = renumber(rows)
    F = inbreeding_tabular(ped)
    f_a = ancestral_inbreeding(ped, F)
    n = len(ped)
    print(json.dumps({
        "pedfile": os.path.basename(a.pedfile),
        "n_animals": n,
        "provenance": "methodology.tex:88 (secondary); Ballou (1997) NOT obtained",
        "f": {str(back[i]): F[i] for i in range(1, n + 1)},
        "f_a": {str(back[i]): f_a[i] for i in range(1, n + 1)},
        "f_a_max": max(f_a[1:]) if n else 0.0,
        "f_a_in_unit_range": all(0.0 <= v <= 1.0 for v in f_a[1:]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
