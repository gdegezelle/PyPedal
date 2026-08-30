#!/usr/bin/env python3
"""
INDEPENDENT oracle for the gene-dropping ancestral inbreeding
coefficient of Suwanlee, Baumung, Solkner & Curik (2007), *Conserv Genet*
8:489-495. Does not import PyPedal.

THE MECHANISM, article p.490, verbatim
--------------------------------------
    "Two unique alleles are assigned to each founder, and the genotypes of all
    descendants along the actual pedigree are generated following Mendelian
    segregation rules. To adapt the approach to the calculation of ancestral
    inbreeding coefficients, one needs to keep track of IBD events in the
    pedigree of an individual. This is done by FLAGGING ALLELES ONCE THEY ARE
    IN IBD STATE FOR THE FIRST TIME... The proportion of already flagged alleles
    out of all loci in an individual genome is considered as its ancestral
    inbreeding coefficient. Alleles which are identical by descent for the first
    time are flagged and contribute to ancestral inbreeding coefficient OF
    OFFSPRING."

So the flag is a persistent property of an ALLELE, inherited by every copy of
it, and an individual's coefficient counts the alleles that arrived ALREADY
flagged. Autozygosity arising in the individual itself flags its alleles for
transmission but does not count towards its own coefficient -- which is the same
convention as Ballou's recursion, where f_a is built from the parents.

HALF-FOUNDERS -- settled by a DIFFERENT source, Baumung et al. (2015)
--------------------------------------------------------------------
Suwanlee (2007) p.490 drops alleles from TWO parents and is **silent** on what a
missing parental side contributes. That silence is a historical fact and nothing
below revises it; it is why this oracle previously refused the case outright.

The rule comes from the GRAIN paper instead -- Baumung R., Farkas J., Boichard
D., Meszaros G., Solkner J. & Curik I. (2015), "GRain: a computer program to
calculate ancestral and partial inbreeding coefficients using a gene dropping
approach", *J Anim Breed Genet* 132:100-108, p.102, "Computational strategy",
verbatim:

    "For each half-founder, that is an animal with just one parent known, a
    dummy founder is created and the unknown second parent is assigned an
    artificially created new identification number and also provided with two
    unique alleles."

Corroborated at p.104, where dummy founders appear in GRAIN's own output files.
This is the same estimator and the same flagging mechanism -- GRAIN cites
Suwanlee (2007) for why the stochastic approach is used, and shares three
authors with it.

Status: **SOURCE-EXPLICIT**. It is a method statement, not a derivation: GRAIN
offers no justification for the rule. That is the same evidential form this
project already accepts from Boichard's appendices, so it clears the same bar,
but it is not a theorem.

THE AMBIGUITY -- and it is the source's, not ours
-------------------------------------------------
"proportion of already flagged alleles out of all loci" does not fix the
denominator. Two readings survive:

    D1  ALLELE-COUNTING   flagged alleles / (2 * loci); an individual with one
                          flagged allele at a locus scores one half there.
    D2  LOCUS-INDICATOR   loci with at least one flagged allele / loci.

They are different estimators, not different roundings: D2 >= D1 always.

THE ADJUDICATION CRITERION IS DECLARED HERE, BEFORE ANY MEASUREMENT
-------------------------------------------------------------------
Suwanlee Figure 1 publishes f_a_g for a 23-animal pedigree, but those values are
Monte-Carlo estimates from 10^6 replicates PRINTED ROUNDED TO THREE DECIMALS.
Choosing whichever candidate lands closer after one run would be fitting noise
and rounding, not adjudicating. So:

1. Both candidates are implemented, neither privileged.
2. Both are evaluated against ALL SIX distinct published values, not a
   favourable subset. Figure 1 prints one row of values per generation, and
   each printed row labels between one and four animals, so the six values
   expand to thirteen animal-level checks:

       f_a_g = 0.000  animals 11, 12
       f_a_g = 0.125  animals 13, 14
       f_a_g = 0.375  animals 15, 16
       f_a_g = 0.585  animals 17, 18, 19, 20
       f_a_g = 0.728  animals 21, 22
       f_a_g = 0.822  animal  23

   The seven extra checks are NOT additional source values. They are the same
   six values applied to co-labelled siblings, so they are independent
   Monte-Carlo draws of the same expectation -- corroborating, not extra
   evidence from the paper. A candidate is judged on the six.
3. ``SEEDS`` independent fixed seeds, ``LOCI`` loci each.
4. Publication rounding and sampling error are NOT combined in quadrature --
   they are different kinds of thing. The rounding interval is a hard half-width
   of ``ROUNDING_HALFWIDTH``; the sampling uncertainties are combined with each
   other and scaled by ``Z``:

       |mean_ours - published| <= ROUNDING_HALFWIDTH
                                  + Z * sqrt(SE_ours^2 + SE_published^2)

   with SE_published the binomial standard error at 10^6 replicates and SE_ours
   the standard error of the mean across seeds.
5. DECISION. Exactly one candidate compatible with every row -> that candidate
   is *empirically adjudicated from the published experiment*, a status below
   source-explicit and labelled as such. Both fit, or neither fits -> the
   denominator remains STILL BLOCKED and the measurement is reported as it came
   out.

These constants are not to be adjusted after seeing results. Changing them is a
change to the experiment, and would have to be justified as such.

NOT A PRODUCTION REPAIR
-----------------------
``pyp_metrics.dropped_ancestral_inbreeding`` is deliberately NOT modified in
this phase. It implements neither candidate: it has no allele flag at all, and
sets a per-locus indicator when an IMMEDIATE PARENT is autozygous, so ancestral
inbreeding more than one generation back is invisible to it. The size of that
gap is measured in ``tests/test_suwanlee_oracle.py`` rather than asserted.

OUTCOME OF THE ADJUDICATION (criterion above unchanged since it was run)
-----------------------------------------------------------------------
    D1  satisfies ALL SIX distinct published values (13/13 animal-level checks)
    D2  fails FIVE of the six (11 of 13 animal-level checks), by 0.086 to
        0.171 -- more than twenty times its tolerance, so this is a difference
        in kind, not a close call. Only the trivial 0.000 row fits.

    An earlier report of this said "fits all 13 published rows". That was a
    REPORTING ERROR: there are six published values, checked at thirteen
    animal-level targets. The criterion, the run and the verdict are unchanged
    -- only the count was misdescribed.

Verdict: **D1**, empirically adjudicated from the published experiment. That is
a status below source-explicit and is labelled so wherever it is relied on.

One residual, recorded rather than glossed: all eleven non-zero D1 deviations
are POSITIVE (+0.0002 to +0.0013). Independent seeds agreeing on the sign points
to a small systematic difference rather than sampling noise -- most plausibly in
exactly when the flag is considered set. Every row is inside the predeclared
tolerance, so the criterion is met as written; the sign pattern is noted because
it is the kind of thing that turns out to matter later.
"""
import argparse
import json
import math
import os

import numpy as np

MISSING = 0

#: Predeclared experiment. See the module docstring.
LOCI = 200_000
SEEDS = (1, 2, 3, 4, 5)
Z = 1.96                        # two-sided 95%
ROUNDING_HALFWIDTH = 5e-4       # published to three decimals
PUBLISHED_REPLICATES = 1_000_000

CANDIDATES = ("D1", "D2")


def read_pedigree(path, sepchar=" ", animal=0, sire=1, dam=2):
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                continue
            parts = [p for p in (line.split(sepchar) if sepchar != " " else line.split())
                     if p != ""]
            rows.append((int(parts[animal]), int(parts[sire]), int(parts[dam])))
    return rows


def gene_drop(rows, loci=LOCI, seed=0, candidate="D1"):
    """
    One replicate block: ``loci`` independent diallelic unlinked loci, dropped
    through the pedigree at once.

    Returns ``{animal_id: f_a}``. Vectorised over loci -- the loci are
    independent, so a locus is a column and the whole genome moves together.
    """
    if candidate not in CANDIDATES:
        raise ValueError(f"candidate must be one of {CANDIDATES}")
    rng = np.random.default_rng(seed)

    allele = {}          # animal -> (labels_1, labels_2), founder-allele identity
    flagged = {}         # animal -> (flags_1, flags_2), persistent IBD flags
    result = {}

    next_label = 0
    for animal_id, sire_id, dam_id in rows:
        if sire_id == MISSING and dam_id == MISSING:
            # Two unique alleles per founder.
            a1 = np.full(loci, next_label, dtype=np.int64)
            a2 = np.full(loci, next_label + 1, dtype=np.int64)
            next_label += 2
            allele[animal_id] = (a1, a2)
            flagged[animal_id] = (np.zeros(loci, dtype=bool),
                                  np.zeros(loci, dtype=bool))
            result[animal_id] = 0.0
            continue

        picks = []
        for parent in (sire_id, dam_id):
            take_first = rng.random(loci) < 0.5
            if parent == MISSING:
                # THE HALF-FOUNDER RULE, Baumung et al. (2015) p.102,
                # "Computational strategy", source-explicit:
                #
                #   "For each half-founder, that is an animal with just one
                #   parent known, a dummy founder is created and the unknown
                #   second parent is assigned an artificially created new
                #   identification number and also provided with two unique
                #   alleles."
                #
                # So the unknown parental side is an ORDINARY FOUNDER: two
                # unique alleles, never flagged, one gamete transmitted by the
                # same Mendelian draw as any other parent.
                #
                # Minted here and deliberately NOT stored in `allele`/`flagged`.
                # A dummy serves exactly one offspring, so it is never looked up
                # again -- and not keying it makes it impossible for two
                # different unknown parental slots to collide, which is the
                # failure this rule exists to avoid. Two half-founders both
                # recording the sentinel 0 are NOT the same individual: the
                # pedigree carries no information that they are, and merging
                # them would invent a relationship.
                #
                # `next_label` is the same allocator the founder branch uses, so
                # a dummy is unrelated to every recorded founder and to every
                # other dummy by construction.
                p_alleles = (np.full(loci, next_label, dtype=np.int64),
                             np.full(loci, next_label + 1, dtype=np.int64))
                p_flags = (np.zeros(loci, dtype=bool),
                           np.zeros(loci, dtype=bool))
                next_label += 2
            else:
                p_alleles, p_flags = allele[parent], flagged[parent]
            picks.append((np.where(take_first, p_alleles[0], p_alleles[1]),
                          np.where(take_first, p_flags[0], p_flags[1])))
        (a1, f1), (a2, f2) = picks

        # The individual's own coefficient counts alleles that arrived ALREADY
        # flagged -- autozygosity arising here contributes to its OFFSPRING.
        if candidate == "D1":
            result[animal_id] = float((f1.sum() + f2.sum()) / (2.0 * loci))
        else:
            result[animal_id] = float(np.count_nonzero(f1 | f2) / float(loci))

        # Alleles identical by descent for the first time are flagged, and the
        # flag then travels with them.
        autozygous = a1 == a2
        allele[animal_id] = (a1, a2)
        flagged[animal_id] = (f1 | autozygous, f2 | autozygous)

    return result


def measure(rows, candidate, loci=LOCI, seeds=SEEDS):
    """Run every seed and return ``{animal_id: (mean, standard error)}``."""
    runs = [gene_drop(rows, loci=loci, seed=seed, candidate=candidate)
            for seed in seeds]
    out = {}
    for animal_id in runs[0]:
        values = [run[animal_id] for run in runs]
        mean = sum(values) / len(values)
        if len(values) > 1:
            variance = sum((v - mean) ** 2 for v in values) / (len(values) - 1)
            sem = math.sqrt(variance / len(values))
        else:
            sem = 0.0
        out[animal_id] = (mean, sem)
    return out


def tolerance(published, sem_ours):
    """The predeclared acceptance half-width for one published value."""
    se_published = math.sqrt(max(published * (1.0 - published), 0.0)
                             / PUBLISHED_REPLICATES)
    return ROUNDING_HALFWIDTH + Z * math.sqrt(sem_ours ** 2 + se_published ** 2)


def adjudicate(rows, published, loci=LOCI, seeds=SEEDS):
    """
    Run the predeclared procedure. ``published`` maps animal id -> f_a_g.

    Returns a report dict; the ``verdict`` is one of ``"D1"``, ``"D2"`` or
    ``"STILL BLOCKED"``.
    """
    report = {"loci": loci, "seeds": list(seeds), "z": Z,
              "rounding_halfwidth": ROUNDING_HALFWIDTH,
              "published_replicates": PUBLISHED_REPLICATES, "candidates": {}}
    fits = []
    for candidate in CANDIDATES:
        measured = measure(rows, candidate, loci=loci, seeds=seeds)
        rows_out, ok = [], True
        for animal_id in sorted(published):
            mean, sem = measured[animal_id]
            want = published[animal_id]
            tol = tolerance(want, sem)
            inside = abs(mean - want) <= tol
            ok = ok and inside
            rows_out.append({"animal": animal_id, "published": want,
                             "mean": mean, "sem": sem, "tolerance": tol,
                             "deviation": mean - want, "inside": inside})
        report["candidates"][candidate] = {"rows": rows_out, "fits_every_row": ok}
        if ok:
            fits.append(candidate)

    if len(fits) == 1:
        report["verdict"] = fits[0]
        report["status"] = ("empirically adjudicated from the published "
                            "experiment -- a status below source-explicit")
    else:
        report["verdict"] = "STILL BLOCKED"
        report["status"] = (
            "both candidates fit" if len(fits) == 2 else "neither candidate fits")
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("pedfile")
    ap.add_argument("--candidate", choices=CANDIDATES, default=None)
    ap.add_argument("--loci", type=int, default=LOCI)
    ap.add_argument("--adjudicate", action="store_true",
                    help="run the predeclared D1-vs-D2 procedure against "
                         "Suwanlee Figure 1")
    args = ap.parse_args()

    rows = read_pedigree(args.pedfile)

    if args.adjudicate:
        # Suwanlee Figure 1, article p.491, the f_a_g column.
        published = {11: 0.000, 12: 0.000, 13: 0.125, 14: 0.125,
                     15: 0.375, 16: 0.375,
                     17: 0.585, 18: 0.585, 19: 0.585, 20: 0.585,
                     21: 0.728, 22: 0.728, 23: 0.822}
        print(json.dumps(adjudicate(rows, published, loci=args.loci),
                         indent=2, sort_keys=True))
        return

    candidate = args.candidate or "D1"
    measured = measure(rows, candidate, loci=args.loci)
    print(json.dumps({
        "pedfile": os.path.basename(args.pedfile),
        "candidate": candidate, "loci": args.loci, "seeds": list(SEEDS),
        "f_a_gene_drop": {str(k): v[0] for k, v in sorted(measured.items())},
        "sem": {str(k): v[1] for k, v in sorted(measured.items())},
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
