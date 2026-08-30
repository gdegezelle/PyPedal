#!/usr/bin/env python3
"""
INDEPENDENT oracle for Boichard, Maignel & Verrier (1997).
Does not import PyPedal.

Implements, verbatim and separately:

  * **Appendix A** (article p.22) -- probabilities of gene origin, and the
    effective number of founders  f_e = 1 / sum(q_k^2)   [eq. 1, p.7]
  * **Appendix B** (article pp.22-23) -- the marginal contributions of the most
    important ancestors, and  f_a = 1 / sum(p_k^2)       [p.8]
  * **pp.9-10** -- the lower and upper bounds f_l and f_u

WHAT IS PAPER-EXPLICIT AND WHAT IS NOT
--------------------------------------
This separation is the point of the file. Everything in ``_appendix_a`` and
``_appendix_b_round`` is transcribed from the appendices. Three things are NOT
in the paper, are switchable here, and must never be presented as the paper's:

``tie_break``     Lowest-ID convention. The paper explicitly declines to specify one:
                  "when two related ancestors have the same marginal
                  contribution, the final result may depend on the chosen one"
                  (p.10). ``lowest_id`` is the PyPedal CONVENTION, chosen for
                  determinism. Boichard's Table I happens to break its round-1
                  tie (7 vs 8, both 0.300) the same way, which is consistent
                  with the convention but is one instance, not a rule.

``half_founder``  Appendix A step 4 halves a half-founder's own
                  contribution. Appendix B never restates it, so whether it
                  applies to f_a is source-silent. All three readings are
                  offered; none is a default that hides the question. Appendix
                  A's halving is NOT switchable -- there the paper is explicit.

                      "none"      no half-founder handling inside Appendix B
                      "halve"     Appendix A step 4 carried in, in place
                      "complete"  READING C -- phantom-founder completion of
                                  the pedigree, then Appendix B unchanged

                  Reading C is what the project adjudicated and implemented.
                  Its evidence class is **mathematically implied /
                  independently supported, NOT source-explicit Appendix-B
                  text**: Appendix B does not contain it. It was selected
                  against invariants the paper does state -- p.7 Sum q = 1,
                  p.8 Sum p = 1 and k <= f, p.17 f_a <= f_e -- which "none"
                  and "halve" each violate.

``phantom_encoding``
                  Which integers stand in for the phantoms under Reading C.
                  The numbers carry NO meaning; the parameter exists so the
                  encoding can be varied and the scientific answer shown to be
                  independent of it, rather than that independence being
                  assumed. ``above_max`` is production's contiguous
                  max(real)+n; ``times_100`` is the adjudication probe's
                  max(real)*100+n.

``reference_pop`` Whether members of the population under study may
                  themselves be selected as ancestors. The paper's figures draw
                  the population as unnumbered nodes, so the question never
                  arises there. ``exclude`` is the long-standing PyPedal and
                  PyPedal-2.0.4 convention and the default.

None of the three is discriminated by the published examples, which is exactly
why they remain open: reproducing Table I and Table II does not validate them.

INVARIANTS -- ASSERTED, NEVER CLAMPED
-------------------------------------
``0 <= a(i) <= 1`` holds structurally: Appendix B step 6 is a *pull*, so a
non-selected animal's a is a convex combination of its parents' settled values,
and selected ancestors are pinned at 1 with their parents deleted. A violation
is an implementation defect, so it raises. Clamping would hide it.

USAGE
-----
    python oracle_boichard.py corpus/boichard_fig1.ped --gen-col 3
    python oracle_boichard.py corpus/boichard_fig2.ped --gen-col 3 --reference 7-14,19,20
"""
import argparse
import json
import os

MISSING = 0

TIE_BREAKS = ("lowest_id", "highest_id")
HALF_FOUNDER_RULES = ("none", "halve", "complete")
REFERENCE_POP_RULES = ("exclude", "include")
PHANTOM_ENCODINGS = ("above_max", "times_100")


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def read_pedigree(path, sepchar=" ", animal=0, sire=1, dam=2, gen_col=None):
    """
    Parse a PyPedal-style .ped file into ordered (id, sire, dam) triples, plus
    a {id: generation} map when ``gen_col`` is given.

    The oracle deliberately does NOT renumber. Its whole purpose is comparison
    with published tables, and renumbering would break the correspondence with
    the papers' own animal labels. Topological order is asserted instead.
    """
    rows, gens = [], {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                continue
            parts = [p for p in (line.split(sepchar) if sepchar != " " else line.split())
                     if p != ""]
            a, s, d = int(parts[animal]), int(parts[sire]), int(parts[dam])
            rows.append((a, s, d))
            if gen_col is not None:
                gens[a] = parts[gen_col]
    return rows, gens


def check_topological(rows):
    """Parents must precede offspring, or every downstream pass is wrong."""
    seen = set()
    for a, s, d in rows:
        for p in (s, d):
            if p != MISSING and p not in seen:
                raise ValueError(
                    f"animal {a} has parent {p} which does not appear earlier; "
                    "this oracle requires a topologically ordered pedigree so "
                    "that the papers' animal labels are preserved")
        seen.add(a)
    return True


# --------------------------------------------------------------------------
# pedigree helpers
# --------------------------------------------------------------------------

def _parents(rows):
    return {a: (s, d) for a, s, d in rows}


def founders(rows):
    """
    Boichard p.7: "A founder is defined as an ancestor with unknown parents.
    Note that when an animal has only one known parent, the animal is
    considered as a founder."

    So the founder set includes half-founders. This is the ``f`` of the bound
    formulas on pp.9-10.
    """
    return [a for a, s, d in rows if s == MISSING or d == MISSING]


def half_founders(rows):
    """Exactly one known parent -- Appendix A step 4's 'half founder'."""
    return [a for a, s, d in rows if (s == MISSING) != (d == MISSING)]


def unknown_parent_slots(rows):
    """
    One entry per unknown parental SIDE of an animal that has at least one
    known parent, as ``[(animal, 'sire'|'dam'), ...]`` in pedigree order.

    A slot, not an animal, is the unit. Two half-founders both recording the
    sentinel are two distinct unknown individuals, and an animal with BOTH
    sides unknown is an ordinary founder and contributes no slot at all.
    """
    slots = []
    for a, s, d in rows:
        if s == MISSING and d == MISSING:
            continue
        if s == MISSING:
            slots.append((a, "sire"))
        if d == MISSING:
            slots.append((a, "dam"))
    return slots


def phantom_complete(rows, encoding="above_max"):
    """
    READING C's pedigree transform, and nothing else.

    Every unknown parental side of an animal that has at least one known parent
    becomes a founder with a record of its own -- Appendix A step 4's own
    words, "This is equivalent to considering the unknown parent as a founder",
    applied to the pedigree rather than to the contribution vector. An animal
    with BOTH parents unknown is already a founder and is left alone.

    Returns ``(completed_rows, {phantom_id: (animal, side)})``.

    Phantoms are emitted FIRST so that parents still precede offspring, which
    Appendix B's single-pass steps 5 and 6 require. Their IDs are strictly
    above every real ID under both encodings, so a real animal wins a
    real-vs-phantom ``lowest_id`` tie and the lowest-ID convention keeps its meaning.

    THE NUMBERS MEAN NOTHING. Two encodings exist so that the independence of
    the scientific result from the encoding can be measured rather than
    asserted from the fact that both sort above the real IDs.
    """
    if encoding not in PHANTOM_ENCODINGS:
        raise ValueError(f"phantom_encoding must be one of {PHANTOM_ENCODINGS}")
    slots = unknown_parent_slots(rows)
    if not slots:
        return list(rows), {}

    top = max(a for a, _, _ in rows)
    base = top if encoding == "above_max" else top * 100
    fill = {slot: base + index + 1 for index, slot in enumerate(slots)}
    phantoms = {pid: slot for slot, pid in fill.items()}

    completed = [(pid, MISSING, MISSING) for pid in sorted(phantoms)]
    for a, s, d in rows:
        if s == MISSING and d == MISSING:
            completed.append((a, s, d))
            continue
        completed.append((a,
                          s if s != MISSING else fill[(a, "sire")],
                          d if d != MISSING else fill[(a, "dam")]))
    return completed, phantoms


def _prepare(rows, half_founder, phantom_encoding):
    """
    Apply the pedigree-level part of the half-founder reading, and report the
    rule that remains for the engine to apply in place.

    Reading C is entirely a transform of the pedigree: once it has run there
    are no half-founders left, so Appendix A step 4 is vacuous and the engine
    below runs with no half-founder handling at all. That is the whole content
    of "then ordinary Appendix B with no special case".
    """
    if half_founder != "complete":
        return list(rows), {}, half_founder
    completed, phantoms = phantom_complete(rows, phantom_encoding)
    return completed, phantoms, "none"


# --------------------------------------------------------------------------
# Appendix A -- probabilities of gene origin
# --------------------------------------------------------------------------

def _appendix_a(rows, reference, parents, apply_half_founder_rule=True):
    """
    Appendix A steps 2-4, verbatim.

    2. initialise q with 1 for animals in the population under study, 0 otherwise
    3. process the pedigree from the YOUNGEST animal to the OLDEST:
           if sire(i) known: q(sire(i)) += 0.5*q(i)
           if dam(i)  known: q(dam(i))  += 0.5*q(i)
    4. if an animal is a half founder, multiply ITS contribution by 0.5.
       Divide the vector q by N.

    The ordering of steps 3 and 4 is load-bearing and is the whole of BL-3(a):
    the known parent receives half of the FULL q(i), and only afterwards is the
    half-founder's own retained contribution halved. Halving first would pass
    0.25*q(i) up instead of 0.5*q(i).
    """
    q = {a: 0.0 for a, _, _ in rows}
    for a in reference:
        q[a] = 1.0

    for a, s, d in reversed(rows):                      # youngest -> oldest
        if s != MISSING:
            q[s] += 0.5 * q[a]
        if d != MISSING:
            q[d] += 0.5 * q[a]

    if apply_half_founder_rule:
        for a in half_founders(rows):
            q[a] *= 0.5

    n = float(len(reference))
    return {a: v / n for a, v in q.items()}


def f_e(rows, reference, half_founder="halve", phantom_encoding="above_max"):
    """
    Effective number of founders, eq. 1 (p.7): f_e = 1 / sum over founders q_k^2.

    Returns (f_e, q_by_founder, q_sum). ``q_sum`` should be 1.0 without being
    normalised to it -- Appendix A step 4 says the division by N is what makes
    founder contributions sum to one, so an unforced 1.0 is a real check and a
    departure from it means the half-founder handling is wrong.

    ``half_founder`` defaults to ``"halve"``, which is Appendix A step 4 as
    written and is what BL-3 is settled on. ``"complete"`` books the same half
    to a phantom node instead of retaining it on the half-founder. Appendix A
    step 4 states the two are equivalent -- "This is equivalent to considering
    the unknown parent as a founder" -- so f_e must come out identical, and
    that identity is a test, not an assumption.
    """
    working, _phantoms, inner = _prepare(rows, half_founder, phantom_encoding)
    parents = _parents(working)
    q = _appendix_a(working, reference, parents,
                    apply_half_founder_rule=(inner != "none"))
    contrib = {a: q[a] for a in founders(working)}
    ssq = sum(v * v for v in contrib.values())
    return ((1.0 / ssq) if ssq > 0 else 0.0), contrib, sum(contrib.values())


# --------------------------------------------------------------------------
# Appendix B -- the shared marginal-contribution engine
# --------------------------------------------------------------------------

def _appendix_b_round(rows, reference, selected, half_founder, reference_pop):
    """
    One round of Appendix B: steps 3-7, returning p for every candidate.

    3. delete the pedigree information of the ancestors already found
    4. q = 1 on the population under study; a = 1 on the already-selected
    5. q processed YOUNGEST -> OLDEST
    6. a processed OLDEST -> YOUNGEST
    7. p(i) = q(i) * (1 - a(i))
    """
    chosen = set(selected)

    # step 3 -- selected ancestors become pseudo founders
    working = [(a, MISSING, MISSING) if a in chosen else (a, s, d)
               for a, s, d in rows]

    # step 4
    q = {a: 0.0 for a, _, _ in working}
    for a in reference:
        q[a] = 1.0
    a_vec = {a: (1.0 if a in chosen else 0.0) for a, _, _ in working}

    # step 5 -- q youngest -> oldest
    for a, s, d in reversed(working):
        if s != MISSING:
            q[s] += 0.5 * q[a]
        if d != MISSING:
            q[d] += 0.5 * q[a]

    # Half-founders in Appendix B. The appendix does not restate Appendix A
    # step 4; both readings are offered and neither is silently adopted.
    if half_founder == "halve":
        for a, s, d in working:
            if (s == MISSING) != (d == MISSING):
                q[a] *= 0.5

    # step 6 -- a oldest -> youngest
    for a, s, d in working:
        if a in chosen:
            continue                       # pinned at 1; parents were deleted
        acc = 0.0
        if s != MISSING:
            acc += 0.5 * a_vec[s]
        if d != MISSING:
            acc += 0.5 * a_vec[d]
        a_vec[a] += acc

    for a, value in a_vec.items():
        if not (-1e-9 <= value <= 1.0 + 1e-9):
            raise AssertionError(
                f"a({a}) = {value!r} is outside [0, 1]. Appendix B step 6 is a "
                "pull, so a is a convex combination of settled parental values "
                "and cannot leave the interval. This is an implementation "
                "defect and must not be clamped away.")

    # step 7
    p = {a: q[a] * (1.0 - a_vec[a]) for a, _, _ in working}

    # May members of the population under study be ancestors of themselves?
    if reference_pop == "exclude":
        for a in reference:
            p[a] = 0.0

    for a in chosen:
        p[a] = 0.0                          # already selected; a=1 gives 0 anyway
    return p, q, a_vec


def marginal_contributions(rows, reference, tie_break="lowest_id",
                           half_founder="none", reference_pop="exclude",
                           tol=1e-12, max_rounds=None,
                           phantom_encoding="above_max"):
    """
    THE SHARED ENGINE. Yields ``(ancestor, marginal_contribution)`` in the
    order Appendix B selects them, contributions already divided by N (step 8).

    Both the exact f_a and the pp.9-10 bounds are truncations of this one
    sequence, so they consume it rather than reimplementing it. That is what
    stops the definite and bounded routines drifting apart.

    Stops when no candidate has a positive contribution: "An exact computation
    of f_a requires the determination of every ancestor with a non-zero
    contribution" (p.8).

    Under ``half_founder="complete"`` a phantom may itself be selected and
    yielded. That is forced by Reading C rather than chosen: a phantom is an
    ordinary founder, it is not in the population under study so the
    reference-population exclusion never reaches it, and withholding its contribution would leave
    Sum p < 1 -- the p.8 invariant Reading C was selected by.
    """
    if tie_break not in TIE_BREAKS:
        raise ValueError(f"tie_break must be one of {TIE_BREAKS}")
    if half_founder not in HALF_FOUNDER_RULES:
        raise ValueError(f"half_founder must be one of {HALF_FOUNDER_RULES}")
    if reference_pop not in REFERENCE_POP_RULES:
        raise ValueError(f"reference_pop must be one of {REFERENCE_POP_RULES}")

    rows, _phantoms, half_founder = _prepare(rows, half_founder, phantom_encoding)

    reference = list(reference)
    if not reference:
        raise ValueError("the population under study is empty")
    n = float(len(reference))
    selected, rounds = [], 0
    limit = max_rounds if max_rounds is not None else len(rows) + 1

    while rounds < limit:
        p, _, _ = _appendix_b_round(rows, reference, selected, half_founder,
                                    reference_pop)
        best = max(p.values()) if p else 0.0
        if best <= tol:
            return
        tied = sorted(a for a, v in p.items() if v >= best - tol)
        # Lowest-ID tie-break is a convention, not a source rule.
        pick = tied[0] if tie_break == "lowest_id" else tied[-1]
        selected.append(pick)
        rounds += 1
        yield pick, p[pick] / n


def f_a(rows, reference, **kw):
    """
    Effective number of ancestors, f_a = 1 / sum(p_k^2) (p.8).

    Returns (f_a, ordered [(ancestor, contribution)], sum of contributions).
    The paper states the marginal contributions sum to one; that is checked by
    the caller rather than enforced here.
    """
    order = list(marginal_contributions(rows, reference, **kw))
    ssq = sum(v * v for _, v in order)
    return ((1.0 / ssq) if ssq > 0 else 0.0), order, sum(v for _, v in order)


# --------------------------------------------------------------------------
# pp.9-10 -- the bounds
# --------------------------------------------------------------------------

def bounds(rows, reference, n_ancestors, **kw):
    """
    f_l and f_u after the first ``n_ancestors`` have been taken (pp.9-10):

        c   = sum of the first n marginal contributions
        f_u = 1 / [ sum(p_i^2) + (1-c)^2 / (f - n) ]
        f_l = 1 / [ sum(p_i^2) + m * p_n^2 ],   m = (1-c)/p_n

    ENDPOINT AND DOMAIN HANDLING. These are algebraic safeguards derived from
    the published formulas, not new semantics:

      * residual mass 1-c zero within tolerance -> the truncation has reached
        the exact answer; return f_l = f_u = f_a and evaluate neither singular
        residual term;
      * never divide by (f - n) when n == f;
      * p_n == 0 with residual mass still positive is an inconsistent internal
        state -- raise rather than invent a bound;
      * the residual-mass test uses an explicit tolerance, never == 0.0, since
        c is a sum of many floating-point contributions.
    """
    tol = 1e-12
    order = list(marginal_contributions(rows, reference, **kw))
    if n_ancestors < 1:
        raise ValueError("n_ancestors must be at least 1")
    taken = order[:n_ancestors]
    if not taken:
        raise ValueError("no ancestor has a positive contribution")

    ssq = sum(v * v for _, v in taken)
    c = sum(v for _, v in taken)
    p_n = taken[-1][1]
    residual = 1.0 - c
    n_taken = len(taken)

    # `f` must be the founder count of the pedigree the engine actually ran on.
    # Under Reading C that is the COMPLETED pedigree, where each phantom is an
    # ordinary founder; counting founders on the caller's rows instead would
    # feed a different f into f_u = 1/(ssq + (1-c)^2/(f-n)) than the sequence
    # was generated from.
    counted, _phantoms, _inner = _prepare(rows,
                                          kw.get("half_founder", "none"),
                                          kw.get("phantom_encoding", "above_max"))
    n_founders = len(founders(counted))

    if residual <= tol:
        exact = (1.0 / ssq) if ssq > 0 else 0.0
        return exact, exact, {"exact": True, "c": c, "n_taken": n_taken}

    if n_taken >= n_founders:
        raise AssertionError(
            f"{n_taken} ancestors taken but the pedigree has {n_founders} "
            f"founders and {residual!r} of the contribution is still "
            "unexplained; f_u would divide by zero. The selection sequence is "
            "inconsistent.")
    if p_n <= tol:
        raise AssertionError(
            f"the {n_taken}th marginal contribution is {p_n!r} while "
            f"{residual!r} of the contribution is still unexplained; m = "
            "(1-c)/p_n is undefined. The selection sequence is inconsistent.")

    f_u = 1.0 / (ssq + residual * residual / float(n_founders - n_taken))
    m = residual / p_n
    f_l = 1.0 / (ssq + m * p_n * p_n)
    return f_l, f_u, {"exact": False, "c": c, "m": m, "n_taken": n_taken,
                      "n_founders": n_founders}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _parse_reference(spec):
    out = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if "-" in chunk:
            lo, hi = chunk.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        elif chunk:
            out.append(int(chunk))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("pedfile")
    ap.add_argument("--animal", type=int, default=0)
    ap.add_argument("--sire", type=int, default=1)
    ap.add_argument("--dam", type=int, default=2)
    ap.add_argument("--gen-col", type=int, default=None,
                    help="column holding the generation; the reference "
                         "population is then the most recent generation")
    ap.add_argument("--reference", default=None,
                    help="explicit population under study, e.g. 7-14,19,20")
    ap.add_argument("--tie-break", choices=TIE_BREAKS, default="lowest_id")
    ap.add_argument("--half-founder", choices=HALF_FOUNDER_RULES, default="none")
    ap.add_argument("--phantom-encoding", choices=PHANTOM_ENCODINGS,
                    default="above_max",
                    help="which integers stand in for phantoms under "
                         "--half-founder complete; the numbers mean nothing")
    ap.add_argument("--reference-pop", choices=REFERENCE_POP_RULES, default="exclude")
    ap.add_argument("--bounds-at", type=int, default=None)
    args = ap.parse_args()

    rows, gens = read_pedigree(args.pedfile, animal=args.animal, sire=args.sire,
                               dam=args.dam, gen_col=args.gen_col)
    check_topological(rows)

    if args.reference:
        reference = _parse_reference(args.reference)
    elif gens:
        newest = max(gens.values(), key=lambda g: float(g))
        reference = [a for a, _, _ in rows if gens[a] == newest]
    else:
        raise SystemExit("give --reference or --gen-col to define the "
                         "population under study")

    kw = dict(tie_break=args.tie_break, half_founder=args.half_founder,
              reference_pop=args.reference_pop,
              phantom_encoding=args.phantom_encoding)
    fe, q, q_sum = f_e(rows, reference,
                       half_founder=("complete" if args.half_founder == "complete"
                                     else "halve"),
                       phantom_encoding=args.phantom_encoding)
    fa, order, p_sum = f_a(rows, reference, **kw)

    counted, phantoms, _inner = _prepare(rows, args.half_founder,
                                         args.phantom_encoding)
    out = {
        "pedfile": os.path.basename(args.pedfile),
        "n_animals": len(rows),
        "reference_population": sorted(reference),
        "n_reference": len(reference),
        "founders": sorted(founders(counted)),
        "n_founders": len(founders(counted)),
        "n_half_founders": len(half_founders(rows)),
        "phantoms": {str(pid): list(slot) for pid, slot in sorted(phantoms.items())},
        "conventions": kw,
        "f_e_boichard": fe,
        "founder_q": {str(k): v for k, v in sorted(q.items())},
        "founder_q_sums_to_one_unforced": abs(q_sum - 1.0) < 1e-9,
        "f_a_boichard": fa,
        "marginal_contributions": [[a, v] for a, v in order],
        "marginal_contributions_sum": p_sum,
        "marginal_sums_to_one_unforced": abs(p_sum - 1.0) < 1e-9,
    }
    if args.bounds_at:
        f_l, f_u, meta = bounds(rows, reference, args.bounds_at, **kw)
        out["bounds"] = {"n": args.bounds_at, "f_l": f_l, "f_u": f_u, **meta}
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
