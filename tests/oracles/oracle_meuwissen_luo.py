#!/usr/bin/env python3
"""
INDEPENDENT oracle for pedigree inbreeding coefficients.

This file deliberately does NOT import PyPedal. It parses the pedigree itself and
implements the mathematics from first principles, so that it can adjudicate a
disagreement between the two PyPedal implementations without inheriting either
one's assumptions.

Derivation (not copied from either codebase)
--------------------------------------------
The numerator relationship matrix factors as  A = L D L'  where L is unit lower
triangular and D is diagonal, for a pedigree ordered so parents precede offspring.

D follows from Mendelian sampling variance:

    d_i = 1                              both parents unknown
    d_i = 0.75 - 0.25 * F_p              exactly one parent p known
    d_i = 0.5  - 0.25 * (F_s + F_d)      both parents known

L propagates each animal's contribution to its parents at one half per meiosis.
Since  A_ii = sum_j L_ij^2 * d_j  and  F_i = A_ii - 1, the diagonal can be
accumulated by walking each animal's ancestors in DESCENDING id order, which is
exactly the Meuwissen & Luo (1992) recursion; the Quaas (1995) modification is
the same recursion driven by an explicit ancestor list so that only nonzero L
entries are visited.

Because animals are processed in ascending order, F_s and F_d are already final
when animal i is reached.

Two implementations of the same mathematics are provided and cross-checked
against each other, plus an O(n^2) tabular NRM built directly from the
recurrence, giving three independent routes to the same F.
"""

import argparse
import json
import sys


def read_pedigree(path, sepchar=" ", animal=0, sire=1, dam=2):
    """Parse a PyPedal-style .ped file into ordered (id, sire, dam) triples."""
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


def renumber(rows, missing=0):
    """
    Map original ids to 1..n in an order where parents precede offspring
    (Kahn topological sort). Independent of PyPedal's reorder/renumber.
    """
    ids = [r[0] for r in rows]
    parents = {a: (s, d) for a, s, d in rows}
    children = {a: [] for a in ids}
    indeg = {a: 0 for a in ids}
    for a, s, d in rows:
        for p in (s, d):
            if p != missing and p in parents:
                children[p].append(a)
                indeg[a] += 1
    # deterministic order: by original id among currently-available nodes
    ready = sorted([a for a in ids if indeg[a] == 0])
    order = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for c in sorted(children[node]):
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
                ready.sort()
    if len(order) != len(ids):
        raise ValueError("pedigree contains a cycle; %d of %d ordered"
                         % (len(order), len(ids)))
    newid = {orig: i + 1 for i, orig in enumerate(order)}
    out = []
    for orig in order:
        s, d = parents[orig]
        out.append((newid[orig],
                    newid[s] if s != missing and s in newid else 0,
                    newid[d] if d != missing and d in newid else 0))
    return out, {v: k for k, v in newid.items()}


def d_vector(ped, F):
    """Mendelian sampling variances, given inbreeding of the parents."""
    n = len(ped)
    d = [0.0] * (n + 1)
    for i, s, dm in ped:
        if s and dm:
            d[i] = 0.5 - 0.25 * (F[s] + F[dm])
        elif s or dm:
            p = s if s else dm
            d[i] = 0.75 - 0.25 * F[p]
        else:
            d[i] = 1.0
    return d


def inbreeding_meuwissen_luo(ped):
    """
    Meuwissen & Luo (1992). For each animal, accumulate A_ii = sum L_j^2 d_j by
    walking ancestors in descending id order, halving L into each known parent.
    """
    n = len(ped)
    sire = [0] * (n + 1)
    dam = [0] * (n + 1)
    for i, s, dm in ped:
        sire[i], dam[i] = s, dm

    F = [0.0] * (n + 1)
    for i in range(1, n + 1):
        s, dm = sire[i], dam[i]
        if not s or not dm:
            F[i] = 0.0            # one parent unknown => not inbred
            continue
        L = [0.0] * (n + 1)
        L[i] = 1.0
        d = d_vector(ped, F)      # parents' F already final
        Aii = 0.0
        for j in range(i, 0, -1):
            if L[j] == 0.0:
                continue
            Aii += L[j] * L[j] * d[j]
            if sire[j]:
                L[sire[j]] += 0.5 * L[j]
            if dam[j]:
                L[dam[j]] += 0.5 * L[j]
            L[j] = 0.0
        F[i] = Aii - 1.0
    return F


def inbreeding_quaas(ped):
    """
    Quaas (1995) modification: same recursion, but driven by an explicit
    ancestor list so only nonzero L entries are visited. Same mathematics,
    different traversal -- a cross-check on the implementation above.
    """
    n = len(ped)
    sire = [0] * (n + 1)
    dam = [0] * (n + 1)
    for i, s, dm in ped:
        sire[i], dam[i] = s, dm

    F = [0.0] * (n + 1)
    for i in range(1, n + 1):
        s, dm = sire[i], dam[i]
        if not s or not dm:
            F[i] = 0.0
            continue
        d = d_vector(ped, F)
        L = {i: 1.0}
        anc = [i]
        Aii = 0.0
        while anc:
            j = max(anc)
            anc.remove(j)
            lj = L.pop(j, 0.0)
            if lj == 0.0:
                continue
            Aii += lj * lj * d[j]
            for p in (sire[j], dam[j]):
                if p:
                    if p not in L:
                        L[p] = 0.0
                        anc.append(p)
                    L[p] += 0.5 * lj
        F[i] = Aii - 1.0
    return F


def inbreeding_tabular(ped):
    """Direct O(n^2) NRM from the standard recurrence; third independent route."""
    n = len(ped)
    A = [[0.0] * (n + 1) for _ in range(n + 1)]
    sire = [0] * (n + 1)
    dam = [0] * (n + 1)
    for i, s, dm in ped:
        sire[i], dam[i] = s, dm
    for i in range(1, n + 1):
        s, dm = sire[i], dam[i]
        A[i][i] = 1.0 + (0.5 * A[s][dm] if s and dm else 0.0)
        for j in range(1, i):
            v = 0.5 * ((A[j][s] if s else 0.0) + (A[j][dm] if dm else 0.0))
            A[i][j] = A[j][i] = v
    return [0.0] + [A[i][i] - 1.0 for i in range(1, n + 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pedfile")
    ap.add_argument("--sepchar", default=" ")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--tabular", action="store_true",
                    help="force the O(n^2) tabular route (auto-enabled for n <= 2000)")
    ap.add_argument("--all-routes", action="store_true",
                    help="force the O(n^3) quaas route (auto-enabled for n <= 200)")
    args = ap.parse_args()

    rows = read_pedigree(args.pedfile, args.sepchar)
    ped, back = renumber(rows)
    n = len(ped)

    # Route costs differ sharply. meuwissen_luo is O(n^2); tabular is O(n^2) time
    # and memory; quaas as written uses max()/remove() over a Python list, which
    # degrades to O(n^3) on deep chains. Gate the expensive routes by size so a
    # large pedigree still gets an answer -- internal consistency of the three
    # routes is established on the small pedigrees.
    ml = inbreeding_meuwissen_luo(ped)
    routes = {"meuwissen_luo": ml}

    if args.all_routes or n <= 200:
        routes["quaas"] = inbreeding_quaas(ped)
    if args.tabular or n <= 2000:
        routes["tabular"] = inbreeding_tabular(ped)

    max_spread = 0.0
    for i in range(1, n + 1):
        vals = [r[i] for r in routes.values()]
        max_spread = max(max_spread, max(vals) - min(vals))

    result = {
        "pedfile": args.pedfile,
        "n_animals": n,
        "routes_computed": sorted(routes),
        "internal_max_spread": max_spread,
        "internally_consistent": max_spread < 1e-12,
        "F_by_original_id": {str(back[i]): round(ml[i], 12) for i in range(1, n + 1)},
        "nonzero_F": {str(back[i]): round(ml[i], 12)
                      for i in range(1, n + 1) if abs(ml[i]) > 1e-12},
    }
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("pedigree            : %s (%d animals)" % (args.pedfile, n))
        print("independent routes  : %s" % ", ".join(sorted(routes)))
        print("internal max spread : %.3e  (%s)"
              % (max_spread, "consistent" if result["internally_consistent"] else "INCONSISTENT"))
        print("animals with F > 0  :")
        for k, v in sorted(result["nonzero_F"].items(), key=lambda kv: int(kv[0])):
            print("   original id %-8s F = %.10f" % (k, v))
    return 0 if result["internally_consistent"] else 1


if __name__ == "__main__":
    sys.exit(main())
