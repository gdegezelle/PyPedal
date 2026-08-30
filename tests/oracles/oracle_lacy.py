#!/usr/bin/env python3
"""
INDEPENDENT oracle for Lacy (1989) effective founder number.
Does not import PyPedal. Reuses the independent NRM built in oracle_meuwissen_luo.

    f_e = 1 / sum_k q_k^2

where q_k is the mean additive relationship of contributor k to the descendant
set. q is a genuine probability vector: it must sum to 1 WITHOUT being
normalised. That unforced sum-to-one is the real validity check -- dividing by
the total, as PyPedal does, makes it true by construction and therefore vacuous.

HALF-FOUNDERS -- SETTLED BY THE PAPER
------------------------------------------------------
Lacy (1989) article p.113, Methods, defines this explicitly:

    "A founder is defined as an animal with no known genetic relationship to
    any other animal in the pedigree except for its own descendants... If one
    parent is known (for example, when a wild-caught female produces an
    offspring sired in the wild), THE UNKNOWN PARENT IS CONSIDERED A FOUNDER of
    the captive population, even though that founder was perhaps never itself
    in the captive population."

That is the ``phantom`` treatment, and it is confirmed twice in the results --
in both cases as SEPARATE founder sources, not one pooled bucket:

* okapi (p.117): one unknown sire, "therefore, that unknown sire must be
  considered a 'founder' in genetic calculations"; 23 real + 1 phantom = 24
  founders, f_e = 16.05 with it and 15.39 without. Figure 1 gives the phantom
  its own bar, labelled "unk".
* Goeldi's monkey (p.118): 19 animals with unknown parents "must be treated as
  'founders'"; f_e = 24.79 with, 20.97 without. Figure 2 shows 19 SEPARATE bars.

``phantom`` is therefore the paper's rule and the default here. The other two
modes are retained ONLY as documented-wrong comparators, so that the change is
explicable; neither reproduces Lacy and neither is scientifically valid:

``mode='lacy'``     contributors are animals with no known parent. A
                    half-founder is a descendant and its unknown side belongs
                    to no contributor, so **q sums to less than 1**.

``mode='half'``     half-founders join the contributor set and leave the
                    descendant set, mirroring the historical
                    ``pyp_metrics.a_effective_founders_lacy(..., half=True)``.
                    **q sums to more than 1**, because the whole animal is
                    credited while the true founders upstream of it are
                    credited as well.

``mode='phantom'``  each unknown parent of a non-founder is named as its own
                    founder. **q sums to exactly 1** on every corpus pedigree,
                    and the mode is a no-op where there are no half-founders.

Measured on the corpus (q_sum, f_e):

    pedigree           lacy                 half                 phantom
    new_lacy.ped       1.000  2.909090909   1.000  2.909090909   1.000  2.909090909
    generations.ped    1.000  4.612612613   1.000  4.612612613   1.000  4.612612613
    boichard2a.ped     1.000  4.000000000   1.000  4.000000000   1.000  4.000000000
    mrode.ped          0.781  3.230283912   1.271  1.850602410   1.000  2.797814208
    hartlandclark.ped  0.643  7.204573215   1.868  1.990101457   1.000  5.831988609

Note the 7.204 in the 'lacy' column against hartlandclark.ped's three founders:
f_e cannot exceed the number of contributors, so that value is impossible. It is
what an ill-conditioned q does to the reciprocal.

The sum-to-one gate was, before the paper was read, only a necessary condition:
'phantom' being the sole mode producing a probability vector was strong evidence
against the other two, not a demonstration that it reproduced Lacy. p.113 closes
that gap -- the gate and the source now agree.

FOUNDER GENOME EQUIVALENTS
--------------------------
``f_g = 1 / sum(p_i^2 / r_i)`` (p.115), with ``r_i`` the expected proportion of
founder i's alleles retained in the descendant population. ``retention()``
implements the closed form from the Table 1 footnote, ``r_i = 1 - .5^x`` for a
founder with x first-generation descendants. The paper is explicit that this
holds for simple pedigrees and that complex ones need a gene-drop simulation, so
``f_g_lacy`` is valid only where that form is. It is not a general estimator and
does not pretend to be.

Validated against Lacy Appendix A (f_e = 2.91, f_g = 2.18 on the pedigree that
is this repository's new_lacy.ped) and against all seven rows of Table 1.
"""
import argparse, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_meuwissen_luo import read_pedigree, renumber

MODES = ("lacy", "half", "phantom")

#: Lacy (1989) p.113 defines the half-founder treatment explicitly, and it is
#: `phantom`. This is the paper's rule, so it is the default. `lacy` and `half`
#: are retained only as documented-wrong comparators -- neither reproduces the
#: paper and neither yields a probability vector on a pedigree that contains a
#: half-founder.
DEFAULT_MODE = "phantom"


def nrm(ped):
    n = len(ped)
    A = [[0.0]*(n+1) for _ in range(n+1)]
    sire=[0]*(n+1); dam=[0]*(n+1)
    for i,s,d in ped: sire[i],dam[i]=s,d
    for i in range(1,n+1):
        s,d = sire[i],dam[i]
        A[i][i] = 1.0 + (0.5*A[s][d] if s and d else 0.0)
        for j in range(1,i):
            A[i][j]=A[j][i]=0.5*((A[j][s] if s else 0.0)+(A[j][d] if d else 0.0))
    return A, sire, dam


def parents(ped):
    n = len(ped)
    sire = {i: 0 for i in range(1, n+1)}
    dam = dict(sire)
    for i, s, d in ped:
        sire[i], dam[i] = s, d
    return sire, dam, n


def half_founders(ped):
    """Animals with exactly one known parent."""
    sire, dam, n = parents(ped)
    return [i for i in range(1, n+1) if bool(sire[i]) != bool(dam[i])]


def add_phantom_founders(ped):
    """
    Name every unknown parent of a non-founder as its own founder.

    Returns ``(extended_pedigree, n_phantoms)``. Phantoms take IDs 1..k and the
    original animals are shifted by k, which keeps the pedigree topologically
    sorted (parents before offspring) because the originals already were.
    """
    sire, dam, n = parents(ped)
    slots = []
    for i in range(1, n+1):
        if sire[i] or dam[i]:                       # not a full founder
            if not sire[i]:
                slots.append((i, "s"))
            if not dam[i]:
                slots.append((i, "d"))
    k = len(slots)
    fill = {slot: j + 1 for j, slot in enumerate(slots)}
    ext = [(j + 1, 0, 0) for j in range(k)]
    for i in range(1, n+1):
        s = sire[i] + k if sire[i] else fill.get((i, "s"), 0)
        d = dam[i] + k if dam[i] else fill.get((i, "d"), 0)
        ext.append((i + k, s, d))
    return ext, k


def partition(sire, dam, n, mode):
    """Split 1..n into (contributors, descendants) according to `mode`."""
    full = [i for i in range(1, n+1) if not sire[i] and not dam[i]]
    if mode in ("lacy", "phantom"):
        # In 'phantom' mode the pedigree has already been completed, so every
        # unknown parent is a real founder and this is the strict definition.
        return full, [i for i in range(1, n+1) if sire[i] or dam[i]]
    if mode == "half":
        halves = [i for i in range(1, n+1) if bool(sire[i]) != bool(dam[i])]
        contributors = sorted(full + halves)
        contrib = set(contributors)
        return contributors, [i for i in range(1, n+1) if i not in contrib]
    raise ValueError("mode must be one of %r, got %r" % (MODES, mode))


def f_e_lacy(ped, mode=DEFAULT_MODE):
    """
    Returns ``(f_e, q, contributors, descendants)``. ``q`` is unnormalised, so
    ``sum(q)`` is meaningful -- see the module docstring.
    """
    if mode not in MODES:
        raise ValueError("mode must be one of %r, got %r" % (MODES, mode))
    if mode == "phantom":
        ped, _k = add_phantom_founders(ped)
    A, sire, dam = nrm(ped)
    n = len(ped)
    contributors, desc = partition(sire, dam, n, mode)
    if not desc:
        return 0.0, [], contributors, desc
    q = [sum(A[c][d] for d in desc)/float(len(desc)) for c in contributors]
    ssq = sum(x*x for x in q)
    return ((1.0/ssq) if ssq else 0.0), q, contributors, desc


def retention(ped, mode=DEFAULT_MODE):
    """
    r_i, the expected proportion of founder i's alleles retained in the
    descendant population.

    Lacy (1989) Table 1 footnote, article p.116: "for a founder with x
    first-generation descendants, r_i = 1 - .5^x". That closed form is what the
    paper uses for its simple pedigrees; for complex ones it says r_i must come
    from a gene-drop simulation, so this is NOT a general implementation and
    says so rather than pretending otherwise.
    """
    sire, dam, n = parents(ped)
    contributors, _desc = partition(sire, dam, n, mode)
    offspring = {c: 0 for c in contributors}
    for i in range(1, n + 1):
        for p in (sire[i], dam[i]):
            if p in offspring:
                offspring[p] += 1
    return {c: 1.0 - 0.5 ** offspring[c] for c in contributors}


def f_g_lacy(ped, mode=DEFAULT_MODE):
    """
    Founder genome equivalents, Lacy (1989) p.115:

        f_g = 1 / sum(p_i^2 / r_i)

    Returns ``(f_g, r)``. Only valid where the closed-form ``r_i`` above is --
    simple pedigrees in which every founder's alleles reach the descendant
    population through its first-generation offspring only.
    """
    f_e, q, contributors, desc = f_e_lacy(ped, mode=mode)
    if not desc:
        return 0.0, {}
    r = retention(ped, mode=mode)
    total = 0.0
    for c, p_i in zip(contributors, q):
        if r[c] <= 0.0:
            return 0.0, r
        total += (p_i * p_i) / r[c]
    return ((1.0 / total) if total else 0.0), r


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("pedfile")
    ap.add_argument("--animal",type=int,default=0); ap.add_argument("--sire",type=int,default=1)
    ap.add_argument("--dam",type=int,default=2)
    ap.add_argument("--mode",choices=MODES,default=DEFAULT_MODE)
    a=ap.parse_args()
    rows=read_pedigree(a.pedfile,animal=a.animal,sire=a.sire,dam=a.dam)
    ped,_back=renumber(rows)
    f_e, q, contributors, desc = f_e_lacy(ped, mode=a.mode)
    print(json.dumps({
        "pedfile": os.path.basename(a.pedfile), "n_animals": len(ped),
        "mode": a.mode,
        "n_contributors": len(contributors), "n_descendants": len(desc),
        "n_half_founders": len(half_founders(ped)),
        "q_sums_to_one_unforced": abs(sum(q)-1.0) < 1e-9,
        "q_sum": sum(q),
        "f_e_lacy": f_e,
        "f_g_lacy": f_g_lacy(ped, mode=a.mode)[0],
    }, indent=2, sort_keys=True))

if __name__=="__main__": main()
