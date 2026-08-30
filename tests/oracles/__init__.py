"""
Thin re-export of independent scientific oracles.

The oracles deliberately do **not** import PyPedal: they parse the pedigree
themselves and implement the mathematics from first principles.

Available names:

    read_pedigree, renumber              pedigree parsing / topological renumber
    inbreeding_meuwissen_luo             Meuwissen & Luo (1992) recursion
    inbreeding_quaas                     Quaas (1995) ancestor-list traversal
    inbreeding_tabular                   direct O(n^2) NRM recurrence
    lacy_nrm, lacy_f_e                   Lacy (1989) effective founder number
    lacy_f_g, lacy_retention             founder genome equivalents and r_i
    lacy_half_founders                   animals with exactly one known parent
    lacy_n_half_founders                 the same, counted, from a pedfile
    LACY_MODES                           the three half-founder treatments
    ballou_f_a                           Ballou ancestral inbreeding, from a pedfile
    ballou_ancestral_inbreeding_oracle   the same, over a renumbered pedigree
    boichard_read                        pedfile -> (rows, generations)
    boichard_f_e                         Appendix A + eq. 1
    boichard_f_a                         Appendix B
    boichard_marginal_contributions      the shared Appendix-B engine
    boichard_bounds                      the pp.9-10 lower/upper bounds
    boichard_founders                    founders INCLUDING half-founders (p.7)
    boichard_unknown_parent_slots        one entry per unknown parental SIDE
    boichard_phantom_complete            Reading C's pedigree transform
    BOICHARD_TIE_BREAKS                  the R1 conventions -- NOT source rules
    BOICHARD_HALF_FOUNDER_RULES          the R2 readings -- 'complete' is the
                                         adjudicated one; source is still silent
    BOICHARD_PHANTOM_ENCODINGS           the two meaningless phantom numberings
    BOICHARD_REFERENCE_POP_RULES         the R3 conventions -- source is silent
    exact_ng, mc_ng                      Boichard eq. 2 N_g -- 
    ng_distribution                      the EXACT distribution of SUM f_k^2
    ng_distribution_bruteforce           the same, by 2**b enumeration -- the
                                         definition the fast path is checked against
    ng_check_bounds                      0.5 <= N_g <= f_e_realised, per outcome
    ng_sum_f2, ng_gene_diversity         the statistic and Lacy's H
    ng_founder_genes, ng_true_founders   the 2f founder genes; both-unknown founders
    NG_HALF_FOUNDER_MODES                'phantom' and 'slot' -- FG-10 shows these
                                         are exactly distributionally equivalent

``exact_ng`` is BOICHARD's N_g and ``lacy_f_g`` is LACY's f_g. They are NOT the
same quantity and must never be substituted for one another -- Lacy's own
Table 1 has them disagreeing by 8-11% on the same pedigrees. See
``the algorithm notes``.

``ballou_f_a`` is **paper-validated**: it implements equation 1 of Ballou (1997)
article p.170, read directly, and reproduces the values printed in that paper's
Figure 1 and in Suwanlee et al. (2007) Figure 1. The earlier secondary grounding
(``the historical manual`` line 88) proved correct -- the recurrence
did not change when the paper arrived, only its status did.

It must still never be compared against ``dropped_ancestral_inbreeding`` as a
correctness check. That is a different estimator, and the divergence is now a
published fact rather than a caution: Suwanlee Figure 1 prints both on one
pedigree and they differ from the third generation onward. Nor may the reverse
ordering become an invariant -- Suwanlee proves no theorem. See .
"""
import os
import sys

_ORACLE_DIR = os.path.abspath(os.path.dirname(__file__))

if _ORACLE_DIR not in sys.path:
    sys.path.insert(0, _ORACLE_DIR)

from oracle_meuwissen_luo import (      # noqa: E402
    read_pedigree,
    renumber,
    inbreeding_meuwissen_luo,
    inbreeding_quaas,
    inbreeding_tabular,
)
from oracle_ballou import (              # noqa: E402
    ancestral_inbreeding as ballou_ancestral_inbreeding_oracle,
    ancestral_inbreeding_from_file as ballou_f_a,
)
from oracle_lacy import (                # noqa: E402
    nrm as lacy_nrm,
    f_e_lacy as _f_e_lacy,
    f_g_lacy as lacy_f_g,
    retention as lacy_retention,
    half_founders as lacy_half_founders,
    MODES as LACY_MODES,
    DEFAULT_MODE as LACY_DEFAULT_MODE,
)
from oracle_boichard import (            # noqa: E402
    read_pedigree as boichard_read,
    check_topological as boichard_check_topological,
    founders as boichard_founders,
    half_founders as boichard_half_founders,
    unknown_parent_slots as boichard_unknown_parent_slots,
    phantom_complete as boichard_phantom_complete,
    f_e as boichard_f_e,
    f_a as boichard_f_a,
    marginal_contributions as boichard_marginal_contributions,
    bounds as boichard_bounds,
    TIE_BREAKS as BOICHARD_TIE_BREAKS,
    HALF_FOUNDER_RULES as BOICHARD_HALF_FOUNDER_RULES,
    PHANTOM_ENCODINGS as BOICHARD_PHANTOM_ENCODINGS,
    REFERENCE_POP_RULES as BOICHARD_REFERENCE_POP_RULES,
)
from oracle_founder_genomes import (     # noqa: E402
    ExactEnumerationTooLarge as NgExactEnumerationTooLarge,
    HALF_FOUNDER_MODES as NG_HALF_FOUNDER_MODES,
    check_bounds as ng_check_bounds,
    distribution as ng_distribution,
    distribution_bruteforce as ng_distribution_bruteforce,
    exact_ng,
    f_e_realised as ng_f_e_realised,
    founder_genes as ng_founder_genes,
    gene_diversity as ng_gene_diversity,
    mc_ng,
    n_binary_choices as ng_n_binary_choices,
    sum_f2 as ng_sum_f2,
    true_founders as ng_true_founders,
)


def lacy_f_e(pedfile, sepchar=" ", animal=0, sire=1, dam=2, mode=LACY_DEFAULT_MODE):
    """
    Lacy (1989) effective founder number, computed independently of PyPedal.

    Returns ``(f_e, q_sums_to_one_unforced)``. The second value is the validity
    gate: q is a genuine probability vector and must sum to 1 *without* being
    normalised. PyPedal divides by the total, which makes sum-to-one true by
    construction and therefore vacuous as a check.

    ``mode`` selects the half-founder treatment. The default is ``'phantom'``,
    which names each unknown parent as its own founder -- **the rule Lacy (1989)
    p.113 states**. ``'lacy'`` (strict) and ``'half'`` (mirrors PyPedal's
    historical ``half=True``) are retained only as documented-wrong comparators.
    On a pedigree with no half-founders all three agree and the gate holds; on
    one with half-founders only ``'phantom'`` satisfies it. See the oracle
    module docstring and audit item .
    """
    rows = read_pedigree(pedfile, sepchar=sepchar, animal=animal, sire=sire, dam=dam)
    ped, _back = renumber(rows)
    f_e, q, _contributors, desc = _f_e_lacy(ped, mode=mode)
    if not desc:
        return 0.0, False
    return f_e, abs(sum(q) - 1.0) < 1e-9


def lacy_n_half_founders(pedfile, sepchar=" ", animal=0, sire=1, dam=2):
    """Count animals with exactly one known parent."""
    rows = read_pedigree(pedfile, sepchar=sepchar, animal=animal, sire=sire, dam=dam)
    ped, _back = renumber(rows)
    return len(lacy_half_founders(ped))


def oracle_inbreeding(pedfile, sepchar=" ", animal=0, sire=1, dam=2):
    """
    Inbreeding coefficients keyed by *original* animal ID, cross-checked across
    all three independent routes. Raises if the routes disagree -- an oracle
    that does not agree with itself must not be used to adjudicate anything.
    """
    rows = read_pedigree(pedfile, sepchar=sepchar, animal=animal, sire=sire, dam=dam)
    ped, back = renumber(rows)
    routes = {
        "meuwissen_luo": inbreeding_meuwissen_luo(ped),
        "quaas": inbreeding_quaas(ped),
        "tabular": inbreeding_tabular(ped),
    }
    ref = routes["meuwissen_luo"]
    for name, vals in routes.items():
        for i in range(1, len(ped) + 1):
            if abs(vals[i] - ref[i]) > 1e-12:
                raise AssertionError(
                    "oracle routes disagree on animal %d: %s=%r meuwissen_luo=%r"
                    % (i, name, vals[i], ref[i])
                )
    return {back[i]: ref[i] for i in range(1, len(ped) + 1)}
