# Lacy and Boichard methods

Lacy (1989) and Boichard, Maignel and Verrier (1997) both measure how
unevenly founders contribute to a later group. They are not the same
calculation. Neither is gene-drop N<sub>g</sub>.

| Method | Typical question | PyPedal functions |
|---|---|---|
| Lacy *f<sub>e</sub>* | Equal-founder equivalent from proportional contributions of descendants | `effective_founders_lacy` (default/scalable); `a_effective_founders_lacy` (dense NRM, small pedigrees only) |
| Boichard *f<sub>e</sub>* | Expected contributions to a named reference population | `a_effective_founders_boichard` |
| Boichard *f<sub>a</sub>* | Marginal ancestors after accounting for bottlenecks | `a_effective_ancestors_definite`, `a_effective_ancestors_indefinite` |
| Boichard N<sub>g</sub> | Simulated founder-gene frequencies | `effective_founder_genomes` |

## When to use which

- You have a small complete pedigree and want Lacy’s published
  Appendix A check: use `effective_founders_lacy` (or
  `a_effective_founders_lacy` on that small file). Expected *f<sub>e</sub>*
  on `new_lacy.ped` is about **2.91**.
  `a_effective_founders_lacy` forms a dense NRM and is not the large-file
  path.
- You are describing a **defined living group** (a birth-year cohort, a
  generation label, or an explicit ID list): use the Boichard routines
  and pass `reference=` or a numeric `gen` column.
- You want segregation, not expected contributions: use gene dropping.

Half-founders (one known parent) are not treated as two-parent founders.
PyPedal 4 uses explicit phantom founders for unknown parental genome
contributions. That is a calculation device, not a new animal in your
file.

## Reference population (Boichard)

Name the animals under study. Do not assume PyPedal will pick “the most
recent generation” from `igen`. See
[Effective ancestors](effective-ancestors.md).

## Files these functions write

Lacy and Boichard metrics write `.dat` summaries next to the pedigree
stem when `output=True` (the default, historical behaviour). Pass
`output=False` to perform the same calculation without writing that
analysis file. Run them in a temporary directory, or set the file tag,
if you do not want those files in your project tree.

[Effective founders](effective-founders.md),
[Effective ancestors](effective-ancestors.md), and
[Gene dropping](gene-dropping.md) split the three questions. This page
exists so the papers stay distinct.
