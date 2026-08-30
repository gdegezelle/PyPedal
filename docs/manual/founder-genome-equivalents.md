# Founder genome equivalents

Gene dropping estimates how many **equally frequent founder genes** would
match the simulated allele frequencies in a later group of animals. The
usual symbol in Boichard, Maignel and Verrier (1997) is **N<sub>g</sub>**.

> N<sub>g</sub> = 1 / (2 × Σ *f<sub>k</sub>*²)

where *f<sub>k</sub>* is the realised frequency of founder gene *k* after
Mendelian segregation (MacCluer et al., 1986). Each replicate produces
one N<sub>g</sub>; PyPedal returns the **arithmetic mean** across
replicates.

**N<sub>g</sub> is not a founder head-count.** A small N<sub>g</sub> on a
large pedigree means the simulated founder-allele frequencies are as
uneven as if few equally represented genomes remained. Many historical
founders can still sit in the file.

## How this differs from related numbers

| Quantity | Why it is different |
|---|---|
| Raw founder count | Head-count of animals with two unknown parents |
| Lacy *f<sub>e</sub>* / Boichard *f<sub>e</sub>* | Expected contributions, not simulated allele frequencies |
| Lacy founder genome equivalents *f<sub>g</sub>* | Lacy (1989) uses allele **retention**. PyPedal 4.0 does **not** implement *f<sub>g</sub>* |

On Lacy’s Appendix A pedigree the published *f<sub>g</sub>* is 2.18; the
N<sub>g</sub> for that example is about 1.84. Do not treat the names as
interchangeable.

The simulation function is
`pyp_metrics.effective_founder_genomes()`. Details of seeds, rounds, and
the Griffon regression number are in [Gene dropping](gene-dropping.md)
and [Large pedigrees](large-pedigrees.md).
