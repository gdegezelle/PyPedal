# Genomic methods

PyPedal 4.0 can attach SNP genotypes to a pedigree and form a **VanRaden
(2008) Method 1** genomic relationship matrix (GRM). That path is
independent of pedigree inbreeding.

## Two different VanRaden methods

| What you call | What it is |
|---|---|
| `pyp_nrm.inbreeding(ped, method="vanraden")` | **VanRaden (1992)** pedigree algorithm for *F* from the recorded pedigree. See [Inbreeding](inbreeding.md). |
| `pyp_snp.form_grm_from_snp` / `compute_genomic_inbreeding_from_grm` | **VanRaden (2008)** Method 1 genomic relationships from SNP dosages. |

They share an author name and nothing else. `method="vanraden"` on
`pyp_nrm.inbreeding` does **not** build a GRM.

## SNP file and load

`pyp_snp.load_snp_file(pedobj, snpfile=...)` reads whitespace-separated
lines:

```
animalID  chip_type  n_snps  genotype
```

`genotype` is a string of per-locus dosages, one character each: `0`,
`1`, or `2` (copies of the reference allele). Animals are matched on
`originalID`. The table is attached as `ped.snp`.

Set `kw['snpfile']` at load to read genotypes **before** IDs are
remapped; `renumber_snp_ids()` then follows a successful renumber. You
can also call `load_snp_file` after load.

Comments (`#`) are skipped. `sepchar` is accepted and ignored; the file
is whitespace-separated.

## Supported genomic APIs

| Function | Role |
|---|---|
| `pyp_snp.load_snp_file` | Parse and attach genotypes |
| `pyp_snp.form_p_matrix_from_snp` | Allele frequencies *p* (one per locus) |
| `pyp_snp.form_m_matrix_from_snp` | Centred marker intermediate used to build *Z* |
| `pyp_snp.form_grm_from_snp` | GRM *G* = *ZZ*' / (2 Σ *p*(1−*p*)) |
| `pyp_snp.compute_genomic_inbreeding_from_grm` | *F*<sub>g</sub> = *G<sub>jj</sub>* − 1 |
| `pyp_snp.compute_genomic_homozygosity_from_snp` | Proportion of typed 0/2 loci, in [0, 1] |

VanRaden (2008) Methods 2 and 3 are **not** implemented.
`form_grm_from_snp(method=...)` other than `1` raises `PyPedalUsageError`.

`form_m_matrix_from_snp(scale_m=False)` is **refused**
(`PyPedalUsageError`). Skipping the centring is not VanRaden’s *Z*.

If every locus is monomorphic, Method 1 is undefined
(`2 Σ p(1−p) = 0`) and `form_grm_from_snp` raises
`PyPedalValidationError`.

### Allele frequencies

`base_frequencies` supplies *p* from the unselected base population, as
VanRaden (2008) p. 4416 asks for. If omitted, frequencies are
**estimated from the genotyped sample**. The paper notes that estimation
can bias genomic inbreeding; the functions log when they take that
fallback.

### Genomic inbreeding from the GRM

```python
import os
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_snp

work = Path(tempfile.mkdtemp())
(work / "ped.ped").write_text("1 0 0\n2 0 0\n3 1 2\n4 1 2\n")
(work / "geno.txt").write_text(
    "1 chip1 4 0120\n"
    "2 chip1 4 1201\n"
    "3 chip1 4 2012\n"
    "4 chip1 4 0000\n"
)
os.chdir(work)

ped = load_pedigree(
    options={
        "pedfile": str(work / "ped.ped"),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
        "snpfile": str(work / "geno.txt"),
    }
)

grm = pyp_snp.form_grm_from_snp(ped)
fg = pyp_snp.compute_genomic_inbreeding_from_grm(ped, grm=grm)
hom = pyp_snp.compute_genomic_homozygosity_from_snp(ped)
print(round(fg[1], 4), round(hom[4], 2))
```

That prints `0.4667` and `1.0` on this four-animal, four-locus toy
panel: genomic *F*<sub>g</sub> for animal 1, and complete homozygosity
for animal 4 (all four loci are `0`).

The keyword is **`grm=`**. There is no supported `g_matrix=` argument.

`compute_genomic_inbreeding_from_grm` returns `{animalID: F_g}`.
*F*<sub>g</sub> is bounded below by −1 and has **no generic finite upper
bound** (an animal homozygous for rare alleles can be arbitrarily large).
Do not apply pedigree-*F* range checks to it. With `store=True` (the
default), values are written onto `NewAnimal.genomic_inbreeding`.

Homozygosity is a direct count of homozygous typed loci. It is **not**
VanRaden *F*<sub>g</sub>. An animal with no typed loci gets
`kw['missing_homozygosity']`, not 0.0.

## AGIL helpers versus `read_agil_pedigree_file`

These file helpers exist in `pyp_snp`:

- `read_agil_chromosome_data`
- `read_agil_genotypes_txt`
- `read_agil_true_frequency`

`read_agil_pedigree_file` is **not a 4.0 function**. Do not invent a
callable stub for it. AGIL *file* helpers above are not a pedigree loader.

## Aguilar / `inbupgf90`

`pyp_nrm.inbreeding(method="aguilar")` is an **external** route: it
shells out to the `inbupgf90` binary. If the binary is missing, the call
raises `FileNotFoundError`.

That is **not** a native, independently validated PyPedal implementation
of Aguilar’s method. It depends on software and files outside this
package. It is not demonstrated in this manual. Status:
[Limitations](limitations.md).

## Related pages

- [Inbreeding](inbreeding.md) — pedigree *F*, including VanRaden (1992)
- [Pedigree objects](object-model.md) — `ped.snp`
- [Recipes](recipes.md) — `new_snp.py` / `new_snp2.py` (`grm=`)
- [Limitations](limitations.md) — AGIL is not a pedigree loader; Aguilar; Method 1 only
