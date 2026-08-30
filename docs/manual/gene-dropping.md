# Gene dropping

Gene dropping assigns two unique alleles to every founder (and to each
unknown parental slot that the implementation treats as a founder
genome), then lets Mendelian segregation carry those alleles down the
pedigree. After many replicates you estimate how concentrated the
founder genes have become.

PyPedal’s `pyp_metrics.effective_founder_genomes()` implements equation 2
of Boichard, Maignel and Verrier (1997):

> N<sub>g</sub> = 1 / (2 × Σ *f<sub>k</sub>*²)

Each replicate produces one N<sub>g</sub>. The function returns their
arithmetic mean.

N<sub>g</sub> is **not** “how many founders are in the file.” See
[Founder genome equivalents](founder-genome-equivalents.md).

## How to call it

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

ng = pyp_metrics.effective_founder_genomes(
    ped, rounds=20, seed=31, output=False
)
print(round(ng, 4))
```

Pass `output=False` unless you want a report file. Pass `seed` for
reproducible replicates. The default `rounds=10` is a demonstration
setting, not a scientific sample size. Assess convergence; use more
replicates for analysis.

`chrometype='autosome'` is the supported domain. Other values, including
`'sex'`, raise `PyPedalUsageError`.

## Population under study

The group whose allele frequencies are scored is the set of animals
whose input `gen` (pedformat `g`) equals the numerically largest
generation label. There is no `reference=` keyword on this function.
Without a `g` column, every animal carries the missing-generation
sentinel and the **whole pedigree** is treated as the study group. That
is why the Mrode example above is valid: it has no `g` column.

## Half-founders

An unknown parent is its own founder-genome source. A half-founder is
not collapsed into the known parent. The simulation does not index the
missing-parent token `0` as if it were an animal record.

## A large-file regression number

On the curated Griffon 2026 export of 98,001 animals, `rounds=3`,
`seed=31`, `chrometype="autosome"`, `output=False` produced
N<sub>g</sub> = **11.018378975785259**. That is a **deterministic
dataset regression**, not a scientific constant and not the number of
historical founders. See [Large pedigrees](large-pedigrees.md).
