# Relationships

Three related quantities are easy to mix up:

| Quantity | Symbol | What it is |
|---|---|---|
| Inbreeding coefficient | *F<sub>i</sub>* | Probability that animal *i*’s two alleles are IBD. See [Inbreeding](inbreeding.md). |
| Additive relationship | *a<sub>ij</sub>* | Numerator-relationship-matrix entry for the pair (*i*, *j*). Diagonal *a<sub>ii</sub>* = 1 + *F<sub>i</sub>*. |
| Wright’s relationship | *r<sub>ij</sub>* | *a<sub>ij</sub>* scaled by the inbreeding of both animals. PyPedal’s pairwise helper does **not** return this scaled form. |

The **numerator relationship matrix** (NRM, sometimes “the A matrix”) is
the square array of all *a<sub>ij</sub>*. IDs in this chapter are current
`animalID` values.

## Pairwise relationship

`pyp_metrics.relationship(anim_a, anim_b, pedobj)` returns the additive
relationship for two **existing** current IDs.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics
from PyPedal.pyp_errors import PyPedalUsageError

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

print(pyp_metrics.relationship(4, 3, ped))
print(pyp_metrics.relationship(1, 2, ped))
```

On Mrode that prints `0.25` then `0.0`.

- Animals 4 and 3 share sire 1 and have no recorded dam in common: a
  half-sib relationship of **0.25**.
- Animals 1 and 2 are unrelated founders: a relationship of **0.0**.

**Zero is a valid result.** It means these two existing animals are
unrelated in this pedigree. It does not mean “an ID could not be
resolved,” and it does not mean “the relationship matrix could not be
formed.” Missing IDs raise `PyPedalUsageError`. A computational or
allocation failure raises `PyPedalError`.

An ID that is missing, non-integral, or an original file ID that is not
also a current ID raises `PyPedalUsageError`.

The desktop app does not ask you to type current IDs. Relationship
searches by name, original ID, or current ID, then calls
`relationship()` with the two selected current IDs. Duplicate names
require an explicit choice; names are not identities. See
[Desktop application](desktop-application.md).

```python
try:
    pyp_metrics.relationship(99999, 3, ped)
except PyPedalUsageError:
    print("unresolved")
```

Pairwise `relationship()` computes an exact selected *a<sub>ij</sub>*
from the ancestors of the requested animals. It does **not** build an
ancestor sub-pedigree numerator relationship matrix for one pair, and it
does not attach a pedigree-wide NRM. On a large studbook that makes a
single pair suitable without forming a dense NRM. Very large explicit
mating *groups* still scale with the number of pairs.

## The full matrix

`pyp_nrm.inbreeding(..., rels=1)` can return relationship summaries when
the method supports them (`tabular` and `vanraden`). Forming the whole
NRM is expensive. For a single pair, use `relationship()`. For a
prospective offspring, use [Test mating](mating.md) — *F* of the
offspring of *i* and *j* is *a<sub>ij</sub>* / 2.

Do not build a dense NRM for a large population. See
[Large pedigrees](large-pedigrees.md).
