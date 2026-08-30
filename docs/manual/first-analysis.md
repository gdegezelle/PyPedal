# Your first analysis

An **inbreeding coefficient** *F* is the probability that the two copies
of a gene an animal carries are identical by descent — they both come
from the same ancestor. *F* lies between 0 and 1.

- *F* = 0 means the pedigree shows no inbreeding.
- *F* = 0.125 means 12.5% inbreeding, the value for offspring of
  half-sibs in a simple pedigree.
- *F* = 1 would mean complete identity by descent.

PyPedal estimates *F* from the **recorded pedigree**, not from genotypes.

## Inbreeding on the Mrode pedigree

Animal 5 in Mrode (2005) Table 2.1 is inbred. The published coefficient
is **0.125**. This script writes the pedigree, loads it, and prints *F*
for animal 5.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_nrm

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("""\
# Pedigree from Mrode (2005) Table 2.1
1 0 0
2 0 0
3 1 2
4 1 0
5 4 3
6 5 2
""")

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
print(len(ped.pedigree))
print(result["fx"][5])
```

That prints `6` and `0.125`.

`result["fx"]` maps **current `animalID` → *F***. After a default load,
animal 5 is the fifth animal in the ordered pedigree.

Always pass `output=False` unless you want PyPedal to write a coefficients
file next to the pedigree.

## Which method to use

`method="tabular"` builds a numerator relationship matrix and reads *F*
from the diagonal. It is appropriate for **small** pedigrees such as this
one.

For **large** pedigrees use `method="meu_luo"` (Meuwissen and Luo, 1992).
That method is linear in the number of animals and does not build a full
matrix.

There is **no automatic switch** at 10,000 animals. If you leave the
default `tabular` method on a file of tens of thousands of animals, you
can exhaust memory. See [Large pedigrees](large-pedigrees.md).

## What the result contains

| Key | Meaning |
|---|---|
| `fx` | *F* for each current `animalID` |
| `metadata` | Summary statistics under `all` (every animal) and `nonzero` (inbred animals only) |

On Mrode, only animal 5 has *F* > 0, so the nonzero count is 1 and the
nonzero mean is 0.125.

## Next

- [IDs and missing parents](ids-and-missing-parents.md) — why `animalID`
  is not always the number in your file
- [Inbreeding](inbreeding.md) — all supported methods
- [Recipes](recipes.md) — copy-and-paste workflows
