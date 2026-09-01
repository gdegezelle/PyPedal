# Inbreeding

An inbreeding coefficient *F* is the probability that the two alleles an
animal carries at a locus are identical by descent. It lies in the closed
interval from 0 to 1.

PyPedal computes *F* from the recorded pedigree. It is not a genomic
estimator. The public function is `pyp_nrm.inbreeding()`.

IDs in this chapter are **current `animalID`** values (1-based after the
default load).

## How to call it

This example is wheel-safe: it writes Mrode’s six-animal pedigree itself.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_nrm

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

result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
print(result["fx"][5])
```

That prints `0.125`. Animal 5 is the inbred individual in Mrode (2005)
Table 2.1.

Pass `output=False` unless you want a coefficients file written next to
the pedigree. The default is `output=True`.

## What is returned

| Key | Content |
|---|---|
| `fx` | Mapping current `animalID` → *F* |
| `metadata` | Summary statistics under `all` and `nonzero` |
| `rel_dict` | Present only when `rels=1` **and** the method can compute relationships |

Look animals up with current IDs, not original file IDs.

The return is an `InbreedingResult`, a `dict` subclass. Existing key access
is unchanged and remains the supported 4.x style: `result["fx"][5]`.
`result.fx` is the same mapping, not a copy. `isinstance(result, dict)`
is still true. `rel_dict` is omitted from the mapping when relationships
were not computed; `result.rel_dict` is then `None`.

## Supported methods

| `method` | What it does | Relationships (`rels=1`) |
|---|---|---|
| `'tabular'` | Default. Builds a numerator relationship matrix and reads *F* from the diagonal as *a<sub>ii</sub>* − 1 | yes |
| `'vanraden'` | VanRaden (1992) pedigree algorithm | yes |
| `'meu_luo'` | Meuwissen and Luo (1992). Coefficients only; linear in the number of animals | no |
| `'mod_meu_luo'` | Meuwissen–Luo as modified by Quaas (1995), following Mrode (2005) Appendix B.2 | no |
| `'aguilar'` | Shells out to the external `inbupgf90` binary | no |

On Mrode, `tabular`, `vanraden`, `meu_luo`, and `mod_meu_luo` all give
animal 5 *F* = 0.125.

`'vanraden'` here is the **1992 pedigree** method. It is not VanRaden
(2008) genomic Method 1.

`'aguilar'` requires the `inbupgf90` program on `PATH`. That is an
external tool, not a native PyPedal implementation.

## Choosing a method

- **Small pedigrees** (hundreds to a few thousand animals): `tabular` is
  simple and also yields relationships when you ask for them.
- **Large pedigrees**: use `meu_luo`. Do not form a dense numerator
  relationship matrix for ~100,000 animals.

There is **no automatic switch** at 10,000 animals. The default remains
`tabular`. See [Large pedigrees](large-pedigrees.md).

## Reading the metadata

`metadata["all"]` summarises every animal. `metadata["nonzero"]`
summarises animals with computed *F* > 0 (the exact stored sign, not a
tolerance). On Mrode the all-animal mean is small because five of six
animals are not inbred; the nonzero mean is 0.125.

Stored coefficients in `fx` and in the metadata dictionary are left
unrounded. On a deep or large pedigree, a few units in the last place of
floating-point residue around zero can occur. Displayed summaries round
those values for readability (six decimal places; a few-ULP negative
residue prints as `0.000000`, not `-0.000000`).

For prospective matings see [Test mating](mating.md). For pairwise
*a<sub>ij</sub>* see [Relationships](relationships.md).
