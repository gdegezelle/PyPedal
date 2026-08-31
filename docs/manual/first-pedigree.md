# Your first pedigree

A pedigree is a list of animals and their parents. In PyPedal each line of
a text file is one animal. The simplest file has three columns:

- **animal** — this animal’s identifier
- **sire** — the father’s identifier
- **dam** — the mother’s identifier

Unknown parents are written as **`0`**. That zero is a missing-parent
token, not an animal in the file.

## The Mrode textbook pedigree

Mrode (2005) Table 2.1 is a six-animal example used throughout PyPedal
tests. Animals 1 and 2 have no recorded parents. Animal 5 is the inbred
individual you will analyse in the next chapter.

This example writes the file itself. It does **not** need
`PyPedal/examples` and therefore works from an installed wheel.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree

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

print(len(ped.pedigree))
for animal in ped.pedigree:
    print(animal.animalID, animal.sireID, animal.damID)
```

That prints `6`, then the six records in oldest-to-youngest order.

## What the columns mean

`pedformat="asd"` means three integer columns: animal, sire, dam, separated
by spaces. Lines that start with `#` are comments.

| Animal | Sire | Dam | In ordinary language |
|---|---|---|---|
| 1 | 0 | 0 | Founder (both parents unknown) |
| 2 | 0 | 0 | Founder |
| 3 | 1 | 2 | Offspring of 1 and 2 |
| 4 | 1 | 0 | Sire known, dam unknown (a half-founder) |
| 5 | 4 | 3 | Offspring of 4 and 3 |
| 6 | 5 | 2 | Offspring of 5 and 2 |

Animal 4 is **not** a founder. A founder in PyPedal has **both** parents
unknown. Animal 4 has one known parent.

## What `load_pedigree` does

`load_pedigree()` builds a `NewPedigree` and calls `load()`. By default it
also **renumbers** the pedigree: animals are ordered so that parents come
before offspring, and each animal receives a sequential `animalID` starting
at 1. The identifier that appeared in the file is kept as `originalID`.

On this example the file already used 1…6, so `animalID` and `originalID`
happen to match. That is a convenience of Mrode’s table, not a general
rule. See [IDs and missing parents](ids-and-missing-parents.md).

Constructing a pedigree opens a logfile next to the pedigree path. The example
uses a temporary directory so that log stays out of your project tree.

`messages="quiet"` and `pedigree_summary=0` keep the console quiet when
the host application has not configured Python logging. The defaults are
verbose.

## Next

[Your first analysis](first-analysis.md) calculates inbreeding for animal 5.
The expected result is *F* = 0.125.
