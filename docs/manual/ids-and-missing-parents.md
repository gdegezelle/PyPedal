# Understanding IDs and missing parents

PyPedal keeps more than one identifier on every animal. Mixing them is
the usual way to look up the wrong parent or print the wrong inbreeding
coefficient.

## Two identity domains

After a default load (`renumber=True`):

| Name | What it is |
|---|---|
| `originalID` | The identifier as it appeared in the file (or the integer hash of a unique string identity) |
| `animalID` | Current internal identity: sequential **1-based** IDs after renumbering |
| `ped.pedigree` index | `animalID - 1` while the list is in renumbered order |
| `ped.idmap` | `originalID → animalID` |
| `ped.backmap` | `animalID → originalID` |

Inbreeding results, relationship calls, and parent slots (`sireID`,
`damID`) use **current `animalID`**. The number in your studbook file is
`originalID`.

On the Mrode example those two numbers happen to be the same. Do not
rely on that. A file that starts at animal 1000, or a string-identity
file, will not match.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree

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

animal = ped.pedigree[4]  # list index 4 is animalID 5
print(animal.animalID, animal.originalID)
print(ped.idmap[animal.originalID], ped.backmap[animal.animalID])
```

## Names are not identities

Integer files use `asd` (animal, sire, dam as numbers). String-identity
files use `ASD`. The uppercase codes are **unique identities**, hashed to
integers for internal use. They are not nicknames.

A display or call name is a separate column (`n`). Two animals may share
a call name. PyPedal will not treat a call name as an ID.

See [Animal identities and names](animal-identities-and-names.md).

## Missing parents

The missing-parent token defaults to integer **`0`**
(`missing_parent`). A sire or dam equal to that token means **the parent
is unknown in this file**. PyPedal does not invent a biological parent.

| Parents in the file | `founder` flag | Ordinary language |
|---|---|---|
| Both unknown | `'y'` | Founder |
| Exactly one unknown | `'n'` | Half-founder |
| Both known | `'n'` | Non-founder |

In the Mrode pedigree, animal 4 has sire 1 and dam 0. It is a
half-founder. `ped.metadata.num_unique_founders` does **not** count
animal 4.

Unknown birth dates are a different concept (`None` on `by` / `bd`). An
unknown date does not mark an unknown parent.

## Renumbering, in user terms

Many calculations need parents to appear before their offspring, with
IDs running 1…*n*. PyPedal does that automatically. You can keep thinking
in file IDs by using `originalID` and `ped.idmap`. You cannot feed a file
ID to `inbreeding()["fx"]` or `mating_coi` unless it is also the current
`animalID`.

If a pedigree contains a cycle (an animal appearing as its own ancestor)
or another structure that cannot be ordered, load **refuses** rather than
guessing an order.

Details of files and format codes are in [Pedigree input](pedigree-input.md).
