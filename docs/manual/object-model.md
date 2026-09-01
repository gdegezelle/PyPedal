# Object model

A loaded pedigree is a `NewPedigree`. `ped.pedigree` is a list of
`NewAnimal` records. Options live on `ped.kw`. Summary counts live on
`ped.metadata`. Identity maps live on the pedigree object.

## Construct, then load

`load_pedigree(options={...})` constructs a `NewPedigree` and calls
`load()`. `loadPedigree` is the same function under the older name.

After a default file load:

1. Records are read and checked against `pedformat`.
2. Animals are ordered so parents precede offspring.
3. Sequential 1-based `animalID` values are assigned (`renumber=True`).
4. `ped.metadata` is a snapshot of counts.
5. Optional flags (`set_generations`, forming an NRM, and similar) run
   only when requested.

Analyses assume that pipeline succeeded.

## Identity on each animal

| Attribute | Meaning |
|---|---|
| `animalID` | Current 1-based ID |
| `originalID` | File identity (or hash of a string identity) |
| `renumberedID` | Mirrors `animalID` after a successful renumber |
| `sireID` / `damID` | Current parent IDs, or `missing_parent` |
| `name` | Display name or the identity string |
| `founder` | `'y'` only when both parents are unknown |
| `by` / `bd` | Recorded year / date, or `None` |
| `gen` | Input generation label |
| `igen` | Inferred depth, after `set_generation` |
| `fa` | Inbreeding coefficient field |
| `age` | Legacy year-offset, not biological age |

## Maps on the pedigree

| Map | Direction |
|---|---|
| `idmap` | `originalID → animalID` |
| `backmap` | `animalID → originalID` |
| `namemap` / `namebackmap` | Unique string identities |

`ped.pedigree[animalID - 1]` is the animal while the list is still in
renumbered order.

## Data state: factual, derived, and cached

A loaded pedigree mixes three kinds of data. The same attribute name can
move from one role to another. PyPedal 4.1.0 does not split those roles
into separate fields.

### Factual / loaded

These values come from the input file (or from an equivalent stream) and
are not computed by an analysis:

- original IDs (`originalID`; or the hash of a unique string identity)
- sire and dam IDs as recorded (`sireID` / `damID`, or `missing_parent`)
- names, sex, breed, user field, and similar optional columns
- recorded birth year / date (`by` / `bd`), including unknown as `None`
- loaded inbreeding `animal.fa` when `pedformat` contains `f`
- loaded genomic inbreeding `animal.fg` when `pedformat` contains `G`
- SNP genotypes when `snpfile` is set (not from a `P` column)

Input generation `gen` (pedformat `g`) is also factual: it is a label
from the file, not the inferred depth `igen`.

### Derived

These values are inferred from pedigree structure or from load-time
options. They can change if the pedigree is reordered, renumbered, or
mutated:

- founder status (`founder == 'y'` only when both parents are unknown)
- half-founder status (exactly one parent unknown; not a stored field)
- implicit parents materialized because they were cited but had no row
- inferred generation `igen` after `set_generation`
- pedigree completeness (`pedcomp`) after that calculation
- estimated chronology ranges after `estimate_birth_dates`
- renumbered / padded identifiers (`animalID`, `paddedID`, `renumberedID`)
- identity maps (`idmap`, `backmap`, `namemap`)

Metadata counts (`num_unique_founders`, `num_unknown_birth_years`, and
similar) are a snapshot of the pedigree at load (or when metadata is
rebuilt). They are derived, not a second copy of the animal records.

### Cached / computed

These values are analysis results stored on the pedigree or on animals.
They are not automatically kept in sync if you later change the pedigree:

- `animal.fa` after `inbreeding()` (see the caveat below)
- `kw["f_computed"]` and `kw["g_computed"]`
- `ped.nrm` when a numerator relationship matrix has been formed

Deleting or merging animals sets `ped.nrm = None`. It does not restore a
previously loaded `fa` column.

### `animal.fa` can change role

`fa` is one field with two lives:

1. After load, if the file had an `f` column, `fa` is the **loaded**
   coefficient and `kw["f_computed"]` is set `True`.
2. After `inbreeding()`, `fa` is overwritten with the **computed**
   coefficient for that animal. PyPedal 4.1.0 keeps a single `fa` field;
it does not store a parallel loaded copy.

If you need the file value after computing inbreeding, keep your own
copy before calling `inbreeding()`.

## Other classes

`PedigreeMetadata` holds counts (records, sires, dams, founders).
`NewAMatrix` holds a numerator relationship matrix when one is formed.
`LightAnimal` is a slimmer record for graph routines. `SimAnimal` is
used internally by pedigree simulation.

## Mutations

`delete_animals` is atomic: it either removes the requested original IDs
or raises. It will not delete a parent still referenced by a surviving
child, and it does not rewrite orphans to `0`. `merge_animals(keep, drop)`
is the explicit redirect path. Either operation sets `ped.nrm = None`.

See [IDs and missing parents](ids-and-missing-parents.md).
