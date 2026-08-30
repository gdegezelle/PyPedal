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
