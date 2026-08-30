# Validation and common data problems

PyPedal checks a pedigree while it loads. The goal is to catch mistakes
that would silently corrupt inbreeding or founder metrics.

## Checks that run at load

- The `pedformat` string must be valid and must match the number of
  columns on a data line.
- Duplicate animal identities in the file are dropped.
- Parents that are mentioned but have no row of their own are added.
- An animal cannot appear as both a sire and a dam.
- Parents must be orderable before their offspring. A cycle (an animal
  as its own ancestor) is refused.
- Sex codes that contradict parent roles are flagged.
- Malformed non-empty birth dates are errors, not quiet `None`.

Failures raise typed exceptions (`PyPedalError` and subclasses such as
`PyPedalUsageError` and `PyPedalPedigreeStructureError`). They do not
exit the process with status 0.

## Problems you will actually see

**Parents younger than offspring.** If birth years are present and a
parent is born after its child, treat it as a data error. PyPedal can
reorder on pedigree links; it will not rewrite your dates.

**Half-founders mistaken for founders.** One unknown parent does not
make an animal a founder. See
[IDs and missing parents](ids-and-missing-parents.md).

**File IDs used as analysis IDs.** After renumbering, `inbreeding`,
`relationship`, and `mating_coi` want current `animalID`. Map with
`ped.idmap` first.

**Dense inbreeding on a large file.** The default `tabular` method
builds a full relationship matrix. For tens of thousands of animals use
`meu_luo`. There is no automatic switch. See
[Large pedigrees](large-pedigrees.md).

**Unknown dates stored as 1800.** In PyPedal 4, unknown is `None`. 1800
and 1900 are real years unless you set a legacy import token.

**Call names used as IDs.** Display names (`n`) are not unique
identities. Use `asd` or `ASD`.

## After load

Print `ped.metadata` (or `ped.metadata.stringme()`) for counts of
records, sires, dams, and founders. If the founder count is far from
what you expect, look for half-founders and missing-parent tokens other
than `0`.
