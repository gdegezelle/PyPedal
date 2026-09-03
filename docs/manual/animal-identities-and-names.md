# Animal identities and names

PyPedal distinguishes **who the animal is** from **what it is called**.

## Integer IDs

With `pedformat='asd'`, the animal, sire, and dam columns are integers.
After the default load:

- `originalID` is the number from the file
- `animalID` is the sequential 1-based ID assigned by renumbering

Look up file IDs with `ped.idmap`. Analysis functions take current
`animalID` values.

## Unique string identities

With `pedformat='ASD'`, the three identity columns are unique strings
(registration numbers, herd-book IDs, and similar). PyPedal hashes each
string to an integer and stores the original string.

`ped.namemap` / `ped.namebackmap` convert between the string identity and
the hashed integer. A hash collision is an error, not a silent merge.

These strings **are** identities. They must be unique.

## Display and call names

Pedformat `n` is a display or call name. It is **not** unique. Two dogs
may both be called “Max”. PyPedal will not use `n` to find an animal in
`mating_coi` or `relationship`.

If you omit `n`, the animal ID is used as the name field.

The desktop Relationship and Mating selectors search call names, original
IDs, and current IDs. They do not use `ped.namemap`. Duplicate names
require an explicit choice.

## Sire and dam names

`S` / `D` (uppercase) are unique string **identities** for sire and dam.
`sireName` / `damName` on the animal object are name fields, not a second
identity system.

## Worked distinction

| You have | Use |
|---|---|
| Studbook number `GB-1234` as the only ID | `ASD` |
| Integer animal numbers | `asd` |
| Integer IDs plus a kennel name | `asd` plus `n` |
| Looking up inbreeding | current `animalID` (or map from `originalID` first) |

See [IDs and missing parents](ids-and-missing-parents.md) for founders
and half-founders.
