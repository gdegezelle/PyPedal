# File formats and pedformat

`pedformat` is one character per input column. Codes are
**case-sensitive**: `a` is not `A`.

A pedigree must identify animal, sire, and dam: `a` or `A`, `s` or `S`,
and `d` or `D`.

The checkout file `PyPedal/PEDIGREE_FORMAT_CODES.txt` lists the same
codes. See [Pedigree format codes](pedigree-format-codes.md) for the
full table.

## Integer identity (`asd`)

| Code | Meaning |
|---|---|
| `a` | Animal ID (integer) |
| `s` | Sire ID (integer) |
| `d` | Dam ID (integer) |

Missing parents use `missing_parent` (default `0`).

## String identity (`ASD`)

`A`, `S`, and `D` are **unique string identities**. PyPedal hashes them
to integers for internal use and keeps the original string. They are not
call names.

| Code | Meaning |
|---|---|
| `A` | Unique animal identity (string) |
| `S` | Unique sire identity (string) |
| `D` | Unique dam identity (string). This is **dam**, not sire |

Do not mix `asd` and `ASD` in one format string as a recipe. Pick integer
IDs or string IDs. Colliding hashes raise `PyPedalStringIDCollisionError`.

A display name is a separate optional column `n`. See
[Animal identities and names](animal-identities-and-names.md).

## Common optional columns

| Code | Meaning |
|---|---|
| `x` | Sex |
| `n` | Display / call name |
| `y` | Birth **year** |
| `b` | Birth **date** |
| `g` | Input generation label |
| `f` | Inbreeding coefficient already present in the file |
| `Z` | Skip this column (repeatable) |

Griffon-style files often use `asdxb`: integer IDs, sex, and a birth date.

Birth year (`y`) and birth date (`b`) are different columns. Unknown
values become `None`. See
[Birth dates and chronology](birth-dates-and-chronology.md).

## Separators and headers

Default `sepchar` is a space. CSV files need `sepchar=","`. Set
`has_header=True` to skip a header line.

You do not need a `% asd` format line inside the pedigree file. If such
a line is present, the loader ignores it.

## What is not a supported load recipe

- GENES stud-file import
- Mixing integer and string identity codes as a documented workflow
- Treating pedformat `p` (Pattie generation coefficient) as a calculation
  PyPedal will perform — the field can be stored; the calculation is not
  supported in PyPedal 4.0
