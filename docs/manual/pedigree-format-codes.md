# Pedigree format codes

`pedformat` is one character per input column. Codes are case-sensitive.
A pedigree must include animal, sire, and dam identification: `a` or `A`,
`s` or `S`, and `d` or `D`.

The checkout file `PyPedal/PEDIGREE_FORMAT_CODES.txt` is the same list.

## Identity

| Code | Meaning |
|---|---|
| `a` | Integer animal ID |
| `s` | Integer sire ID |
| `d` | Integer dam ID |
| `A` | Unique string animal identity (hashed; not a call name) |
| `S` | Unique string sire identity |
| `D` | Unique string dam identity (**dam**, not sire) |

Do not mix `asd` and `ASD` as a documented recipe.

## Optional fields

| Code | Meaning | Notes |
|---|---|---|
| `x` | Sex | Missing default `u` |
| `n` | Display / call name | Not a unique identity |
| `y` | Birth year | Unknown → `None` |
| `b` | Birth date | Exact dates become `datetime.date`; a four-digit year sets `by` and leaves `bd` as `None` |
| `g` | Input generation | Not `igen` |
| `f` | Input inbreeding | Sets `f_computed` |
| `r` | Breed | |
| `u` | User field | Default `'Unknown'` |
| `l` | Alive (`1`) / dead (`0`) | |
| `e` | Demographic year-offset (`animal.age`) | Not biological age |
| `L` | Alleles | Split by `alleles_sepchar` (default `/`) |
| `Z` | Skip column | Repeatable |
| `H` | Herd as a string | Use `H`, not `h` |
| `G` | Genomic inbreeding column | Presence is not an estimator |
| `Y` | Genomic homozygosity column | |
| `p` | Pattie coefficient **storage** | Calculation is not supported in PyPedal 4.0 |

`T` (traits) and `P` (SNP string) are listed in the code table but are
not assigned by the file loader. SNP genotypes load from `snpfile`.

See [File formats](pedigree-formats.md) and
[Animal identities and names](animal-identities-and-names.md).
