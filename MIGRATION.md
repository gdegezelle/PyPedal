# Migrating from PyPedal 2.0.4 to PyPedal 4

This guide is for users of PyPedal 2.0.4 who are moving to PyPedal 4.0.
It describes user-facing differences, not an engineering history. Callers
already on PyPedal 4.0.x should skip to
[PyPedal 4.0.x to PyPedal 4.1.0](#pypedal-40x-to-pypedal-410).

PyPedal 4 is a Python 3 reimplementation of Cole’s 2.0.4 library. Some
results match 2.0.4 exactly. Some differ because the 4.0 behaviour is the
intended scientific or data-integrity rule. Where a number can change,
this page says so in ordinary terms.

## Runtime and installation

- **Python 3.12 or newer** is required. Python 2.7 is not supported.
- Install from a local checkout or a wheel you built yourself. This
  candidate **has not been published to PyPI**.
- Use an isolated environment. A conda env named `pypedal3_env` is the
  repository convention:

  ```bash
  conda create -n pypedal3_env python=3.12
  conda activate pypedal3_env
  python -m pip install -e ".[gui,graphics,reports,test]"
  ```

- There is no `setup.py install` path. Metadata lives in `pyproject.toml`.
- Optional extras: `gui`, `graphics`, `graphviz-extra`, `reports`, `test`,
  `docs`. Core analysis needs NumPy, pandas, SciPy, and NetworkX.

## Imports and aliases

Module names are still `PyPedal.pyp_*`. Many functions have snake_case
names with camelCase aliases where 2.x callers depended on them
(`loadPedigree` / `load_pedigree`, `connectToDatabase` /
`connect_to_database`, and similar). Prefer snake_case in new code.

`from PyPedal.pyp_newclasses import load_pedigree` is the usual entry
point.

## Identifiers after renumbering

Default load still **renumbers**. Animals are ordered oldest-first and
given sequential **1-based** `animalID` values. The identifier from the
file is `originalID`.

```
pedigree list index == animalID - 1
```

Inbreeding dictionaries, `relationship`, `mating_coi`, and parent slots
use **current `animalID`**. The number in your studbook file is
`originalID`. Use `ped.idmap` (`originalID → animalID`) and `ped.backmap`
the other way.

On the six-animal Mrode example those two numbers happen to coincide.
Do not rely on that.

## String identities versus names

`asd` means integer IDs. `ASD` means **unique string identities**, hashed
to integers internally. A display or call name is a separate column
(`n`). Two animals may share a call name; PyPedal will not treat a call
name as an ID.

Do not mix `asd` and `ASD` as a documented recipe.

## Missing parents

The missing-parent token still defaults to integer `0`. That means the
parent is unknown in this file. PyPedal does not invent a biological
parent.

| Parents in the file | Ordinary language |
|---|---|
| Both unknown | Founder |
| Exactly one unknown | Half-founder |
| Both known | Non-founder |

Half-founders are not counted as two-parent founders in metadata.

## Chronology: `None`, not a fake year

Unknown birth year (`by`) and birth date (`bd`) are **`None`**.

`1800` and `1900` in a file are ordinary recorded years unless you set
`legacy_missing_byear_token` to treat a chosen token as unknown. That
legacy token is opt-in.

`animal.age` remains a year-offset from 1800 when a year is known. It is
not current biological age and does not fall back to inferred generation
`igen`.

Old pickle files that stored 1800 or 1900 as “unknown” cannot be
auto-disambiguated. Re-load from text if you need true unknowns.

## `gen` versus `igen`

- `gen` is an **input annotation** (pedformat `g`). Load does not compute it.
- `igen` is **inferred pedigree depth**, assigned only by `set_generation`
  (or `"set_generations": True` at load). Founders are 1.

Do not copy `igen` into `gen`. Boichard metrics and gene dropping that
read a generation label use `gen`.

## Inbreeding

`pyp_nrm.inbreeding(ped, method=...)` is the dispatcher.

- Small pedigrees: `method="tabular"` (the default).
- Large pedigrees: `method="meu_luo"`. There is **no automatic switch**
  at 10,000 animals. A dense numerator relationship matrix for ~98,000
  animals needs on the order of 80 GB; do not form one.

Mrode animal 5 is still *F* = 0.125 under tabular.

Return shape is `{'metadata': ..., 'fx': {animalID: F}, ...}`. Coefficients
are keyed by current `animalID`. Pass `output=False` unless you want
side-effect files in the working directory.

## Test mating

`mating_coi(sire, dam, ped)` and `mating_coi_group(pairs, ped)` compute
the inbreeding of a **prospective** offspring without adding that animal.
For distinct parents, *F* = *A<sub>ij</sub>* / 2. The pedigree is not
mutated. IDs are current `animalID` values. `gens` may be `0` or `-1`
only.

## Deletion and merge

`delete_animals` is **atomic**. If any requested ID is invalid, or if a
surviving child would still name a deleted parent, the call raises and
the pedigree is unchanged.

`merge_animals` is an explicit API for combining records. It is not an
automatic duplicate scavenger.

## Reports, GUI, database

- PDF pedigree reports need the `reports` extra (ReportLab). They are
  library functions, not GUI menus.
- The desktop app is CustomTkinter (`python -m PyPedal` or `pypedal`).
  wxPython is gone.
- SQLite uses stdlib `sqlite3`. SQLAlchemy, pyDAL, and ADOdb are not
  PyPedal 4 options.

## Removed components

- **GENES** stud-file import/export.
- **ADOdb** (not distributed; not a project license).
- **wxPython** drawing windows.
- **Pattie (1965) generation-coefficient calculation**. Pedformat `p` can
  *store* a supplied value; `gen_coeff=True` is refused.

## Genomics

Supported domain: VanRaden (2008) **Method 1** GRM from SNP dosages
attached to a pedigree. `inbreeding(..., method="vanraden")` is the
**1992 pedigree** algorithm, not a GRM. AGIL pedigree loaders are not
part of this product.

## Founder and ancestor metrics

Lacy *f<sub>e</sub>* on `new_lacy.ped` is still about **2.91**. PyPedal 4
uses Lacy’s phantom-founder treatment of unknown parents (the default
`mode="phantom"`). Historical `half=True` / `half=False` and the
provisional `mode="absorb"` / `mode="strict"` variants are **not
supported**: they do not form valid founder probability partitions when
half-founders are present, and they raise `PyPedalUsageError`. Use the
default.

Boichard *f<sub>a</sub>* needs a named **reference population**
(`reference=` of current IDs, or the input `gen` column). It does not
fall back to `igen`.

Unknown parental slots are represented with analysis-local phantom
founders for Boichard Appendix B. That is a calculation device, not a
new animal in your file.

Gene-drop N<sub>g</sub> (`effective_founder_genomes`) is the mean of
replicate founder-gene frequencies (Boichard eq. 2). It is not “how many
founders are in the file.” Sex-chromosome gene drop is not supported.

Some 2.0.4 numerical quirks (duplicate materialisation of a string
missing-parent placeholder, non-deterministic reorder ties) are **not**
reproduced. If you compared PyPedal 2 and 4 on the same file and a count
or order differs, check IDs, half-founders, and input order before
assuming a scientific regression.

## Large pedigrees

The repository ships one Griffon sample,
`PyPedal/examples/griffonbruxellois_2026_pyp.ped` (checkout/sdist only;
not in the wheel). It is a 2026 export with recorded births through 2025.
Load with `pedformat="asdxb"` and `sepchar=","`. Observed load: 98,001
records, 6,689 founders, 915 half-founders, 3,997 unknown chronology
dates, `igen` 1…70. Use `meu_luo` for inbreeding. See the
[large pedigrees](docs/manual/large-pedigrees.md) chapter.

## PyPedal 4.0.x to PyPedal 4.1.0

This section is for callers already on PyPedal 4.0.x. Scientific values
on successful inputs are unchanged. No migration step is required if you
only pass valid arguments and keep using existing result keys.

### Dict-compatible result objects

Successful `inbreeding`, Lacy effective-founder, and `mating_coi_group`
results remain mappings. Existing key access still works:

```python
result["fx"]
```

The new types (`InbreedingResult`, `EffectiveFoundersResult`,
`MatingCoIGroupResult`) are `dict` subclasses. Named properties are
convenience accessors; they do not replace the mapping.

### Optional analysis files

Additional analysis writers accept `output=False`. The default remains
`output=True`, so historical `.dat` writing is unchanged unless you opt
out.

### Stricter analysis arguments

Invalid direct analysis arguments that 4.0.x coerced may now raise
`PyPedalUsageError`. Configuration and presentation fallbacks (INI
values, paper size, `NewPedigree` defaults) still warn and coerce.

### Computation and resource failures

Scientific computation or resource failures that previously returned
plausible sentinels (`0.0`, `{}`, a zero matrix, `False`) now raise
typed PyPedal errors. A genuine unrelated pair still has relationship
`0.0`. An empty `a_coefficients` mapping still means nobody has *F* > 0.

### Theoretical Ne

`theoretical_ne_from_metadata` returns the calculated Ne float, not
`True`/`False`. Compare against the number, not `is True`.

### Implicit Lacy renumbering

`effective_founders_lacy` still auto-renumbers an unnumbered pedigree in
4.1. When that automatic renumbering actually happens it emits
`DeprecationWarning`. Renumber explicitly to avoid the warning.
Automatic Lacy renumbering has not been removed.

### Optional progress callbacks

Selected long-running operations accept an optional `progress=`
callback. The default is `None` and preserves 4.0.x results. Not every
function reports progress. There is no cancellation API.

## Where to read next

- [User manual](docs/manual/index.md)
- [Changelog](CHANGELOG.md)
- [Installation](docs/manual/installation.md)
