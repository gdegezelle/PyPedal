# Configuration and options

PyPedal options are a flat dictionary `ped.kw` on `NewPedigree`. Defaults
are established only in `NewPedigree.__init__`. Pass a dict to
`load_pedigree(options={...})` or an INI file (`options_file=`).

If you construct `NewPedigree` with no options dict, it reads
**`pypedal.ini` in the current working directory**, not a file inside the
installed package.

`PyPedal/PyPedal.ini` in a checkout is a template with the same defaults
as the code.

## Load

| Option | Default | Notes |
|---|---|---|
| `pedfile` | required for file load | Input path; also the stem for logs |
| `pedformat` | `'asd'` | One character per column |
| `sepchar` | space | Column delimiter |
| `has_header` | `False` | Skip line 1 when true |
| `missing_parent` | `0` | Unknown sire/dam token |
| `renumber` | `True` | Sequential 1-based `animalID` |
| `messages` | `'verbose'` | Use `'quiet'` in scripts |
| `pedigree_summary` | `1` | Set `0` to skip the dump |

## Missing chronology

| Option | Default | Notes |
|---|---|---|
| `legacy_missing_byear_token` | `None` | Map that integer year to unknown |
| `legacy_missing_bdate_token` | `None` | Map that date token to unknown |
| `estimate_birth_dates` | `False` | Off unless a vital-rate profile is also given |

1800 and 1900 are real years unless a legacy token says otherwise.

## Analysis flags at load

| Option | Default | Notes |
|---|---|---|
| `set_generations` | `False` | Assigns `igen` after metadata exists |
| `set_ancestors` | `False` | |
| `set_sexes` | `False` | |
| `form_nrm` | `False` | Builds `ped.nrm` |
| `matrix_type` | `'sparse'` | `'sparse'` or `'dense'` |
| `gen_coeff` | `False` | Leave false; `True` is refused. Pattie calculation is not supported |
| `validate` | `True` | Runtime checks; do not turn off |

Inbreeding **method** is an argument to `pyp_nrm.inbreeding()`, not a
`kw` default.

## Output and logging

Constructing a pedigree **reconfigures process-wide logging** and
truncates `{filetag}.log`. Analysis functions such as `inbreeding` take
their own `output=False` so they do not write `.dat` files.

`filetag` is derived from the `pedfile` stem and prefixes generated
names.

## INI files

Section headers are flattened. `0`/`1` become integers, so
`renumber = 0` is false. A tab as `sepchar` cannot be stored reliably in
INI; use a visible delimiter.

SQLAlchemy, pyDAL, and ADOdb are not options in PyPedal 4.0.
SQLite uses stdlib `sqlite3` (`database_file`, `database_table`).
