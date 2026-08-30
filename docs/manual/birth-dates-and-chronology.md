# Birth dates and chronology

PyPedal stores recorded chronology as:

- `by` — birth **year**, an `int` or `None`
- `bd` — birth **date**, a `datetime.date` or `None`

**Unknown is `None`.** That is the PyPedal 4 rule.

## Recorded versus estimated

**Recorded** chronology is what the file said. A missing year or date
becomes `None`. A four-digit year such as 1800 or 1900 is a **real
year** unless you configure an importer token that maps that integer to
unknown (`legacy_missing_byear_token` / `legacy_missing_bdate_token`).

**Estimated** chronology is a separate, optional load-time step
(`estimate_birth_dates`, which also needs a `vital_rate_profile`). It is
off by default. There is no built-in species preset. Do not treat an
estimated date as if it had been recorded.

## How unknown values are written in files

| Input | Meaning |
|---|---|
| `.` | Modern missing chronology token → `None` |
| Historical year `0` | Accepted as unknown → `None` |
| Empty / omitted year or date | Unknown → `None` |
| `1800`, `1900` | Ordinary recorded years |
| Malformed non-empty value | Error — the load refuses |

A four-digit year in the `b` (date) column sets `by` and leaves `bd` as
`None`. PyPedal does not invent 1 January from a year-only value.

## Age is not biological age

`animal.age` is a **legacy demographic year-offset** (`by` minus 1800
when a year is known). It is not current biological age in years. It
does not fall back to the inferred generation `igen` when the year is
unknown.

## Why this matters

Generation-interval calculations skip animals whose recorded year is
`None`. Grouping inbreeding “by year” must skip missing years rather
than treating them as 1800. Old pickle files that stored 1800 or 1900
as a stand-in for “unknown” cannot be auto-corrected; PyPedal cannot
tell a genuine year from a historical placeholder.

See [Generation intervals](generation-intervals.md) for `gen` versus
`igen`.
