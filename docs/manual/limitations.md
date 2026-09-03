# Limitations

This page is the **PyPedal 4.0 public contract** for genuine bounds of
supported features and for historical domains that are not on this product.
It is not a backlog of unfinished advertised APIs.

A capability listed as **supported** in this manual works inside the bounds
below. Bounds are not a substitute for a missing feature.

| Status | Meaning |
|---|---|
| **Supported, bounded** | Works inside a documented contract. A bound is a limited domain, not a stub. |
| **Outside 4.0 / removed** | Historical functionality is not a current workflow. |
| **External** | Optional shell-out; not a fully native implementation. |
| **Ambiguous signal** | A return value does not distinguish success from failure. |

“Unsupported” on this page means a domain that is not a current 4.0
workflow, not an advertised feature that still needs writing.

## Supported-feature bounds

### Test mating (`mating_coi` / `mating_coi_group`)

Prospective offspring inbreeding is a **supported** calculation.
See [Test mating](mating.md).

- IDs are current / renumbered `animalID` values. `originalID` is not
  translated. The pedigree must already be renumbered.
- `gens` may be `0` or `-1` (both: full available pedigree, read-only).
  Other values raise `PyPedalUsageError`.
- The functions evaluate the pair or list of pairs you supply. They are
  not a mate-selection engine and do not form a Cartesian product.
- They do not automatically form a dense NRM.
- Each pair uses the same exact selected relationship calculation as
  `relationship()`; they do not build an ancestor sub-NRM for one pair.

### PDF pedigree reports

Headless PDF reports are **supported** with the `reports` extra (ReportLab).
See [PDF reports](pdf-reports.md).

- ReportLab is optional. Calling a PDF API without it raises
  `PyPedalDependencyError`.
- The APIs are headless. The PySide6 desktop can also export those PDFs
  from **File → Export** when ReportLab is installed.
- Three-generation subjects are current `animalID` values. Layout depth is
  the historical 15-slot pedigree (proband + three ancestral generations).

### `delete_animals` / `merge_animals`

Atomic mutation, not a general database.
See [IDs, missing parents, and half-founders](ids-and-missing-parents.md).

- Requests use `idmap` keys (`originalID` after a normal load).
- Unrenumbered state raises `PyPedalUsageError`.
- No cascade deletion and no orphan rewrite to `missing_parent`.
- `merge_animals(keep, drop)` is the only offspring-redirect path.
- Structural mutation sets `ped.nrm = None`.

### Unknown birth chronology

Unknown recorded birth year and date are **`None`**.
See [Generation and demography](birth-dates-and-chronology.md).

- 1800 and 1900 are ordinary recorded years unless
  `legacy_missing_byear_token` / `legacy_missing_bdate_token` is set.
- Old pickles that stored 1800 or 1900 cannot be auto-disambiguated.
- Optional birth-date range estimation is off by default and has **no**
  built-in species vital-rate preset.
- `animal.age` is a year-offset from 1800, not biological age, and never
  falls back to `igen`.

### Tabular inbreeding on large pedigrees

Default `pyp_nrm.inbreeding(method='tabular')` builds a full NRM.
There is **no automatic switch** at 10,000 animals. For large files use
`method="meu_luo"`. See [Inbreeding](inbreeding.md).

### Progress reporting

Optional `progress(done, total)` callbacks exist on selected long-running
operations (Meuwissen–Luo, gene-drop rounds, Boichard ancestor
selection, pedigree record reading). They are not present on every
function. There is no cancellation token; a callback that raises
propagates that exception. See [Inbreeding](inbreeding.md).

### Effective founder genomes

`chrometype='autosome'` is the supported domain. Other values, including
`'sex'`, raise `PyPedalUsageError`. There is no `reference=` keyword; the
population under study is the highest numeric `gen` label (or the whole
pedigree). See [Gene dropping](gene-dropping.md).

### Genomic Method 1 only

VanRaden (2008) Method 1 is the supported GRM. Methods 2 and 3 are outside
the 4.0 domain (`PyPedalUsageError`). `scale_m=False` is refused.
See [Genomic methods](genomics.md).

### `renumber=False`

Load may succeed. Many analyses assume renumbered 1-based IDs.
See [Reordering and renumbering](ids-and-missing-parents.md).

### Pedformat `p` / Pattie generation coefficient

`p` **stores** a supplied generation coefficient. PyPedal does **not**
compute Pattie (1965) coefficients. `kw['gen_coeff']=True` raises
`PyPedalUsageError`. See [Pedigree format codes](pedigree-format-codes.md).

### Pickle

Supported Python serialisation for a session, not a stable archive format.
See [Saving, pickle, and SQLite](saving-and-exporting.md).

### Deprecated NRM helpers

`pyp_nrm.a_matrix` is marked deprecated. Prefer `fast_a_matrix` or
`inbreeding`. Do not form a dense NRM on a large pedigree merely to obtain
*F*.

### Matplotlib plot-to-file helpers

`plot_line_xy` and related wrappers are not a reliable PNG export path.
Supported pedigree drawing is `draw_pedigree`. See [Graphics](graphics.md).

### Logging

Importing PyPedal does not configure the root logger. Constructing a
`NewPedigree` attaches a PyPedal-owned logfile handler to the `PyPedal`
package logger and truncates `<filetag>.log`. Host-application handlers
are not removed. `messages="quiet"` does not override logging the host
has already configured.

### `a_coefficients` / `fast_a_coefficients`

A returned `{}` means nobody in the pedigree has *F* > 0. Matrix
construction failure raises `PyPedalError`. It does not return `{}`.

### Direct API versus configuration

Invalid arguments to analysis functions raise `PyPedalUsageError`.
Configuration and presentation values (INI coercion, paper size,
`NewPedigree` defaults, unknown pedformat codes on save) may still warn
and fall back. That distinction is intentional.

### `pyp_nrm.inbreeding(method="aguilar")`

Shells out to external **`inbupgf90`**. Missing binary →
`FileNotFoundError`. Not a fully native, independently validated
implementation.

## Historical domains outside PyPedal 4.0

These are **not** current 4.0 workflows. They are listed once here rather
than on every user page.

### GENES import/export

GENES import/export is **not supported** in PyPedal 4.0. There is no
current import or export path. GEDCOM and text save are not a GENES
replacement. See [Saving, pickle, and SQLite](saving-and-exporting.md).

### Historical wxPython GUI

The 4.2 application is PySide6.
See [Desktop application](desktop-application.md).

### `pyp_network.dyad_census`

Outside the 4.0 domain (`PyPedalUsageError`). NetworkX 3 has no equivalent.
Active examples do not call it.

### AGIL pedigree reader

`read_agil_pedigree_file` is not a 4.0 function. Related AGIL *file*
helpers exist; they are not a pedigree loader.
See [Genomic methods](genomics.md).

## Related pages

- [Inbreeding](inbreeding.md)
- [PDF reports](pdf-reports.md)
- [Genomic methods](genomics.md)
- [Graphics](graphics.md)
- [Desktop application](desktop-application.md)
- [Recipes](recipes.md)
- [Options and configuration](configuration.md)
- [References](references.md)
