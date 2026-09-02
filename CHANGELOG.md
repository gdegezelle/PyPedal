# Changelog

All notable changes to PyPedal are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
where practical. Version identifiers follow [Semantic Versioning](https://semver.org/)
and PEP 440 (`4.1.0`).

## [Unreleased]

### Added

- Added a Qt-free application/session layer in preparation for the PySide6
  desktop replacement.
- Initial PySide6 desktop shell and scalable pedigree browser
  (`pypedal-qt`, `python -m PyPedal.desktop`). CustomTkinter remains the
  `pypedal` / `pypedal-gui` application during this development cycle.

## [4.1.0] — 2026-09-01

Minor release. It has not been published to PyPI.

### Added

- Successful `inbreeding`, Lacy effective-founder, and `mating_coi_group`
  results are dict subclasses (`InbreedingResult`,
  `EffectiveFoundersResult`, `MatingCoIGroupResult`). Existing key access
  is unchanged; named properties are convenience accessors.
- Selected long-running operations accept an optional `progress(done,
  total)` callback (Meuwissen–Luo, gene-drop rounds, Boichard ancestor
  selection, pedigree record reading). The default `progress=None` keeps
  existing results. Callback exceptions propagate unchanged. There is no
  cancellation API, and not every function reports progress.
- CI runs the non-integration suite on macOS and Windows with Python 3.12,
  in addition to the Ubuntu 3.12/3.13/3.14 matrix.
- `tests/README.md` classifies the suite (oracles, regression, product,
  docs, packaging, static policy, integration, archaeology). Optional
  `docs` and `packaging` pytest markers are selectors only.

### Changed

- Analysis metrics that historically always wrote `.dat` files now accept
  `output=False` to compute the same result without those analysis files.
  `output=True` remains the default and keeps historical file contents.
- `theoretical_ne_from_metadata` returns the calculated Ne float. It no
  longer returns `True`/`False`. `output=True` still writes the historical
  `.dat` file; `output=False` still writes nothing. The Ne formula is
  unchanged.
- Pedigree file parsing still goes through `NewPedigree.preprocess` /
  `load()`. Column mapping, record iteration, and implicit-parent
  detection now live in private helpers with no public parser API.
- `inbreeding` method names and dense/sparse matrix storage are documented
  as `Literal` aliases. Callers still pass strings; invalid values still
  raise `PyPedalUsageError`.
- The object-model manual page documents factual, derived, and cached
  pedigree data, including that `animal.fa` can hold a loaded coefficient
  and later be overwritten by a computed one.

### Deprecated

- `effective_founders_lacy` still auto-renumbers an unnumbered pedigree
  in 4.1, but emits `DeprecationWarning` when that automatic renumbering
  actually happens. Implicit Lacy renumbering has not been removed.

### Fixed

- Computational failures in relationship, inbreeding, `a_coefficients`,
  NRM construction, `summary_inbreeding`, and `min_max_f` raise a typed
  PyPedal error instead of returning a plausible scientific sentinel
  (`0.0`, `{}`, a 1×1 zero matrix, `"0"`, or `False`).
- Exhausted memory while forming a float64 relationship matrix raises
  instead of silently falling back to a float32 memory-mapped file.
- Direct Python analysis arguments that used to be silently coerced now
  raise `PyPedalUsageError` (`inbreeding` method/`gens`, `foundercoi`,
  `dropped_ancestral_inbreeding` rounds/loci/seed, `fast_a_matrix` method,
  `inbreeding_aguilar` amethod, `find_ancestors_g` gens). Configuration
  and presentation fallbacks (INI values, paper size, `NewPedigree`
  defaults) are unchanged.
- Improved Windows compatibility for pedigree resource cleanup and
  canonical dataset checkout.

## [4.0.1] — 2026-08-31

Patch release. It has not been published to PyPI.

### Fixed

- Derive `filetag` from `os.path.splitext` so `./mrode.ped` and directories
  that contain dots keep a usable prefix beside the pedigree.
- Pedigree logging uses the `PyPedal` package logger with a PyPedal-owned
  FileHandler. Importing PyPedal does not configure the root logger.
- Missing SNP data is logged at DEBUG, not ERROR.
- Inbreeding summaries round for display and never print `-0.000000`.
  Stored coefficients are unchanged.
- The desktop app installs a loaded pedigree on the UI thread, keeps the
  previous pedigree active after a failed load, caps large previews, shows
  a busy state, and uses a single neutral zero-inbreeding sentence.

### Changed

- The `dev` extra includes setuptools and wheel so packaging tests run
  after `pip install -e ".[dev]"`.

## [4.0.0] — 2026-08-30

First release of the Python 3 PyPedal 4 line. It has not been published to
PyPI.

### Changed

- Python 3.12 or newer, modern packaging, tests, and a task-oriented user
  manual.
- Pedigree validation and mutation safety, scientific algorithm validation,
  string identity support, unknown chronology handling, mating analysis,
  reporting, and desktop support.
- Improved Meuwissen-Luo inbreeding and Lacy effective-founder calculations
  on large pedigrees without changing validated coefficients.
- Lacy effective-founder calculations use the validated phantom-founder
  treatment. Historical strict/absorb/half variants are refused instead of
  returning non-probability effective-founder values.
- The 2021 Griffon Bruxellois sample was replaced by the 2026 export
  canonical dataset (`PyPedal/examples/griffonbruxellois_2026_pyp.ped`;
  checkout/sdist only). The file is a 2026 export with recorded births
  through 2025.

### Fixed

- The desktop application now normalizes pedigree separators so a comma
  (with or without surrounding spaces) loads comma-separated files, and
  empty/whitespace separators still mean space. Column-count load errors
  name the separator in use.

## [4.0.0-rc8] — 2026-08-30

Large-pedigree work on the tagged 4.0.0-rc7 tree. Not published to PyPI,
not tagged in this implementation pass, and not a final 4.0.0 release.

### Changed

- Improved the performance of the Meuwissen-Luo inbreeding calculation on
  large pedigrees without changing calculated coefficients.
- Reworked Lacy effective-founder contribution propagation for large
  pedigrees while preserving validated effective-founder results.
- Restricted Lacy effective-founder calculations to the validated
  phantom-founder treatment. Historical strict/absorb/half variants are
  now refused instead of returning non-probability effective-founder
  values.
- The 2021 Griffon Bruxellois sample was replaced by the 2026 export
  canonical dataset (`PyPedal/examples/griffonbruxellois_2026_pyp.ped`;
  checkout/sdist only). The file is a 2026 export with recorded births
  through 2025.

### Fixed

- The desktop application now normalizes pedigree separators so a comma
  (with or without surrounding spaces) loads comma-separated files, and
  empty/whitespace separators still mean space. Column-count load errors
  name the separator in use.

## [4.0.0-rc7] — 2026-08-28

Publication polish of the tagged 4.0.0-rc6 tree. Not published to PyPI, and
not a final 4.0.0 release.

### Fixed

- `tools/release/build_pristine_wheel.sh` now derives the expected wheel
  version from `pyproject.toml` instead of asserting a hard-coded release
  identifier.

### Changed

- mypy is configured to type-check production `PyPedal` modules only
  (not `PyPedal/examples/`). It remains reporting-only for 4.0.

### Removed

- Obsolete `PyPedal/examples/profile/` scratch (mixed tabs and spaces;
  not a supported example).
- Leftover Latin-1 `newfoundland.ped.py`. The example already loads
  `newfoundland.ped`.

## [4.0.0-rc6] — 2026-08-28

Local publication-candidate tree. Not tagged, not published to PyPI, and
not a final 4.0.0 release.

### Added

- Task-oriented user manual under `docs/manual/` (MkDocs).
- `CHANGELOG.md`, `MIGRATION.md`, `CONTRIBUTING.md`, and `THIRD-PARTY.md`.
- Read-only test mating: `mating_coi` and `mating_coi_group` return the
  inbreeding of a prospective offspring (`A_ij / 2`) without inserting a
  child.
- Atomic `delete_animals`: either the named animals are removed or the
  pedigree is left unchanged. Deletion is refused when a surviving animal
  would still name a deleted parent.
- Explicit `merge_animals` for combining duplicate records.
- Unknown recorded chronology stored as `None` on birth year and birth
  date. Optional explicit chronology estimation is a separate step.
- Headless PDF pedigree reports (`reports` extra, ReportLab).
- CustomTkinter desktop application (`gui` extra).
- VanRaden (2008) Method 1 genomic relationship matrix from SNP dosages.
- Conservative LGPL modification notices on original modules that were
  changed for the Python 3 release.

### Changed

- Python 3 port of PyPedal 2.0.4. Requires **Python 3.12 or newer**.
- Packaging is `pyproject.toml` / setuptools. There is no `setup.py` install
  path.
- Pedigree reordering is deterministic: founders first, then a stable
  tie-break that preserves input order among equally eligible animals.
- After default `renumber=True`, calculations use sequential 1-based
  `animalID`. The file identifier is `originalID`.
- Unique string identities (`ASD`) are hashed identities, not call names.
- Validation failures raise domain exceptions rather than returning
  sentinel codes in ordinary analysis paths.
- Gene-drop founder genomes (`effective_founder_genomes`) keep simulation
  state local to the call and return the arithmetic mean of per-replicate
  N<sub>g</sub>. An unknown parent is its own founder-genome source.
- Boichard ancestor metrics complete unknown parental slots with analysis-
  local phantom founders, then run Appendix B on the completed pedigree.
- Dropped ancestral inbreeding follows the GRAIN dummy-founder rule for
  half-founders.
- Graphics use Graphviz hex colours rather than nearest CSS names.
- The repository ships **one** curated Griffon Bruxellois sample
  (`PyPedal/examples/griffonbruxellois_2021_pyp.ped`; checkout/sdist only).

### Removed

- Python 2 runtime.
- GENES stud-file import and export.
- Bundled ADOdb. SQLite uses the standard library `sqlite3` module.
- wxPython GUI. The supported desktop app is CustomTkinter.
- Automatic dense NRM at a 10,000-animal cutoff. Large pedigrees must
  choose `inbreeding(..., method="meu_luo")` explicitly.
- Example scripts from the **wheel**. They remain in a checkout and sdist.

### Fixed

- Reliable string-identity maps (`idmap` / `backmap` / name maps).
- `oldsave` missing-parent handling for unknown sires and dams.
- Fast reordering when the missing-parent token is present.
- `set_age` treats unknown chronology as missing rather than a fake year.
- `set_generation` on load assigns inferred `igen` only when requested,
  and does not copy `igen` into the input `gen` field.

## [2.0.4] — 2010-09-29

Last Python 2 release by John B. Cole. See `MIGRATION.md` for differences
relative to PyPedal 4.
