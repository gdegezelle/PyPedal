# Changelog

All notable changes to PyPedal are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
where practical. Version identifiers follow [Semantic Versioning](https://semver.org/)
and PEP 440 (`4.0.0`).

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
