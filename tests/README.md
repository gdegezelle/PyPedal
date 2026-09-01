# PyPedal tests

Classification for discoverability. The tree is not reorganized by this
document. Independent oracles stay in `tests/oracles/` and must remain
physically separate from production code.

Run the default suite with `pytest tests/ -m "not integration"`. CI also
runs `pytest tests/ -m integration`. Optional selectors `-m docs` and
`-m packaging` exist only where enough tests share a command-line use;
they are not extra release gates.

## Independent scientific oracles

`tests/oracles/` contains standalone implementations of published methods
(Meuwissen–Luo, Lacy, Boichard, Ballou, Suwanlee, VanRaden, founder
genomes). Production tests compare against these oracles.

Independent oracles must **never** import or reuse production algorithms.
An oracle that called `pyp_nrm.inbreeding_meuwissen_luo` would not be an
independent check. They exist so a defect in PyPedal cannot hide by
agreeing with itself.

## Scientific regression

Production results pinned against oracles, published worked examples, and
dataset regressions (Mrode *F*(5) = 0.125, Lacy Appendix A, Boichard
controls, Griffon Meuwissen–Luo summaries). Examples:
`test_meuwissen_luo.py`, `test_lacy_oracle.py`, `test_boichard_oracle.py`,
`test_griffon_sample_contract.py`.

These are **release gates**.

## Product behavior

Public API behaviour that is not a published scientific number: IDs,
mutation, logging, GUI helpers, configuration, output files, error
types. Examples: `test_core.py`, `test_error_contract.py`,
`test_gui_ux.py`, `test_progress.py`.

These are **release gates**.

## Documentation

User-manual existence, nav, and example contracts. Marked `docs`.
Examples: `test_user_manual.py`, `test_manual_pages.py`,
`test_manual_analyses.py`.

These are **release gates** (they run in the default suite). The `docs`
marker is only a selector.

## Packaging / release

Version identity, license text, wheel helper, `verify_setup.py`.
Marked `packaging`. Examples: `test_release_version.py`,
`test_release_license.py`, `test_pristine_wheel_helper.py`.

These are **release gates**. The `packaging` marker is only a selector.

## Static policy

Lint and typing *policy* tests (the blocking Ruff subset, mypy scope).
They do not type-check the library by themselves. Example:
`test_lint_gates.py`, `test_mypy_scope.py`.

These are **release gates**.

## Integration

Example scripts and large-pedigree jobs. Marked `integration`. Slow;
deselected by `pytest tests/ -m "not integration"` and run explicitly in
Linux CI. Examples: `test_examples_integration.py`, the Griffon load in
`test_griffon_sample_contract.py`.

These are **release gates** on the Ubuntu job. macOS/Windows CI runs the
non-integration suite only.

## Archaeology / legacy compatibility

Pinned historical signatures, 2.0.4 numeric comparisons, and known
divergences. Examples: `test_public_signatures.py`,
`test_legacy_numeric_compatibility.py`,
`test_known_scientific_divergences.py`.

Do **not** expand these casually. They record the past; they are not a
place to add new product behaviour. Do not delete them to make a cleanup
look cleaner.

These remain **release gates** so accidental signature or numeric drift
is visible.
