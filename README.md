# PyPedal

PyPedal is a Python 3 pedigree-analysis toolkit and desktop application,
modernized from PyPedal 2.0.4 originally developed by John B. Cole.

It provides tools for pedigree management, population-genetics analysis,
relationship and mating calculations, and interactive analysis of large
pedigrees.

## What PyPedal can do

### Pedigree analysis

- Calculate individual inbreeding coefficients
- Summarize inbreeding by year
- Calculate pairwise additive relationships
- Evaluate prospective offspring inbreeding for mating pairs
- Estimate the effective number of founders
- Calculate population-level metrics such as theoretical effective population size

### Desktop application

- Open and browse large pedigree files
- Search animals by name, original pedigree ID, or current PyPedal ID
- Distinguish animals with duplicate names using pedigree metadata
- Run Relationship and Mating analyses interactively
- Export analysis results with source pedigree IDs, names, and explicit current IDs
- Work with pedigrees containing approximately 100,000 animals without
  constructing a dense relationship matrix for pairwise analyses

## Installation

Install from a local checkout. There is no PyPI package yet.

```bash
conda create -n pypedal3_env python=3.12
conda activate pypedal3_env
python -m pip install -e ".[gui]"
```

A standard virtual environment works as well.

For development, documentation, testing, and optional dependencies, see [Installation](docs/manual/installation.md).

## Minimal example

Example files are included in a checkout or source distribution, but not in
the wheel.

This example constructs Mrode's textbook pedigree inline:

```
from pathlib import Path
from tempfile import TemporaryDirectory

from PyPedal import pyp_nrm
from PyPedal.pyp_newclasses import load_pedigree

with TemporaryDirectory() as work:
    pedfile = Path(work) / "mrode.ped"

    pedfile.write_text(
        "1 0 0\n"
        "2 0 0\n"
        "3 1 2\n"
        "4 1 0\n"
        "5 4 3\n"
        "6 5 2\n"
    )

    ped = load_pedigree(
        options={
            "pedfile": str(pedfile),
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )

    result = pyp_nrm.inbreeding(
        ped,
        method="tabular",
        output=False,
    )

    print(result["fx"][5])  # 0.125
```

**Large pedigrees:** use `method="meu_luo"`.

PyPedal does not automatically switch methods at a particular pedigree size.

## Example datasets
The source distribution includes two curated Griffon Bruxellois pedigree
datasets with 98,001 animals.

- `griffonbruxellois_2026_pyp.ped`
  - scientific dataset
  - comma-separated
  - pedigree format `asdxb`
- `griffonbruxellois_2026_named_pyp.ped`
  - same genealogy as the scientific dataset
  - includes display names
  - pedigree format `asdxbn`

The named dataset is intended for desktop and animal-selection workflows. 
Names are labels and are not guaranteed to be unique identities.

Both datasets are included in a checkout or source distribution, but not in
the wheel.

## Desktop application

Install the GUI dependencies with:

```bash
python -m pip install -e ".[gui]"
```

Start the desktop application with:

```bash
python -m PyPedal
```

or:

```bash
pypedal
```

The desktop application provides interactive access to pedigree browsing,
inbreeding analysis, effective founder analysis, relationship calculations,
mating analysis, and population-level summaries.

## Analysis exports

Animal-level analysis exports preserve pedigree provenance explicitly.

Depending on the analysis, exported rows may contain:

- `original_id` — the identifier from the source pedigree
- `name` — the stored display name, when available
- `current_id` — the current internal PyPedal ID after renumbering
- raw scientific coefficients
- supplementary percentage values where useful

Names are not used as unique identities.

Scientific CSV exports are locale-independent:

- comma-delimited
- decimal point
- UTF-8
- standard CSV quoting

PyPedal does not switch to locale-dependent decimal commas.

Spreadsheet applications configured for decimal-comma locales may require
 explicit CSV import settings.

## Large pedigrees

PyPedal is designed to work with large pedigrees without requiring a dense
 relationship matrix for every analysis.

Pairwise Relationship and Mating calculations use scalable selected-pair
 methods and do not construct a pedigree-wide dense numerator relationship
 matrix.

For the canonical 98,001-animal Griffon Bruxellois pedigree, the desktop
 application can browse and search the pedigree directly while keeping the
 scientific calculations in the Python reference implementation.

Some full-pedigree calculations, such as Meuwissen-Luo inbreeding, remain
 computationally more expensive and may take significantly longer than
 pairwise analyses.

## Verification

Verify the installation with:

```bash
python verify_setup.py
```

Run the default test suite with:

```bash
python -m pytest tests/ -m "not integration" -q
```

Run integration tests with:

```bash
python -m pytest tests/ -m integration -q
```


## Documentation

- [User manual](docs/manual/index.md)
- [Installation](docs/manual/installation.md)
- [Changelog](CHANGELOG.md)
- [Migration guide](MIGRATION.md)
- [Contributing](CONTRIBUTING.md)
- [Third-party notices](THIRD-PARTY.md)

Build the user manual with:

```bash
python -m pip install -e ".[docs]"
mkdocs build --strict -d /tmp/pypedal-user-manual
```

## Project history

PyPedal was originally developed by John B. Cole.

The current codebase modernizes PyPedal 2.0.4 for Python 3 and continues its
 development as a pedigree-analysis toolkit and desktop application.

Copyright 2001–2025 John B. Cole.

Current maintainer: Geert Degezelle.

## License

Licensed under the GNU Lesser General Public License, version 2.1 or later.