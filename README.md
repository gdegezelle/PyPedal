# PyPedal 4.2.0

PyPedal is a Python package for **pedigree analysis**. It loads a recorded
animal pedigree, checks common data problems, and computes inbreeding,
additive relationships, founder and ancestor contributions, and related
summaries.

**Original author:** John B. Cole
**Current maintainer:** Geert Degezelle
**License:** GNU LGPL 2.1 or later (`LGPL-2.1-or-later`) — see [`LICENSE`](LICENSE)
**Python:** Python 3.12 or newer (3.12, 3.13, and 3.14)

This repository is a standalone Python 3 port of Cole’s PyPedal 2.0.4.
This is **PyPedal 4.2.0**. It **has not been published** to PyPI.

## What you can compute

- Inbreeding coefficients (*F*) and additive (numerator) relationships
- Read-only test mating (`mating_coi` / `mating_coi_group`)
- Effective founders, effective ancestors, and founder genomes
- Gene dropping
- Optional Graphviz / matplotlib drawings (`graphics` extra)
- Optional headless PDF reports (`reports` extra)
- Optional PySide6 desktop app (`gui` extra)

Bounds of those capabilities are in the [user manual](docs/manual/index.md).

## Installation

Install from a local checkout. There is no PyPI package yet.

```bash
conda create -n pypedal3_env python=3.12
conda activate pypedal3_env
python -m pip install -e ".[gui,graphics,reports,test]"
```

A standard venv works as well. See [Installation](docs/manual/installation.md).

## Minimal example

Example files are in `PyPedal/examples/` in a **checkout or sdist**, not
in the wheel. This snippet constructs Mrode’s textbook pedigree inline:

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_nrm

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
print(result["fx"][5])  # 0.125
```

For large files use `method="meu_luo"`. There is no automatic switch at
10,000 animals. The repository ships one curated Griffon sample
(`PyPedal/examples/griffonbruxellois_2026_pyp.ped`; checkout/sdist only;
comma-separated `asdxb`, 2026 export with recorded births through 2025).

Desktop app: `python -m PyPedal` or `pypedal` (requires `PyPedal[gui]`).

```bash
python verify_setup.py
python -m pytest tests/ -m "not integration" -q
```

## Documentation

- [User manual](docs/manual/index.md)
- [Changelog](CHANGELOG.md)
- [Migration from 2.0.4, 4.0.x, and 4.1.x](MIGRATION.md)
- [Contributing](CONTRIBUTING.md)
- [Third-party notices](THIRD-PARTY.md)

```bash
python -m pip install -e ".[docs]"
mkdocs build --strict -d /tmp/pypedal-user-manual
```

Copyright 2001–2025 John B. Cole. Maintained by Geert Degezelle.
Licensed under the GNU LGPL 2.1 or later.
