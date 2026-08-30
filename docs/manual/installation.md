# Installation

PyPedal 4.0 requires **Python 3.12 or newer**. The package metadata names
3.12, 3.13, and 3.14.

This project has **not** been published to PyPI. There is no conda-forge
package. Install from a local checkout or from a wheel you built yourself.

## Isolated environment

Use a virtual environment. A conda environment named `pypedal3_env` is the
convention used in this repository.

```bash
conda create -n pypedal3_env python=3.12
conda activate pypedal3_env
```

A standard venv works as well:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

## Install from a checkout

From the repository root:

```bash
python -m pip install -e ".[gui,graphics,reports,test]"
```

That is an editable install of the library plus optional extras. Core
analysis needs NumPy, pandas, SciPy, and NetworkX; those install without
extras.

| Extra | What it adds |
|---|---|
| `gui` | CustomTkinter desktop app (`python -m PyPedal` or `pypedal`) |
| `graphics` | matplotlib, pydot, Graphviz Python bindings, Pillow |
| `graphviz-extra` | pygraphviz (needed by some drawing functions) |
| `reports` | ReportLab; headless PDF pedigree reports |
| `test` | pytest and pypdf (tests only) |
| `docs` | MkDocs, for building this manual |
| `dev` | pytest, ruff, mypy, pre-commit |
| `all` | feature extras (`gui`, `graphics`, `reports`), not the toolchain |

`docs` and `dev` are contributor extras. They are not part of `all`.

There is no `setup.py` install path in PyPedal 4.0.

## Confirm the install

```python
import PyPedal
print(PyPedal.__version__)
```

On this line that prints `4.0.0`.

From a repository checkout you can also run:

```bash
python verify_setup.py
```

`verify_setup.py` is a smoke check. It is not installed as a console
script.

## Build this manual

```bash
python -m pip install -e ".[docs]"
mkdocs build --strict -d /tmp/pypedal-user-manual
```

Write the HTML **outside** the repository. A default `mkdocs build` would
create a `site/` directory in the working tree.

## Example files versus the wheel

Scripts and `.ped` files under `PyPedal/examples/` ship in a **repository
checkout** and in the **source distribution**. They are **not** installed
with the wheel. Beginner examples in this manual write a small pedigree
in a temporary directory so they work after `pip install` of a wheel.
The large Griffon sample is checkout/sdist only; see
[Large pedigrees](large-pedigrees.md).

## What this chapter does not cover

- Python 2
- `python setup.py install`
- SourceForge or Enthought installers
- Numarray
- `pip install PyPedal` from PyPI (there is no published package yet)
