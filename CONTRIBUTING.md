# Contributing to PyPedal

This is a **developer** checkout guide. Users who only want to analyse a
pedigree should start at [`docs/manual/installation.md`](docs/manual/installation.md).

## Supported Python

Python **3.12, 3.13, and 3.14**. Use an isolated environment. A conda
environment named `pypedal3_env` is the convention in this repository.

```bash
conda create -n pypedal3_env python=3.12
conda activate pypedal3_env
python -m pip install -e ".[gui,graphics,reports,test,docs,dev]"
# Optional Qt desktop during 4.2-B:
# python -m pip install -e ".[desktop-test]"
```

A venv works the same way with `python3 -m venv .venv`.

Do not activate a leftover Python 2 virtualenv.

## Verify the install

```bash
python verify_setup.py
python -c "import PyPedal; print(PyPedal.__version__)"
```

## Tests

pytest collects only `tests/` (`testpaths` in `pyproject.toml`).
See [`tests/README.md`](tests/README.md) for categories (oracles,
regression, product, docs, packaging, integration, archaeology).

```bash
# Default suite (skips slow example/large-pedigree jobs)
python -m pytest tests/ -m "not integration" -q

# Integration: example scripts and large pedigrees
python -m pytest tests/ -m integration -q

# One module
python -m pytest tests/test_core.py -q
```

Example scripts must be run from their own directory:

```bash
cd PyPedal/examples && python new_methods.py
```

The suite must leave the repository byte-identical. Confine generated
output to a temporary directory (`tests/_pedhelpers.py`: `load_corpus`,
`load_example`, `chdir_tmp`). Do not add generated `.dat` / `.log` /
`.ped` side effects to the tree or to `.gitignore` as a workaround.

## Documentation

The public user manual is `docs/manual/`. MkDocs `docs_dir` is that
directory.

```bash
python -m pip install -e ".[docs]"
mkdocs build --strict -d /tmp/pypedal-user-manual
```

## Application and desktop layers

`PyPedal.application` is a Qt-free session and load layer between the
PySide6 desktop and the scientific library. `PyPedal.desktop` is the
PySide6 UI. `pypedal`, `pypedal-gui`, and `python -m PyPedal` launch it.

Permanent dependency direction:

```
PyPedal.desktop  →  PyPedal.application  →  scientific pyp_* modules
```

Scientific modules must not import application or desktop. Application
must not import GUI toolkits (`PySide6`, `tkinter`, `customtkinter`).
Desktop should reach science through application adapters.

Architecture tests under `tests/test_application/` and
`tests/test_desktop/` enforce that direction.

## Lint

Three Ruff rules are blocking:

```bash
ruff check --select F821,F811,F601 .
```

Everything else Ruff reports, and all mypy output on legacy `pyp_*`
modules, is reporting-only.
Do not run a whole-tree format on scientific modules. New
`PyPedal/application` and `PyPedal/desktop` code, and their tests, must
pass `ruff format --check` and `ruff check` under the repository rules.
Targeted mypy for `PyPedal.application` and `PyPedal.desktop` must stay
clean.

## Packaging smoke

```bash
python -m pip install build
python -m build
```

Wheels must not contain `PyPedal/examples`. Source distributions may.
Do not commit `dist/`, `build/`, or `*.egg-info`.

`tools/release/build_pristine_wheel.sh` is optional packaging hygiene.

## What belongs in a pull request

- Behaviour, tests, or public documentation for the library.
- New options belong as `kw.setdefault(...)` in `NewPedigree.__init__`,
  not as ad-hoc parameters halfway down a call chain.

Do not add generated analysis files, verification dossiers, or a second
user manual.

Parser internals used by `NewPedigree.preprocess` are private. Do not
re-export them from `PyPedal.__init__` or document them as a public
parser API.
