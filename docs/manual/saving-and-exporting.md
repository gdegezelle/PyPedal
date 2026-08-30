# Saving and exporting

Work in a temporary directory. Load and save both write files next to the
pedigree stem.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_io

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
```

## Text save

`NewPedigree.save()` writes columns described by a `pedformat` string
(default `'asd'`), separated by `sepchar` (default space).

- `originalID=False` (default) writes **current** `animalID` values.
- `originalID=True` writes **file** IDs for animal, sire, and dam.

```python
out = work / "saved.ped"
ped.save(str(out), pedformat="asd", originalID=True)
print(out.read_text())
```

If `filename` is omitted, the file is `{filetag}_saved.ped`. Prefer
`save()` over the older `oldsave()` writer.

## Pickle

`pyp_io.pickle_pedigree` / `pyp_io.unpickle_pedigree` serialise the whole
`NewPedigree` object. Pass a **stem without** `.pkl`; both functions
append `.pkl`.

Pickle is Python object serialisation. It is not a stable public pedigree
format. Only unpickle data you trust.

Unknown `by` / `bd` pickle as `None`. Older pickles that stored 1800 or
1900 keep those integers.

## SQLite

`NewPedigree.savedb()` writes a SQLite table using the standard library
`sqlite3` module. ADOdb and other SQL backends are not supported. See
[SQLite](sqlite.md).

## GEDCOM

`savegedcom` writes GEDCOM 5.5. Import via `pedsource='gedcomfile'` is
available for files PyPedal can parse. Human genealogy files vary; treat
GEDCOM as a convenience, not a lossless round-trip for every exporter.

## What is not supported

GENES stud-file import and export are not part of PyPedal 4.0.
