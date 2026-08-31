# Pedigree input

PyPedal reads a pedigree from a text file into a `NewPedigree`. The usual
entry point is `load_pedigree()`, which constructs the object and calls
`load()`.

## Default load

One animal per line. Lines starting with `#` are comments. Columns are
split on `sepchar` (default: a space). Column meaning comes from
`pedformat` (default: `'asd'`).

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree

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

print(len(ped.pedigree))
print(ped.metadata.num_records)
```

`NewPedigree({...})` followed by `ped.load()` is the same workflow.
`load_pedigree()` is the one-step wrapper used by examples and the
desktop app.

`pedsource` defaults to `'file'`. Other sources exist (`db`, `graph`,
`graphfile`, `null`, `animallist`, `gedcomfile`, `textstream`). This
chapter documents the file path.

## Options that matter at load

| Option | Default | Role |
|---|---|---|
| `pedfile` | required for file load | Input path |
| `pedformat` | `'asd'` | One character per column |
| `sepchar` | `' '` | Column separator |
| `has_header` | `False` | Skip a header line when `True` |
| `missing_parent` | `0` | Token for an unknown sire or dam |
| `renumber` | `True` | Assign 1-based sequential `animalID` |
| `messages` | `'verbose'` | Use `'quiet'` to suppress console chatter |
| `pedigree_summary` | `1` | Set `0` to skip the metadata dump |

Constructing a pedigree opens a logfile named from the pedigree path
(`os.path.splitext(pedfile)[0]`, keeping the directory prefix).
Load from a working copy, not from a file you cannot overwrite beside.

Example files live in a **repository checkout** or **sdist** under
`PyPedal/examples/`. They are not installed with the wheel.

Comma-separated files need `sepchar=","`. In an `.ini` options file,
quote the comma: `sepchar = ","`.

## After load

Animals are stored in `ped.pedigree`, a Python list. Metadata such as
counts of sires, dams, and founders lives on `ped.metadata`.

Parents that appear as sire or dam but have no row of their own are
added as records. Duplicate animal IDs are removed. An animal cannot
appear as both a sire and a dam.

Next: [File formats](pedigree-formats.md).
