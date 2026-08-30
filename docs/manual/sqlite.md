# SQLite

`NewPedigree.savedb()` writes a SQLite table using the Python standard
library `sqlite3` module. ADOdb and other SQL backends are not supported.

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

ped.kw["database_file"] = str(work / "mrode.db")
ped.kw["database_table"] = "animals"
assert ped.savedb(drop=True)
```

`savedb()` stores four columns in ASDx style: `animalName`, `sireName`,
`damName`, `sex`. Missing parents are written as the missing-parent
token. `drop=True` recreates the table.

Lower-level helpers in `pyp_db` (`connect_to_database`,
`create_pedigree_table`, `populate_pedigree_table`) can write a wider
animal table. Unknown recorded birth years are SQL `NULL`.

Reload-from-database (`pedsource='db'`) exists as a specialised path
(`pedformat` is forced to `'ASDx'`). Prefer a text file or pickle if you
only need to round-trip a session.

See [Saving and exporting](saving-and-exporting.md) for text save,
pickle, and GEDCOM.
