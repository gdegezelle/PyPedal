# Historical example scripts

These programs are **not** part of the active PyPedal 4.0 example suite.
They remain in the repository as evidence of APIs that 4.0 does not
support.

They must not be run as current product examples.

| Script | Why it is historical |
|---|---|
| `new_db_sqa.py` | SQLAlchemy / `getCursorSQA` / pyDAL-style database access. Current 4.0 SQLite support is stdlib `sqlite3` (`new_db.py`). |

The user manual is in `docs/manual/`. Active examples are the
top-level `*.py` files in this directory that the integration harness
collects.
