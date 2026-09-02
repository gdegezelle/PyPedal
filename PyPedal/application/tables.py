"""GUI-independent browse contract for a loaded pedigree.

``PedigreeTableSource`` reads cells from ``ped.pedigree`` in place. It
does not copy animals into dicts, row tuples, or display strings. A
future desktop table model can wrap this source; GUI toolkit types do
not belong here.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from PyPedal.pyp_newclasses import NewPedigree


@dataclass(frozen=True, slots=True)
class TableColumn:
    """One browse column.

    ``key`` is the stable application identifier. ``attribute`` is the
    ``NewAnimal`` field the column reads. Titles are labels, not data.
    """

    key: str
    title: str
    attribute: str


# Restrained initial browse schema. Parent columns use sireID/damID
# because those are the NewAnimal fields; the application keys stay
# ``sire`` / ``dam`` / ``year`` for a future table model.
BROWSE_COLUMNS: tuple[TableColumn, ...] = (
    TableColumn("originalID", "Original ID", "originalID"),
    TableColumn("animalID", "Animal ID", "animalID"),
    TableColumn("sire", "Sire", "sireID"),
    TableColumn("dam", "Dam", "damID"),
    TableColumn("year", "Year", "by"),
    TableColumn("sex", "Sex", "sex"),
    TableColumn("name", "Name", "name"),
    TableColumn("fa", "F", "fa"),
)


class PedigreeTableSource:
    """Random-access view of ``pedigree.pedigree``.

    Row count is ``len(ped.pedigree)``. Cell access is ``O(1)`` list
    indexing, so a 98k-row pedigree is not duplicated for browsing.
    Sorting and filtering stay out of this type; presentation sorting
    belongs in the desktop layer later.
    """

    def __init__(self, pedigree: NewPedigree) -> None:
        self._pedigree = pedigree

    def row_count(self) -> int:
        return len(self._pedigree.pedigree)

    def column_count(self) -> int:
        return len(BROWSE_COLUMNS)

    def columns(self) -> Sequence[TableColumn]:
        return BROWSE_COLUMNS

    def column(self, index: int) -> TableColumn:
        if index < 0 or index >= len(BROWSE_COLUMNS):
            raise IndexError(f"column {index} is out of range")
        return BROWSE_COLUMNS[index]

    def value(self, row: int, column: int) -> object:
        """Return the raw domain value at ``(row, column)``.

        Missing parents, unknown years, and default sex/inbreeding tokens
        are whatever ``NewAnimal`` stored. This method does not invent
        display strings such as ``"unknown"``.
        """
        if row < 0 or row >= self.row_count():
            raise IndexError(f"row {row} is out of range")
        if column < 0 or column >= self.column_count():
            raise IndexError(f"column {column} is out of range")
        animal = self._pedigree.pedigree[row]
        return getattr(animal, BROWSE_COLUMNS[column].attribute)
