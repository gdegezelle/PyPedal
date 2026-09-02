"""Qt-free pedigree open options and load service.

The desktop (and the temporary CustomTkinter app) pass user-facing open
fields through ``PedigreeOpenOptions``. Loading a file into a
``PedigreeSession`` happens here so GUI code does not own the
replace-on-success / retain-on-failure contract.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PyPedal.pyp_newclasses import NewPedigree
from PyPedal.pyp_results import ProgressCallback

from .session import PedigreeSession


def normalize_sepchar(raw: str | None) -> str:
    """Map a separator field to a ``sepchar``.

    An empty field means a space. A comma with leftover spaces from that
    default -- ``', '`` or ``' ,'`` -- must still be a comma, or a CSV
    with no space after the delimiter is read as one column. A tab is
    kept. Only spaces are trimmed, so a tab is not stripped to empty and
    then replaced by a space.
    """
    if raw is None:
        return " "
    text = str(raw)
    if text == "\t":
        return "\t"
    stripped = text.strip(" ")
    if stripped == "":
        return " "
    return stripped


@dataclass(frozen=True, slots=True)
class PedigreeOpenOptions:
    """User-facing pedigree-open fields for the application layer.

    ``messages="quiet"`` and ``pedigree_summary=0`` are the intended
    application load configuration: no console chatter and no metadata
    dump. The temporary CustomTkinter wrapper historically set
    ``messages="quiet"`` without ``pedigree_summary=0``; that CTk dict
    is unchanged. This dataclass is the contract future PySide6 loading
    should use.

    ``separator`` is stored after ``normalize_sepchar`` when constructed
    through :meth:`normalized`. Callers may pass a raw field value and
    then call :meth:`normalized`.
    """

    pedformat: str = "asd"
    separator: str | None = None
    renumber: bool = True
    messages: str = "quiet"
    pedigree_summary: int = 0

    def normalized(self) -> PedigreeOpenOptions:
        """Return a copy with pedformat and separator canonicalised."""
        fmt = (self.pedformat or "").strip() or "asd"
        return PedigreeOpenOptions(
            pedformat=fmt,
            separator=normalize_sepchar(self.separator),
            renumber=bool(self.renumber),
            messages=self.messages,
            pedigree_summary=self.pedigree_summary,
        )

    def to_library_options(self, source: Path) -> dict[str, str | bool | int]:
        """Keyword dict ``NewPedigree`` accepts for a file load."""
        options = self.normalized()
        separator = options.separator
        if separator is None:
            separator = " "
        return {
            "pedfile": str(source),
            "pedname": source.name,
            "pedformat": options.pedformat,
            "sepchar": separator,
            "renumber": options.renumber,
            "messages": options.messages,
            "pedigree_summary": options.pedigree_summary,
        }


def resolve_source_path(source: Path) -> Path:
    """Return an explicit filesystem path for ``source``.

    The application layer does not assume a repository checkout, conda
    prefix, or process working directory beyond what ``Path.resolve``
    uses for a relative path the caller supplied.
    """
    return source.expanduser().resolve()


def load_into_session(
    session: PedigreeSession,
    source: Path,
    options: PedigreeOpenOptions | None = None,
    progress: ProgressCallback | None = None,
) -> NewPedigree:
    """Load ``source`` into ``session`` after a successful library load.

    The previous pedigree, path, options, and cached analysis results
    remain in place if construction or ``load()`` raises. The session is
    replaced only after ``NewPedigree.load`` returns. Unexpected
    exceptions propagate unchanged; PyPedal errors stay typed.
    """
    path = resolve_source_path(source)
    open_options = (options or PedigreeOpenOptions()).normalized()
    pedigree = NewPedigree(open_options.to_library_options(path))
    pedigree.load(progress=progress)
    session.replace_pedigree(pedigree, path, open_options)
    return pedigree
