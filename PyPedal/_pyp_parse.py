"""Private pedigree parse helpers.

This module is not a public API. Callers should continue to use
``NewPedigree.preprocess`` / ``load_pedigree``. It exists so pedformat
mapping, record iteration, and implicit-parent detection can be tested
without constructing a full ``NewPedigree``.
"""
from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from typing import Any, Literal

PEDFORMAT_ABSENT = -999

_OPTIONAL_FIELDS = (
    ("generation", "g", "generation"),
    ("gencoeff", "p", "generation coefficient"),
    ("sex", "x", "sex"),
    ("birthyear", "y", "birth date (YYYY)"),
    ("inbreeding", "f", "coeffcient of inbreeding"),
    ("genomic_inbreeding", "G", "coeffcient of genomic inbreeding"),
    ("homozygosity", "Y", None),
    ("breed", "r", "breed"),
    ("name", "n", "name"),
    ("birthdate", "b", "birth date (MMDDYYYY)"),
    ("alive", "l", "alive/dead"),
    ("age", "e", "age"),
    ("userfield", "u", "user-defined field"),
)


def canonicalize_pedformat(
    pedformat: str,
    pedformat_codes: Sequence[str],
) -> tuple[list[str], list[tuple[str, str]]]:
    """Return the character list used for column mapping, plus log events.

    ``Z`` becomes ``'.'``. Unknown codes become ``'.'``. Codes are
    case-sensitive; this function does not rename or case-fold them.
    Events are ``('Z', char)`` or ``('invalid', char)`` in input order.
    """
    codes = set(pedformat_codes)
    canonical: list[str] = []
    events: list[tuple[str, str]] = []
    for char in pedformat:
        if char in codes and char != "Z":
            canonical.append(char)
        elif char in codes and char == "Z":
            canonical.append(".")
            events.append(("Z", char))
        else:
            canonical.append(".")
            events.append(("invalid", char))
    return canonical, events


def build_pedformat_locations(
    canonical: Sequence[str],
    *,
    alleles_sepchar: str,
    sepchar: str,
    pedformat: str | None = None,
) -> tuple[dict[str, int], int, bool, list[str]]:
    """Map a canonical pedformat onto column indices.

    Built once. The historical loop rebuilt the same dict for every
    format character; the resulting locations for a valid string are
    identical.

    The ``h`` then ``H`` lookup is preserved: a missing ``H`` overwrites
    the ``h`` column even when ``h`` was present.

    Returns ``(locations, critical_count, alleles_sepchar_collision,
    missing_optional_debug_messages)``.
    """
    locations: dict[str, int] = {}
    critical_count = 0
    debug_messages: list[str] = []
    pedformat_string = pedformat if pedformat is not None else "".join(canonical)

    try:
        locations["animal"] = canonical.index("a")
    except ValueError:
        try:
            locations["animal"] = canonical.index("A")
        except ValueError:
            critical_count += 1

    try:
        locations["sire"] = canonical.index("s")
    except ValueError:
        try:
            locations["sire"] = canonical.index("S")
        except ValueError:
            critical_count += 1

    try:
        locations["dam"] = canonical.index("d")
    except ValueError:
        try:
            locations["dam"] = canonical.index("D")
        except ValueError:
            critical_count += 1

    for key, code, label in _OPTIONAL_FIELDS:
        try:
            locations[key] = canonical.index(code)
        except ValueError:
            locations[key] = PEDFORMAT_ABSENT
            if label is not None:
                if key == "sex":
                    debug_messages.append(
                        f"[DEBUG]: No {label} code was specified in the pedigree "
                        f"format string {pedformat_string}. This program  will continue."
                    )
                elif key == "birthyear":
                    debug_messages.append(
                        f"[DEBUG]: No {label} code was specified in the pedigree "
                        f"format string {pedformat_string}.  This program will continue."
                    )
                elif key in {"gencoeff", "userfield"}:
                    debug_messages.append(
                        f"[DEBUG]: No {label} was specified in the pedigree format "
                        f"string {pedformat_string}. This program will continue."
                    )
                else:
                    debug_messages.append(
                        f"[DEBUG]: No {label} code was specified in the pedigree "
                        f"format string {pedformat_string}. This program will continue."
                    )

    alleles_collision = False
    try:
        locations["alleles"] = canonical.index("L")
        if alleles_sepchar == sepchar:
            alleles_collision = True
            locations["alleles"] = PEDFORMAT_ABSENT
    except ValueError:
        locations["alleles"] = PEDFORMAT_ABSENT
        debug_messages.append(
            f"[DEBUG]: No alleles code was specified in the pedigree format "
            f"string {pedformat_string}. This program will continue."
        )

    # h then H: a missing H overwrites a found h. Do not "fix" this.
    try:
        locations["herd"] = canonical.index("h")
    except ValueError:
        locations["herd"] = PEDFORMAT_ABSENT
        debug_messages.append(
            f"[DEBUG]: No herd code was specified in the pedigree format "
            f"string {pedformat_string}. This program will continue."
        )
    try:
        locations["herd"] = canonical.index("H")
    except ValueError:
        locations["herd"] = PEDFORMAT_ABSENT
        debug_messages.append(
            f"[DEBUG]: No herd code was specified in the pedigree format "
            f"string {pedformat_string}. This program will continue."
        )

    return locations, critical_count, alleles_collision, debug_messages


def implicit_parent_locations(locations: Mapping[str, int]) -> dict[str, int]:
    """Column map used to construct a parent that has no input record."""
    null_locations = {key: PEDFORMAT_ABSENT for key in locations}
    null_locations["animal"] = 0
    null_locations["sire"] = 1
    null_locations["dam"] = 2
    return null_locations


def iter_implicit_parent_tokens(
    sires: Mapping[Any, Any],
    dams: Mapping[Any, Any],
    idmap: Mapping[Any, Any],
    pedformat: str,
    missing_parent: Any,
    missing_name: Any,
) -> Iterator[tuple[Literal["sire", "dam"], Any]]:
    """Yield ``(role, source_token)`` for parents that need their own record.

    Sires are considered first, then dams. A token already materialized as a
    sire is not created again as a dam. The ``idmap`` probe uses the raw
    source token, which for ``A``/``S``/``D`` pedigrees is the string
    identity rather than the hashed integer.
    """
    materialized: set[str] = set()
    for token in sires:
        if str(token) in materialized:
            continue
        try:
            idmap[token]
        except KeyError:
            if (
                ("S" in pedformat and str(token) != str(missing_name))
                or ("s" in pedformat and str(token) != str(missing_parent))
            ):
                materialized.add(str(token))
                yield "sire", token
    for token in dams:
        if str(token) in materialized:
            continue
        try:
            idmap[token]
        except KeyError:
            if (
                ("D" in pedformat and str(token) != str(missing_name))
                or ("d" in pedformat and str(token) != str(missing_parent))
            ):
                materialized.add(str(token))
                yield "dam", token


class PedigreeRecordSource:
    """Normalize file, text-stream, and database records into raw lines.

    Database records are joined with a comma, matching the historical
    loader. The caller still splits on ``kw['sepchar']``.
    """

    def __init__(
        self,
        pedfile: str,
        textstream: str = "",
        dbstream: Any = "",
    ) -> None:
        self.pedfile = pedfile
        self._file = None
        self._lines: list[str] | None = None
        self._db = None
        self._db_index = 0
        self._known_total: int | None = None
        if textstream == "" and dbstream == "":
            self.kind = "file"
            self._file = open(pedfile, "r", encoding="utf-8-sig")
        elif dbstream == "":
            self.kind = "text"
            self._lines = textstream.split("\n")[:-1]
            self._known_total = len(self._lines)
        else:
            self.kind = "db"
            self._db = dbstream
            try:
                self._known_total = len(dbstream)
            except TypeError:
                self._known_total = None

    @property
    def known_total(self) -> int | None:
        """Cheap record count when already in memory; ``None`` for files.

        File sources do not scan twice merely to report a total.
        """
        return self._known_total

    def readline(self, line_counter: int, logger: Any) -> str | Literal[False]:
        """Return the next raw record, or False when a stream is exhausted."""
        if self.kind == "file":
            assert self._file is not None
            return self._file.readline()
        if self.kind == "text":
            assert self._lines is not None
            try:
                return self._lines.pop(0)
            except IndexError:
                logger.warning(
                    "Reached the end of the textstream after reading %s records.",
                    line_counter,
                )
                return False
        assert self._db is not None
        try:
            dbline = self._db[self._db_index]
            line = ",".join(map(str, dbline))
            self._db_index += 1
            return line
        except IndexError:
            logger.info(
                "Reached the end of the dbstream after reading %s records.",
                line_counter,
            )
            return False

    def close(self) -> None:
        """Close an owned file handle. Text and DB sources own nothing.

        Idempotent. Does not close ``dbstream`` or other caller-owned
        objects; those were never opened by this class.
        """
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self) -> PedigreeRecordSource:
        return self

    def __exit__(self, exc_type, exc, tb) -> Literal[False]:
        self.close()
        return False
