"""Qt-free animal lookup for desktop analyses.

Breeders locate animals by display/call name, original ID, or current
``animalID``. Scientific relationship and mating APIs stay ID-based.
This index does **not** use ``pedobj.namemap``: that map is the unique
string-identity table for ``ASD`` formats, not a call-name directory.

Names are not identities. Duplicate names return every matching animal
and never silently prefer the first hit.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from PyPedal.pyp_newclasses import NewPedigree

DEFAULT_RESULT_LIMIT = 50
_MISSING_NAME_FOLDED = "unknown_name"
_OMITTED_SEX = frozenset({"", "u", "unknown", "none"})
_SEX_SYMBOLS = {"f": "♀", "female": "♀", "m": "♂", "male": "♂"}

RANK_CURRENT_ID = 0
RANK_ORIGINAL_ID = 1
RANK_EXACT_NAME = 2
RANK_PREFIX_NAME = 3
RANK_SUBSTRING_NAME = 4


def normalize_query(text: str) -> str:
    """Case-fold and collapse whitespace without changing stored data."""
    return " ".join(str(text).split()).casefold()


def _as_int(value: object) -> int | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if text.isdigit():
        return int(text)
    return None


def _birth_year(value: object) -> int | None:
    if value is None or value == "" or value == 0 or value == -999:
        return None
    year = _as_int(value)
    if year is None or year <= 0:
        return None
    return year


def _sex_label(sex: object) -> str | None:
    if sex is None:
        return None
    token = str(sex).strip()
    if not token:
        return None
    folded = token.casefold()
    if folded in _OMITTED_SEX:
        return None
    return _SEX_SYMBOLS.get(folded, token)


def _display_name(
    raw: object,
    *,
    original_id: int | str,
    animal_id: int,
    missing_name: str,
) -> str:
    if raw is None:
        return ""
    name = str(raw).strip()
    if not name:
        return ""
    folded = name.casefold()
    if folded == _MISSING_NAME_FOLDED or folded == missing_name.strip().casefold():
        return ""
    if name == str(original_id) or name == str(animal_id):
        return ""
    return name


def format_animal_label(
    *,
    name: str,
    original_id: int | str,
    sex: object,
    birth_year: int | None,
    animal_id: int,
) -> str:
    """Human-facing row text. Current ID is never the primary identity."""
    parts: list[str] = []
    display = name.strip()
    if display:
        parts.append(display)
    parts.append(str(original_id))
    sex_label = _sex_label(sex)
    if sex_label is not None:
        parts.append(sex_label)
    if birth_year is not None:
        parts.append(str(birth_year))
    parts.append(f"ID {animal_id}")
    return " — ".join(parts)


@dataclass(frozen=True, slots=True)
class AnimalLookupHit:
    """One animal in search results. Lightweight; not a ``NewAnimal`` copy."""

    animal_id: int
    original_id: int | str
    name: str
    name_normalized: str
    sex: str
    birth_year: int | None
    label: str


@dataclass(frozen=True, slots=True)
class AnimalLookupResult:
    """Ranked, bounded search hits plus whether more matches exist."""

    hits: tuple[AnimalLookupHit, ...]
    truncated: bool
    total: int


class AnimalLookupIndex:
    """O(n) search index over a loaded pedigree.

    Stores integer IDs, a display name, a normalized name, sex, and birth
    year per animal. It does not clone ``NewAnimal`` objects.
    """

    def __init__(self, hits: Sequence[AnimalLookupHit]) -> None:
        ordered = tuple(hits)
        self._hits = ordered
        by_animal: dict[int, AnimalLookupHit] = {}
        by_original: dict[int, list[AnimalLookupHit]] = {}
        by_original_text: dict[str, list[AnimalLookupHit]] = {}
        named: list[AnimalLookupHit] = []
        for hit in ordered:
            by_animal[hit.animal_id] = hit
            original_int = _as_int(hit.original_id)
            if original_int is not None:
                by_original.setdefault(original_int, []).append(hit)
            by_original_text.setdefault(str(hit.original_id), []).append(hit)
            if hit.name_normalized:
                named.append(hit)
        self._by_animal_id = by_animal
        self._by_original_id = by_original
        self._by_original_text = by_original_text
        self._named = tuple(named)

    def __len__(self) -> int:
        return len(self._hits)

    @property
    def named_count(self) -> int:
        return len(self._named)

    def hit_for_animal_id(self, animal_id: int) -> AnimalLookupHit | None:
        return self._by_animal_id.get(animal_id)

    @classmethod
    def from_pedigree(cls, pedigree: NewPedigree) -> AnimalLookupIndex:
        missing_name = "Unknown_Name"
        kw = getattr(pedigree, "kw", None)
        if isinstance(kw, Mapping):
            raw_missing = kw.get("missing_name", missing_name)
            if raw_missing is not None:
                missing_name = str(raw_missing)
        animals = getattr(pedigree, "pedigree", ())
        return cls(_hits_from_animals(animals, missing_name=missing_name))

    def search(self, query: str, *, limit: int = DEFAULT_RESULT_LIMIT) -> AnimalLookupResult:
        """Return deterministically ranked matches for ``query``.

        Ranking: exact current ID, exact original ID, exact case-insensitive
        name, name prefix, name substring. Ties break by original ID then
        current ``animalID``.
        """
        if limit <= 0:
            raise ValueError("limit must be a positive integer")
        needle = normalize_query(query)
        if not needle:
            return AnimalLookupResult(hits=(), truncated=False, total=0)

        best: dict[int, tuple[int, AnimalLookupHit]] = {}

        def consider(hit: AnimalLookupHit, rank: int) -> None:
            previous = best.get(hit.animal_id)
            if previous is None or rank < previous[0]:
                best[hit.animal_id] = (rank, hit)

        as_id = _as_int(needle)
        if as_id is not None:
            current = self._by_animal_id.get(as_id)
            if current is not None:
                consider(current, RANK_CURRENT_ID)
            for hit in self._by_original_id.get(as_id, ()):
                consider(hit, RANK_ORIGINAL_ID)

        for hit in self._by_original_text.get(needle, ()):
            consider(hit, RANK_ORIGINAL_ID)

        for hit in self._named:
            name = hit.name_normalized
            if name == needle:
                consider(hit, RANK_EXACT_NAME)
            elif name.startswith(needle):
                consider(hit, RANK_PREFIX_NAME)
            elif needle in name:
                consider(hit, RANK_SUBSTRING_NAME)

        ranked = sorted(best.values(), key=_rank_sort_key)
        total = len(ranked)
        truncated = total > limit
        hits = tuple(hit for _rank, hit in ranked[:limit])
        return AnimalLookupResult(hits=hits, truncated=truncated, total=total)


def _rank_sort_key(item: tuple[int, AnimalLookupHit]) -> tuple[int, int, int | str, int]:
    rank, hit = item
    original = hit.original_id
    if isinstance(original, int):
        return (rank, 0, original, hit.animal_id)
    original_int = _as_int(original)
    if original_int is not None:
        return (rank, 0, original_int, hit.animal_id)
    return (rank, 1, str(original), hit.animal_id)


def _hits_from_animals(animals: Iterable[object], *, missing_name: str) -> list[AnimalLookupHit]:
    hits: list[AnimalLookupHit] = []
    for animal in animals:
        animal_id = _as_int(getattr(animal, "animalID", None))
        if animal_id is None:
            continue
        original_raw = getattr(animal, "originalID", animal_id)
        original_int = _as_int(original_raw)
        original_id: int | str = original_int if original_int is not None else str(original_raw)
        name = _display_name(
            getattr(animal, "name", ""),
            original_id=original_id,
            animal_id=animal_id,
            missing_name=missing_name,
        )
        sex_raw = getattr(animal, "sex", "")
        sex = "" if sex_raw is None else str(sex_raw).strip()
        birth_year = _birth_year(getattr(animal, "by", None))
        hits.append(
            AnimalLookupHit(
                animal_id=animal_id,
                original_id=original_id,
                name=name,
                name_normalized=normalize_query(name) if name else "",
                sex=sex,
                birth_year=birth_year,
                label=format_animal_label(
                    name=name,
                    original_id=original_id,
                    sex=sex,
                    birth_year=birth_year,
                    animal_id=animal_id,
                ),
            )
        )
    return hits
