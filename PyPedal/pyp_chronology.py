###############################################################################
# NAME: pyp_chronology.py
# VERSION: see PyPedal.__version__
# LICENSE: LGPL
# Written for PyPedal 4.0; this module is not part of the original 2.0.4 sources.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
"""
Recorded and estimated birth chronology for PyPedal 4.0.

Recorded facts
--------------
``bd`` is ``datetime.date`` or ``None``. ``by`` is ``int`` or ``None``.
Unknown is never a numeric or date sentinel (not 1800, 1900, ``01011800``,
0 on output, …).

Estimated chronology
--------------------
``BirthDateEstimate`` lives on ``animal.birth_date_estimate`` and must never
overwrite ``bd`` / ``by``. Default load does not estimate. There is no
built-in species vital-rate preset.
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass
from typing import Any, Iterable, Optional

from .pyp_errors import PyPedalUsageError, PyPedalValidationError

TEXT_MISSING_TOKEN = "."
_SOURCE_VITAL_RATE = "vital-rate-profile"


def is_text_missing_token(token: Any) -> bool:
    """True for absent, blank, or the documented chronology missing token."""
    if token is None:
        return True
    if isinstance(token, datetime.date):
        return False
    text = str(token).strip()
    return text == "" or text == TEXT_MISSING_TOKEN


def recorded_year(animal: Any) -> Optional[int]:
    """Known recorded year: ``bd.year`` if an exact date exists, else ``by``."""
    bd = getattr(animal, "bd", None)
    if isinstance(bd, datetime.date):
        return bd.year
    by = getattr(animal, "by", None)
    if by is None:
        return None
    try:
        return int(by)
    except (TypeError, ValueError):
        return None


def _ident(animal: Any) -> str:
    return str(getattr(animal, "originalID", getattr(animal, "animalID", "?")))


def parse_recorded_year(
    token: Any,
    *,
    legacy_missing_byear_token: Any = None,
    animal_id: Any = None,
) -> Optional[int]:
    """
    Parse a birth-year token.

    Missing (blank, whitespace, ``.``) and historical input year ``0`` become
    ``None``. ``legacy_missing_byear_token``, when an integer, maps only that
    year to ``None``. Any other integer, including 1800 and 1900, is a real
    year. Non-empty malformed tokens raise ``PyPedalValidationError``.
    """
    if is_text_missing_token(token):
        return None
    if isinstance(token, bool):
        raise PyPedalValidationError(
            _malformed_year_message(token, animal_id)
        )
    if isinstance(token, int):
        year = token
    else:
        text = str(token).strip()
        try:
            year = int(text)
        except (TypeError, ValueError) as exc:
            raise PyPedalValidationError(
                _malformed_year_message(token, animal_id)
            ) from exc
    if year == 0:
        return None
    if (
        legacy_missing_byear_token is not None
        and legacy_missing_byear_token != ""
        and year == int(legacy_missing_byear_token)
    ):
        return None
    return year


def parse_recorded_date(
    token: Any,
    *,
    legacy_missing_bdate_token: Any = None,
    legacy_missing_byear_token: Any = None,
    animal_id: Any = None,
) -> tuple[Optional[datetime.date], Optional[int]]:
    """
    Parse a pedformat ``b`` token.

    Returns ``(bd, by)``:

    * exact MMDDYYYY -> ``(date, year)``
    * seven-digit MMDYYYY (unpadded day 1-9) -> ``(date, year)``
    * four-digit year -> ``(None, year)``  (year-only; no invented day)
    * missing token / year 0 -> ``(None, None)``

    ``legacy_missing_bdate_token`` compares the stripped input string and maps
    only that token to unknown. It is not inferred from the year option.
    """
    if isinstance(token, datetime.datetime):
        token = token.date()
    if isinstance(token, datetime.date):
        return token, token.year
    if is_text_missing_token(token):
        return None, None
    text = str(token).strip()
    if (
        legacy_missing_bdate_token is not None
        and str(legacy_missing_bdate_token).strip() != ""
        and text == str(legacy_missing_bdate_token).strip()
    ):
        return None, None
    if text.isdigit() and int(text) == 0:
        return None, None
    if text.isdigit() and len(text) == 4:
        return None, parse_recorded_year(
            text,
            legacy_missing_byear_token=legacy_missing_byear_token,
            animal_id=animal_id,
        )
    if text.isdigit() and len(text) == 7:
        month = int(text[0:2])
        day = int(text[2:3])
        year = int(text[3:7])
        if year == 0:
            return None, None
        try:
            parsed = datetime.date(year, month, day)
        except ValueError as exc:
            raise PyPedalValidationError(
                _malformed_date_message(token, animal_id)
            ) from exc
        year_out = parse_recorded_year(
            year,
            legacy_missing_byear_token=legacy_missing_byear_token,
            animal_id=animal_id,
        )
        if year_out is None:
            return None, None
        return parsed, year_out
    if text.isdigit() and len(text) == 8:
        month = int(text[0:2])
        day = int(text[2:4])
        year = int(text[4:8])
        if year == 0:
            return None, None
        try:
            parsed = datetime.date(year, month, day)
        except ValueError as exc:
            raise PyPedalValidationError(
                _malformed_date_message(token, animal_id)
            ) from exc
        year_out = parse_recorded_year(
            year,
            legacy_missing_byear_token=legacy_missing_byear_token,
            animal_id=animal_id,
        )
        if year_out is None:
            return None, None
        return parsed, year_out
    raise PyPedalValidationError(_malformed_date_message(token, animal_id))


def reconcile_recorded_chronology(
    bd: Optional[datetime.date],
    by: Optional[int],
    *,
    animal_id: Any = None,
) -> tuple[Optional[datetime.date], Optional[int]]:
    """Combine year and date columns. Conflicting exact values refuse."""
    if bd is not None and by is not None and bd.year != by:
        ident = "" if animal_id is None else f" for animal {animal_id}"
        raise PyPedalValidationError(
            f"Recorded birth date {bd.isoformat()} and birth year {by} "
            f"disagree{ident}."
        )
    if bd is not None and by is None:
        return bd, bd.year
    return bd, by


def parse_animal_chronology(
    year_token: Any,
    date_token: Any,
    *,
    has_year_column: bool,
    has_date_column: bool,
    legacy_missing_byear_token: Any = None,
    legacy_missing_bdate_token: Any = None,
    animal_id: Any = None,
) -> tuple[Optional[datetime.date], Optional[int]]:
    """Parse recorded chronology from the year and/or date columns present."""
    bd = None
    by = None
    if has_date_column:
        bd, by_from_date = parse_recorded_date(
            date_token,
            legacy_missing_bdate_token=legacy_missing_bdate_token,
            legacy_missing_byear_token=legacy_missing_byear_token,
            animal_id=animal_id,
        )
        by = by_from_date
    if has_year_column:
        by_from_year = parse_recorded_year(
            year_token,
            legacy_missing_byear_token=legacy_missing_byear_token,
            animal_id=animal_id,
        )
        if has_date_column:
            if (
                by is not None
                and by_from_year is not None
                and by != by_from_year
            ):
                ident = "" if animal_id is None else f" for animal {animal_id}"
                raise PyPedalValidationError(
                    f"Recorded birth date/year {by} and birth year "
                    f"{by_from_year} disagree{ident}."
                )
            bd, by = reconcile_recorded_chronology(
                bd, by_from_year if by_from_year is not None else by,
                animal_id=animal_id,
            )
        else:
            by = by_from_year
    return bd, by


def format_year_token(by: Any) -> str:
    """Canonical text for a ``y`` field. Unknown is ``.``, never 0/1800/1900."""
    if by is None:
        return TEXT_MISSING_TOKEN
    return str(int(by))


def format_date_token(
    bd: Any,
    by: Any = None,
    *,
    allow_year_only: bool = True,
) -> str:
    """
    Canonical text for a ``b`` field.

    Exact dates are MMDDYYYY. Year-only (``bd is None``, ``by`` known) writes
    a four-digit year when ``allow_year_only`` is true. Unknown writes ``.``.
    """
    if isinstance(bd, datetime.datetime):
        bd = bd.date()
    if isinstance(bd, datetime.date):
        return f"{bd.month:02d}{bd.day:02d}{bd.year:04d}"
    if allow_year_only and by is not None:
        return str(int(by))
    return TEXT_MISSING_TOKEN


def format_pedigree_field(code: str, animal: Any, pedformat: str = "") -> Any:
    """Format one save-column value. Chronology codes never emit sentinels."""
    if code == "y":
        return format_year_token(getattr(animal, "by", None))
    if code == "b":
        allow_year_only = "y" not in pedformat
        return format_date_token(
            getattr(animal, "bd", None),
            getattr(animal, "by", None),
            allow_year_only=allow_year_only,
        )
    return None


def padded_identity(animal_id: Any, by: Any = None) -> str:
    """
    Deterministic padded ID. Unknown year uses identity zero-padding, never
    ``None`` / ``1800`` / ``1900`` as a fake year.
    """
    aid = str(animal_id)
    if by is None:
        return aid.zfill(15)
    length = len(aid)
    pad = max(0, 15 - length - len(str(by)))
    return f"{by}{'0' * pad}{aid}"


def light_padded_id(by: Any, animal_id: Any) -> str:
    """LightAnimal pad_id formula; unknown year uses identity padding."""
    if by is None:
        return str(animal_id).zfill(15)
    length = len(str(animal_id))
    pad = 15 - length - 1
    if pad > 0:
        return f"{by}{'0' * pad}{animal_id}{length}"
    return f"{by}{animal_id}{length}"


def _malformed_year_message(token: Any, animal_id: Any) -> str:
    ident = "" if animal_id is None else f" for animal {animal_id}"
    return f"Malformed birth year {token!r}{ident}."


def _malformed_date_message(token: Any, animal_id: Any) -> str:
    ident = "" if animal_id is None else f" for animal {animal_id}"
    return f"Malformed or impossible birth date {token!r}{ident}."


def validate_recorded_chronology(pedobj: Any) -> None:
    """
    Refuse detectable parent/offspring chronology impossibilities.

    Unknown chronology is valid. Same calendar year at year resolution is
    valid. Exact same-day parent and child is not.
    """
    missing = pedobj.kw["missing_parent"]
    by_id = {animal.animalID: animal for animal in pedobj.pedigree}
    for child in pedobj.pedigree:
        for role, parent_id in (("sire", child.sireID), ("dam", child.damID)):
            if parent_id == missing or str(parent_id) == str(missing):
                continue
            parent = by_id.get(parent_id)
            if parent is None:
                continue
            _check_parent_child_order(parent, child, role)


def _check_parent_child_order(parent: Any, child: Any, role: str) -> None:
    parent_bd = parent.bd if isinstance(getattr(parent, "bd", None), datetime.date) else None
    child_bd = child.bd if isinstance(getattr(child, "bd", None), datetime.date) else None
    if parent_bd is not None and child_bd is not None:
        if parent_bd >= child_bd:
            raise PyPedalValidationError(
                f"{role.capitalize()} {_ident(parent)} birth date "
                f"{parent_bd.isoformat()} is not strictly earlier than child "
                f"{_ident(child)} birth date {child_bd.isoformat()}."
            )
        return
    parent_year = recorded_year(parent)
    child_year = recorded_year(child)
    if parent_year is None or child_year is None:
        return
    if parent_year > child_year:
        raise PyPedalValidationError(
            f"{role.capitalize()} {_ident(parent)} birth year {parent_year} "
            f"is later than child {_ident(child)} birth year {child_year}."
        )


def ranges_are_impossible(
    parent_min: datetime.date,
    parent_max: datetime.date,
    child_min: datetime.date,
    child_max: datetime.date,
) -> bool:
    """Inclusive estimated ranges: impossible iff parent_min >= child_max."""
    return parent_min >= child_max


@dataclass
class BirthDateEstimate:
    """Inferred chronology. Never written into ``bd`` / ``by``."""

    earliest: Optional[datetime.date] = None
    latest: Optional[datetime.date] = None
    typical: Optional[datetime.date] = None
    source: Optional[str] = None
    profile: Optional[str] = None

    def is_empty(self) -> bool:
        return (
            self.earliest is None
            and self.latest is None
            and self.typical is None
        )


@dataclass
class VitalRateProfile:
    """
    Caller-supplied vital-rate bounds. No built-in species preset exists.

    Day counts are nonnegative integers. Role-specific range inference
    requires gestation plus that role's min and max conception ages.
    Incomplete roles are skipped rather than filled with hidden defaults.
    """

    name: Optional[str] = None
    gestation_days: Optional[int] = None
    sire_min_age_at_conception_days: Optional[int] = None
    sire_max_age_at_conception_days: Optional[int] = None
    dam_min_age_at_conception_days: Optional[int] = None
    dam_max_age_at_conception_days: Optional[int] = None
    founder_typical_age_at_progeny_days: Optional[int] = None

    def validate(self) -> None:
        _require_nonnegative_int("gestation_days", self.gestation_days)
        _require_nonnegative_int(
            "sire_min_age_at_conception_days",
            self.sire_min_age_at_conception_days,
        )
        _require_nonnegative_int(
            "sire_max_age_at_conception_days",
            self.sire_max_age_at_conception_days,
        )
        _require_nonnegative_int(
            "dam_min_age_at_conception_days",
            self.dam_min_age_at_conception_days,
        )
        _require_nonnegative_int(
            "dam_max_age_at_conception_days",
            self.dam_max_age_at_conception_days,
        )
        _require_nonnegative_int(
            "founder_typical_age_at_progeny_days",
            self.founder_typical_age_at_progeny_days,
        )
        _require_min_le_max(
            "sire",
            self.sire_min_age_at_conception_days,
            self.sire_max_age_at_conception_days,
        )
        _require_min_le_max(
            "dam",
            self.dam_min_age_at_conception_days,
            self.dam_max_age_at_conception_days,
        )

    def sire_range_ready(self) -> bool:
        return (
            self.gestation_days is not None
            and self.sire_min_age_at_conception_days is not None
            and self.sire_max_age_at_conception_days is not None
        )

    def dam_range_ready(self) -> bool:
        return (
            self.gestation_days is not None
            and self.dam_min_age_at_conception_days is not None
            and self.dam_max_age_at_conception_days is not None
        )


def _require_nonnegative_int(name: str, value: Any) -> None:
    if value is None:
        return
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise PyPedalUsageError(
            f"VitalRateProfile.{name} must be a nonnegative integer or None, "
            f"not {value!r}."
        )


def _require_min_le_max(role: str, minimum: Any, maximum: Any) -> None:
    if minimum is None or maximum is None:
        return
    if minimum > maximum:
        raise PyPedalUsageError(
            f"VitalRateProfile {role} min age at conception ({minimum}) "
            f"exceeds max ({maximum})."
        )


def _is_founder(animal: Any, missing: Any) -> bool:
    return (
        animal.sireID == missing or str(animal.sireID) == str(missing)
    ) and (animal.damID == missing or str(animal.damID) == str(missing))


def _is_half_founder(animal: Any, missing: Any) -> bool:
    sire_missing = animal.sireID == missing or str(animal.sireID) == str(missing)
    dam_missing = animal.damID == missing or str(animal.damID) == str(missing)
    return sire_missing != dam_missing


def _has_recorded_chronology(animal: Any) -> bool:
    return recorded_year(animal) is not None or isinstance(
        getattr(animal, "bd", None), datetime.date
    )


def estimate_birth_date_ranges(pedobj: Any, profile: Optional[VitalRateProfile] = None) -> None:
    """
    Fill ``birth_date_estimate`` from an explicit ``VitalRateProfile``.

    Exact offspring ``bd`` drives range estimation. Year-only offspring are
    not converted into invented dates. Recorded ``bd`` / ``by`` are never
    overwritten. Default load does not call this.
    """
    if profile is None:
        profile = pedobj.kw.get("vital_rate_profile")
    if profile is None:
        return
    if not isinstance(profile, VitalRateProfile):
        raise PyPedalUsageError(
            "estimate_birth_date_ranges() requires a VitalRateProfile instance."
        )
    profile.validate()

    missing = pedobj.kw["missing_parent"]
    by_id = {animal.animalID: animal for animal in pedobj.pedigree}
    children: dict[Any, list[tuple[str, Any]]] = {}
    for child in pedobj.pedigree:
        for role, parent_id in (("sire", child.sireID), ("dam", child.damID)):
            if parent_id == missing or str(parent_id) == str(missing):
                continue
            children.setdefault(parent_id, []).append((role, child))

    for animal in pedobj.pedigree:
        estimate = getattr(animal, "birth_date_estimate", None)
        if estimate is None:
            animal.birth_date_estimate = BirthDateEstimate()
            estimate = animal.birth_date_estimate
        if _has_recorded_chronology(animal):
            continue
        _apply_offspring_range(
            animal, children.get(animal.animalID, []), profile, by_id
        )
        _apply_founder_typical(animal, children.get(animal.animalID, []), profile, missing)


def _apply_offspring_range(
    parent: Any,
    offspring_roles: Iterable[tuple[str, Any]],
    profile: VitalRateProfile,
    by_id: dict,
) -> None:
    earliest_candidates = []
    latest_candidates = []
    for role, child in offspring_roles:
        child_bd = getattr(child, "bd", None)
        if not isinstance(child_bd, datetime.date):
            continue
        window = _parent_window_from_offspring(role, child_bd, profile)
        if window is None:
            continue
        earliest_candidates.append(window[0])
        latest_candidates.append(window[1])
    if not earliest_candidates:
        return
    earliest = max(earliest_candidates)
    latest = min(latest_candidates)
    if earliest > latest:
        raise PyPedalValidationError(
            f"Vital-rate profile {profile.name!r} is inconsistent with dated "
            f"offspring of animal {_ident(parent)}: earliest {earliest.isoformat()} "
            f"is after latest {latest.isoformat()}."
        )
    estimate = parent.birth_date_estimate
    estimate.earliest = earliest
    estimate.latest = latest
    estimate.source = _SOURCE_VITAL_RATE
    estimate.profile = profile.name


def _parent_window_from_offspring(
    role: str,
    offspring_bd: datetime.date,
    profile: VitalRateProfile,
) -> Optional[tuple[datetime.date, datetime.date]]:
    if role == "sire":
        if not profile.sire_range_ready():
            return None
        min_age = profile.sire_min_age_at_conception_days
        max_age = profile.sire_max_age_at_conception_days
    elif role == "dam":
        if not profile.dam_range_ready():
            return None
        min_age = profile.dam_min_age_at_conception_days
        max_age = profile.dam_max_age_at_conception_days
    else:
        return None
    conception = offspring_bd - datetime.timedelta(days=int(profile.gestation_days))
    earliest = conception - datetime.timedelta(days=int(max_age))
    latest = conception - datetime.timedelta(days=int(min_age))
    return earliest, latest


def _apply_founder_typical(
    animal: Any,
    offspring_roles: Iterable[tuple[str, Any]],
    profile: VitalRateProfile,
    missing: Any,
) -> None:
    if profile.founder_typical_age_at_progeny_days is None:
        return
    if not _is_founder(animal, missing) or _is_half_founder(animal, missing):
        return
    if _has_recorded_chronology(animal):
        return
    dated = [
        child.bd
        for _role, child in offspring_roles
        if isinstance(getattr(child, "bd", None), datetime.date)
    ]
    if not dated:
        return
    typical = min(dated) - datetime.timedelta(
        days=int(profile.founder_typical_age_at_progeny_days)
    )
    estimate = animal.birth_date_estimate
    if estimate.earliest is not None and typical < estimate.earliest:
        return
    if estimate.latest is not None and typical > estimate.latest:
        return
    estimate.typical = typical
    estimate.source = _SOURCE_VITAL_RATE
    estimate.profile = profile.name
