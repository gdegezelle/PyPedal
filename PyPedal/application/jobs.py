"""Application adapters for scientific analysis jobs.

This module does not reimplement formulas. Each function calls an existing
library entry point with ``output=False`` so the desktop never writes
analysis ``.dat`` files implicitly.

Mutation
--------
- ``run_inbreeding``: mutates ``animal.fa`` (existing ``inbreeding`` API).
- ``run_effective_founders``: mutates the pedigree only if Lacy's
  compatibility path actually renumbers.
- ``run_relationship``, ``run_mating_coi``, ``run_mating_coi_group``,
  ``run_theoretical_ne``, ``save_pedigree``, PDF exports: read-only with
  respect to coefficients (save/PDF write a file the user chose).

Progress
--------
Only Meuwissen-Luo forwards ``progress``. Lacy, relationship, mating, and
theoretical Ne have no library progress callback; callers must not invent
fake percentages.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from PyPedal.application.inbreeding import InbreedingByYearRow, group_inbreeding_by_year
from PyPedal.application.session import PedigreeSession
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_metrics import (
    effective_founders_lacy,
    mating_coi,
    mating_coi_group,
    relationship,
    theoretical_ne_from_metadata,
)
from PyPedal.pyp_newclasses import NewPedigree
from PyPedal.pyp_nrm import inbreeding
from PyPedal.pyp_reports import pdf_pedigree_metadata, pdf_three_gen_ped
from PyPedal.pyp_results import (
    EffectiveFoundersResult,
    InbreedingResult,
    MatingCoIGroupResult,
    ProgressCallback,
)

_LACY_RENUMBER_MARK = "automatically renumbers the pedigree"


@dataclass(frozen=True, slots=True)
class PairwiseResult:
    """One pairwise coefficient (relationship or mating CoI)."""

    animal_a: int
    animal_b: int
    coefficient: float


@dataclass(frozen=True, slots=True)
class YearInbreedingOutcome:
    """Year grouping, plus whether Meuwissen-Luo ran on this call."""

    rows: tuple[InbreedingByYearRow, ...]
    computed_inbreeding: bool


@dataclass(frozen=True, slots=True)
class FoundersOutcome:
    """Lacy result plus whether implicit renumbering actually happened."""

    result: EffectiveFoundersResult
    implicit_renumber: bool


def require_pedigree(session: PedigreeSession) -> NewPedigree:
    """Return the loaded pedigree or raise a typed usage error."""
    pedigree = session.pedigree
    if pedigree is None:
        raise PyPedalUsageError("A pedigree must be open before running this analysis.")
    return pedigree


def parse_animal_id(raw: str, *, label: str = "Animal ID") -> int:
    """Parse a current/renumbered animal ID. No silent coercion."""
    text = str(raw).strip()
    if not text:
        raise PyPedalUsageError(f"{label} is required.")
    try:
        value = int(text)
    except (TypeError, ValueError) as exc:
        raise PyPedalUsageError(f"{label} must be an integer; got {raw!r}.") from exc
    return value


def run_inbreeding(
    session: PedigreeSession,
    progress: ProgressCallback | None = None,
) -> InbreedingResult:
    """Meuwissen-Luo inbreeding. Mutates ``animal.fa``. Caches the result."""
    pedigree = require_pedigree(session)
    result = inbreeding(
        pedigree,
        method="meu_luo",
        output=False,
        progress=progress,
    )
    session.inbreeding_result = result
    return result


def ensure_inbreeding(
    session: PedigreeSession,
    progress: ProgressCallback | None = None,
) -> tuple[InbreedingResult, bool]:
    """Return the cached inbreeding result, computing it once if needed."""
    cached = session.inbreeding_result
    if cached is not None:
        return cached, False
    return run_inbreeding(session, progress=progress), True


def run_inbreeding_by_year(
    session: PedigreeSession,
    progress: ProgressCallback | None = None,
) -> YearInbreedingOutcome:
    """Group existing (or newly computed) *F* by birth year.

    Does not run Meuwissen-Luo a second time when a cache is present.
    """
    pedigree = require_pedigree(session)
    result, computed = ensure_inbreeding(session, progress=progress)
    rows = group_inbreeding_by_year(pedigree, result.fx)
    return YearInbreedingOutcome(rows=rows, computed_inbreeding=computed)


def run_effective_founders(session: PedigreeSession) -> FoundersOutcome:
    """Lacy effective founders with ``output=False``.

    If the library actually auto-renumbers, pedigree-dependent caches other
    than this result are cleared because animal IDs changed.
    """
    pedigree = require_pedigree(session)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        result = effective_founders_lacy(pedigree, output=False)
    implicit = any(
        issubclass(item.category, DeprecationWarning) and _LACY_RENUMBER_MARK in str(item.message)
        for item in caught
    )
    if implicit:
        session.clear_analysis_cache()
    session.effective_founders_result = result
    return FoundersOutcome(result=result, implicit_renumber=implicit)


def run_relationship(
    session: PedigreeSession,
    animal_a: int,
    animal_b: int,
) -> PairwiseResult:
    """Pairwise numerator relationship. Read-only."""
    pedigree = require_pedigree(session)
    coefficient = float(relationship(animal_a, animal_b, pedigree))
    outcome = PairwiseResult(animal_a=animal_a, animal_b=animal_b, coefficient=coefficient)
    session.relationship_result = outcome
    return outcome


def run_mating_coi(
    session: PedigreeSession,
    animal_a: int,
    animal_b: int,
) -> PairwiseResult:
    """Prospective offspring inbreeding for one explicit pair. Read-only."""
    pedigree = require_pedigree(session)
    coefficient = float(mating_coi(animal_a, animal_b, pedigree))
    outcome = PairwiseResult(animal_a=animal_a, animal_b=animal_b, coefficient=coefficient)
    session.mating_pair_result = outcome
    return outcome


def run_mating_coi_group(
    session: PedigreeSession,
    pairs: Sequence[tuple[int, int]],
) -> MatingCoIGroupResult:
    """Prospective offspring inbreeding for an explicit list of pairs.

    The caller supplies the pairs. This never forms a Cartesian product.
    """
    pedigree = require_pedigree(session)
    if not pairs:
        raise PyPedalUsageError("At least one mating pair is required.")
    result = mating_coi_group(list(pairs), pedigree, names=0, gens=0)
    session.mating_group_result = result
    return result


def run_theoretical_ne(session: PedigreeSession) -> float:
    """Theoretical Ne from sire/dam metadata. Not a census Ne estimate."""
    pedigree = require_pedigree(session)
    value = float(theoretical_ne_from_metadata(pedigree, output=False))
    session.theoretical_ne = value
    return value


def _destination(path: Path, *, overwrite: bool) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.exists() and not overwrite:
        raise PyPedalUsageError(f"{resolved} already exists.")
    return resolved


def export_metadata_pdf(
    session: PedigreeSession,
    destination: Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write the existing pedigree-metadata PDF to an explicit path."""
    pedigree = require_pedigree(session)
    path = _destination(destination, overwrite=overwrite)
    pdf_pedigree_metadata(pedigree, reportfile=str(path))
    return path


def export_three_gen_pdf(
    session: PedigreeSession,
    animal_id: int,
    destination: Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write the existing three-generation pedigree PDF for one animal ID."""
    pedigree = require_pedigree(session)
    path = _destination(destination, overwrite=overwrite)
    pdf_three_gen_ped(animal_id, pedigree, reportfile=str(path))
    return path


def save_pedigree(
    session: PedigreeSession,
    destination: Path,
    *,
    overwrite: bool = False,
    pedformat: str | None = None,
    sepchar: str | None = None,
) -> Path:
    """Write the loaded pedigree to ``destination`` via ``NewPedigree.save``."""
    pedigree = require_pedigree(session)
    path = _destination(destination, overwrite=overwrite)
    fmt = pedformat
    if fmt is None and session.load_options is not None:
        fmt = session.load_options.pedformat
    if not fmt:
        fmt = "asd"
    separator = sepchar
    if separator is None and session.load_options is not None:
        separator = session.load_options.separator
    if separator is None:
        separator = " "
    ok = pedigree.save(filename=str(path), pedformat=fmt, sepchar=separator)
    if not ok:
        raise PyPedalUsageError(f"Could not save the pedigree to {path}.")
    return path
