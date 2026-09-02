"""Application-owned pedigree session.

A ``PedigreeSession`` is the explicit owner of the currently loaded
pedigree, the source path it came from, the options that produced a
successful load, and a small set of cached analysis results.

There is no process-wide singleton. Desktop code constructs a session
and passes it to application operations. Scientific state stays on
``NewPedigree``; this object only holds references.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from PyPedal.application.jobs import PairwiseResult
    from PyPedal.application.load import PedigreeOpenOptions
    from PyPedal.pyp_newclasses import NewPedigree
    from PyPedal.pyp_results import (
        EffectiveFoundersResult,
        InbreedingResult,
        MatingCoIGroupResult,
    )

SessionLifecycle = Literal["empty", "loaded"]


class PedigreeSession:
    """Explicit current-pedigree owner for a desktop or console application.

    Lifecycle
    ---------
    empty
        No pedigree. ``source_path`` and ``load_options`` are ``None``.
    loaded
        A pedigree that loaded successfully. Failed loads never enter
        this state and never replace an existing pedigree.
    replaced
        A successful load overwrites the previous pedigree and clears
        cached analysis results.
    cleared
        ``clear()`` returns the session to empty and drops the cache.
    """

    def __init__(self) -> None:
        self.pedigree: NewPedigree | None = None
        self.source_path: Path | None = None
        self.load_options: PedigreeOpenOptions | None = None
        self.inbreeding_result: InbreedingResult | None = None
        self.effective_founders_result: EffectiveFoundersResult | None = None
        self.relationship_result: PairwiseResult | None = None
        self.mating_pair_result: PairwiseResult | None = None
        self.mating_group_result: MatingCoIGroupResult | None = None
        self.theoretical_ne: float | None = None

    @property
    def is_empty(self) -> bool:
        return self.pedigree is None

    @property
    def state(self) -> SessionLifecycle:
        return "empty" if self.pedigree is None else "loaded"

    def replace_pedigree(
        self,
        pedigree: NewPedigree,
        source_path: Path,
        options: PedigreeOpenOptions,
    ) -> None:
        """Install a successfully loaded pedigree and drop cached results."""
        self.pedigree = pedigree
        self.source_path = source_path
        self.load_options = options
        self.clear_analysis_cache()

    def clear_analysis_cache(self) -> None:
        """Drop pedigree-dependent analysis results. Pedigree stays loaded."""
        self.inbreeding_result = None
        self.effective_founders_result = None
        self.relationship_result = None
        self.mating_pair_result = None
        self.mating_group_result = None
        self.theoretical_ne = None

    def clear(self) -> None:
        """Drop the current pedigree, path, options, and cached results."""
        self.pedigree = None
        self.source_path = None
        self.load_options = None
        self.clear_analysis_cache()
