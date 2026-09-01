###############################################################################
# NAME: pyp_results.py
# VERSION: see PyPedal.__version__
# LICENSE: LGPL
# Written for PyPedal 4.1; this module is not part of the original 2.0.4 sources.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
"""
Dict-compatible analysis result types.

These are thin ``dict`` subclasses. Existing 4.0.x callers that use
``result["key"]``, ``result.get(...)``, ``isinstance(result, dict)``,
``dict(result)``, iteration, and equality with a plain dict keep working.

Convenience properties read the same stored keys. They do not cache,
copy, mutate, or invent values. Import them from this module::

    from PyPedal.pyp_results import InbreedingResult, ProgressCallback

This module imports nothing else from PyPedal.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

# Current animalID after a default (renumbered) load is an int. String IDs
# remain possible on unrestored or non-renumbered records, so the public
# alias is the union actually stored as mapping keys.
AnimalId = int | str

# Optional long-running-operation reporter. Callers pass ``progress=None``
# (the default) to keep existing behaviour. ``done`` is completed units;
# ``total`` is the known unit count, or ``None`` when the length is not
# known cheaply. Callback exceptions propagate unchanged; they are not
# translated into ``PyPedalError``. There is no cancellation API.
ProgressCallback = Callable[[int, int | None], None]


class InbreedingResult(dict):
    """Successful ``pyp_nrm.inbreeding`` result.

    Keys are exactly ``fx``, ``metadata``, and optionally ``rel_dict``.
    ``rel_dict`` is omitted when relationships were not requested or the
    method cannot supply them. The ``rel_dict`` property then returns
    ``None`` rather than fabricating an empty mapping.
    """

    @property
    def fx(self) -> dict[AnimalId, float]:
        return self["fx"]

    @property
    def metadata(self) -> dict[str, Any]:
        return self["metadata"]

    @property
    def rel_dict(self) -> dict[str, Any] | None:
        return self.get("rel_dict")


class EffectiveFoundersResult(dict):
    """Successful Lacy effective-founder result.

    Shared by ``effective_founders_lacy`` and ``a_effective_founders_lacy``.
    Keys are exactly ``fa_animal_count``, ``fa_founder_count``,
    ``fa_descendant_count``, and ``fa_effective_founders``.
    """

    @property
    def fa_animal_count(self) -> int:
        return self["fa_animal_count"]

    @property
    def fa_founder_count(self) -> int:
        return self["fa_founder_count"]

    @property
    def fa_descendant_count(self) -> int:
        return self["fa_descendant_count"]

    @property
    def fa_effective_founders(self) -> float:
        return self["fa_effective_founders"]


class MatingCoIGroupResult(dict):
    """Successful ``mating_coi_group`` result.

    Keys are exactly ``matings`` and ``metadata``. Nested mating keys remain
    ``(animal_a, animal_b)`` tuples of current animal IDs.
    """

    @property
    def matings(self) -> dict[tuple[AnimalId, AnimalId], float]:
        return self["matings"]

    @property
    def metadata(self) -> dict[str, Any]:
        return self["metadata"]
