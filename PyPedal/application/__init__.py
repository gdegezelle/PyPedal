"""Qt-free application layer between a future desktop and the library.

``PyPedal.application`` owns session/load orchestration, browse access,
GUI-independent error presentation, and analysis job adapters. It does
not implement scientific formulas and it does not import GUI toolkits.

The intended dependency direction is::

    PyPedal.desktop
            ↓
    PyPedal.application
            ↓
    PyPedal scientific modules

Scientific ``pyp_*.py`` modules must not import this package or
``PyPedal.desktop``. Application adapters may call scientific APIs such
as ``pyp_nrm`` and ``pyp_metrics``. The reverse is not allowed.
"""

from __future__ import annotations

from PyPedal.application.errors import (
    EXIT_STATUS,
    ApplicationErrorInfo,
    describe_exception,
    exit_status_for,
)
from PyPedal.application.export import (
    export_inbreeding_csv,
    export_mating_group_csv,
    export_year_inbreeding_csv,
    write_csv,
    write_text,
)
from PyPedal.application.inbreeding import (
    UNKNOWN_YEAR_TOKENS,
    InbreedingByYearRow,
    group_inbreeding_by_year,
)
from PyPedal.application.jobs import (
    FoundersOutcome,
    PairwiseResult,
    YearInbreedingOutcome,
    ensure_inbreeding,
    export_metadata_pdf,
    export_three_gen_pdf,
    parse_animal_id,
    require_pedigree,
    run_effective_founders,
    run_inbreeding,
    run_inbreeding_by_year,
    run_mating_coi,
    run_mating_coi_group,
    run_relationship,
    run_theoretical_ne,
    save_pedigree,
)
from PyPedal.application.load import (
    PedigreeOpenOptions,
    load_into_session,
    normalize_sepchar,
    resolve_source_path,
)
from PyPedal.application.session import PedigreeSession
from PyPedal.application.tables import (
    BROWSE_COLUMNS,
    PedigreeTableSource,
    TableColumn,
)
from PyPedal.pyp_results import (
    EffectiveFoundersResult,
    InbreedingResult,
    MatingCoIGroupResult,
)

__all__ = [
    "BROWSE_COLUMNS",
    "EXIT_STATUS",
    "UNKNOWN_YEAR_TOKENS",
    "ApplicationErrorInfo",
    "EffectiveFoundersResult",
    "FoundersOutcome",
    "InbreedingByYearRow",
    "InbreedingResult",
    "MatingCoIGroupResult",
    "PairwiseResult",
    "PedigreeOpenOptions",
    "PedigreeSession",
    "PedigreeTableSource",
    "TableColumn",
    "YearInbreedingOutcome",
    "describe_exception",
    "ensure_inbreeding",
    "exit_status_for",
    "export_inbreeding_csv",
    "export_mating_group_csv",
    "export_metadata_pdf",
    "export_three_gen_pdf",
    "export_year_inbreeding_csv",
    "group_inbreeding_by_year",
    "load_into_session",
    "normalize_sepchar",
    "parse_animal_id",
    "require_pedigree",
    "resolve_source_path",
    "run_effective_founders",
    "run_inbreeding",
    "run_inbreeding_by_year",
    "run_mating_coi",
    "run_mating_coi_group",
    "run_relationship",
    "run_theoretical_ne",
    "save_pedigree",
    "write_csv",
    "write_text",
]
