# API overview

This is a map of the functions most users call. It is not a dump of every
internal helper.

## Load

| Call | Module | Role |
|---|---|---|
| `load_pedigree` | `pyp_newclasses` | Construct and load; optional `progress=` on record reading |
| `loadPedigree` | `pyp_newclasses` | Same function, older name |
| `NewPedigree.save` | `pyp_newclasses` | Write a text pedigree |
| `NewPedigree.savedb` | `pyp_newclasses` | Write SQLite |

## Inbreeding and relationships

| Call | Module | Role |
|---|---|---|
| `inbreeding` | `pyp_nrm` | *F* by method (`tabular`, `vanraden`, `meu_luo`, …); returns an `InbreedingResult` dict; `meu_luo` / `mod_meu_luo` accept `progress=` |
| `relationship` | `pyp_metrics` | Pairwise *a<sub>ij</sub>* |
| `mating_coi` | `pyp_metrics` | One prospective offspring *F* |
| `mating_coi_group` | `pyp_metrics` | Explicit list of pairs; returns a `MatingCoIGroupResult` dict |

## Founders, ancestors, gene drop

| Call | Module | Role |
|---|---|---|
| `effective_founders_lacy` / `a_effective_founders_lacy` | `pyp_metrics` | Lacy *f<sub>e</sub>* (scalable default / dense-NRM small-pedigree form); both return `EffectiveFoundersResult` |
| `a_effective_founders_boichard` | `pyp_metrics` | Boichard *f<sub>e</sub>* |
| `a_effective_ancestors_definite` | `pyp_metrics` | Boichard *f<sub>a</sub>*; optional `progress=` |
| `a_effective_ancestors_indefinite` | `pyp_metrics` | Approximate *f<sub>a</sub>*; optional `progress=` |
| `effective_founder_genomes` | `pyp_metrics` | Gene-drop N<sub>g</sub>; optional `progress=` after each round |
| `dropped_ancestral_inbreeding` | `pyp_metrics` | Ancestral *F* by gene dropping; optional `progress=` after each round |

## Generations and utilities

| Call | Module | Role |
|---|---|---|
| `set_generation` | `pyp_utils` | Assign `igen` |
| `set_age` | `pyp_utils` | Demographic year-offset |
| `generation_intervals` | `pyp_metrics` | Parent age at oldest son/daughter |
| `reorder` / `renumber` | `pyp_utils` | Order and 1-based IDs (also run from `load`) |

## Output

| Call | Module | Role |
|---|---|---|
| `pdf_pedigree_metadata` | `pyp_reports` | Metadata PDF (`reports` extra) |
| `pdf_three_gen_ped` | `pyp_reports` | 15-slot pedigree PDF |
| `draw_pedigree` | `pyp_graphics` | Graphviz drawing (`graphics` extra) |
| `pickle_pedigree` | `pyp_io` | Serialise a `NewPedigree` |

## Desktop

`python -m PyPedal`, `pypedal`, and `pypedal-gui` launch the PySide6
desktop application (`gui` extra). PyPedal 4.2.0 does not provide a
command-oriented CLI and has no analysis subcommands.

The type alias `ProgressCallback` lives in `pyp_results`. It is
`Callable[[int, int | None], None]`. Not every analysis accepts
`progress=`. There is no cancellation API.

## Genomic note

`pyp_snp` implements VanRaden (2008) **Method 1** genomic relationship
matrices. Methods 2 and 3 are outside PyPedal 4.0. That path is not the
1992 pedigree VanRaden inbreeding method.

## Errors

Public failures raise `PyPedalError` subclasses (`PyPedalUsageError`,
`PyPedalDependencyError`, `PyPedalPedigreeStructureError`, and others).
They do not exit 0.
