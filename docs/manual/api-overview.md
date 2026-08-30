# API overview

This is a map of the functions most users call. It is not a dump of every
internal helper.

## Load

| Call | Module | Role |
|---|---|---|
| `load_pedigree` | `pyp_newclasses` | Construct and load |
| `loadPedigree` | `pyp_newclasses` | Same function, older name |
| `NewPedigree.save` | `pyp_newclasses` | Write a text pedigree |
| `NewPedigree.savedb` | `pyp_newclasses` | Write SQLite |

## Inbreeding and relationships

| Call | Module | Role |
|---|---|---|
| `inbreeding` | `pyp_nrm` | *F* by method (`tabular`, `vanraden`, `meu_luo`, …) |
| `relationship` | `pyp_metrics` | Pairwise *a<sub>ij</sub>* |
| `mating_coi` | `pyp_metrics` | One prospective offspring *F* |
| `mating_coi_group` | `pyp_metrics` | Explicit list of pairs |

## Founders, ancestors, gene drop

| Call | Module | Role |
|---|---|---|
| `effective_founders_lacy` / `a_effective_founders_lacy` | `pyp_metrics` | Lacy *f<sub>e</sub>* |
| `a_effective_founders_boichard` | `pyp_metrics` | Boichard *f<sub>e</sub>* |
| `a_effective_ancestors_definite` | `pyp_metrics` | Boichard *f<sub>a</sub>* |
| `a_effective_ancestors_indefinite` | `pyp_metrics` | Approximate *f<sub>a</sub>* |
| `effective_founder_genomes` | `pyp_metrics` | Gene-drop N<sub>g</sub> |

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

`python -m PyPedal` or `pypedal` launches `PyPedal.pyp_app:main`
(`gui` extra).

## Genomic note

`pyp_snp` implements VanRaden (2008) **Method 1** genomic relationship
matrices. Methods 2 and 3 are outside PyPedal 4.0. That path is not the
1992 pedigree VanRaden inbreeding method.

## Errors

Public failures raise `PyPedalError` subclasses (`PyPedalUsageError`,
`PyPedalDependencyError`, `PyPedalPedigreeStructureError`, and others).
They do not exit 0.
