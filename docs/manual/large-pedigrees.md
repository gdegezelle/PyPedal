# Large pedigrees

Small textbook files (Mrode, Lacy) can use the default tabular inbreeding
method. A population of tens of thousands of animals cannot.

## Algorithm selection

| Task | Small pedigree | Large pedigree |
|---|---|---|
| Inbreeding *F* | `method="tabular"` | `method="meu_luo"` |
| One pairwise relationship | `relationship()` | `relationship()` (do not form a full NRM) |
| Test mating | `mating_coi` | `mating_coi` (still read-only; still no dense NRM) |
| Lacy *f<sub>e</sub>* | `effective_founders_lacy` | same function; work in a temp directory. Do not use `a_effective_founders_lacy` (dense NRM) |
| Gene-drop N<sub>g</sub> | modest `rounds` | more RAM/time; do not also form a dense NRM |

There is **no automatic switch** at 10,000 animals. If you leave
`inbreeding(method="tabular")` on a large file, PyPedal will try to build
a full numerator relationship matrix.

## Dense NRM warning

A dense float64 matrix for *n* animals needs about *n*² × 8 bytes.

For the curated Griffon sample, *n* = 98,001:

> 98,001 × 98,001 × 8 bytes ≈ **77 GB**

before Python and library overhead — on the order of **80 GB**. Do **not**
form a dense NRM for this sample. Use `meu_luo` for inbreeding.

If a required float64 allocation fails, PyPedal raises `PyPedalError`.
It does not silently continue in float32 or write a `*.bin` memmap in
the current directory.

## Meuwissen–Luo

`pyp_nrm.inbreeding(ped, method="meu_luo", output=False)` is linear in
the number of animals and returns coefficients only (no relationship
summaries). `mod_meu_luo` is the Quaas (1995) / Mrode Appendix B.2
variant. The desktop app uses `meu_luo` for its inbreeding button.

Pass `progress=` on `meu_luo` / `mod_meu_luo` if you want per-animal
completion events. The calculation with `progress=None` is unchanged.
See [Inbreeding](inbreeding.md).

The method requires contiguous 1-based `animalID` values with parents
before offspring. Default `renumber=True` establishes that. Calling
`inbreeding_meuwissen_luo` (or `inbreeding(..., method="meu_luo")`) on a
pedigree that does not satisfy the numbering raises `PyPedalUsageError`
rather than returning silently wrong coefficients.

## The Griffon sample

The repository ships **one** Griffon Bruxellois pedigree file:

`PyPedal/examples/griffonbruxellois_2026_pyp.ped`

It is a **Griffon Bruxellois 2026 export** (sample updated/exported in
2026, with recorded births through 2025). It is **checkout / sdist
only**. Wheels do not install `PyPedal/examples`. It is **curated project
data**, not an independently authoritative registry or studbook dump.

Load it as comma-separated `asdxb` with no padded delimiters and no
name column:

```python
ped = load_pedigree(options={
    "pedfile": "griffonbruxellois_2026_pyp.ped",
    "pedformat": "asdxb",
    "sepchar": ",",
    "messages": "quiet",
    "pedigree_summary": 0,
})
```

Observed load on that file (PyPedal 4 behaviour):

| Quantity | Value |
|---|---|
| Records | 98,001 |
| Founders (both parents unknown) | 6,689 |
| Half-founders (exactly one parent unknown) | 915 |
| Unknown recorded chronology (`None`) | 3,997 |
| Known birth years | 1870 … 2025 |
| Inferred generation `igen` after `set_generation` | 1 … 70 |

A gene-drop regression with `rounds=3`, `seed=31`, `chrometype="autosome"`,
`output=False` gave

> N<sub>g</sub> = 11.018378975785259

Lacy effective founders on the same load gave

> *f<sub>e</sub>* = 193.31434658869796

Those numbers are **deterministic dataset regressions**. They are not
scientific constants and they are not “the number of historical
founders.” *rounds=3* was a smoke setting for timing, not a recommended
analysis sample size.

Do not commit derived Griffon CSV or subset files. Tests that need a
small extract write a temporary subset from this file.

## Observed runtime (not a guarantee)

The following times were measured on a MacBook Pro M2 Pro, late 2023,
32 GB RAM, loading the 2026 canonical file with `asdxb` and
`sepchar=","`. They are **observed hardware-specific characterization**,
not a guaranteed benchmark, not a test threshold, and not a scientific
result.

| Analysis | Before the RC8 optimization | After |
|---|---|---|
| Meuwissen–Luo inbreeding | approximately 16 minutes | approximately 80 seconds |
| `effective_founders_lacy` | more than 20 minutes without completing | approximately 11 seconds |

## Performance notes

- Load from a working copy; constructing a pedigree writes a logfile on
  the PyPedal logger (not the process root logger).
- Use `messages="quiet"` and `pedigree_summary=0` if you do not need
  console dumps. Quiet does not override logging configured by the host.
- Pass `output=False` on inbreeding, gene dropping, Lacy and Boichard
  metrics, and coefficient dumps unless you want side-effect files.
- Drawings of 98,000 nodes are not a practical Graphviz workflow.
- The GUI inbreeding path already uses `meu_luo`.
- Load time is dominated by topological reordering when many offspring
  appear **before** their parents in the file. A smaller file that is
  badly ordered can take longer to load than a larger file that is
  already parent-before-offspring. Writing animals oldest-first is the
  cheapest preparation step you can do outside PyPedal.

See [Recipes](recipes.md) for a load-and-inspect sketch of this file.
