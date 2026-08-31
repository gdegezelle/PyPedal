# Effective ancestors

The **effective number of ancestors** *f<sub>a</sub>* accounts for
bottlenecks and overlapping ancestors. It asks: after the major
contributors are taken in turn, how many **marginal** ancestors explain
the gene pool of a defined group?

This is not a founder head-count, not Lacy’s *f<sub>e</sub>*, and not
gene-drop *N<sub>g</sub>*.

## Reference population

Boichard et al. (1997) require you to **define the population under
study** — the animals whose gene pool you are describing.

PyPedal supports two ways to name that group on the Boichard routines:

1. **`reference=`** — an iterable of current `animalID` integers. Order
   does not matter. Original or string IDs are not translated. Empty
   sets, duplicates, unknown IDs, and the missing-parent token are
   refused.
2. **`gen=` / default** — animals whose input `gen` field (pedformat
   `g`) equals a chosen generation. If `gen` is omitted, the
   **numerically** largest generation label is used. Non-numeric labels
   are refused rather than sorted as strings.

`reference=` and `gen=` cannot be combined.

These routines do **not** fall back to `igen`, birth year, or “animals
without offspring.” Mrode (`pedformat='asd'`) has no generation column;
analyse it with an explicit `reference=` list.

## Functions

- `pyp_metrics.a_effective_ancestors_definite` — Boichard’s “definite”
  ancestor algorithm
- `pyp_metrics.a_effective_ancestors_indefinite` — the approximate
  form that stops after the largest remaining contributions

Both write a `.dat` summary when `output=True` (the default). Pass
`output=False` to compute *f<sub>a</sub>* without that analysis file.

If an approximate contribution sequence is internally inconsistent,
PyPedal **refuses** the calculation rather than quietly repairing the
numbers. That is a data or algorithm-domain problem, not a number you
should round away.

Worked API notes and established checks are in
[Lacy and Boichard](lacy-and-boichard.md).
