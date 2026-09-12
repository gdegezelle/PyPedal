# Equivalent Complete Generations

**Equivalent Complete Generations** (ECG) measures how much of an
animal's ancestry is known. A fully recorded generation contributes 1.
Incomplete generations contribute in proportion to the fraction of
ancestors that are known. Unknown ancestry contributes 0.

This is the completeness metric defined by Maignel, Boichard and Verrier
(1996), *Interbull Bulletin* 14:49–54, p.50: the number of complete
generation equivalents is the sum of the proportion of known ancestors
over all generations traced.

## What it reports

For each known ancestral **pedigree slot** at generation *g* from the
animal (parents *g* = 1, grandparents *g* = 2, and so on), add
1 / 2<sup>*g*</sup>. The animal itself is not a slot and contributes 0.

Consequences:

- A founder (both parents unknown) has ECG = 0.
- An animal with two known founder parents has ECG = 1.
- Two complete recorded generations contribute 2.
- An unknown parent adds nothing. There is no phantom founder for ECG.
- The same ancestor occupying several slots counts once **per slot**,
  not once per identity.
- There is no depth parameter. The entire known pedigree is used.

ECG is not the historical `pedigree_completeness` function.
`pedigree_completeness` remains available as a **legacy** metric with
different arithmetic (unique identities on each parental side, a chosen
depth, and a slot denominator). Prefer ECG.

## Function

```python
from PyPedal import pyp_metrics

ecg = pyp_metrics.equivalent_complete_generations(ped, output=False)
print(ecg[5], ped.pedigree[4].ecg)
```

The return value is a mapping of current `animalID` to ECG. The same
number is stored on `animal.ecg`. Parents must already precede offspring
in the pedigree (the usual load path does this). The call does not
reorder the pedigree.

Pass `output=False` to calculate without writing `{filetag}_ecg_.dat`.
The default `output=True` writes that analysis file, like neighbouring
metrics.

## Reference

Maignel, L., Boichard, D., Verrier, E. 1996. Genetic variability of
French dairy breeds estimated from pedigree information. *Interbull
Bulletin* 14:49–54.
