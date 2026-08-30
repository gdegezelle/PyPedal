# Generation intervals

PyPedal keeps two generation-like fields. They are **not** aliases.

| Field | What it is |
|---|---|
| `gen` | An **input annotation** from pedformat `g`. Never computed by `set_generation`. |
| `igen` | An **inferred pedigree depth**: founders at 1, each animal one more than the deepest known parent. Assigned only after `set_generation`. |

Do not copy `igen` into `gen`. Boichard metrics and gene dropping that
read a generation label use `gen`, not `igen`.

## Assigning inferred generations

Default **load does not assign `igen`**. After an ordinary load every
animal still has the missing-generation sentinel (`-999.0`).

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_utils

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

print([a.igen for a in ped.pedigree])
pyp_utils.set_generation(ped)
print([a.igen for a in ped.pedigree])
print([a.gen for a in ped.pedigree])
```

That prints the sentinel six times, then `igen` `[1.0, 1.0, 2.0, 2.0, 3.0, 4.0]`,
while `gen` stays at the sentinel.

You can also pass `"set_generations": True` in the `load_pedigree`
options dict.

## Generation intervals from ages

`pyp_metrics.generation_intervals(pedobj, units='y')` estimates the
average age of parents at the birth of their **oldest son and oldest
daughter**, along four paths (sire–son, sire–daughter, dam–son,
dam–daughter) plus an overall mean of those path means. That is the
definition Cole documented for this function: selection is treated as
happening at the first (oldest) offspring of each parent.

Animals whose recorded birth year is `None` are skipped. That is why
unknown chronology must stay `None` rather than a fake year. See
[Birth dates and chronology](birth-dates-and-chronology.md).

`generation_intervals_all` is a related path-mean routine. It still
records **one** offspring per parent per path (the last qualifying child
encountered), not a true average over every recorded offspring. Use
`generation_intervals` when you want the documented oldest-offspring
definition.

These calculations need known sexes and recorded years. They are
demographic summaries, not *N<sub>e</sub>* from inbreeding change.

## Pattie coefficients

Pedformat `p` can **store** a supplied generation coefficient. PyPedal 4
does **not** compute Pattie (1965) coefficients. Setting
`gen_coeff=True` is outside the supported domain.
