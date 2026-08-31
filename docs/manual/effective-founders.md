# Effective founders

A **founder** is an animal with both parents unknown in this pedigree.
The raw founder count is `ped.metadata.num_unique_founders`. That count
says nothing about whether those founders contributed evenly to later
generations.

The **effective number of founders** *f<sub>e</sub>* answers: if founders
had contributed equally, how many would produce the same imbalance of
gene origin that we observe? When contributions are equal, *f<sub>e</sub>*
equals the founder count. When a few founders dominate, *f<sub>e</sub>* is
smaller.

This is **not** effective population size *N<sub>e</sub>*, **not** the
effective number of ancestors *f<sub>a</sub>*, and **not** the effective
number of founder genomes *N<sub>g</sub>*.

## Lacy’s *f<sub>e</sub>*

Lacy (1989) defines the founder equivalent as

> 1 / Σ *p<sub>k</sub>*²

where *p<sub>k</sub>* is founder *k*’s proportional contribution to the
descendants (animals that are not themselves founders).

PyPedal implements that with explicit phantom founders for unknown
parental genome contributions. Two entry points exist:

- `pyp_metrics.effective_founders_lacy(pedobj)` — the scalable default
  supported path. It walks founder-row contributions and does not form a
  dense numerator relationship matrix.
- `pyp_metrics.a_effective_founders_lacy(pedobj)` — the dense-NRM form,
  for small pedigrees and published-reference checks only. Do not use it
  on a large pedigree.

On Lacy’s Appendix A pedigree (`new_lacy.ped` in a **checkout or sdist**):

```python
import shutil
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

examples = Path("PyPedal/examples")  # checkout / sdist only
work = Path(tempfile.mkdtemp())
pedfile = work / "new_lacy.ped"
shutil.copy(examples / "new_lacy.ped", pedfile)

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

got = pyp_metrics.effective_founders_lacy(ped)
print(got["fa_founder_count"], round(got["fa_effective_founders"], 2))
```

That prints `3` and `2.91`. Three founders do not contribute equally, so
the effective number is a little under 3. Extra floating-point digits are
not biological precision.

These functions write a `*_fe_lacy*.dat` file next to the pedigree stem
when `output=True` (the default). Pass `output=False` to run the same
calculation without that analysis file. Work in a temporary directory.

## Boichard’s *f<sub>e</sub>*

Boichard, Maignel and Verrier (1997) estimate a related but distinct
*f<sub>e</sub>* from expected contributions to a **reference population**
you name. That is not automatically “the last inferred generation.” See
[Lacy and Boichard](lacy-and-boichard.md).

Use the function that matches the paper you mean.
