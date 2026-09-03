# Recipes

Copy these scripts as they are. Each one writes a small pedigree in a
temporary directory unless it says **checkout / sdist only**.

Common mistakes to avoid in every recipe:

- Using a file ID where a current `animalID` is required
- Leaving `output=True` (the default) and scattering `.dat` files
- Using `tabular` inbreeding on a ~100,000-animal file

---

## 1. Calculate inbreeding

**Input.** A three-column integer pedigree. **Code.** `pyp_nrm.inbreeding`.
**Output.** *F* for each current ID. **Interpretation.** On Mrode, animal 5
is 0.125. **Mistake.** Looking up `result["fx"][originalID]` after
renumbering has changed IDs.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_nrm

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
print(result["fx"][5])
```

Prints `0.125`. See [Inbreeding](inbreeding.md).

---

## 2. Calculate the relationship between two animals

**Input.** Same pedigree. **Code.** `pyp_metrics.relationship`.
**Output.** Additive *a<sub>ij</sub>*. **Interpretation.** 4 and 3 are
half-sibs (0.25); 1 and 2 are unrelated founders (0.0). **Mistake.**
Treating 0.0 as “ID not found.” Unresolved IDs raise an error.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
print(pyp_metrics.relationship(4, 3, ped))
print(pyp_metrics.relationship(1, 2, ped))
```

Prints `0.25` then `0.0`. See [Relationships](relationships.md).

---

## 3. Compare one prospective mating

**Input.** Two current IDs. **Code.** `mating_coi`. **Output.** Offspring
*F* = *A<sub>ij</sub>* / 2. **Interpretation.** Mating 4 × 3 gives 0.125.
**Mistake.** Expecting a new animal to appear in `ped.pedigree`. The
call is read-only.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
print(pyp_metrics.mating_coi(4, 3, ped))
```

Prints `0.125`. See [Test mating](mating.md).

---

## 4. Compare several prospective matings

**Input.** A list of pairs. **Code.** `mating_coi_group`. **Output.** A
dict of pair → *F*. **Interpretation.** Compare candidates you named;
PyPedal does not form every possible pair. **Mistake.** Passing
underscore strings as the primary API. Prefer `[(4, 3), (1, 2)]`.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
got = pyp_metrics.mating_coi_group([(4, 3), (1, 2)], ped)
print(got["matings"][(4, 3)], got["matings"][(1, 2)])
```

Prints `0.125 0.0`.

---

## 5. Calculate effective founder metrics

**Input.** Lacy’s Appendix A file from a **checkout or sdist**.
**Code.** `effective_founders_lacy`. **Output.** Founder count 3,
*f<sub>e</sub>* ≈ 2.91. **Interpretation.** Three founders contribute
unequally. **Mistake.** Treating *f<sub>e</sub>* as *N<sub>e</sub>* or as
a head-count of living animals.

```python
import shutil
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

examples = Path("PyPedal/examples")
work = Path(tempfile.mkdtemp())
pedfile = work / "new_lacy.ped"
shutil.copy(examples / "new_lacy.ped", pedfile)
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
got = pyp_metrics.effective_founders_lacy(ped)
print(got["fa_founder_count"], round(got["fa_effective_founders"], 2))
```

See [Effective founders](effective-founders.md) and
[Lacy and Boichard](lacy-and-boichard.md).

---

## 6. Run gene dropping

**Input.** Any loaded pedigree. **Code.** `effective_founder_genomes`.
**Output.** Mean N<sub>g</sub> across replicates. **Interpretation.** Not
a founder head-count. **Mistake.** Using `rounds=3` as if it were a
published precision.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
ng = pyp_metrics.effective_founder_genomes(
    ped, rounds=20, seed=31, output=False
)
print(round(ng, 4))
```

See [Gene dropping](gene-dropping.md).

---

## 7. Calculate generation intervals

**Input.** A pedigree with recorded birth years. Mrode has none, so
intervals that need years will skip those animals. **Code.**
`set_generation` then `generation_intervals`. **Output.** Mean parent
ages on four paths. **Interpretation.** `igen` is pedigree depth; `gen`
is an input label. **Mistake.** Copying `igen` into `gen` for Boichard
metrics.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_utils

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
pyp_utils.set_generation(ped)
print([int(a.igen) for a in ped.pedigree])
```

Prints `[1, 1, 2, 2, 3, 4]`. See [Generation intervals](generation-intervals.md).

---

## 8. Generate a metadata PDF

**Input.** A loaded pedigree and the `reports` extra. **Code.**
`pdf_pedigree_metadata`. **Output.** A PDF path. **Interpretation.**
Summary counts, not a drawing. **Mistake.** Expecting a GUI menu.

Requires ReportLab (`pip install -e ".[reports]"`).

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_reports

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
print(pyp_reports.pdf_pedigree_metadata(
    ped, reportfile=str(work / "metadata.pdf"),
))
```

See [PDF reports](pdf-reports.md).

---

## 9. Generate a three-generation pedigree PDF

**Input.** A current `animalID` (here, 5). **Code.** `pdf_three_gen_ped`.
**Output.** A 15-slot PDF page. **Interpretation.** Unknown parents show
as `(Unknown Parent)`. **Mistake.** Passing `originalID` after
renumbering has changed it.

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_reports

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")
ped = load_pedigree(options={
    "pedfile": str(pedfile), "pedformat": "asd",
    "messages": "quiet", "pedigree_summary": 0,
})
print(pyp_reports.pdf_three_gen_ped(
    5, ped, reportfile=str(work / "three_gen.pdf"),
))
```

---

## 10. Load and inspect the large Griffon sample

**Input.** `PyPedal/examples/griffonbruxellois_2026_pyp.ped` —
**checkout / sdist only**. Comma-separated `asdxb` (no name column).
The named desktop companion is `griffonbruxellois_2026_named_pyp.ped`
(`asdxbn`); genealogy matches this file.
**Code.** Load with `asdxb` and `sepchar=","`, then print counts.
**Output.** 98,001 records, 6,689 founders, 915 half-founders, 3,997
unknown dates. **Interpretation.** Griffon Bruxellois 2026 export
(recorded births through 2025); curated project sample, not a registry
extract. **Mistake.** Running `inbreeding(method="tabular")` on this
file (~80 GB dense matrix). Use `meu_luo`.

```python
import shutil
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree

examples = Path("PyPedal/examples")
work = Path(tempfile.mkdtemp())
pedfile = work / "griffonbruxellois_2026_pyp.ped"
shutil.copy(examples / "griffonbruxellois_2026_pyp.ped", pedfile)

ped = load_pedigree(options={
    "pedfile": str(pedfile),
    "pedformat": "asdxb",
    "sepchar": ",",
    "messages": "quiet",
    "pedigree_summary": 0,
})
print(len(ped.pedigree))
print(ped.metadata.num_unique_founders)
```

See [Large pedigrees](large-pedigrees.md).
