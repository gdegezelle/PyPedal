# Test mating

A **test mating** answers: if these two animals produced an offspring,
what would that offspring’s inbreeding coefficient be?

PyPedal computes that number **without adding a child to the pedigree**.
No phantom record is inserted. The calculation is read-only.

## The formula

For two **different** prospective parents *i* and *j*:

> *F*<sub>offspring</sub> = *A<sub>ij</sub>* / 2

that is, half the additive relationship of the parents.

For a **self** mating (*i* with *i*):

> *F*<sub>offspring</sub> = (1 + *F<sub>i</sub>*) / 2

IDs are current `animalID` values. The pedigree must already be
renumbered (the default load). The desktop Mating page selects those
IDs through the same name / original ID / current ID search as
Relationship. Duplicate names must be disambiguated explicitly. Sex is
shown in the match list; PyPedal does not swap the two animals.

## One pair: `mating_coi`

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics

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

print(pyp_metrics.mating_coi(4, 3, ped))
print(pyp_metrics.mating_coi(1, 2, ped))
print(pyp_metrics.mating_coi(5, 5, ped))
```

On Mrode that prints `0.125`, `0.0`, and `0.5625`.

- 4 × 3: half-sibs; offspring *F* = 0.25 / 2 = **0.125**.
- 1 × 2: unrelated founders; offspring *F* = **0.0**.
- 5 × 5: selfing; (1 + 0.125) / 2 = **0.5625**.

`gens` may be `0` or `-1` (both: use the full available pedigree). Other
values raise `PyPedalUsageError`. Truncated-generation approximations
are not implemented.

## Several pairs: `mating_coi_group`

Pass an explicit list of pairs. PyPedal does **not** form a Cartesian
product and does not choose mates for you.

```python
got = pyp_metrics.mating_coi_group([(4, 3), (1, 2)], ped)
print(got["matings"][(4, 3)])
print(got["matings"][(1, 2)])
```

That prints `0.125` and `0.0`. Result keys are tuples of current
animal IDs. The return is a `MatingCoIGroupResult` (`dict` subclass);
`got["matings"]` remains the supported 4.x access, and `got.matings` is
the same mapping.

Underscore strings such as `["4_3", "1_2"]` are accepted as
**compatibility** input when parsing is unambiguous. Prefer the pair
list. `names=1` is deprecated: it looks up unique string identities, not
call names.

## What this is not

- Not a mate-selection engine
- Not a write to the pedigree
- Not a reason to form a dense NRM for a large file

`mating_coi` and `mating_coi_group` use the same exact selected
relationship calculation as `relationship()` (half of *A<sub>ij</sub>*
for a distinct pair). They do not build an ancestor sub-NRM per pair.
A single test mating is suitable on large pedigrees; an explicit group
still costs one calculation per supplied pair.

See [Relationships](relationships.md) for *A<sub>ij</sub>* itself.
