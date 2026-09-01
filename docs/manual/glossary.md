# Glossary

**age**
: Legacy demographic year-offset (`by` minus 1800 when a year is known). Not current biological age.

**additive relationship (*a<sub>ij</sub>*)**
: Entry of the numerator relationship matrix. Diagonal *a<sub>ii</sub>* = 1 + *F<sub>i</sub>*.

**animalID**
: Current 1-based identifier after renumbering. Analysis functions use this domain.

**call name / display name**
: Pedformat `n`. Not a unique identity.

**coefficient of inbreeding (*F*)**
: Probability that two alleles drawn from an individual are identical by descent.

**effective ancestor number (*f<sub>a</sub>*)**
: Number of equally contributing ancestors, not necessarily founders, needed to match the gene-origin structure of a defined group (Boichard et al., 1997).

**effective founder number (*f<sub>e</sub>*)**
: Number of equally contributing founders needed to match observed contribution imbalance (Lacy, 1989; related Boichard estimator exists).

**fa**
: Animal inbreeding field. After load this may be a file column (`f`); after `inbreeding()` it is the computed coefficient. One field, not a loaded/computed pair. See [Object model](object-model.md).

**founder**
: Animal with both parents unknown (`founder == 'y'`).

**founder genome equivalent / N<sub>g</sub>**
: Effective number of equally frequent founder genes after gene dropping. Not a head-count of founder animals.

**gen**
: Input generation label (pedformat `g`). Not computed by `set_generation`.

**half-founder**
: Exactly one parent unknown. Not counted as a founder.

**igen**
: Inferred pedigree depth (founders at 1). Assigned only after `set_generation`.

**missing parent**
: Sentinel (default `0`) meaning the parent is unknown in this file.

**numerator relationship matrix (NRM, A)**
: Matrix of additive relationships among animals in the pedigree.

**originalID**
: Identifier as it appeared in the input file (or the integer hash of a unique string identity).

**pedformat**
: One-character-per-column description of the input file.

**reference population**
: The group of animals whose gene pool a Boichard metric describes. Must be named; not silently taken from `igen`.

**renumbering**
: Ordering parents before offspring and assigning sequential `animalID` values starting at 1.

**reordering**
: Arranging animals so parents appear before offspring, without necessarily assigning new IDs.

**string identity**
: Unique `A`/`S`/`D` codes. Not a call name.

**test mating**
: Inbreeding of a prospective offspring, computed without adding that animal to the pedigree.

**unknown chronology**
: Recorded `by` / `bd` stored as `None`.
