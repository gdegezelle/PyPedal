"""Canonical Griffon file integrity for the adopted N=97,002 reference.

This is reference-data validation. It does not detect fuzzy duplicates.
The named file is the identity authority. The scientific file is the same
genealogy with the name column dropped.
"""
from _pedhelpers import canonical_griffon_path, named_griffon_path

EXPECTED_N = 97002


def _read_rows(path, nfields):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split(",", nfields - 1)
            assert len(fields) == nfields, line
            rows.append(fields)
    return rows


def test_canonical_griffon_named_and_scientific_are_graph_equivalent():
    scientific = _read_rows(canonical_griffon_path(), 5)
    named = _read_rows(named_griffon_path(), 6)
    assert len(scientific) == EXPECTED_N
    assert len(named) == EXPECTED_N

    sci_ids = [int(row[0]) for row in scientific]
    named_ids = [int(row[0]) for row in named]
    assert sci_ids == named_ids
    assert len(set(sci_ids)) == EXPECTED_N
    present = set(sci_ids)

    dangling = []
    self_parent = []
    for sci_row, named_row in zip(scientific, named):
        assert sci_row[:5] == named_row[:5]
        animal, sire, dam, _sex, _bdate = sci_row
        if sire == animal or dam == animal:
            self_parent.append(animal)
        for parent in (sire, dam):
            if parent != "0" and int(parent) not in present:
                dangling.append((animal, parent))
    assert self_parent == []
    assert dangling == []
    assert all(name.strip() for *_fields, name in named)
