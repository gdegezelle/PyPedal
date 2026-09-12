"""Canonical Griffon file integrity after the confirmed identity merges.

This is reference-data validation. It does not detect fuzzy duplicates.
"""
from _pedhelpers import canonical_griffon_path, named_griffon_path

EXPECTED_N = 97999
SURVIVING = {37481, 54586}
REMOVED = {37482, 54587}


def _read_rows(path, nfields):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split(",")
            assert len(fields) == nfields, line
            rows.append(fields)
    return rows


def test_canonical_griffon_merged_identities_and_parent_closure():
    scientific = _read_rows(canonical_griffon_path(), 5)
    named = _read_rows(named_griffon_path(), 6)
    assert len(scientific) == EXPECTED_N
    assert len(named) == EXPECTED_N

    sci_ids = [int(row[0]) for row in scientific]
    named_ids = [int(row[0]) for row in named]
    assert sci_ids == named_ids
    assert len(set(sci_ids)) == EXPECTED_N
    present = set(sci_ids)
    assert SURVIVING <= present
    assert present.isdisjoint(REMOVED)

    dangling = []
    stale_parents = []
    for sci_row, named_row in zip(scientific, named):
        assert sci_row[:5] == named_row[:5]
        animal, sire, dam, _sex, _bdate = sci_row
        for parent in (sire, dam):
            if parent in {"37482", "54587"}:
                stale_parents.append((animal, parent))
            elif parent != "0" and int(parent) not in present:
                dangling.append((animal, parent))
    assert stale_parents == []
    assert dangling == []
