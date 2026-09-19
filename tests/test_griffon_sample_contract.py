"""Checked-in Griffon Bruxellois pedigree datasets.

The repository ships two Griffon pedigree files with identical genealogy:

    PyPedal/examples/griffonbruxellois_2026_pyp.ped       (asdxb, science)
    PyPedal/examples/griffonbruxellois_2026_named_pyp.ped  (asdxbn, desktop names)

Scientific tests that historically used smaller extracts derive those
extracts into a temporary directory from the scientific file. Those
subsets are not committed.
"""
import hashlib
import os
import struct
import tarfile

import pytest
from _pedhelpers import (
    CANONICAL_GRIFFON_PED,
    GRIFFON_1871_1890_IDS,
    GRIFFON_TEST_SMALL_IDS,
    NAMED_GRIFFON_PED,
    REPO,
    canonical_griffon_path,
    load_canonical_griffon,
    load_griffon_1871_1890,
    named_griffon_path,
    write_canonical_griffon_subset,
)

from PyPedal import pyp_chronology, pyp_metrics, pyp_nrm, pyp_utils


def _griffon_data_files(root):
    hits = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [
            d for d in dirnames
            if d not in {".git", "__pycache__", ".pytest_cache", "build", "dist"}
            and not d.endswith(".egg-info")
        ]
        for name in filenames:
            lowered = name.lower()
            if "griffon" not in lowered:
                continue
            ext = os.path.splitext(lowered)[1]
            if ext in {".ped", ".csv", ".dat", ".txt"}:
                hits.append(os.path.relpath(os.path.join(dirpath, name), root))
    return sorted(hits)


def test_repository_has_exactly_two_griffon_pedigree_datasets():
    hits = _griffon_data_files(REPO)
    expected = sorted(
        [
            os.path.join("PyPedal", "examples", CANONICAL_GRIFFON_PED),
            os.path.join("PyPedal", "examples", NAMED_GRIFFON_PED),
        ]
    )
    assert hits == expected, hits
    examples = os.path.join(REPO, "PyPedal", "examples")
    assert not os.path.exists(os.path.join(examples, "test_descendants.txt"))
    ini_hits = []
    for dirpath, _, filenames in os.walk(os.path.join(REPO, "PyPedal", "examples")):
        for name in filenames:
            if "griffon" in name.lower() and name.lower().endswith(".ini"):
                ini_hits.append(name)
    assert ini_hits == [], ini_hits


def _parse_asdxb(path):
    records = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            animal, sire, dam, sex, bdate = line.strip().split(",")
            records[str(int(animal))] = (sire, dam, sex, bdate)
    return records


def test_1871_1890_subset_matches_canonical_records(tmp_path):
    dest = tmp_path / "griffon_1871_1890.ped"
    write_canonical_griffon_subset(dest, GRIFFON_1871_1890_IDS)
    derived = _parse_asdxb(dest)
    canonical = _parse_asdxb(canonical_griffon_path())
    assert set(derived) == {str(i) for i in GRIFFON_1871_1890_IDS}
    assert len(derived) == 164
    for animal_id, fields in derived.items():
        assert fields == canonical[animal_id]
        sire, dam, _sex, bdate = fields
        _bd, by = pyp_chronology.parse_recorded_date(bdate)
        if by is not None:
            assert by <= 1890
        if sire != "0":
            assert sire in canonical
        if dam != "0":
            assert dam in canonical
    dangling = [
        sire
        for sire, dam, _sex, _bdate in derived.values()
        for parent in (sire, dam)
        if parent != "0" and parent not in derived
    ]
    assert dangling == ["32240"]


def test_test_small_subset_matches_canonical_records(tmp_path):
    dest = tmp_path / "griffon_test_small.ped"
    write_canonical_griffon_subset(dest, GRIFFON_TEST_SMALL_IDS)
    derived = _parse_asdxb(dest)
    canonical = _parse_asdxb(canonical_griffon_path())
    assert set(derived) == {str(i) for i in GRIFFON_TEST_SMALL_IDS}
    assert len(derived) == 18
    for animal_id, fields in derived.items():
        assert fields == canonical[animal_id]


def test_1871_1890_load_still_materializes_one_implicit_parent():
    ped = load_griffon_1871_1890()
    assert len(ped.pedigree) == 165
    assert ped.metadata.num_implicit_parents == 1
    assert pyp_utils.set_generation(ped)
    distribution = {}
    for animal in ped.pedigree:
        distribution[animal.igen] = distribution.get(animal.igen, 0) + 1
    assert distribution == {1: 122, 2: 20, 3: 8, 4: 12, 5: 3}


def test_sdist_contains_both_griffon_pedigree_datasets(tmp_path):
    from test_product_surface import _build_setuptools_artifact

    sdist = _build_setuptools_artifact(tmp_path / "dist", "sdist")
    with tarfile.open(sdist) as archive:
        names = archive.getnames()
    griffon_peds = sorted(
        name for name in names
        if "griffon" in name.lower() and name.endswith(".ped")
    )
    assert len(griffon_peds) == 2, griffon_peds
    assert any(name.endswith(CANONICAL_GRIFFON_PED) for name in griffon_peds)
    assert any(name.endswith(NAMED_GRIFFON_PED) for name in griffon_peds)
    leftover = [
        name for name in names
        if "griffon" in name.lower()
        and name.lower().endswith((".csv", ".dat", ".ini", ".txt"))
    ]
    assert leftover == [], leftover


# Exact LF bytes are part of the scientific data contract. Git checks this
# file out with eol=lf (.gitattributes) so Windows autocrlf cannot change it.
CANONICAL_SHA256 = "520f9d626384119aee6a93d279a9a0af9a3805161adbfc43243eb6f4f32bfbc8"
NAMED_SHA256 = "6bd07d40c4dbaf5d8dab95f0de59cf106515a7f80d95fdaf2e21f2a78d5d536e"
CANONICAL_N = 97002


def test_canonical_griffon_is_comma_asdxb_without_padded_delimiters():
    path = canonical_griffon_path()
    digest = hashlib.sha256()
    records = 0
    with open(path, "rb") as handle:
        for line in handle:
            digest.update(line)
            text = line.decode("utf-8").rstrip("\n")
            if not text.strip():
                continue
            records += 1
            assert ", " not in text
            assert " ," not in text
            fields = text.split(",")
            assert len(fields) == 5, text
    assert digest.hexdigest() == CANONICAL_SHA256
    assert records == CANONICAL_N


def _parse_named_asdxbn(path):
    records = {}
    names = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            animal, sire, dam, sex, bdate, name = line.rstrip("\n").split(",", 5)
            records[str(int(animal))] = (sire, dam, sex, bdate)
            names.append(name)
    return records, names


def test_named_griffon_matches_scientific_genealogy_and_adds_names():
    scientific = _parse_asdxb(canonical_griffon_path())
    named, names = _parse_named_asdxbn(named_griffon_path())
    assert len(scientific) == CANONICAL_N
    assert len(named) == CANONICAL_N
    assert set(scientific) == set(named)
    for animal_id, fields in scientific.items():
        assert named[animal_id] == fields
    nonempty = [name for name in names if name.strip()]
    assert len(nonempty) == CANONICAL_N
    unique = set(nonempty)
    assert len(unique) == 97000
    counts: dict[str, int] = {}
    for name in nonempty:
        counts[name] = counts.get(name, 0) + 1
    duplicated = {name for name, count in counts.items() if count > 1}
    assert duplicated == {"Colette", "Stella"}
    assert max(counts.values()) == 2
    assert named["98685"] == scientific["98685"]
    assert named["98667"] == scientific["98667"]
    with open(named_griffon_path(), "rb") as handle:
        assert hashlib.sha256(handle.read()).hexdigest() == NAMED_SHA256
    with open(named_griffon_path(), encoding="utf-8") as handle:
        by_id = {}
        for line in handle:
            if not line.strip():
                continue
            animal, _sire, _dam, _sex, _bdate, name = line.rstrip("\n").split(",", 5)
            by_id[str(int(animal))] = name
    assert by_id["98685"] == "Hierners Heartbreaker"
    assert by_id["98667"] == "Morning Bell Virgine"


@pytest.mark.integration
def test_canonical_griffon_dataset_regression_metrics():
    """Dataset regressions on the 2026 export. Not scientific constants."""
    ped = load_canonical_griffon()
    assert len(ped.pedigree) == CANONICAL_N
    assert pyp_utils.set_generation(ped)
    igens = [animal.igen for animal in ped.pedigree]
    assert min(igens) == 1
    assert max(igens) == 70
    ng = pyp_metrics.effective_founder_genomes(
        ped, rounds=3, seed=31, chrometype="autosome", output=False, quiet=True
    )
    assert ng == 12.421689363554368
    lacy = pyp_metrics.effective_founders_lacy(ped)
    assert lacy["fa_effective_founders"] == 193.46506304667966
    assert lacy["fa_animal_count"] == CANONICAL_N
    assert lacy["fa_founder_count"] == 7574
    assert lacy["fa_descendant_count"] == 90343
    ne = pyp_metrics.theoretical_ne_from_metadata(ped, output=False)
    assert ne == 28538.554380711772
    result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
    fx = result["fx"]
    values = list(fx.values())
    assert len(values) == CANONICAL_N
    assert sum(1 for value in values if value == 0.0) == 10864
    assert sum(1 for value in values if value > 0.0) == 83514
    assert min(values) == -5.551115123125783e-15
    assert max(values) == 0.546875
    mean = sum(values) / len(values)
    assert abs(mean - 0.09328960441457022) < 1e-12
    positional = [fx[i] for i in range(1, CANONICAL_N + 1)]
    digest = hashlib.sha256()
    digest.update(struct.pack(">Q", len(positional)))
    for value in positional:
        digest.update(struct.pack(">d", value))
    assert digest.hexdigest() == (
        "68f5b89fef8b15e051417bd996201f9db3b0a4d11122f10c10ea716d18fc2a35"
    )
    by_oid = {int(animal.originalID): fx[animal.animalID] for animal in ped.pedigree}
    assert 37481 not in by_oid and 54586 not in by_oid
    assert by_oid[37482] == 0.3060716643589345
    assert by_oid[54587] == 0.07958395780637417
    assert by_oid[37476] == 0.1450196736803535
    assert by_oid[98685] == 0.08597606312833217
    assert by_oid[98667] == 0.08683932119141691
