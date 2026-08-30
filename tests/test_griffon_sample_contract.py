"""Single checked-in Griffon Bruxellois pedigree dataset.

The repository ships exactly one Griffon pedigree file:

    PyPedal/examples/griffonbruxellois_2026_pyp.ped

Scientific tests that historically used smaller extracts derive those
extracts into a temporary directory from that file. Those subsets are not
committed.
"""
import hashlib
import os
import tarfile

import pytest

from PyPedal import pyp_chronology, pyp_metrics, pyp_nrm, pyp_utils
from _pedhelpers import (
    CANONICAL_GRIFFON_PED,
    GRIFFON_1871_1890_IDS,
    GRIFFON_TEST_SMALL_IDS,
    REPO,
    canonical_griffon_path,
    load_canonical_griffon,
    load_griffon_1871_1890,
    write_canonical_griffon_subset,
)


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


def test_repository_has_exactly_one_griffon_pedigree_dataset():
    hits = _griffon_data_files(REPO)
    assert hits == [os.path.join("PyPedal", "examples", CANONICAL_GRIFFON_PED)], hits
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
    assert len(derived) == 166
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
    assert len(ped.pedigree) == 167
    assert ped.metadata.num_implicit_parents == 1
    assert pyp_utils.set_generation(ped)
    distribution = {}
    for animal in ped.pedigree:
        distribution[animal.igen] = distribution.get(animal.igen, 0) + 1
    assert distribution == {1: 124, 2: 20, 3: 8, 4: 12, 5: 3}


def test_sdist_contains_exactly_one_griffon_pedigree_dataset(tmp_path):
    from test_product_surface import _build_setuptools_artifact

    sdist = _build_setuptools_artifact(tmp_path / "dist", "sdist")
    with tarfile.open(sdist) as archive:
        names = archive.getnames()
    griffon_peds = [
        name for name in names
        if "griffon" in name.lower() and name.endswith(".ped")
    ]
    assert len(griffon_peds) == 1, griffon_peds
    assert griffon_peds[0].endswith("griffonbruxellois_2026_pyp.ped")
    leftover = [
        name for name in names
        if "griffon" in name.lower()
        and name.lower().endswith((".csv", ".dat", ".ini", ".txt"))
    ]
    assert leftover == [], leftover


CANONICAL_SHA256 = "f288e1ab00eb710e8cfdb2df6175e8513c4d15df57ea16220d8397bc91389443"


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
    assert records == 98001


@pytest.mark.integration
def test_canonical_griffon_dataset_regression_metrics():
    """Dataset regressions on the 2026 export. Not scientific constants."""
    ped = load_canonical_griffon()
    assert len(ped.pedigree) == 98001
    assert pyp_utils.set_generation(ped)
    igens = [animal.igen for animal in ped.pedigree]
    assert min(igens) == 1
    assert max(igens) == 70
    ng = pyp_metrics.effective_founder_genomes(
        ped, rounds=3, seed=31, chrometype="autosome", output=False, quiet=True
    )
    assert ng == 11.018378975785259
    lacy = pyp_metrics.effective_founders_lacy(ped)
    assert lacy["fa_effective_founders"] == 193.31434658869796
    result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
    fx = result["fx"]
    values = list(fx.values())
    assert len(values) == 98001
    assert sum(1 for value in values if value > 0.0) == 84442
    assert max(values) == 0.546875
    mean = sum(values) / len(values)
    assert abs(mean - 0.09313044278029989) < 1e-12
