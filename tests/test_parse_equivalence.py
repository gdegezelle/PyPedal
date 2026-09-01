"""Pin pedigree parse/load results before the 4.1-C helper extraction.

These snapshots cover preprocess() and a default load(). They do not run
inbreeding, NRM formation, or other scientific calculations.
"""
import hashlib
import json
import logging
import os
import shutil
import tempfile

import pytest

from _pedhelpers import (
    CORPUS,
    chdir_tmp,
    load_canonical_griffon,
    load_corpus,
    load_corpus_from_path,
)
from PyPedal.pyp_errors import (
    PyPedalConfigurationError,
    PyPedalPedigreeFormatError,
)
from PyPedal.pyp_newclasses import NewPedigree

CORPUS_FORMATS = (
    ("mrode.ped", "asd", " "),
    ("new_lacy.ped", "asd", " "),
    ("boichard2a.ped", "asdg", " "),
    ("doug.ped", "ASDx", " "),
    ("new_ids.ped", "ASD", " "),
    ("horse.ped", "ASD", ","),
    ("generations.ped", "asdbx", " "),
    ("hartlandclark.ped", "asdb", " "),
)

PREPROCESS_DIGESTS = {
    "mrode.ped": "8cb85e23c954f9108161ae626bb73bc9a81b8278a25af1c0017621cb8f1de815",
    "new_lacy.ped": "1351d9bfe5aca40fe106ec6f5735eb64ba89e232c01de2c684ebe15cf89fe0e3",
    "boichard2a.ped": "f660a588f36a7d28846563930e0a89415e626a5c976d5f05b24dc12bab3be7b7",
    "doug.ped": "9f5958d687cfe7dea9674879ccfa9096da76e2b824872ef61b008b0374a01bcf",
    "new_ids.ped": "56deaf19ea4c207bc0beb0ee07484d2684de2d43124417c4865270462bce788d",
    "horse.ped": "ea5d9ebd3e56fad385a6ea7af96fdd3f8ccdf96861fbd56425b1ef72cc4fa8ec",
    "generations.ped": "458e6ec0a7c8cdd5c09bb06de02f75fd331083af3f93694c74c960df23ae088f",
    "hartlandclark.ped": "b03335236dfd151b456db9c40018a2c33eacbf9b82264875b576bdee419708fd",
}

LOAD_DIGESTS = {
    "mrode.ped": "a1b482bed4804af4c3506d28254dc15cb5b867aa3eff759b04b5d5803a71cd3d",
    "new_lacy.ped": "95e6496136331d678d3f3c8939b86f7dfd5a05b386967101eaae5a9a90864869",
    "boichard2a.ped": "00cd7d443456a7d1378aaed99c6546391da33c006f4a95a29b1a63d5af8ccfe5",
    "doug.ped": "12694463573067a3e2fd4f152e037614260408ce890da902b7ceaf332d561b23",
    "new_ids.ped": "f1c3aa5636a508efaf601efed07660adb39a9299df9293fbcf5ff4647dc7dbc2",
    "horse.ped": "380df3eb7138ac4024b54a412b295cea2c8c802eb75dcc5c4cdc497b9a69393f",
    "generations.ped": "a9dd74538117a8b94ee49b855abc33d1fcb31214697139b324df2297509d6c26",
    "hartlandclark.ped": "4fbcdd90c9e6b860fccda0c310170b53b374f0e9c29a4bffa39cdea4882f88ef",
}

GRIFFON_ANIMAL_DIGEST = (
    "528cdf79f12d137e8c9a193701743ed3a62c813483ee9e13870cac7212dc3315"
)
GRIFFON_IDMAP_DIGEST = (
    "9f9c369e2771fde1fa99ae1bdb06ded2ed66d4cfc68086ce3c7bf6b6c6131234"
)
GRIFFON_BACKMAP_DIGEST = (
    "9039d8ea926e1c26cf02680a44e0ac79d6b9e3a137a61d9c636daf2bdcad993d"
)

MRODE_PREPROCESS_ANIMALS = (
    (1, 1, 0, 0, "u", "1", None, None, 0.0, -999.0, -999.0, "y",
     "Unknown_Name", "Unknown_Name"),
    (2, 2, 0, 0, "u", "2", None, None, 0.0, -999.0, -999.0, "y",
     "Unknown_Name", "Unknown_Name"),
    (3, 3, 1, 2, "u", "3", None, None, 0.0, -999.0, -999.0, "n", "1", "2"),
    (4, 4, 1, 0, "u", "4", None, None, 0.0, -999.0, -999.0, "n",
     "1", "Unknown_Name"),
    (5, 5, 4, 3, "u", "5", None, None, 0.0, -999.0, -999.0, "n", "4", "3"),
    (6, 6, 5, 2, "u", "6", None, None, 0.0, -999.0, -999.0, "n", "5", "2"),
)


def animal_snapshot(animal):
    bd = animal.bd.isoformat() if getattr(animal, "bd", None) is not None else None
    return {
        "originalID": animal.originalID,
        "animalID": animal.animalID,
        "sireID": animal.sireID,
        "damID": animal.damID,
        "sex": animal.sex,
        "name": animal.name,
        "by": animal.by,
        "bd": bd,
        "fa": animal.fa,
        "gen": animal.gen,
        "igen": animal.igen,
        "founder": animal.founder,
        "sireName": animal.sireName,
        "damName": animal.damName,
    }


def _map_snapshot(mapping):
    return {str(key): mapping[key] for key in sorted(mapping, key=lambda item: str(item))}


def pedigree_snapshot(ped):
    return {
        "n": len(ped.pedigree),
        "animals": [animal_snapshot(animal) for animal in ped.pedigree],
        "idmap": _map_snapshot(ped.idmap),
        "backmap": _map_snapshot(ped.backmap),
        "namemap": _map_snapshot(ped.namemap),
        "namebackmap": _map_snapshot(ped.namebackmap),
        "implicit_parent_count": len(ped._implicit_parents),
        "implicit_parent_ids": list(ped._implicit_parents),
        "f_computed": ped.kw.get("f_computed"),
        "g_computed": ped.kw.get("g_computed"),
        "renumber": ped.kw.get("renumber"),
        "pedigree_is_renumbered": ped.kw.get("pedigree_is_renumbered"),
        "pedformat": ped.kw.get("pedformat"),
    }


def snapshot_digest(obj):
    blob = json.dumps(obj, sort_keys=True, default=str, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def preprocess_only(src, pedformat, sepchar=" "):
    tmp = tempfile.mkdtemp(prefix="pypedal_parse_")
    local = os.path.join(tmp, os.path.basename(src))
    shutil.copy(src, local)
    cwd = os.getcwd()
    os.chdir(tmp)
    try:
        ped = NewPedigree({
            "pedfile": local,
            "pedformat": pedformat,
            "sepchar": sepchar,
            "messages": "quiet",
            "renumber": False,
            "pedigree_summary": 0,
        })
        ped.preprocess()
        return ped
    finally:
        os.chdir(cwd)
        logging.disable(logging.CRITICAL)


def load_rows(rows, pedformat, sepchar=" ", **overrides):
    tmp = tempfile.mkdtemp(prefix="pypedal_parse_")
    path = os.path.join(tmp, "rows.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, sepchar, **overrides)


def _animal_tuple(animal):
    snap = animal_snapshot(animal)
    return (
        snap["originalID"], snap["animalID"], snap["sireID"], snap["damID"],
        snap["sex"], snap["name"], snap["by"], snap["bd"], snap["fa"],
        snap["gen"], snap["igen"], snap["founder"], snap["sireName"],
        snap["damName"],
    )


@pytest.mark.parametrize("name,pedformat,sepchar", CORPUS_FORMATS)
def test_preprocess_snapshot_digest(name, pedformat, sepchar):
    ped = preprocess_only(os.path.join(CORPUS, name), pedformat, sepchar)
    assert snapshot_digest(pedigree_snapshot(ped)) == PREPROCESS_DIGESTS[name]


@pytest.mark.parametrize("name,pedformat,sepchar", CORPUS_FORMATS)
def test_load_snapshot_digest(name, pedformat, sepchar):
    ped = load_corpus(name, pedformat, sepchar=sepchar)
    assert snapshot_digest(pedigree_snapshot(ped)) == LOAD_DIGESTS[name]


def test_mrode_preprocess_animal_fields():
    ped = preprocess_only(os.path.join(CORPUS, "mrode.ped"), "asd")
    assert len(ped.pedigree) == 6
    assert [_animal_tuple(animal) for animal in ped.pedigree] == list(
        MRODE_PREPROCESS_ANIMALS
    )
    assert ped.idmap == {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    assert ped.backmap == {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    assert ped._implicit_parents == []
    assert ped.kw["f_computed"] is False
    assert ped.kw["g_computed"] is False


def test_file_textstream_and_dbstream_agree_on_mrode():
    src = os.path.join(CORPUS, "mrode.ped")
    file_ped = preprocess_only(src, "asd")
    file_animals = [_animal_tuple(animal) for animal in file_ped.pedigree]
    body = "\n".join(
        line.strip() for line in open(src, encoding="utf-8")
        if line.strip() and not line.lstrip().startswith("#")
    ) + "\n"
    with chdir_tmp() as tmp:
        dummy = os.path.join(tmp, "stream.ped")
        open(dummy, "w", encoding="utf-8").write("")
        text_ped = NewPedigree({
            "pedfile": dummy,
            "pedformat": "asd",
            "messages": "quiet",
            "renumber": False,
            "pedigree_summary": 0,
        })
        text_ped.preprocess(textstream=body)
        db_ped = NewPedigree({
            "pedfile": dummy,
            "pedformat": "asd",
            "sepchar": ",",
            "messages": "quiet",
            "renumber": False,
            "pedigree_summary": 0,
        })
        db_rows = [tuple(line.split()) for line in body.strip().split("\n")]
        db_ped.preprocess(dbstream=db_rows)
    assert [_animal_tuple(animal) for animal in text_ped.pedigree] == file_animals
    assert [_animal_tuple(animal) for animal in db_ped.pedigree] == file_animals


def test_blank_line_stops_file_parse():
    ped = load_rows(["1 0 0", "", "2 0 0"], "asd", renumber=False)
    assert [animal.originalID for animal in ped.pedigree] == [1]


def test_textstream_without_trailing_newline_drops_last_record():
    with chdir_tmp() as tmp:
        dummy = os.path.join(tmp, "stream.ped")
        open(dummy, "w", encoding="utf-8").write("")
        ped = NewPedigree({
            "pedfile": dummy,
            "pedformat": "asd",
            "messages": "quiet",
            "renumber": False,
            "pedigree_summary": 0,
        })
        ped.preprocess(textstream="1 0 0\n2 0 0")
    assert [animal.originalID for animal in ped.pedigree] == [1]


def test_sire_only_implicit_parent():
    ped = load_rows(["2 1 0"], "asd", renumber=False)
    assert len(ped.pedigree) == 2
    assert ped.metadata.num_implicit_parents == 1
    by_oid = {animal.originalID: animal for animal in ped.pedigree}
    assert by_oid[1].founder == "y"
    assert by_oid[1].sireID == 0
    assert by_oid[1].damID == 0
    assert by_oid[1].sex == "u"
    assert by_oid[2].sireID == 1
    assert by_oid[2].damID == 0
    assert by_oid[2].founder == "n"


def test_dam_only_implicit_parent():
    ped = load_rows(["2 0 3"], "asd", renumber=False)
    assert len(ped.pedigree) == 2
    assert ped.metadata.num_implicit_parents == 1
    by_oid = {animal.originalID: animal for animal in ped.pedigree}
    assert by_oid[3].founder == "y"
    assert by_oid[3].sireID == 0
    assert by_oid[3].damID == 0
    assert by_oid[2].sireID == 0
    assert by_oid[2].damID == 3


def test_both_parents_implicit():
    ped = load_rows(["3 1 2"], "asd", renumber=False)
    assert len(ped.pedigree) == 3
    assert ped.metadata.num_implicit_parents == 2
    by_oid = {animal.originalID: animal for animal in ped.pedigree}
    assert by_oid[1].founder == "y"
    assert by_oid[2].founder == "y"
    assert by_oid[3].sireID == 1
    assert by_oid[3].damID == 2
    assert ped.metadata.implicit_parent_ids == [1, 2]


def test_shared_implicit_parent_is_materialized_once():
    ped = load_rows(["3 1 0", "4 1 0"], "asd", renumber=False)
    assert len(ped.pedigree) == 3
    assert ped.metadata.num_implicit_parents == 1
    assert ped.metadata.implicit_parent_ids == [1]
    offspring = [animal for animal in ped.pedigree if animal.originalID in {3, 4}]
    assert {animal.sireID for animal in offspring} == {1}


def test_string_id_sire_only_implicit_parent():
    ped = load_rows(["child,sireX,0"], "ASD", sepchar=",", renumber=False)
    assert len(ped.pedigree) == 2
    assert ped.metadata.num_implicit_parents == 1
    by_name = {animal.name: animal for animal in ped.pedigree}
    assert set(by_name) == {"child", "sireX"}
    assert by_name["sireX"].founder == "y"
    assert by_name["child"].sireID == by_name["sireX"].animalID
    assert by_name["child"].damID == 0


def test_string_id_shared_implicit_parent():
    ped = load_rows(["a,sx,0", "b,sx,0"], "ASD", sepchar=",", renumber=False)
    assert len(ped.pedigree) == 3
    assert ped.metadata.num_implicit_parents == 1
    by_name = {animal.name: animal for animal in ped.pedigree}
    assert by_name["sx"].founder == "y"
    assert by_name["a"].sireID == by_name["sx"].animalID
    assert by_name["b"].sireID == by_name["sx"].animalID


def test_skip_column_z_does_not_shift_ids():
    ped = load_rows(["1 extra 0 0"], "aZsd")
    assert len(ped.pedigree) == 1
    animal = ped.pedigree[0]
    assert animal.originalID == 1
    assert animal.sireID == 0
    assert animal.damID == 0


def test_loaded_fa_sets_f_computed():
    ped = load_rows(["1 0 0 0.125"], "asdf")
    assert ped.kw["f_computed"] is True
    assert ped.pedigree[0].fa == 0.125


def test_herd_code_h_is_overwritten_absent_without_H():
    ped = load_rows(["1 0 0 7"], "asdh")
    animal = ped.pedigree[0]
    assert animal.herd == "Unknown_Herd"
    assert animal.originalHerd == "Unknown_Herd"


def test_herd_code_H_populates_herd():
    ped = load_rows(["1 0 0 7"], "asdH")
    animal = ped.pedigree[0]
    assert animal.herd == 7
    assert animal.originalHerd == "7"


def test_case_sensitive_integer_versus_string_identity():
    integer_ped = load_rows(["10 0 0"], "asd", renumber=False)
    string_ped = load_rows(["10,0,0"], "ASD", sepchar=",", renumber=False)
    assert integer_ped.pedigree[0].originalID == 10
    assert string_ped.pedigree[0].name == "10"
    assert string_ped.pedigree[0].originalID != 10


def test_missing_animal_code_raises_configuration_error():
    with pytest.raises(PyPedalConfigurationError, match="required field") as raised:
        load_rows(["1 0 0"], "sd")
    assert "sd" in str(raised.value)


def test_column_count_mismatch_raises_format_error():
    with pytest.raises(PyPedalPedigreeFormatError, match="2 columns"):
        load_rows(["1 0"], "asd")


def test_canonical_griffon_load_counts_and_digest():
    ped = load_canonical_griffon()
    missing = ped.kw["missing_parent"]
    present = {animal.animalID for animal in ped.pedigree}
    dangling = [
        (animal.originalID, parent)
        for animal in ped.pedigree
        for parent in (animal.sireID, animal.damID)
        if str(parent) != str(missing) and parent not in present
    ]
    half = sum(
        1
        for animal in ped.pedigree
        if animal.founder == "n"
        and (
            (animal.sireID == missing and animal.damID != missing)
            or (animal.sireID != missing and animal.damID == missing)
        )
    )
    animals = [animal_snapshot(animal) for animal in ped.pedigree]
    animals_sorted = sorted(
        animals, key=lambda row: (str(row["originalID"]), str(row["animalID"]))
    )
    assert len(ped.pedigree) == 98001
    assert len({animal.originalID for animal in ped.pedigree}) == 98001
    assert dangling == []
    assert ped.metadata.num_implicit_parents == 0
    assert ped.metadata.num_unique_founders == 6689
    assert half == 915
    assert ped.metadata.num_unknown_birth_years == 3997
    assert snapshot_digest(animals_sorted) == GRIFFON_ANIMAL_DIGEST
    assert snapshot_digest(_map_snapshot(ped.idmap)) == GRIFFON_IDMAP_DIGEST
    assert snapshot_digest(_map_snapshot(ped.backmap)) == GRIFFON_BACKMAP_DIGEST
