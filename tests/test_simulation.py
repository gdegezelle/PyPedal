"""
an earlier revision workstream B: simulated pedigrees must complete the normal load lifecycle.

The Python 3 port of NewPedigree.load() returned immediately after
simulate(), skipping reorder/renumber/PedigreeMetadata. Both 2.0.4 and
pypedal3 fall through. These tests pin the restored contract without
asserting incidental unseeded individual IDs.
"""
import numpy as np

from PyPedal import pyp_nrm
from PyPedal.pyp_newclasses import PedigreeMetadata, NewPedigree, load_pedigree

from _pedhelpers import chdir_tmp, load_corpus

SEED = 42


def _sim_opts(**overrides):
    opts = {
        "simulate_pedigree": True,
        "simulate_n": 20,
        "simulate_g": 3,
        "simulate_ns": 4,
        "simulate_nd": 4,
        "simulate_ir": 0.0,
        "simulate_sr": 0.5,
        "renumber": True,
        "messages": "quiet",
        "pedigree_summary": 0,
        "pedformat": "asdxg",
        "pedigree_save": False,
    }
    opts.update(overrides)
    return opts


def _load_simulated(seed=SEED, **overrides):
    np.random.seed(seed)
    with chdir_tmp():
        return load_pedigree(_sim_opts(**overrides))


def _raw_simulated_graph(seed=SEED, **overrides):
    """Animal records as simulate() emits them, before load finalization."""
    np.random.seed(seed)
    with chdir_tmp():
        ped = NewPedigree(_sim_opts(**overrides))
        ped.simulate()
        return [
            (int(a.animalID), int(a.sireID), int(a.damID), a.sex, int(a.gen))
            for a in ped.pedigree
        ]


def _finalized_original_graph(ped):
    missing = str(ped.kw["missing_parent"])
    rows = []
    for animal in ped.pedigree:
        sire = (
            int(ped.kw["missing_parent"])
            if str(animal.sireID) == missing
            else int(ped.backmap[animal.sireID])
        )
        dam = (
            int(ped.kw["missing_parent"])
            if str(animal.damID) == missing
            else int(ped.backmap[animal.damID])
        )
        rows.append(
            (int(animal.originalID), sire, dam, animal.sex, int(animal.gen))
        )
    return rows


def test_sim1_supported_simulation_load_completes():
    ped = _load_simulated()
    assert len(ped.pedigree) > 0


def test_sim2_metadata_is_pedigree_metadata_not_empty_dict():
    ped = _load_simulated()
    assert isinstance(ped.metadata, PedigreeMetadata)
    assert ped.metadata is not {}
    assert ped.metadata.num_records == len(ped.pedigree)
    assert ped.metadata.unique_gen_list
    assert len(ped.metadata.unique_gen_list) >= 1


def test_sim3_renumber_true_assigns_sequential_animal_ids():
    ped = _load_simulated(renumber=True)
    n = len(ped.pedigree)
    assert [a.animalID for a in ped.pedigree] == list(range(1, n + 1))
    assert ped.kw["pedigree_is_renumbered"] is True


def test_sim4_list_index_invariant():
    ped = _load_simulated(renumber=True)
    for i, animal in enumerate(ped.pedigree):
        assert animal.animalID == i + 1


def test_sim5_idmap_backmap_agree():
    ped = _load_simulated(renumber=True)
    assert ped.idmap
    assert ped.backmap
    for original, current in ped.idmap.items():
        assert ped.backmap[current] == original
    for current, original in ped.backmap.items():
        assert ped.idmap[original] == current


def test_sim6_parent_references_resolve():
    ped = _load_simulated()
    missing = str(ped.kw["missing_parent"])
    present = {a.animalID for a in ped.pedigree}
    unresolved = []
    for animal in ped.pedigree:
        for parent in (animal.sireID, animal.damID):
            if str(parent) != missing and parent not in present:
                unresolved.append((animal.animalID, parent))
    assert unresolved == []


def test_sim7_parent_before_child():
    ped = _load_simulated(renumber=True)
    missing = str(ped.kw["missing_parent"])
    index = {a.animalID: i for i, a in enumerate(ped.pedigree)}
    for animal in ped.pedigree:
        child_i = index[animal.animalID]
        for parent in (animal.sireID, animal.damID):
            if str(parent) != missing:
                assert index[parent] < child_i


def test_sim8_downstream_inbreeding_vanraden_and_tabular():
    ped = _load_simulated()
    van = pyp_nrm.inbreeding(ped, method="vanraden", output=False)
    tab = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    assert len(van["fx"]) == len(ped.pedigree)
    assert len(tab["fx"]) == len(ped.pedigree)
    for fx in (van["fx"], tab["fx"]):
        for value in fx.values():
            assert 0.0 <= float(value) <= 1.0


def test_sim11_renumber_false_still_builds_metadata_without_forcing_ids():
    ped = _load_simulated(renumber=False)
    assert isinstance(ped.metadata, PedigreeMetadata)
    assert ped.kw["pedigree_is_renumbered"] is not True
    ids = [a.animalID for a in ped.pedigree]
    assert ids != list(range(1, len(ids) + 1))
    assert len(set(ids)) == len(ids)


def test_sim12_repeated_simulation_does_not_share_state():
    first = _load_simulated(seed=1)
    second = _load_simulated(seed=2)
    assert first is not second
    assert first.pedigree is not second.pedigree
    first.pedigree[0].animalID = 999999
    assert second.pedigree[0].animalID != 999999


def test_sim_founder_flags_match_missing_parents():
    ped = _load_simulated()
    missing = str(ped.kw["missing_parent"])
    for animal in ped.pedigree:
        both_missing = (
            str(animal.sireID) == missing and str(animal.damID) == missing
        )
        if both_missing:
            assert animal.founder == "y"
        else:
            assert animal.founder == "n"


def test_sim_raw_generated_graph_is_preserved_through_finalization():
    raw = _raw_simulated_graph()
    ped = _load_simulated()
    finalized = _finalized_original_graph(ped)
    assert sorted(finalized) == sorted(raw)
    assert [row[0] for row in finalized] != [a.animalID for a in ped.pedigree]


def test_sim_ordinary_file_load_still_builds_metadata():
    ped = load_corpus("mrode.ped")
    assert isinstance(ped.metadata, PedigreeMetadata)
    result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    fx = {int(k): float(v) for k, v in result["fx"].items()}
    assert abs(fx[5] - 0.125) < 1e-12
