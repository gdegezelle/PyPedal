"""Named Griffon lookup scale: 98,001 animals, no Qt widgets."""

from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

import pytest
from _pedhelpers import (
    close_owned_pypedal_log_handlers,
    named_griffon_path,
)

from PyPedal.application import (
    PedigreeOpenOptions,
    PedigreeSession,
    load_into_session,
    run_mating_coi,
    run_relationship,
)
from PyPedal.application.lookup import (
    DEFAULT_RESULT_LIMIT,
    AnimalLookupHit,
    AnimalLookupIndex,
)

EXPECTED_N = 98_001
A_EXPECTED = 0.20191301769610437
F_EXPECTED = 0.10095650884805218
OID_A = 98685
OID_B = 98667
CURRENT_A = 98001
CURRENT_B = 97984
NAME_A = "Hierners Heartbreaker"
NAME_B = "Morning Bell Virgine"


@pytest.mark.integration
def test_named_griffon_lookup_scale_and_benchmark_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Correctness and structural contracts for named-Griffon lookup.

    Wall-clock samples are printed for characterization. They are not a
    correctness gate: prefix search over 98k names is O(n) and hosted-runner
    load routinely exceeds a 250 ms microbenchmark.
    """
    session = PedigreeSession()
    # Logs and side-effect files follow pedfile; copy into pytest-owned tmp.
    local = tmp_path / "griffonbruxellois_2026_named_pyp.ped"
    shutil.copy(named_griffon_path(), local)
    monkeypatch.chdir(tmp_path)
    try:
        pedigree = load_into_session(
            session,
            local,
            PedigreeOpenOptions(pedformat="asdxbn", separator=",", renumber=True),
        )
        assert len(pedigree.pedigree) == EXPECTED_N
        index = session.animal_lookup
        assert index is not None
        assert isinstance(index, AnimalLookupIndex)
        assert len(index) == EXPECTED_N
        assert index.named_count == EXPECTED_N

        builds: list[int] = []

        def forbid_rebuild(cls: type[AnimalLookupIndex], pedigree_obj: object) -> AnimalLookupIndex:
            builds.append(1)
            raise AssertionError("AnimalLookupIndex rebuilt during lookup query")

        monkeypatch.setattr(AnimalLookupIndex, "from_pedigree", classmethod(forbid_rebuild))

        t0 = time.perf_counter()
        prefix = index.search("a")
        prefix_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        by_name_a = index.search(NAME_A)
        name_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        by_original = index.search(str(OID_A))
        original_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        by_current = index.search(str(CURRENT_A))
        current_s = time.perf_counter() - t0

        assert builds == []
        assert prefix.truncated is True
        assert prefix.total > DEFAULT_RESULT_LIMIT
        assert len(prefix.hits) == DEFAULT_RESULT_LIMIT
        assert all(isinstance(hit, AnimalLookupHit) for hit in prefix.hits)
        assert not hasattr(prefix.hits[0], "animalID")
        assert by_name_a.total >= 1
        assert by_name_a.hits[0].animal_id == CURRENT_A
        assert by_name_a.hits[0].name == NAME_A
        assert by_original.hits[0].animal_id == CURRENT_A
        assert by_current.hits[0].animal_id == CURRENT_A
        assert by_current.hits[0].original_id == OID_A

        by_name_b = index.search(NAME_B)
        assert by_name_b.hits[0].animal_id == CURRENT_B
        assert by_name_b.hits[0].original_id == OID_B

        colettes = index.search("Colette")
        assert colettes.hits[0].name == "Colette"
        assert colettes.hits[1].name == "Colette"
        assert {hit.original_id for hit in colettes.hits[:2]} == {20196, 20209}
        assert colettes.total >= 2
        assert session.animal_lookup is index

        rough_bytes = sys.getsizeof(index) + EXPECTED_N * 256
        assert rough_bytes < 80 * 1024 * 1024

        t0 = time.perf_counter()
        related = run_relationship(session, CURRENT_A, CURRENT_B)
        rel_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        mated = run_mating_coi(session, CURRENT_A, CURRENT_B)
        mate_s = time.perf_counter() - t0
        assert abs(related.coefficient - A_EXPECTED) < 1e-12
        assert abs(mated.coefficient - F_EXPECTED) < 1e-12
        print(
            "named-griffon-lookup "
            f"prefix={prefix_s:.4f}s/{prefix.total} "
            f"name={name_s:.4f}s original={original_s:.4f}s "
            f"current={current_s:.4f}s rel={rel_s:.4f}s mate={mate_s:.4f}s"
        )
    finally:
        close_owned_pypedal_log_handlers()
