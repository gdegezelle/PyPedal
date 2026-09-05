"""Selected-pair additive relationship via shared Meuwissen-Luo core."""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import time
import unittest

import pytest

from _pedhelpers import (
    chdir_tmp,
    load_corpus,
    load_corpus_from_path,
    write_temp_pedigree,
)
from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_newclasses import NewAMatrix


def _load_asd(rows, **overrides):
    path = write_temp_pedigree(rows)
    try:
        return load_corpus_from_path(path, "asd", **overrides)
    finally:
        os.remove(path)


def _dense_additive_matrix(ped):
    opts = dict(ped.kw)
    opts["messages"] = "quiet"
    return pyp_nrm.fast_a_matrix(ped.pedigree, opts, method="dense")


def _assert_pairwise_matches_dense(ped, places=12):
    matrix = _dense_additive_matrix(ped)
    n = len(ped.pedigree)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            want = float(matrix[i - 1, j - 1])
            got = float(pyp_nrm._pairwise_additive_relationship(ped, i, j))
            assert got == pytest.approx(want, abs=10 ** (-places)), (i, j, got, want)
            if i == j:
                assert pyp_metrics.relationship(i, j, ped) == 1.0
            else:
                pub = float(pyp_metrics.relationship(i, j, ped))
                assert pub == pytest.approx(want, abs=10 ** (-places)), (i, j, pub, want)


def _fibonacci_chain_rows(n):
    rows = ["1 0 0", "2 0 0"]
    for i in range(3, n + 1):
        rows.append("%d %d %d" % (i, i - 1, i - 2))
    return rows


def _ped_snapshot(ped):
    return {
        "n": len(ped.pedigree),
        "animal_ids": [int(a.animalID) for a in ped.pedigree],
        "sire_ids": [a.sireID for a in ped.pedigree],
        "dam_ids": [a.damID for a in ped.pedigree],
        "fa": [a.fa for a in ped.pedigree],
        "nrm": getattr(ped, "nrm", None),
        "meta": ped.metadata.num_records,
        "renumbered": ped.kw.get("pedigree_is_renumbered"),
    }


class TestMrodeScientificControls(unittest.TestCase):
    def test_unrelated_founders(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertEqual(0.0, pyp_metrics.relationship(1, 2, ped))

    def test_parent_offspring(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertAlmostEqual(0.5, pyp_metrics.relationship(1, 3, ped), places=12)

    def test_half_sibs(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertAlmostEqual(0.25, pyp_metrics.relationship(4, 3, ped), places=12)

    def test_full_sibs(self):
        ped = _load_asd(["1 0 0", "2 0 0", "3 1 2", "4 1 2"])
        with chdir_tmp():
            self.assertAlmostEqual(0.5, pyp_metrics.relationship(3, 4, ped), places=12)

    def test_inbred_parent_pair(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertAlmostEqual(0.625, pyp_metrics.relationship(3, 5, ped), places=12)

    def test_half_founder_unrelated_side(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertEqual(0.0, pyp_metrics.relationship(2, 4, ped))

    def test_all_mrode_pairs_match_dense_nrm(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            _assert_pairwise_matches_dense(ped)


class TestSelfCompatibility(unittest.TestCase):
    def test_public_self_shortcut_and_true_diagonal(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            self.assertEqual(1.0, pyp_metrics.relationship(5, 5, ped))
            true = pyp_nrm._pairwise_additive_relationship(ped, 5, 5)
            self.assertAlmostEqual(1.125, true, places=12)
            self.assertAlmostEqual(0.5625, pyp_metrics.mating_coi(5, 5, ped), places=12)


class TestFibonacciExplosionRegression(unittest.TestCase):
    def test_chain_n40_completes_quickly(self):
        ped = _load_asd(_fibonacci_chain_rows(40))
        with chdir_tmp():
            t0 = time.perf_counter()
            got = pyp_metrics.relationship(40, 39, ped)
            elapsed = time.perf_counter() - t0
        self.assertLess(elapsed, 5.0)
        matrix = _dense_additive_matrix(ped)
        self.assertAlmostEqual(float(matrix[39, 38]), got, places=10)


class TestStructuralGuards(unittest.TestCase):
    def test_readonly_repeated_calls(self):
        ped = load_corpus("mrode.ped")
        before = _ped_snapshot(ped)
        with chdir_tmp():
            for _ in range(3):
                pyp_metrics.relationship(4, 3, ped)
                pyp_metrics.mating_coi(4, 3, ped)
        after = _ped_snapshot(ped)
        self.assertEqual(before, after)


def test_relationship_avoids_matrix_and_recurse(monkeypatch):
    ped = load_corpus("mrode.ped")

    def _fail(*_args, **_kwargs):
        raise AssertionError("expensive pedigree matrix path must not run")

    monkeypatch.setattr(pyp_nrm, "fast_a_matrix", _fail)
    monkeypatch.setattr(pyp_nrm, "fast_a_matrix_r", _fail)
    monkeypatch.setattr(pyp_nrm, "recurse_pedigree", _fail)
    with chdir_tmp():
        assert pyp_metrics.relationship(4, 3, ped) == pytest.approx(0.25)


class TestRandomSmallDags(unittest.TestCase):
    def test_random_pedigrees_match_dense_nrm(self):
        rng = random.Random(20260303)
        for _ in range(12):
            rows = ["1 0 0", "2 0 0", "3 0 0"]
            next_id = 4
            known = {1, 2, 3}
            for _ in range(rng.randint(3, 10)):
                sire = rng.choice(sorted(known))
                dam = rng.choice(sorted(known))
                rows.append("%d %d %d" % (next_id, sire, dam))
                known.add(next_id)
                next_id += 1
            ped = _load_asd(rows)
            with chdir_tmp():
                _assert_pairwise_matches_dense(ped, places=10)


class TestAttachedNrmFastPath(unittest.TestCase):
    def test_form_nrm_lookup_still_used(self):
        ped = load_corpus("mrode.ped")
        ped.nrm = NewAMatrix(dict(ped.kw))
        ped.nrm.kw["messages"] = "quiet"
        ped.kw["form_nrm"] = True
        ped.nrm.form_a_matrix(ped.pedigree)
        with chdir_tmp():
            self.assertAlmostEqual(0.25, pyp_metrics.relationship(4, 3, ped), places=12)


@pytest.mark.integration
def test_griffon_observed_pair_completes_under_timeout(tmp_path):
    """Canonical Griffon pair originalID 98685 x 98667 -> current 98001 x 97984."""
    import shutil

    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    tmp = str(tmp_path)
    ped_src = os.path.join(repo, "PyPedal", "examples", "griffonbruxellois_2026_pyp.ped")
    ped_copy = os.path.join(tmp, "griffon.ped")
    shutil.copy2(ped_src, ped_copy)
    worker = """
import json, os, sys, time
sys.path.insert(0, sys.argv[1])
os.chdir(sys.argv[2])
from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_metrics
ped = load_pedigree({
    "pedfile": "griffon.ped",
    "pedformat": "asdxb",
    "sepchar": ",",
    "messages": "quiet",
    "pedigree_summary": 0,
    "renumber": True,
    "form_nrm": False,
})
t0 = time.perf_counter()
val = pyp_metrics.relationship(98001, 97984, ped)
print(json.dumps({"ok": True, "value": float(val), "seconds": time.perf_counter() - t0}))
"""
    script = os.path.join(tmp, "worker.py")
    with open(script, "w", encoding="utf-8") as handle:
        handle.write(worker)
    proc = subprocess.run(
        [sys.executable, script, repo, tmp],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stderr[-500:]
    payload = json.loads(proc.stdout.strip().splitlines()[-1])
    assert payload["ok"] is True
    assert isinstance(payload["value"], float)
    assert payload["seconds"] < 30.0
