"""Optional progress callbacks: science-neutral, exceptions propagate."""

import pytest
from _pedhelpers import chdir_tmp, load_corpus

from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_errors import PyPedalError
from PyPedal.pyp_newclasses import NewPedigree, load_pedigree
from PyPedal.pyp_results import ProgressCallback


class Recorder:
    def __init__(self):
        self.events = []

    def __call__(self, done, total):
        self.events.append((done, total))


def _assert_monotonic_known_total(events, total):
    assert events, "expected at least a final progress event"
    dones = [done for done, _total in events]
    assert dones == sorted(dones)
    assert dones[0] >= 1
    assert all(done <= total for done, _total in events)
    assert all(_total == total for _done, _total in events)
    assert events[-1] == (total, total)
    assert dones.count(total) == 1


def test_progress_callback_type_is_public():
    assert callable(ProgressCallback) or ProgressCallback is not None

    def sample(done: int, total: int | None) -> None:
        return None

    typed: ProgressCallback = sample
    typed(1, 1)


def test_meuwissen_luo_progress_none_matches_recording_callback():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        none = pyp_nrm.inbreeding_meuwissen_luo(ped)
        recorder = Recorder()
        recorded = pyp_nrm.inbreeding_meuwissen_luo(ped, progress=recorder)
        via_dispatch = pyp_nrm.inbreeding(ped, method="meu_luo", output=False, progress=Recorder())
        assert none == recorded
        assert {int(k): float(v) for k, v in via_dispatch["fx"].items()} == {
            int(k): float(v) for k, v in none.items()
        }
        assert abs(none[5] - 0.125) < 1e-12
        n = len(ped.pedigree)
        _assert_monotonic_known_total(recorder.events, n)


def test_modified_meuwissen_luo_progress_none_matches_recording_callback():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        none = pyp_nrm.inbreeding_modified_meuwissen_luo(ped)
        recorder = Recorder()
        recorded = pyp_nrm.inbreeding_modified_meuwissen_luo(ped, progress=recorder)
        assert none == recorded
        assert abs(none[5] - 0.125) < 1e-12
        _assert_monotonic_known_total(recorder.events, len(ped.pedigree))


def test_meuwissen_luo_callback_exception_propagates_and_stops():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        seen = []

        def boom(done, total):
            seen.append((done, total))
            raise RuntimeError("abort inbreeding")

        with pytest.raises(RuntimeError, match="abort inbreeding") as caught:
            pyp_nrm.inbreeding_meuwissen_luo(ped, progress=boom)
        assert not isinstance(caught.value, PyPedalError)
        assert seen == [(1, len(ped.pedigree))]


def test_dropped_ancestral_inbreeding_same_seed_with_and_without_callback():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        none = pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=5, loci=8, seed=31)
        recorder = Recorder()
        recorded = pyp_metrics.dropped_ancestral_inbreeding(
            ped, rounds=5, loci=8, seed=31, progress=recorder
        )
        assert none == recorded
        _assert_monotonic_known_total(recorder.events, 5)


def test_effective_founder_genomes_same_seed_with_and_without_callback():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        none = pyp_metrics.effective_founder_genomes(
            ped, rounds=4, seed=31, output=False, quiet=True
        )
        recorder = Recorder()
        recorded = pyp_metrics.effective_founder_genomes(
            ped, rounds=4, seed=31, output=False, quiet=True, progress=recorder
        )
        assert none == recorded
        _assert_monotonic_known_total(recorder.events, 4)


def test_gene_drop_callback_exception_propagates():
    with chdir_tmp():
        ped = load_corpus("mrode.ped")
        seen = []

        def boom(done, total):
            seen.append((done, total))
            raise ValueError("stop drop")

        with pytest.raises(ValueError, match="stop drop") as caught:
            pyp_metrics.dropped_ancestral_inbreeding(ped, rounds=4, loci=4, seed=7, progress=boom)
        assert not isinstance(caught.value, PyPedalError)
        assert seen == [(1, 4)]


def test_boichard_progress_is_monotonic_with_unknown_total():
    with chdir_tmp():
        ped = load_corpus("boichard2a.ped")
        none = pyp_metrics.a_effective_ancestors_definite(ped, output=False)
        recorder = Recorder()
        recorded = pyp_metrics.a_effective_ancestors_definite(ped, output=False, progress=recorder)
        assert none == recorded
        assert recorder.events
        dones = [done for done, total in recorder.events]
        assert dones == sorted(dones)
        assert dones[0] >= 1
        assert all(total is None for _done, total in recorder.events)
        assert dones[-1] == len(recorder.events)


def test_boichard_indefinite_progress_matches_none():
    with chdir_tmp():
        ped = load_corpus("boichard2a.ped")
        none = pyp_metrics.a_effective_ancestors_indefinite(ped, n=2, output=False)
        recorder = Recorder()
        recorded = pyp_metrics.a_effective_ancestors_indefinite(
            ped, n=2, output=False, progress=recorder
        )
        assert none == recorded
        assert recorder.events
        assert all(total is None for _done, total in recorder.events)


def test_preprocess_file_progress_has_unknown_total(tmp_path):
    pedfile = tmp_path / "tiny.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    recorder = Recorder()
    with chdir_tmp():
        ped = NewPedigree(
            {
                "pedfile": str(pedfile),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
                "renumber": False,
            }
        )
        ped.preprocess(progress=recorder)
    assert recorder.events
    assert all(total is None for _done, total in recorder.events)
    dones = [done for done, _total in recorder.events]
    assert dones == sorted(dones)
    assert dones[0] >= 1
    assert dones[-1] == 3


def test_preprocess_textstream_progress_has_known_total():
    recorder = Recorder()
    with chdir_tmp():
        ped = NewPedigree(
            {
                "pedfile": "textstream.ped",
                "pedformat": "ASD",
                "sepchar": ",",
                "messages": "quiet",
                "pedigree_summary": 0,
                "renumber": False,
            }
        )
        ped.preprocess(textstream="1,0,0\n2,0,0\n3,1,2\n", progress=recorder)
    _assert_monotonic_known_total(recorder.events, 3)


def test_preprocess_dbstream_progress_has_known_total():
    recorder = Recorder()
    rows = [("1", "0", "0"), ("2", "0", "0")]
    with chdir_tmp():
        ped = NewPedigree(
            {
                "pedfile": "db.ped",
                "pedformat": "asd",
                "sepchar": ",",
                "messages": "quiet",
                "pedigree_summary": 0,
                "renumber": False,
            }
        )
        ped.preprocess(dbstream=rows, progress=recorder)
    _assert_monotonic_known_total(recorder.events, 2)


def test_load_callback_exception_is_not_translated(tmp_path):
    pedfile = tmp_path / "tiny.ped"
    pedfile.write_text("1 0 0\n2 0 0\n", encoding="utf-8")
    seen = []

    def boom(done, total):
        seen.append((done, total))
        raise RuntimeError("abort load")

    with pytest.raises(RuntimeError, match="abort load") as caught:
        load_pedigree(
            options={
                "pedfile": str(pedfile),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
            },
            progress=boom,
        )
    assert not isinstance(caught.value, PyPedalError)
    assert seen
    assert seen[0][0] >= 1
