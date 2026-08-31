"""Failed GUI loads must not leave a previous pedigree silently active."""
from types import SimpleNamespace

from PyPedal.pyp_app import apply_pedigree_load


def _session(pedigree=None, filename=""):
    return SimpleNamespace(pedigree=pedigree, filename=filename)


def test_successful_load_installs_the_new_pedigree():
    session = _session()
    pedigree = object()
    outcome = apply_pedigree_load(session, "/tmp/a.ped", pedigree=pedigree)
    assert outcome["ok"] is True
    assert session.pedigree is pedigree
    assert session.filename == "/tmp/a.ped"
    assert "a.ped" in outcome["status"]


def test_failed_load_keeps_previous_pedigree_and_says_so():
    previous = object()
    session = _session(pedigree=previous, filename="/tmp/a.ped")
    outcome = apply_pedigree_load(
        session,
        "/tmp/b.ped",
        error="PyPedalInputError\n\nbad file",
    )
    assert outcome["ok"] is False
    assert session.pedigree is previous
    assert session.filename == "/tmp/a.ped"
    assert "Load of b.ped failed" in outcome["output"]
    assert "a.ped remains the active pedigree" in outcome["output"]
    assert "bad file" in outcome["output"]
    assert "a.ped remains active" in outcome["status"]


def test_failed_load_with_no_previous_pedigree_stays_empty():
    session = _session()
    outcome = apply_pedigree_load(session, "/tmp/b.ped", error="boom")
    assert outcome["ok"] is False
    assert session.pedigree is None
    assert session.filename == ""
    assert "Load of b.ped failed" in outcome["output"]
    assert "No pedigree is active" in outcome["output"]
    assert "no pedigree is active" in outcome["status"]
