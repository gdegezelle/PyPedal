"""Canonical Griffon scale: 98,001 rows without duplicating the pedigree."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from _pedhelpers import canonical_griffon_path, close_owned_pypedal_log_handlers

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from PyPedal.application import (
    PedigreeOpenOptions,
    PedigreeSession,
    PedigreeTableSource,
    load_into_session,
)
from PyPedal.desktop.models.pedigree_table import PedigreeFilterProxy, PedigreeTableModel

if QApplication.instance() is None:
    QApplication(["pypedal-desktop-tests"])

EXPECTED_N = 98_001


@pytest.mark.integration
def test_griffon_table_model_addresses_all_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    session = PedigreeSession()
    local = tmp_path / "griffonbruxellois_2026_pyp.ped"
    shutil.copy(canonical_griffon_path(), local)
    monkeypatch.chdir(tmp_path)
    try:
        pedigree = load_into_session(
            session,
            local,
            PedigreeOpenOptions(pedformat="asdxb", separator=",", renumber=True),
        )
        assert len(pedigree.pedigree) == EXPECTED_N
        table = PedigreeTableSource(pedigree)
        model = PedigreeTableModel(table)
        assert model.rowCount() == EXPECTED_N
        assert model.columnCount() == 8
        sample = model.data(model.index(0, 0), Qt.ItemDataRole.DisplayRole)
        assert sample not in {None, ""}
        mid = EXPECTED_N // 2
        assert model.data(model.index(mid, 1), Qt.ItemDataRole.UserRole) is not None
        last = model.data(model.index(EXPECTED_N - 1, 0), Qt.ItemDataRole.DisplayRole)
        assert last not in {None, ""}
        proxy = PedigreeFilterProxy()
        proxy.setSourceModel(model)
        first_id = str(model.data(model.index(0, 0), Qt.ItemDataRole.UserRole))
        proxy.set_query(first_id)
        assert proxy.rowCount() >= 1
        proxy.set_query("zzzz-not-a-griffon-name")
        assert proxy.rowCount() == 0
        proxy.set_query("")
        assert proxy.rowCount() == EXPECTED_N
    finally:
        close_owned_pypedal_log_handlers()
