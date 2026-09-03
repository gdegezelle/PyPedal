"""macOS bundle helpers derive identity from project metadata."""

from __future__ import annotations

import importlib.util
import sys
import tomllib
from pathlib import Path

from _pedhelpers import REPO

MACOS = Path(REPO) / "tools" / "macos"


def _bundle_config():
    path = MACOS / "bundle_config.py"
    spec = importlib.util.spec_from_file_location("pypedal_bundle_config", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_version_comes_from_pyproject() -> None:
    config = _bundle_config()
    with (Path(REPO) / "pyproject.toml").open("rb") as handle:
        expected = tomllib.load(handle)["project"]["version"]
    assert config.read_project_version(Path(REPO)) == expected


def test_spec_does_not_hardcode_the_release_version() -> None:
    spec = (MACOS / "PyPedal.spec").read_text(encoding="utf-8")
    assert "read_project_version" in spec
    assert 'CFBundleShortVersionString": "4.' not in spec
    assert "4.1.0" not in spec


def test_submodule_filter_drops_ctk_and_tests() -> None:
    config = _bundle_config()
    assert config.keep_submodule("PyPedal.application.jobs") is True
    assert config.keep_submodule("PyPedal.desktop.pages.inbreeding") is True
    assert config.keep_submodule("PyPedal.pyp_app") is True
    assert config.keep_submodule("PyPedal.pyp_tests") is True
    assert config.keep_submodule("PyPedal.examples.new_methods") is False


def test_qt_webengine_and_griffon_are_dropped() -> None:
    config = _bundle_config()
    assert config.drop_collected_path("PySide6/QtWebEngineCore") is True
    assert config.drop_collected_path("PySide6/QtWidgets") is False
    assert config.drop_collected_path("plugins/sqldrivers/libqsqlite.dylib") is True
    assert config.drop_collected_path("griffonbruxellois_2026_pyp.ped") is True
    assert config.drop_collected_path("griffonbruxellois_2026_named_pyp.ped") is True
    assert config.drop_collected_path("customtkinter/windows.py") is True


def test_info_plist_identity() -> None:
    config = _bundle_config()
    plist = config.info_plist("4.1.0")
    assert plist["CFBundleName"] == "PyPedal"
    assert plist["CFBundleDisplayName"] == "PyPedal"
    assert plist["CFBundleIdentifier"] == "org.pypedal.PyPedal"
    assert plist["CFBundlePackageType"] == "APPL"
    assert plist["CFBundleShortVersionString"] == "4.1.0"
    assert plist["CFBundleExecutable"] == "PyPedal"


def test_icns_is_optional_until_an_asset_is_supplied() -> None:
    config = _bundle_config()
    assert config.icns_path(Path(REPO)) is None
    assert not (MACOS / "PyPedal.icns").exists()
