"""LIC-3..LIC-6: project license identity and packaging hygiene.

LIC-8/LIC-9 live in this file as well (README and authors/maintainers).
"""

import tomllib
from pathlib import Path

import pytest
from _pedhelpers import REPO

from PyPedal.__version__ import version as PACKAGE_VERSION

pytestmark = pytest.mark.packaging

ROOT = Path(REPO)
LICENSE = ROOT / "LICENSE"
PYPROJECT = ROOT / "pyproject.toml"
MANIFEST = ROOT / "MANIFEST.in"


def _project():
    with open(PYPROJECT, "rb") as fh:
        return tomllib.load(fh)


def test_lic3_pyproject_license_is_lgpl_2_1_or_later():
    project = _project()["project"]
    license_table = project["license"]
    assert license_table == {"text": "LGPL-2.1-or-later"}
    classifiers = project.get("classifiers", [])
    assert not any(c.startswith("License ::") for c in classifiers)


def test_lic3_setuptools_ships_only_the_project_license_file():
    setuptools = _project()["tool"]["setuptools"]
    assert setuptools["license-files"] == ["LICENSE"]


def test_lic4_root_license_is_official_lgpl_2_1_text():
    text = LICENSE.read_text(encoding="utf-8")
    stripped = text.lstrip()
    assert stripped.startswith("GNU LESSER GENERAL PUBLIC LICENSE")
    assert "Version 2.1, February 1999" in text
    assert "Version 2, June 1991" not in stripped[:400]
    assert "This license, the Lesser General Public License" in text
    assert "version 2.1 of the License, or (at your option) any later version" in text


def test_lic5_root_gpl_v2_licence_file_is_gone():
    assert not (ROOT / "LICENCE.txt").exists()
    assert not (ROOT / "LICENSE.txt").exists()
    assert LICENSE.is_file()


def test_lic6_adodb_is_not_the_project_license():
    assert not (ROOT / "PyPedal" / "LICENSE.txt").exists()
    assert not (ROOT / "ADODB-LICENSE.txt").exists()
    manifest = MANIFEST.read_text(encoding="utf-8")
    assert "include LICENSE" in manifest
    assert "LICENSE.txt" not in manifest
    assert "ADODB" not in manifest


def test_lic8_readme_is_current_candidate_and_not_published_final():
    project = _project()["project"]
    assert project["version"] == PACKAGE_VERSION
    assert project["requires-python"] == ">=3.12"
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{PACKAGE_VERSION}]" in changelog

    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert text.startswith("# PyPedal\n")
    assert not text.startswith("# PyPedal ")
    lowered = text.lower()
    assert "lgpl-2.1-or-later" in lowered or (
        "lesser general public license" in lowered and "2.1" in lowered and "later" in lowered
    )
    assert "3.12" in text
    assert "John B. Cole" in text
    assert "Geert Degezelle" in text
    assert (
        "has not been published" in lowered
        or "not published" in lowered
        or "no pypi package" in lowered
    )
    assert "available on pypi" not in lowered
    assert "v4.0.0-rc3 tag already exists" not in lowered
    assert "docs/RELEASE-4.0.0-FINAL.md" not in text


def test_lic9_authors_and_maintainers_are_separated():
    project = _project()["project"]
    authors = project["authors"]
    maintainers = project["maintainers"]
    assert authors == [{"name": "John B. Cole, PhD", "email": "john.b.cole@gmail.com"}]
    assert maintainers == [{"name": "Geert Degezelle", "email": "geertdegezelle@telenet.be"}]
    assert authors[0]["name"] != maintainers[0]["name"]
