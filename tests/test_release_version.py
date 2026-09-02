"""Release identity: package metadata and runtime version must agree.

The PEP 440 identity of this release is ``4.2.0``. Human-readable docs
use the same spelling.
"""
import tomllib
from importlib.metadata import version as distribution_version
from pathlib import Path

import pytest

import PyPedal
from PyPedal.__version__ import version as module_version

from _pedhelpers import REPO

pytestmark = pytest.mark.packaging

EXPECTED = "4.2.0"


def _project():
    with open(Path(REPO) / "pyproject.toml", "rb") as fh:
        return tomllib.load(fh)["project"]


def test_pyproject_version_is_current_release():
    assert _project()["version"] == EXPECTED


def test_runtime_package_version_is_current_release():
    assert PyPedal.__version__ == EXPECTED


def test_version_module_agrees_with_package():
    assert module_version == EXPECTED
    assert module_version == PyPedal.__version__


def test_installed_distribution_version_is_current_release():
    assert distribution_version("PyPedal") == EXPECTED


def test_requires_python_is_3_12_or_newer():
    assert _project()["requires-python"] == ">=3.12"


def test_classifiers_name_3_12_through_3_14_only():
    classifiers = _project()["classifiers"]
    python_classifiers = [
        c for c in classifiers if c.startswith("Programming Language :: Python :: 3.")
    ]
    assert python_classifiers == [
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ]
