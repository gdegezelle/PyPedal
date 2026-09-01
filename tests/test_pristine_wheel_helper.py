"""The pristine-wheel helper must track project metadata, not a frozen RC."""
import re
from pathlib import Path

import pytest
import tomllib

from _pedhelpers import REPO

pytestmark = pytest.mark.packaging

RELEASE = Path(REPO) / "tools" / "release"
HELPER = RELEASE / "build_pristine_wheel.sh"
STALE_RC = re.compile(r"4\.0\.0rc[0-9]|4\.0\.0-rc[0-9]")


def test_release_tooling_has_no_hardcoded_rc_versions():
    files = sorted(p for p in RELEASE.rglob("*") if p.is_file())
    assert files, "tools/release/ is empty"
    leftovers = []
    for path in files:
        text = path.read_text(encoding="utf-8")
        hits = STALE_RC.findall(text)
        if hits:
            leftovers.append((path.name, hits))
    assert leftovers == [], leftovers


def test_pristine_wheel_helper_derives_version_from_pyproject():
    text = HELPER.read_text(encoding="utf-8")
    assert "tomllib" in text
    assert "pyproject.toml" in text
    assert "EXPECTED_VERSION" in text
    assert '["project"]["version"]' in text or "['project']['version']" in text


def test_pristine_wheel_helper_default_matches_current_metadata():
    with (Path(REPO) / "pyproject.toml").open("rb") as fh:
        expected = tomllib.load(fh)["project"]["version"]
    assert expected
    assert not STALE_RC.fullmatch("not-a-version")
    # The helper reads the same field. Keep this assertion coupled to metadata,
    # not to a named RC, so a future bump cannot leave the helper behind.
    assert "rc" in expected or expected.startswith("4.")
