"""
Session-scoped guard: the test suite must leave the repository as it found it.

The invariant is a **delta**, not global cleanliness::

    repository state after the run == repository state before the run

A developer with legitimate uncommitted work must be able to run the tests
without spurious failures, so pre-existing modifications are deliberately
invisible here: they appear identically in both snapshots and cancel out. What
must not happen is the *suite* creating, modifying or deleting anything.

Why a filesystem manifest rather than ``git status``
----------------------------------------------------
``git status`` does not see ignored files, and ignored files are exactly where
this repository's leaks hid. ``pyp_utils.renumber(io='yes')`` with no filetag
writes ``_renumbered__renum.ped`` and ``_renumbered__id_map.map`` into the
working directory; both are matched by the ``_renumbered_*`` line in
``.gitignore``, so they were rewritten on every single test run and ``git
status`` stayed silent about it for as long as anyone had been looking. A guard
built on ``git status`` would have reported a clean run forever.

So the manifest walks the filesystem and hashes contents, ignoring nothing
except infrastructure that is not part of the repository's meaningful state
(``.git``, caches, ``build/``, egg-info).

Known blind spot
----------------
A byte-identical rewrite of a file that already existed before the run is not
detected, because the repository state genuinely did not change. On a clean
checkout -- which is what CI runs -- such a file does not exist beforehand, so
it is caught as a creation. Every writer known at the time of writing has been
fixed to use a temporary directory, so there is nothing left for this blind spot
to hide.
"""
import hashlib
import os

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Directories that are not part of the repository's meaningful state.
PRUNE_DIRS = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".hypothesis",
    "htmlcov",
    "build",
    "dist",
    ".tox",
    ".venv",
    "venv",
}

# Files written by the test *tooling* rather than by the code under test.
#
# Excluded explicitly rather than by accident. `pytest --cov` writes `.coverage`
# after session fixtures have torn down, so it currently escapes the snapshot by
# timing alone; naming it here means the exclusion survives a change in that
# ordering instead of silently becoming a false pass.
PRUNE_FILES = {".coverage", "coverage.xml"}


def _is_tooling_artifact(relpath):
    name = os.path.basename(relpath)
    return name in PRUNE_FILES or name.startswith(".coverage.")


def _manifest(root):
    """Map repo-relative path -> sha256 of contents, for every file under root."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [
            d for d in dirnames if d not in PRUNE_DIRS and not d.endswith(".egg-info")
        ]
        for name in filenames:
            path = os.path.join(dirpath, name)
            rel = os.path.relpath(path, root)
            if _is_tooling_artifact(rel):
                continue
            try:
                with open(path, "rb") as handle:
                    out[rel] = hashlib.sha256(handle.read()).hexdigest()
            except OSError:
                # A file that cannot be read is still a file that exists; record
                # its presence so appearance and disappearance are both visible.
                out[rel] = "<unreadable>"
    return out


def _describe(created, modified, deleted):
    lines = [
        "The test suite changed the repository. Tests must confine their output "
        "to temporary directories.",
        "",
        "Use tests/_pedhelpers.py: load_example() repoints pedfile at a temporary "
        "copy so <filetag>_*.dat output follows it, and chdir_tmp() handles "
        "routines that write relative to the working directory "
        "(fast_reorder(io='yes'), renumber(io='yes')).",
        "",
    ]
    for label, items in (
        ("created", created),
        ("modified", modified),
        ("deleted", deleted),
    ):
        if items:
            lines.append(f"{len(items)} {label}:")
            lines.extend(f"    {item}" for item in sorted(items))
            lines.append("")
    return "\n".join(lines)


@pytest.fixture(scope="session", autouse=True)
def repository_delta_guard(request):
    """Fail the session if the suite created, modified or deleted any file."""
    if request.config.getoption("--no-repo-guard"):
        yield
        return

    before = _manifest(REPO)
    yield
    after = _manifest(REPO)

    created = set(after) - set(before)
    deleted = set(before) - set(after)
    modified = {k for k in set(before) & set(after) if before[k] != after[k]}

    if created or modified or deleted:
        pytest.fail(_describe(created, modified, deleted), pytrace=False)


def pytest_addoption(parser):
    parser.addoption(
        "--no-repo-guard",
        action="store_true",
        default=False,
        help="Skip the session-scoped repository-delta guard.",
    )
