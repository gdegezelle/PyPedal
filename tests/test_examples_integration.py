"""
End-to-end example programs for the current PyPedal 4.0 product surface.

Each collected top-level script runs as a subprocess in a throwaway copy of
``PyPedal/examples``. The copy is required because the scripts resolve ``.ini``
and pedigree files by relative path; copying is what stops them writing
``.log`` / ``.dat`` / graphics into the repository.

Classification
--------------
* **SUPPORTED** active examples are collected and must exit 0.
* **Optional-dependency** examples are skipped with a reason when the extra
  is absent (never a false pass from a silent no-op).
* **Historical** scripts live under subdirectories (``historical/``,
  ``merge/``, ``sets/``, ``horse/``) and are **not** collected.
  ``profile/`` was removed in 4.0.0-rc7 and must stay absent.
* **NOT_EXAMPLES** are top-level ``.py`` files that are not product demos.
* **KNOWN_FAILING** is only for remaining current 4.0rc3 production bounds
  that an otherwise supported script still hits. Every reason must describe
  present behaviour. An unexpected pass fails the test.

This file is not a 2.0.4 ``maketest.sh`` transcription.
"""
import importlib.util
import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXAMPLES = os.path.join(REPO, "PyPedal", "examples")

# Per-script timeout. The examples are smoke tests over small pedigrees; one
# that runs longer than this has changed character and should be looked at.
TIMEOUT = 300

# Top-level .py files that are not PyPedal 4 product examples.
NOT_EXAMPLES = {
    # Matches python_files = ["test_*.py"]. Harmless today because testpaths
    # restricts collection to tests/, but it would be collected by a bare
    # `pytest PyPedal/`.
    "test_reordering.py",
    # A one-off CSV fixer, not a demonstration of the library.
    "fix_csv.py",
    # Standalone pandas CSV utility; not a PyPedal API example.
    "animal_analysis.py",
    # Loads new_amatrix.ini / new_amatrix.ped (asdgb). That year column is
    # real chronology; dam 2047 is later than named offspring. RC4 refuses
    # the file. Kept as evidence, not an active successful example.
    "new_amatrix.py",
}

# Subdirectories that may contain .py files but are not the active suite.
# Discovery is top-level only; this map is the explicit reason those trees
# stay out of the harness.
HISTORICAL_SUBDIRS = {
    "historical": "removed APIs archived from the active 4.0 suite",
    "merge": "historical merge-pedigree experiments; not a current 4.0 workflow",
    "sets": "historical set-operation experiments",
    "horse": "historical horse-pedigree scratch, not an active example",
}

# Capability -> (probe module, extra that provides it).
CAPABILITIES = {
    "graphics": ("pydot", "graphics"),
    "graphviz": ("pygraphviz", "graphviz-extra"),
    "reports": ("reportlab", "reports"),
}

# Example -> capability it needs to do its advertised work.
#
# Determined empirically, not by reading imports: importing pyp_graphics does
# not require matplotlib or pydot (both are imported lazily), so grepping the
# source over-reports. What is listed here is what actually failed, or silently
# no-opped, in a controlled environment.
REQUIRES = {
    # ModuleNotFoundError: pydot without [graphics]
    "new_classes.py": "graphics",
    "new_doug.py": "graphics",
    "new_graphics.py": "graphics",
    "new_graphics2.py": "graphics",
    "new_graphics4.py": "graphics",
    "new_hartl.py": "graphics",
    "new_jbc.py": "graphics",
    "new_reporting.py": "reports",
    "new_simulate.py": "graphics",
    # exits 0 having drawn nothing without [graphviz-extra], because
    # pyp_graphics.new_draw_pedigree catches the ImportError and returns
    "new_graphics3.py": "graphviz",
}


def _capability_available(capability):
    module, _extra = CAPABILITIES[capability]
    return importlib.util.find_spec(module) is not None


def _skip_reason(script):
    """Reason string if `script` cannot run here, else None."""
    capability = REQUIRES.get(script)
    if capability is None or _capability_available(capability):
        return None
    module, extra = CAPABILITIES[capability]
    return (f"{script} needs the '{capability}' capability ({module}), which is "
            f"not installed. Install with: pip install -e \".[{extra}]\"")


# Remaining failures that are current 4.0rc3 production bounds, not example
# drift. Do not add repaired scripts here. Do not list unsupported historical
# APIs that have been archived or rewritten as explicit demonstrations.
KNOWN_FAILING = {}


def _example_scripts():
    if not os.path.isdir(EXAMPLES):
        return []
    return sorted(
        name for name in os.listdir(EXAMPLES)
        if name.endswith(".py") and name not in NOT_EXAMPLES
    )


SCRIPTS = _example_scripts()


@pytest.fixture()
def examples_copy(tmp_path):
    """
    A throwaway copy of PyPedal/examples, ONE PER SCRIPT.

    Per-script, not shared, and that is not incidental. An earlier revision used
    a module-scoped copy and reported 15 failures where a per-script copy
    reports 20: scripts write output that later scripts then consume, so a
    shared directory makes the result depend on collection order and on which
    subset was selected. An integration suite whose answer changes with -k is
    not measuring anything.

    Copying costs a few hundred kilobytes per script. That is the price of a
    reproducible number. The copy lives under pytest ``tmp_path``.
    """
    workdir = tmp_path / "examples"
    shutil.copytree(EXAMPLES, workdir, ignore=_ignore_bulk)
    yield str(workdir)


# PyPedal/examples is ~305 MB, but 269 MB of that is one gitignored pickle
# (boichard2_pedigree.bin) and most of the rest is previously generated output.
# Copying all of it 37 times would be ~11 GB of I/O for no benefit. Inputs --
# .py, .ini, .ped, .csv -- are about 19 MB.
_BULK_SUFFIXES = (".bin", ".dat", ".log", ".dot", ".png", ".jpg", ".jpeg",
                  ".pdf", ".ps", ".eps", ".pkl", ".db", ".sqlite")


def _ignore_bulk(_directory, names):
    return {name for name in names
            if name.lower().endswith(_BULK_SUFFIXES)}


def _run(script, workdir):
    env = dict(os.environ)
    # Import PyPedal from the working tree, not from any installed copy.
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    # Non-interactive: several examples build plots.
    env["MPLBACKEND"] = "Agg"
    return subprocess.run(
        [sys.executable, script],
        cwd=workdir, env=env, capture_output=True, text=True, timeout=TIMEOUT,
    )


@pytest.mark.integration
@pytest.mark.parametrize("script", SCRIPTS)
def test_example_runs(script, examples_copy):
    # Gate BEFORE running. Skipping after the fact would still let the script
    # exit 0 having done nothing and be recorded as a pass.
    reason = _skip_reason(script)
    if reason:
        pytest.skip(reason)

    proc = _run(script, examples_copy)
    expected_failure = KNOWN_FAILING.get(script)

    if expected_failure:
        assert proc.returncode != 0, (
            f"{script} is on the known-failing list but SUCCEEDED.\n"
            f"Recorded reason: {expected_failure}\n"
            "If the underlying gap has been closed, remove it from "
            "KNOWN_FAILING -- the list is meant to shrink."
        )
        return

    assert proc.returncode == 0, (
        f"{script} exited {proc.returncode}.\n"
        f"--- stdout (tail) ---\n{proc.stdout[-2000:]}\n"
        f"--- stderr (tail) ---\n{proc.stderr[-2000:]}"
    )


@pytest.mark.integration
def test_the_examples_directory_is_not_touched(tmp_path):
    """
    The failure mode the original maketest.sh had: running the examples in
    place wrote their output into the source tree, which is where the 21
    tracked .dat files came from.
    """
    before = sorted(os.listdir(EXAMPLES))
    workdir = tmp_path / "examples"
    shutil.copytree(EXAMPLES, workdir, ignore=_ignore_bulk)
    _run("new_lacy.py", str(workdir))
    assert before == sorted(os.listdir(EXAMPLES))


def test_every_known_failure_still_exists_as_a_script():
    """
    Guard on the list itself: a KNOWN_FAILING entry naming a script that has
    been renamed or deleted would sit there forever, exempting nothing.

    Not marked integration -- it is a cheap consistency check on the list and
    should run in the default suite.
    """
    missing = [name for name in KNOWN_FAILING
               if not os.path.exists(os.path.join(EXAMPLES, name))]
    assert missing == [], f"KNOWN_FAILING names scripts that do not exist: {missing}"


def test_the_suite_actually_discovers_examples():
    """A discovery bug would turn this whole file into a silent no-op."""
    assert len(SCRIPTS) > 30, f"only discovered {len(SCRIPTS)} example scripts"


def test_discovery_is_top_level_only():
    """Historical subdirectory scripts must not enter the active suite."""
    assert all(os.sep not in name for name in SCRIPTS)
    for subdir, reason in HISTORICAL_SUBDIRS.items():
        path = os.path.join(EXAMPLES, subdir)
        assert os.path.isdir(path), f"{subdir} ({reason}) is missing"
    for name in SCRIPTS:
        assert os.path.isfile(os.path.join(EXAMPLES, name)), name
        assert not os.path.isdir(os.path.join(EXAMPLES, name))
    assert not os.path.isdir(os.path.join(EXAMPLES, "profile"))


def test_known_failing_reasons_are_current():
    """A KNOWN_FAILING entry without a current 4.0 reason is a silent skip."""
    stale_tokens = (
        "dropped in the port",
        "Stage C",
        "not investigated",
    )
    for name, reason in KNOWN_FAILING.items():
        assert reason.strip(), f"{name} has an empty KNOWN_FAILING reason"
        lowered = reason.lower()
        for token in stale_tokens:
            assert token.lower() not in lowered, (
                f"{name} still cites a stale reason ({token!r}): {reason}"
            )


def test_every_capability_names_a_real_extra():
    """
    A capability pointing at a non-existent extra prints uninstallable advice
    in its skip reason.

    Not marked integration -- a cheap consistency check that belongs in the
    default suite.
    """
    import tomllib

    with open(os.path.join(REPO, "pyproject.toml"), "rb") as handle:
        extras = set(tomllib.load(handle)["project"]["optional-dependencies"])
    for capability, (_module, extra) in CAPABILITIES.items():
        assert extra in extras, f"capability {capability!r} names unknown extra {extra!r}"


def test_every_gated_script_exists_and_is_collected():
    """A REQUIRES entry for a renamed script would gate nothing."""
    unknown = [s for s in REQUIRES if s not in SCRIPTS]
    assert unknown == [], f"REQUIRES names scripts that are not collected: {unknown}"


def test_the_gate_takes_precedence_over_the_expected_failure_list():
    """
    The two mechanisms answer different questions -- "can this run here?" and
    "is this known broken?" -- and they compose, in that order.

    A script may legitimately be both gated and known-failing: the gate
    answers "can this run here?" and KNOWN_FAILING answers "is this known
    broken?" They compose, in that order. After an earlier revision, no active script
    is on KNOWN_FAILING.

    The property asserted here is the ordering that makes the composition safe:
    the gate is consulted first, so a script that cannot run is never judged
    against an expectation about running it.
    """
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(test_example_runs))
    body = tree.body[0].body
    statements = [n for n in body if not isinstance(n, ast.Expr)]
    first = statements[0]
    assert isinstance(first, ast.Assign), "expected the gate check to come first"
    assert "_skip_reason" in ast.dump(first), (
        "test_example_runs must consult _skip_reason() before running the "
        "script or reading KNOWN_FAILING")


def test_gate_reports_available_capabilities_as_runnable():
    """When a capability is present, nothing is skipped for it."""
    for script, capability in REQUIRES.items():
        if _capability_available(capability):
            assert _skip_reason(script) is None


def test_gate_reason_names_the_extra_to_install():
    """A skip must tell the reader how to make it not skip."""
    for capability, (module, extra) in CAPABILITIES.items():
        if not _capability_available(capability):
            script = next(s for s, c in REQUIRES.items() if c == capability)
            reason = _skip_reason(script)
            assert reason is not None
            assert module in reason and extra in reason
