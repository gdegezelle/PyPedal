"""
The blocking Ruff correctness subset must stay clean.

Three rules are gated, and they are correctness rather than style:

===== =============================== =========================================
F821  undefined name                  a name that is never bound
F811  redefinition of unused name     a definition silently replaced by a later
                                      one
F601  repeated dict key literal       a key whose earlier value is discarded
===== =============================== =========================================

Each has found a real defect here: F821 the ``pyp_graphics`` and ``pyp_jbc``
``NameError``s repaired in an earlier revision, F811 the duplicate ``logging`` import in
``pyp_io``, and F601 the ``STUD_ID`` collision that became Findings 28 and 29.
The remaining ~390 Ruff findings are stylistic and stay reporting-only.

WHY THIS IS A TEST AND NOT ONLY A CI STEP
-----------------------------------------
The gate also runs in ``.github/workflows/ci.yml``. But this repository is
private and has never been pushed, so that job has never executed -- the local
matrix is what actually runs. A gate that only exists in CI is a gate that does
not exist yet.

E999 IS DELIBERATELY NOT IN THE SUBSET
--------------------------------------
Ruff removed ``E999`` as a selectable rule in 0.16; syntax errors are now always
reported. Passing ``--select E999`` is a hard error ("Rule `E999` was removed and
cannot be selected"), which would make the gate fail for a reason unrelated to
the code. ``test_the_subset_is_selectable`` guards against that returning.
"""
import os
import shutil
import subprocess
import unittest

import pytest

from _pedhelpers import REPO

# The blocking subset. Kept in one place; CI names the same three.
BLOCKING = ["F821", "F811", "F601"]

ruff = shutil.which("ruff")


@pytest.mark.skipif(ruff is None, reason="ruff is not installed (pip install -e '.[dev]')")
class TestBlockingRuffSubset(unittest.TestCase):

    def _run(self, *args):
        return subprocess.run(
            [ruff, "check", "--select", ",".join(BLOCKING), *args],
            cwd=REPO, capture_output=True, text=True)

    def test_the_subset_is_clean(self):
        result = self._run("--output-format", "concise", ".")
        self.assertEqual(
            0, result.returncode,
            "The blocking Ruff correctness subset is not clean:\n"
            + result.stdout + result.stderr)

    def test_the_subset_is_selectable(self):
        """
        A removed or renamed rule makes Ruff exit 2 with an error on stderr,
        which is not the same thing as a clean run and must not be mistaken for
        one. This is how E999 would have broken the gate silently if the check
        above were the only one.
        """
        result = self._run("--output-format", "concise", ".")
        self.assertNotIn("cannot be selected", result.stderr)
        self.assertNotEqual(2, result.returncode, result.stderr)

    def test_ci_gates_the_same_rules_this_test_does(self):
        """
        The subset is written down in two places. If they drift, the local gate
        and the CI gate stop meaning the same thing -- so the drift is a test
        failure rather than a surprise on someone else's machine.
        """
        workflow = os.path.join(REPO, ".github", "workflows", "ci.yml")
        text = open(workflow).read()
        expected = "--select " + ",".join(BLOCKING)
        self.assertIn(
            expected, text,
            "ci.yml no longer gates %s; update BLOCKING here and the workflow "
            "together" % expected)
        # And it must not be neutered with `|| true`, which is how the
        # reporting-only steps in the same job are written.
        for line in text.splitlines():
            if expected in line:
                self.assertNotIn("|| true", line, "the gate is not blocking")

    def test_ci_keeps_linux_matrix_and_adds_cost_conscious_platforms(self):
        workflow = os.path.join(REPO, ".github", "workflows", "ci.yml")
        text = open(workflow).read()
        self.assertIn('python-version: ["3.12", "3.13", "3.14"]', text)
        self.assertIn("ubuntu-latest", text)
        self.assertIn("macos-latest", text)
        self.assertIn("windows-latest", text)
        self.assertIn('pytest tests/ -m "not integration"', text)
        self.assertIn("python-version: \"3.12\"", text)
        self.assertIn('pip install -e ".[dev]"', text)

    def test_ci_adds_dedicated_desktop_jobs_without_converting_library_jobs(self):
        workflow = os.path.join(REPO, ".github", "workflows", "ci.yml")
        text = open(workflow).read()
        self.assertIn("desktop-test", text)
        self.assertIn("QT_QPA_PLATFORM", text)
        self.assertIn("tests/test_desktop/", text)
        install_lines = [
            line.strip()
            for line in text.splitlines()
            if "pip install -e" in line
        ]
        self.assertTrue(
            any(line == 'pip install -e ".[dev]"' for line in install_lines),
            install_lines,
        )
        self.assertTrue(
            any("desktop-test" in line for line in install_lines),
            install_lines,
        )


if __name__ == "__main__":
    unittest.main()
