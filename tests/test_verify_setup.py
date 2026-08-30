"""
``verify_setup.py`` must agree with the dependency model in ``pyproject.toml``.

an earlier revision stage E1 moved ``matplotlib`` out of the runtime dependencies and into
the ``graphics`` extra, because ``pyp_graphics`` imports it lazily inside each
plotting function and the library is fully usable without it.
``verify_setup.py`` was not updated, so it went on listing matplotlib as
*required* and declared a perfectly good ``pip install -e ".[dev]"`` broken.

That was invisible in a development environment with every extra installed, and
only surfaced when the 3.12/3.13/3.14 matrix was run in clean venvs.

These tests pin the contract in both directions: a core install must verify, and
an optional dependency must be *reported* rather than *required*.
"""
import importlib.util
import os
import subprocess
import sys
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VERIFY = os.path.join(REPO, "verify_setup.py")

_spec = importlib.util.spec_from_file_location("verify_setup_under_test", VERIFY)
verify_setup = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(verify_setup)


class TestDependencyModelMatchesPackaging(unittest.TestCase):

    def test_core_dependencies_are_exactly_the_runtime_ones(self):
        """
        The core list must track ``[project] dependencies`` in pyproject.toml.
        If they drift, one of them is lying about what PyPedal needs.
        """
        import tomllib

        with open(os.path.join(REPO, "pyproject.toml"), "rb") as handle:
            declared = tomllib.load(handle)["project"]["dependencies"]

        # "numpy>=1.26" -> "numpy"
        names = {d.split(">")[0].split("=")[0].split("[")[0].strip().lower()
                 for d in declared}
        core = {m.lower() for m in verify_setup.CORE_DEPENDENCIES}
        self.assertEqual(
            names, core,
            "verify_setup.CORE_DEPENDENCIES and pyproject [project].dependencies "
            "disagree about what is required")

    def test_matplotlib_is_optional_not_core(self):
        """The specific regression: it belongs to the graphics extra."""
        self.assertNotIn("matplotlib", verify_setup.CORE_DEPENDENCIES)
        self.assertIn("matplotlib",
                      verify_setup.OPTIONAL_DEPENDENCIES["graphics"])

    def test_every_optional_group_names_a_real_extra(self):
        """A group naming a non-existent extra prints uninstallable advice."""
        import tomllib

        with open(os.path.join(REPO, "pyproject.toml"), "rb") as handle:
            extras = set(tomllib.load(handle)["project"]["optional-dependencies"])
        for group in verify_setup.OPTIONAL_DEPENDENCIES:
            with self.subTest(extra=group):
                self.assertIn(group, extras)


class TestVerificationOutcome(unittest.TestCase):

    def test_passes_with_the_current_environment(self):
        self.assertTrue(verify_setup.check_dependencies())

    def test_passes_when_an_optional_dependency_is_absent(self):
        """
        The case that was broken. Simulated by making matplotlib unimportable,
        which is what a plain ``pip install -e ".[dev]"`` looks like.
        """
        real_import = verify_setup._importable

        def without_matplotlib(module):
            return False if module == "matplotlib" else real_import(module)

        verify_setup._importable = without_matplotlib
        try:
            self.assertTrue(
                verify_setup.check_dependencies(),
                "a missing optional dependency must not fail core verification")
        finally:
            verify_setup._importable = real_import

    def test_fails_when_a_core_dependency_is_absent(self):
        """The check must still have teeth."""
        real_import = verify_setup._importable

        def without_numpy(module):
            return False if module == "numpy" else real_import(module)

        verify_setup._importable = without_numpy
        try:
            self.assertFalse(verify_setup.check_dependencies())
        finally:
            verify_setup._importable = real_import

    def test_script_exits_zero_end_to_end(self):
        proc = subprocess.run([sys.executable, VERIFY], cwd=REPO,
                              capture_output=True, text=True, timeout=300)
        self.assertEqual(0, proc.returncode,
                         f"verify_setup.py failed:\n{proc.stdout[-1500:]}")


if __name__ == "__main__":
    unittest.main()
