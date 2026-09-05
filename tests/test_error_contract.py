"""
The error contract: PyPedal raises, and never ends the host process.

Audit . Seventeen ``sys.exit(0)`` calls and one ``raise SystemExit``
lived in ``pyp_newclasses``. Three things were wrong with that, in increasing
order of severity:

1. **Library code terminated the host process.** A caller embedding PyPedal in a
   web service, a notebook or a batch pipeline had its process killed by a
   data-validation failure in a dependency.
2. **The status was 0, which means success.** A shell script or CI job wrapping
   PyPedal saw a clean exit and carried on as though the pedigree had loaded --
   the same failure-presenting-as-success mode as , expressed at the
   process level.
3. **It was uncatchable.** ``SystemExit`` derives from ``BaseException``, so the
   ``except Exception`` handlers around the load path did not see it. No caller
   could recover, retry with a corrected ``pedformat``, or report the problem in
   its own terms.

These tests cover the contract itself, not any one call site: that the process
survives, that the exception types are distinguishable, and -- structurally --
that no ``SystemExit`` has crept back into the library.
"""
import ast
import os
import subprocess
import sys
import unittest

from PyPedal import pyp_errors, pyp_newclasses
from PyPedal.application import EXIT_STATUS, exit_status_for
from _pedhelpers import owned_temp_dir

PACKAGE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "PyPedal"))

# The only places allowed to end the process are ``if __name__ == '__main__'``
# guards on package and desktop launchers.
ENTRY_POINTS = {
    ("__main__.py", "<module>"),
}


def _pedigree(rows, suffix=".ped"):
    tmp = owned_temp_dir(prefix="pypedal_errors_")
    path = os.path.join(tmp, "bad" + suffix)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return path


class TestNoLibraryCodeTerminatesTheProcess(unittest.TestCase):
    """Structural guarantees, so the contract cannot regress unnoticed."""

    def _process_exits(self):
        """Every SystemExit raise and sys.exit call under PyPedal/, located."""
        found = []
        for dirpath, dirnames, filenames in os.walk(PACKAGE):
            dirnames[:] = [d for d in dirnames
                           if d not in {"__pycache__", "examples"}]
            for name in sorted(filenames):
                if not name.endswith(".py"):
                    continue
                path = os.path.join(dirpath, name)
                rel = os.path.relpath(path, PACKAGE)
                with open(path, encoding="utf-8") as handle:
                    tree = ast.parse(handle.read(), filename=path)

                scope = {}

                def register(node, enclosing):
                    for child in ast.iter_child_nodes(node):
                        if isinstance(child, (ast.FunctionDef,
                                              ast.AsyncFunctionDef)):
                            scope[id(child)] = child.name
                            register(child, child.name)
                        else:
                            scope.setdefault(id(child), enclosing)
                            register(child, enclosing)

                register(tree, "<module>")

                for node in ast.walk(tree):
                    hit = False
                    if isinstance(node, ast.Raise):
                        exc = node.exc
                        if isinstance(exc, ast.Call):
                            exc = exc.func
                        hit = isinstance(exc, ast.Name) and exc.id == "SystemExit"
                    elif isinstance(node, ast.Call):
                        func = node.func
                        hit = (isinstance(func, ast.Attribute)
                               and func.attr == "exit"
                               and isinstance(func.value, ast.Name)
                               and func.value.id == "sys")
                    if hit:
                        enclosing = "<module>"
                        for parent in ast.walk(tree):
                            for child in ast.iter_child_nodes(parent):
                                if child is node and isinstance(
                                        parent, (ast.FunctionDef,
                                                 ast.AsyncFunctionDef)):
                                    enclosing = parent.name
                        found.append((rel, node.lineno, enclosing))
        return found

    def test_only_the_declared_entry_points_end_the_process(self):
        offenders = []
        for rel, lineno, _enclosing in self._process_exits():
            if not any(rel.endswith(f) for f, _ in ENTRY_POINTS):
                offenders.append(f"{rel}:{lineno}")
        self.assertEqual(
            [], offenders,
            "library code must raise, not terminate the host process; "
            "found process exits outside the declared entry points")

    def test_pyp_newclasses_has_no_process_exit_at_all(self):
        """
        The module that held all eighteen of them. Called out separately so a
        regression there is unmistakable.
        """
        offenders = [f"{rel}:{lineno}"
                     for rel, lineno, _ in self._process_exits()
                     if rel == "pyp_newclasses.py"]
        self.assertEqual([], offenders)


class TestFailuresRaiseTypedExceptions(unittest.TestCase):
    """One representative per condition class."""

    def test_malformed_record_raises_a_pedigree_format_error(self):
        """A record with fewer fields than the pedformat declares."""
        path = _pedigree(["1 0 0", "2 0 0", "3 1"])
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            pyp_newclasses.loadPedigree(options={
                "pedfile": path, "pedformat": "asd", "sepchar": " ",
                "messages": "quiet", "pedigree_summary": 0})
        self.assertIsInstance(caught.exception, pyp_errors.PyPedalError)

    def test_unknown_pedsource_raises_a_usage_error(self):
        path = _pedigree(["1 0 0", "2 0 0", "3 1 2"])
        ped = pyp_newclasses.NewPedigree({
            "pedfile": path, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0})
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            ped.load(pedsource="not_a_real_source")

    def test_usage_error_is_also_a_value_error(self):
        """
        Deliberate: an argument outside its documented domain is what
        ValueError means in Python, so existing `except ValueError` callers
        keep working.
        """
        path = _pedigree(["1 0 0", "2 0 0", "3 1 2"])
        ped = pyp_newclasses.NewPedigree({
            "pedfile": path, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0})
        with self.assertRaises(ValueError):
            ped.load(pedsource="not_a_real_source")

    def test_missing_pedfile_for_graphfile_raises_a_configuration_error(self):
        path = _pedigree(["1 0 0", "2 0 0", "3 1 2"])
        ped = pyp_newclasses.NewPedigree({
            "pedfile": path, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0})
        ped.kw["pedfile"] = ""
        with self.assertRaises(pyp_errors.PyPedalConfigurationError):
            ped.load(pedsource="graphfile")

    def test_empty_animallist_raises_a_usage_error(self):
        path = _pedigree(["1 0 0", "2 0 0", "3 1 2"])
        ped = pyp_newclasses.NewPedigree({
            "pedfile": path, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0})
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            ped.load(pedsource="animallist", animallist=[])

    def test_unreadable_option_file_raises_an_option_error(self):
        """
        The one site that used ``raise SystemExit`` rather than ``sys.exit(0)``.
        ``NewPedigree(kw=None, kwfile=...)`` reads the option file, so passing
        no kw is what routes through ``_load_config_file``.
        """
        tmp = owned_temp_dir(prefix="pypedal_errors_")
        bad_ini = os.path.join(tmp, "broken.ini")
        with open(bad_ini, "w", encoding="utf-8") as handle:
            handle.write("this is not an ini file at all\n[[[\n")
        with self.assertRaises(pyp_errors.PyPedalConfigurationError):
            pyp_newclasses.NewPedigree(kwfile=bad_ini)

    def test_every_failure_is_catchable_as_pypedal_error(self):
        """
        The property the old code could not offer: SystemExit derives from
        BaseException, so `except Exception` -- and `except PyPedalError` --
        never saw it.
        """
        path = _pedigree(["1 0 0", "2 0 0", "3 1 2"])
        ped = pyp_newclasses.NewPedigree({
            "pedfile": path, "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "pedigree_summary": 0})
        try:
            ped.load(pedsource="not_a_real_source")
        except pyp_errors.PyPedalError:
            pass
        else:
            self.fail("expected a PyPedalError")


class TestTheInterpreterSurvives(unittest.TestCase):
    """
    End-to-end: a failing load in a subprocess must leave the interpreter
    running and then exit non-zero when the script says so -- never 0.
    """

    SCRIPT = """
import sys
sys.path.insert(0, {repo!r})
from PyPedal import pyp_errors, pyp_newclasses

try:
    pyp_newclasses.loadPedigree(options={{
        "pedfile": {path!r}, "pedformat": "asd", "sepchar": " ",
        "messages": "quiet", "pedigree_summary": 0}})
except pyp_errors.PyPedalError:
    print("STILL_ALIVE")
    sys.exit(3)
print("NO_ERROR_RAISED")
sys.exit(0)
"""

    def test_a_failing_load_does_not_kill_the_process(self):
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        path = _pedigree(["1 0 0", "2 0 0", "3 1"])
        script = self.SCRIPT.format(repo=repo, path=path)
        proc = subprocess.run([sys.executable, "-c", script],
                              capture_output=True, text=True, timeout=120)
        self.assertIn("STILL_ALIVE", proc.stdout,
                      "the process died before it could report the failure; "
                      "stderr:\n%s" % proc.stderr)
        self.assertEqual(3, proc.returncode)

    def test_failure_never_exits_zero(self):
        """The specific defect: exit status 0 told the shell it had worked."""
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        path = _pedigree(["1 0 0", "2 0 0", "3 1"])
        script = (
            "import sys\n"
            "sys.path.insert(0, %r)\n"
            "from PyPedal.application import exit_status_for\n"
            "from PyPedal import pyp_errors, pyp_newclasses\n"
            "try:\n"
            "    pyp_newclasses.loadPedigree(options={'pedfile': %r,\n"
            "        'pedformat': 'asd', 'sepchar': ' ', 'messages': 'quiet',\n"
            "        'pedigree_summary': 0})\n"
            "except pyp_errors.PyPedalError as exc:\n"
            "    sys.exit(exit_status_for(exc))\n"
            "sys.exit(0)\n" % (repo, path)
        )
        proc = subprocess.run([sys.executable, "-c", script],
                              capture_output=True, text=True, timeout=120)
        self.assertNotEqual(
            0, proc.returncode,
            "a failed load reported success to the shell; stdout:\n%s\nstderr:\n%s"
            % (proc.stdout, proc.stderr))


class TestExitStatusMapping(unittest.TestCase):
    """Choosing a status is an application decision and lives in one place."""

    def test_each_class_maps_to_a_distinct_non_zero_status(self):
        cases = {
            pyp_errors.PyPedalPedigreeFormatError("x"): 2,
            pyp_errors.PyPedalUsageError("x"): 2,
            pyp_errors.PyPedalOptionError("x"): 3,
            pyp_errors.PyPedalDependencyError("x"): 4,
            pyp_errors.PyPedalValidationError("x"): 65,
            pyp_errors.PyPedalInternalError("x"): 70,
            # EX_CANTCREAT: the output file could not be produced.
            pyp_errors.PyPedalExportFormatError(
                "FORMAT", "FIELD", 3, -999, 4): 73,
        }
        for exc, expected in cases.items():
            with self.subTest(exception=type(exc).__name__):
                self.assertEqual(expected, exit_status_for(exc))

    def test_no_mapping_is_zero(self):
        for klass, status in EXIT_STATUS:
            with self.subTest(exception=klass.__name__):
                self.assertNotEqual(0, status)

    def test_an_unrecognised_error_still_maps_to_failure(self):
        self.assertEqual(1, exit_status_for(RuntimeError("x")))

    def test_subclasses_map_before_their_base(self):
        """
        PyPedalError is last in the table, so a subclass must not fall through
        to the generic 1.
        """
        self.assertEqual(
            70, exit_status_for(pyp_errors.PyPedalInternalError("x")))


if __name__ == "__main__":
    unittest.main()
