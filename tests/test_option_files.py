"""
Option-file (.ini) loading -- audit .

Four independent defects made option files unusable, and ~28 of the 39 example
scripts could not run:

  1. section-less files were REJECTED, though 17 of the 26 example option files
     use exactly that format;
  2. the loader returned ``{section: {key: value}}``, a shape nothing consumes
     -- every option is read as ``kw['renumber']``;
  3. every value stayed a string, so ``renumber = 0`` arrived as ``'0'``, which
     is true in Python: turning an option off turned it on;
  4. ``loadPedigree()`` normalised ``options=None`` to ``{}`` before passing it
     on, and the file is read only when ``kw is None`` -- so
     ``loadPedigree(optionsfile=...)`` ignored the file entirely and then failed
     for want of a pedfile.

These are integration tests over the real example files rather than synthetic
ones, because the formats in the repository are the specification here.
"""
import os
import shutil
import unittest

from PyPedal import pyp_io
from PyPedal.pyp_newclasses import loadPedigree

from _pedhelpers import owned_temp_dir, EXAMPLES, corpus


class TestValueCoercion(unittest.TestCase):
    """
    Coercion rules are pinned exactly. Type inference that is merely
    approximately right is worse than none, because it fails silently.
    """

    def test_integers_booleans_and_floats(self):
        cases = {
            "0": 0, "1": 1, "100": 100, "-999": -999,
            "true": True, "True": True, "yes": True, "on": True,
            "false": False, "False": False, "no": False, "off": False,
            "0.05": 0.05, "-1.5": -1.5,
            "none": None,
            "asd": "asd", "quiet": "quiet",
        }
        for raw, expected in cases.items():
            with self.subTest(raw=raw):
                self.assertEqual(expected, pyp_io.coerce_ini_value(raw))

    def test_booleans_are_not_coerced_to_int(self):
        """`True` must not fall through to the int branch and become 1."""
        self.assertIs(True, pyp_io.coerce_ini_value("true"))
        self.assertIs(False, pyp_io.coerce_ini_value("false"))

    def test_integers_stay_integers(self):
        self.assertIsInstance(pyp_io.coerce_ini_value("100"), int)

    def test_surrounding_quotes_are_stripped(self):
        """`messages = 'verbose'` must compare equal to `verbose`."""
        self.assertEqual("verbose", pyp_io.coerce_ini_value("'verbose'"))
        self.assertEqual("verbose", pyp_io.coerce_ini_value('"verbose"'))
        self.assertEqual(",", pyp_io.coerce_ini_value('","'))

    def test_quoted_digits_stay_strings(self):
        """Quoting is how a caller says 'this is text, not a number'."""
        self.assertEqual("0", pyp_io.coerce_ini_value("'0'"))


class TestRealOptionFiles(unittest.TestCase):

    def test_section_less_legacy_file_loads(self):
        options = pyp_io.read_ini_file(os.path.join(EXAMPLES, "new_lacy.ini"))
        self.assertEqual("new_lacy.ped", options["pedfile"])
        self.assertEqual("asd", options["pedformat"])
        self.assertEqual("quiet", options["messages"])
        self.assertEqual(1, options["renumber"])

    def test_sectioned_file_is_flattened(self):
        options = pyp_io.read_ini_file(os.path.join(EXAMPLES, "new_format.ini"))
        # flat, not {section: {...}}
        self.assertEqual("boichard2a.ped", options["pedfile"])
        for value in options.values():
            self.assertNotIsInstance(value, dict)

    def test_every_example_option_file_loads(self):
        """
        The whole point of the fix. One file legitimately holds several
        alternative configurations and needs a section chosen; every other one
        must load unaided.
        """
        multi_config = {"new_inbreeding2multiple.ini"}
        failures = []
        checked = 0
        for name in sorted(os.listdir(EXAMPLES)):
            if not name.endswith(".ini") or name in multi_config:
                continue
            checked += 1
            try:
                pyp_io.read_ini_file(os.path.join(EXAMPLES, name))
            except Exception as exc:
                failures.append("%s: %s: %s" % (name, type(exc).__name__, exc))
        self.assertGreater(checked, 20)
        self.assertEqual([], failures)

    def test_multi_configuration_file_requires_a_section(self):
        """
        new_inbreeding2multiple.ini defines a [noinbreeding] pedigree and a
        [horse] pedigree. Flattening would merge two different pedigrees into
        one, so it must raise rather than guess -- and the message must say how
        to resolve it.
        """
        path = os.path.join(EXAMPLES, "new_inbreeding2multiple.ini")
        with self.assertRaises(ValueError) as caught:
            pyp_io.read_ini_file(path)
        self.assertIn("section=", str(caught.exception))

        horse = pyp_io.read_ini_file(path, section="horse")
        self.assertEqual("horse.ped", horse["pedfile"])
        self.assertEqual("ASD", horse["pedformat"])
        self.assertEqual(",", horse["sepchar"])

    def test_unknown_section_names_the_available_ones(self):
        path = os.path.join(EXAMPLES, "new_inbreeding2multiple.ini")
        with self.assertRaises(ValueError) as caught:
            pyp_io.read_ini_file(path, section="nonexistent")
        self.assertIn("horse", str(caught.exception))


class TestOptionFileEndToEnd(unittest.TestCase):

    def _write(self, body):
        tmp = owned_temp_dir(prefix="pypedal_ini_")
        shutil.copy(corpus("new_lacy.ped"), os.path.join(tmp, "new_lacy.ped"))
        path = os.path.join(tmp, "options.ini")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(body)
        return tmp, path

    def _load(self, body):
        tmp, path = self._write(body)
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            return loadPedigree(optionsfile=path)
        finally:
            os.chdir(cwd)

    def test_options_file_is_actually_read(self):
        ped = self._load(
            "messages = quiet\npedigree_summary = 0\n"
            "pedfile = new_lacy.ped\npedformat = asd\nrenumber = 1\n")
        self.assertEqual(7, len(ped.pedigree))

    def test_renumber_zero_actually_disables_renumbering(self):
        """
        The sharpest consequence of defect 3. As a string, '0' is true, so
        `renumber = 0` in an option file switched renumbering ON.
        """
        ped = self._load(
            "messages = quiet\npedigree_summary = 0\n"
            "pedfile = new_lacy.ped\npedformat = asd\nrenumber = 0\n")
        self.assertEqual(0, ped.kw["renumber"])
        self.assertFalse(ped.kw["renumber"])

    def test_explicit_options_still_take_precedence(self):
        """Passing an options dict must not start reading the file instead."""
        tmp, path = self._write(
            "messages = quiet\npedfile = new_lacy.ped\npedformat = asd\n")
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            ped = loadPedigree(
                options={"pedfile": os.path.join(tmp, "new_lacy.ped"),
                         "pedformat": "asd", "sepchar": " ",
                         "messages": "quiet", "pedigree_summary": 0,
                         "pedname": "explicit"},
                optionsfile=path)
        finally:
            os.chdir(cwd)
        self.assertEqual("explicit", ped.kw["pedname"])



class TestOutputDisciplineAndImplicitParents(unittest.TestCase):
    """
    Findings 16 and 19.

    16: ``messages='quiet'`` did not suppress the ``[NOTE]:`` lines emitted
    while adding implicit parents. A historical large Griffon extract printed
    7,680 of them, which breaks any machine-readable use of stdout and dominated
    the load time.

    19: PyPedal silently adds an animal for any sire or dam named as a parent
    without a record of its own, so the post-load count exceeds the input row
    count. That was discoverable only by reading those same NOTE lines, and
    there was no programmatic way to reconcile the two.

    horse.ped is the fixture: 10 data rows, 16 animals after load.
    """

    def _load(self, **overrides):
        import contextlib
        import io

        from _pedhelpers import load_corpus

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            ped = load_corpus("horse.ped", "ASD", sepchar=",", **overrides)
        notes = [line for line in buffer.getvalue().splitlines()
                 if "[NOTE]" in line]
        return ped, notes

    def test_quiet_suppresses_implicit_parent_notes(self):
        _ped, notes = self._load(messages="quiet")
        self.assertEqual([], notes)

    def test_verbose_still_reports_them(self):
        """Suppression must be a setting, not a removal of the diagnostic."""
        _ped, notes = self._load(messages="verbose")
        self.assertGreater(len(notes), 0)

    def test_implicit_parents_are_counted_in_metadata(self):
        ped, _notes = self._load(messages="quiet")
        self.assertEqual(6, ped.metadata.num_implicit_parents)
        self.assertEqual(6, len(ped.metadata.implicit_parent_ids))

    def test_input_rows_reconcile_with_the_record_count(self):
        """The property  existed to make checkable."""
        ped, _notes = self._load(messages="quiet")
        with open(corpus("horse.ped"), "r", encoding="utf-8") as handle:
            data_rows = [line for line in handle
                         if line.strip() and not line.strip().startswith("#")]
        self.assertEqual(len(data_rows) + ped.metadata.num_implicit_parents,
                         ped.metadata.num_records)

    def test_pedigree_without_implicit_parents_reports_zero(self):
        from _pedhelpers import load_corpus

        ped = load_corpus("mrode.ped", "asd")
        self.assertEqual(0, ped.metadata.num_implicit_parents)
        self.assertEqual([], ped.metadata.implicit_parent_ids)


if __name__ == "__main__":
    unittest.main()
