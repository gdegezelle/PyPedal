"""Desktop-app pedigree open options, including CSV separator padding."""

import os
import tempfile
import unittest
from pathlib import Path

from _pedhelpers import close_owned_pypedal_log_handlers

from PyPedal.application import PedigreeOpenOptions, normalize_sepchar
from PyPedal.pyp_errors import PyPedalPedigreeFormatError
from PyPedal.pyp_newclasses import load_pedigree


class TestNormalizeSepchar(unittest.TestCase):
    def test_empty_and_spaces_mean_a_space(self):
        self.assertEqual(" ", normalize_sepchar(""))
        self.assertEqual(" ", normalize_sepchar(" "))
        self.assertEqual(" ", normalize_sepchar("  "))
        self.assertEqual(" ", normalize_sepchar(None))

    def test_comma_with_leftover_spaces_is_a_comma(self):
        """The GUI default was a space; typing a comma left ', '."""
        self.assertEqual(",", normalize_sepchar(","))
        self.assertEqual(",", normalize_sepchar(", "))
        self.assertEqual(",", normalize_sepchar(" ,"))
        self.assertEqual(",", normalize_sepchar(" , "))

    def test_tab_is_kept(self):
        self.assertEqual("\t", normalize_sepchar("\t"))


class TestPedigreeOpenOptions(unittest.TestCase):
    def test_csv_open_options_use_a_bare_comma(self):
        opts = PedigreeOpenOptions(
            pedformat="asdxbn", separator=", ", renumber=True
        ).to_library_options(Path("/tmp/dogs.ped"))
        self.assertEqual(",", opts["sepchar"])
        self.assertEqual("asdxbn", opts["pedformat"])
        self.assertTrue(opts["renumber"])
        self.assertEqual("dogs.ped", opts["pedname"])
        self.assertEqual("quiet", opts["messages"])
        self.assertEqual(0, opts["pedigree_summary"])


class TestAsdxbnCsvLoad(unittest.TestCase):
    """The Griffon export shape: comma CSV, MMDDYYYY, names with spaces."""

    ROWS = (
        "1,30497,52843,f,03132018,A Day Before Sunrise de Mar&Mar\n"
        "2,12401,68419,f,01111995,A Galaxie Mii Jimijo Nubegin\n"
    )

    def _load(self, sepchar):
        fd, path = tempfile.mkstemp(suffix=".ped")
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(self.ROWS)
            return load_pedigree(
                options={
                    "pedfile": path,
                    "pedformat": "asdxbn",
                    "sepchar": sepchar,
                    "renumber": True,
                    "messages": "quiet",
                    "pedigree_summary": 0,
                }
            )
        finally:
            close_owned_pypedal_log_handlers()
            os.remove(path)
            log = path[:-4] + ".log"
            if os.path.exists(log):
                os.remove(log)

    def test_bare_comma_loads_names_with_spaces(self):
        ped = self._load(",")
        named = [a for a in ped.pedigree if a.originalID == 1]
        self.assertEqual(1, len(named))
        self.assertEqual("A Day Before Sunrise de Mar&Mar", named[0].name)

    def test_comma_space_is_the_one_column_failure(self):
        with self.assertRaises(PyPedalPedigreeFormatError) as raised:
            self._load(", ")
        message = str(raised.exception)
        self.assertIn("1 columns", message)
        self.assertIn("separator ', '", message)

    def test_gui_options_load_the_padded_comma(self):
        fd, path = tempfile.mkstemp(suffix=".ped")
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(self.ROWS)
            opts = PedigreeOpenOptions(
                pedformat="asdxbn", separator=", ", renumber=True
            ).to_library_options(Path(path))
            ped = load_pedigree(options=opts)
            named = [a for a in ped.pedigree if a.originalID == 1]
            self.assertEqual("A Day Before Sunrise de Mar&Mar", named[0].name)
        finally:
            close_owned_pypedal_log_handlers()
            os.remove(path)
            log = path[:-4] + ".log"
            if os.path.exists(log):
                os.remove(log)


class TestSimpleSeparatorLoads(unittest.TestCase):
    """Space, tab, and bare-comma loads without a display server."""

    def _load(self, text, sepchar):
        fd, path = tempfile.mkstemp(suffix=".ped")
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            return load_pedigree(
                options={
                    "pedfile": path,
                    "pedformat": "asd",
                    "sepchar": sepchar,
                    "renumber": True,
                    "messages": "quiet",
                    "pedigree_summary": 0,
                }
            )
        finally:
            close_owned_pypedal_log_handlers()
            os.remove(path)
            log = path[:-4] + ".log"
            if os.path.exists(log):
                os.remove(log)

    def test_space_separated_loads(self):
        ped = self._load("1 0 0\n2 0 0\n3 1 2\n", " ")
        self.assertEqual(3, len(ped.pedigree))
        ids = {int(a.originalID) for a in ped.pedigree}
        self.assertEqual({1, 2, 3}, ids)

    def test_tab_separated_loads(self):
        ped = self._load("1\t0\t0\n2\t0\t0\n3\t1\t2\n", "\t")
        self.assertEqual(3, len(ped.pedigree))

    def test_bare_comma_without_spaces_loads(self):
        ped = self._load("1,0,0\n2,0,0\n3,1,2\n", ",")
        self.assertEqual(3, len(ped.pedigree))

    def test_gui_empty_separator_loads_space_file(self):
        fd, path = tempfile.mkstemp(suffix=".ped")
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("1 0 0\n2 0 0\n3 1 2\n")
            opts = PedigreeOpenOptions(
                pedformat="asd", separator="", renumber=True
            ).to_library_options(Path(path))
            self.assertEqual(" ", opts["sepchar"])
            ped = load_pedigree(options=opts)
            self.assertEqual(3, len(ped.pedigree))
        finally:
            close_owned_pypedal_log_handlers()
            os.remove(path)
            log = path[:-4] + ".log"
            if os.path.exists(log):
                os.remove(log)
