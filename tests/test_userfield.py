"""
``NewAnimal.userField`` and the ``u`` pedformat code ().

PyPedal 2.0.4 assigns the attribute unconditionally in ``NewAnimal.__init__``
(``pyp_newclasses.py:3234-3238``): the pedigree column when ``locations['userfield']``
is not the ``-999`` sentinel, otherwise ``kw['missing_userfield']``, stripped. The
Python 3 port kept the whole surrounding machinery -- the ``u`` code is accepted by
the format validator, ``load()`` computes the column position, and the default option
is declared -- and dropped only the assignment. Loading any pedigree with ``u`` raised
``AttributeError`` from ``PedigreeMetadata.nufields()``, the first consumer reached.

WHY THIS FILE IS LARGE
----------------------
The repair is five lines. Its consequence is that ten direct and five indirect
consumer sites go live for the first time in the port, none of which had any test.
One class below per consumer group, so a regression names the consumer it broke.

THE DEFAULT
-----------
``kw['missing_userfield']`` is ``'Unknown'``, not ``''``. PyPedal 2.0.4 sets that key
twice in ``NewPedigree.__init__`` -- ``'Unknown'`` at ``:181-182``, which runs first
and wins, and ``''`` at ``:283-284``, which is guarded by ``not in kw.keys()`` and is
therefore dead code. The port had transcribed the dead branch. Since no other corpus
pedigree carries a ``u`` column, every animal in the differential harness takes this
default, so the choice is measured directly against 2.0.4 on all nine pedigrees.
"""
import os
import sqlite3
import unittest

import pytest

from PyPedal import pyp_db, pyp_io, pyp_utils
from PyPedal.pyp_newclasses import load_pedigree

from _pedhelpers import CORPUS, chdir_tmp, load_corpus

# The fixture, restated here so a test failure says what it expected rather than
# making the reader open the .ped file.
EXPECTED = {1: "hatch", 2: "hatch", 3: "other", 4: "other", 5: "third"}
UNIQUE = {"hatch", "other", "third"}


def _userfield_ped(**overrides):
    return load_corpus("userfield.ped", "asdu", **overrides)


def _plain_ped(**overrides):
    """A pedigree with no 'u' column, so every animal takes the default."""
    return load_corpus("mrode.ped", "asd", **overrides)


class TestLoadsAtAll(unittest.TestCase):
    """The literal reproduction from the algorithm notes."""

    def test_a_pedigree_with_a_u_column_loads(self):
        ped = _userfield_ped()
        self.assertEqual(5, len(ped.pedigree))

    def test_every_animal_has_the_attribute(self):
        for animal in _userfield_ped().pedigree:
            self.assertTrue(
                hasattr(animal, "userField"),
                f"animal {animal.originalID} has no 'userField'; PyPedal 2.0.4 "
                "sets it at pyp_newclasses.py:3235")


class TestTheValueItself(unittest.TestCase):

    def test_userfield_matches_the_column(self):
        for animal in _userfield_ped().pedigree:
            self.assertEqual(EXPECTED[animal.originalID], animal.userField)

    def test_the_value_is_stripped(self):
        """
        2.0.4 applies string.strip() to the column; the port's equivalent is the
        local safe_strip() helper, which every sibling field already uses.

        Comma-separated, because whitespace padding cannot be expressed in a
        space-separated file -- the loader would read it as extra columns and
        reject the record against the pedformat.
        """
        with chdir_tmp() as tmp:
            path = os.path.join(tmp, "padded.ped")
            with open(path, "w") as handle:
                handle.write("1,0,0,  spaced  \n2,0,0,  spaced  \n3,1,2,tight\n")
            ped = load_pedigree({
                "pedfile": path, "pedformat": "asdu", "sepchar": ",",
                "messages": "quiet", "pedigree_summary": 0, "renumber": True,
            })
        self.assertEqual(
            ["spaced", "spaced", "tight"],
            [a.userField for a in ped.pedigree])


class TestTheDefault(unittest.TestCase):

    def test_a_pedigree_without_u_gets_the_missing_value(self):
        ped = _plain_ped()
        self.assertEqual("Unknown", ped.kw["missing_userfield"])
        for animal in ped.pedigree:
            self.assertEqual("Unknown", animal.userField)

    def test_an_explicit_override_is_honoured(self):
        ped = _plain_ped(missing_userfield="NOTHING")
        for animal in ped.pedigree:
            self.assertEqual("NOTHING", animal.userField)

    def test_the_default_is_not_applied_when_the_column_is_present(self):
        """The override must not leak into pedigrees that really have the column."""
        ped = _userfield_ped(missing_userfield="NOTHING")
        for animal in ped.pedigree:
            self.assertEqual(EXPECTED[animal.originalID], animal.userField)


class TestConsumerMetadataNufields(unittest.TestCase):
    """Consumer 1: PedigreeMetadata.nufields(), pyp_newclasses.py:3053."""

    def test_nufields_counts_unique_values(self):
        """
        nufields() is a construction-time method: PedigreeMetadata.__init__ ends
        with ``self.myped = []`` ("Detaching pedigree...", pyp_newclasses.py:2983),
        so calling it afterwards sees an empty pedigree and returns (0, set()).
        That detach is inherited verbatim from PyPedal 2.0.4 (:3835-3836) and
        applies to every nu*() sibling equally, so it is faithful behaviour rather
        than a port defect -- but it means the method has to be handed a pedigree
        to be exercised directly.
        """
        ped = _userfield_ped()
        md = ped.metadata

        self.assertEqual((0, set()), md.nufields(), "the detach no longer happens")

        md.myped = ped.pedigree
        try:
            count, fields = md.nufields()
        finally:
            md.myped = []
        self.assertEqual(UNIQUE, set(fields))
        self.assertEqual(3, count)

    def test_metadata_attributes_are_populated(self):
        md = _userfield_ped().metadata
        self.assertEqual(3, md.num_unique_fields)
        self.assertEqual(UNIQUE, set(md.unique_field_list))

    def test_stringme_reports_the_count(self):
        md = _userfield_ped().metadata
        self.assertIn("Unique userFields", md.stringme())


class TestConsumerWriters(unittest.TestCase):
    """Consumers 2, 5 and 6: oldsave(), save(), save_newanimals_to_file()."""

    def test_oldsave_writes_the_real_value(self):
        ped = _userfield_ped()
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "old.ped")
            self.assertTrue(ped.oldsave(filename=out))
            body = [ln for ln in open(out) if not ln.startswith("#")]
        self.assertEqual(5, len(body))
        for line in body:
            self.assertTrue(
                any(line.rstrip().endswith(v) for v in UNIQUE),
                f"oldsave wrote no userField on: {line!r}")

    def test_save_writes_the_real_value_not_an_empty_column(self):
        """
        save() reads userField indirectly, through new_animal_attr['u'], with a
        getattr default of ''. Before the repair it silently wrote an empty
        column rather than raising -- the failure mode this asserts against.
        """
        ped = _userfield_ped()
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "saved.ped")
            self.assertTrue(ped.save(filename=out, pedformat="asdu"))
            body = [ln for ln in open(out) if not ln.startswith("#")]
        written = [ln.split()[3] for ln in body if ln.strip()]
        self.assertEqual([EXPECTED[i] for i in range(1, 6)], written)

    def test_save_newanimals_to_file_writes_the_real_value(self):
        ped = _userfield_ped()
        kw = dict(ped.kw)
        kw["pedformat"] = "asdu"
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "animals.ped")
            self.assertTrue(pyp_io.save_newanimals_to_file(
                ped.pedigree, out, kw, ped.new_animal_attr))
            lines = [ln for ln in open(out) if ln.strip()]
        written = [ln.split()[3] for ln in lines]
        self.assertEqual([EXPECTED[i] for i in range(1, 6)], written)


class TestConsumerMatchRules(unittest.TestCase):
    """
    Consumers 3 and 4: __sub__ and intersection, pyp_newclasses.py:518 and :631.

    Both resolve 'u' through new_animal_attr and call getattr with **no**
    default, so a match_rule containing 'u' raised AttributeError outright.
    """

    def test_subtraction_with_a_u_match_rule_completes(self):
        """
        What is under test is that the 'u' branch resolves at all. __sub__ returns
        a freshly loaded NewPedigree on success and False on failure, so a
        NewPedigree back is the evidence the getattr no longer raises.
        """
        from PyPedal.pyp_newclasses import NewPedigree

        with chdir_tmp() as tmp:
            a = _userfield_ped(match_rule="au")
            b = _userfield_ped(match_rule="au")
            result = a.__sub__(b, filename=os.path.join(tmp, "diff.ped"))
        self.assertIsInstance(result, NewPedigree)

    def test_intersection_with_a_u_match_rule_completes(self):
        with chdir_tmp():
            a = _userfield_ped(match_rule="au")
            b = _userfield_ped(match_rule="au")
            result = a.intersection(b, newpedname="intersected")
        self.assertIsNotNone(result)
        self.assertEqual(5, len(result.pedigree))


class TestConsumerGedcom(unittest.TestCase):
    """Consumer 7: pyp_io.save_to_gedcom, pyp_io.py:891."""

    def test_gedcom_export_emits_the_userfield_as_the_name(self):
        ped = _userfield_ped()
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "out.ged")
            pyp_io.save_to_gedcom(ped, out)
            text = open(out).read()
        for value in UNIQUE:
            self.assertIn(f"1 NAME {value}", text)


class TestConsumerGuessPedformat(unittest.TestCase):
    """
    Consumer 8: pyp_utils.guess_pedformat, pyp_utils.py:1411.

    This is the site that pins the 'Unknown' default end-to-end: it decides
    whether to emit a 'u' code by comparing userField against
    kw['missing_userfield'], so a mismatched default silently mislabels every
    pedigree as carrying a user field.
    """

    def test_a_pedigree_with_u_is_guessed_with_u(self):
        ped = _userfield_ped()
        guessed = pyp_utils.guess_pedformat(ped.pedigree[0], ped.kw)
        self.assertIn("u", guessed)

    def test_a_pedigree_without_u_is_guessed_without_u(self):
        ped = _plain_ped()
        guessed = pyp_utils.guess_pedformat(ped.pedigree[0], ped.kw)
        self.assertNotIn("u", guessed)


class TestConsumerSqlite(unittest.TestCase):
    """
    Consumer 9: pyp_db.populate_pedigree_table, pyp_db.py:218.

    The same INSERT reads record.ancestor at :214, which is  -- this
    round-trip is the end-to-end test that needed both repairs.
    """

    def test_pedigree_round_trips_through_sqlite(self):
        with chdir_tmp() as tmp:
            dbfile = os.path.join(tmp, "userfield.db")
            ped = _userfield_ped(database_file=dbfile, database_table="uf")

            conn = pyp_db.connect_to_database(ped)
            self.assertIsNotNone(conn)
            try:
                self.assertTrue(pyp_db.create_pedigree_table(ped, conn=conn, drop=True))
                self.assertTrue(pyp_db.populate_pedigree_table(ped, conn=conn))
                rows = conn.execute(
                    "SELECT animalID, userField, ancestor FROM uf ORDER BY animalID"
                ).fetchall()
            finally:
                conn.close()

        self.assertEqual(5, len(rows))
        for animal_id, user_field, ancestor in rows:
            self.assertEqual(EXPECTED[animal_id], user_field)
            # `ancestor` is declared TEXT in the schema (pyp_db.py:63), so SQLite
            # returns the integer back as a string.
            self.assertIn(ancestor, ("0", "1"))
        # Animals 1-4 are parents in this fixture, 5 is not. Asserting the shape
        # rather than just "not null" is what makes this a  test too.
        self.assertEqual("0", dict((r[0], r[2]) for r in rows)[5])


class TestConsumerGraphics(unittest.TestCase):
    """Consumer 10: pyp_graphics colorByUser palette, pyp_graphics.py:1310/:1314."""

    def test_color_by_user_reaches_the_palette_lookup(self):
        pytest.importorskip("matplotlib")
        pytest.importorskip("pydot")
        # No webcolors gate: the colorByUser path used to funnel every node
        # colour through get_colour_name(), which needed a dependency no extra
        # declared and an API no installable version still provides. It now
        # emits the colormap RGB directly, so the declared `graphics` extra is
        # sufficient. See TestColorByUserEmitsColormapRGB in
        # tests/test_graphics_colorbyuser.py.
        from PyPedal import pyp_graphics

        ped = _userfield_ped()
        with chdir_tmp() as tmp:
            # What matters is that the userField-keyed lookup resolves rather
            # than raising; whether a PNG is produced depends on Graphviz being
            # installed, which is a separate capability.
            pyp_graphics.new_draw_pedigree(
                ped, gfilename=os.path.join(tmp, "cbu"), gtitle="colorByUser",
                colorByUser=True)


if __name__ == "__main__":
    unittest.main()
