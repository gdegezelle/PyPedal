"""
 /  -- ``NewPedigree.oldsave()`` must not invent a parent.

WHAT THIS FILE IS ABOUT
-----------------------
``oldsave(..., idformat='o')`` resolves each parent as
``pedigree[int(sireID) - 1]`` / ``pedigree[int(damID) - 1]`` without consulting
``kw['missing_parent']``. The shipped sentinel is ``0``, so a missing parent
becomes ``pedigree[-1]`` -- the last real animal -- and is written under that
animal's ``originalID``. A custom sentinel typically ``IndexError``s; the
method swallows it and returns ``False``.

``save(originalID=True)`` already guards the same lookup. ``oldsave`` is a
legacy format (renamed off ``save()`` on 09/20/2010) and is not rewritten as
``save()``. See ``the algorithm notes``.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with every test
asserting the repaired contract marked ``xfail(strict=True)``. Those markers
were removed in the commit that implemented the repair. The characterisation
tests that pinned the *broken* bytes were kept rather than deleted, restated
as ``TestTheReproducerNoLongerReproduces``.

``self.subTest`` was kept out of the tests that were xfailed.

WHAT IS *NOT* CLAIMED HERE
--------------------------
Nothing here changes ``save()``, reorder, renumber, or  birth-year
semantics. ``oldsave`` always appends a birth-year column; Python 2 writes
``1900`` and Python 3 writes ``1800``. That offset is . Parent
columns are 1 and 2 (0-based) regardless.
"""
import os
import tempfile
import unittest

from PyPedal import pyp_errors
from _pedhelpers import chdir_tmp, load_corpus_from_path


# Original IDs are deliberately not 1..n so Finding-37 renumbering cannot be
# mistaken for the file contract. 999 is last in memory after a normal load.
FOUNDER = "10"
HALF_SIRE = "20"
HALF_DAM = "30"
ORDINARY = "40"
LAST = "999"

DEFAULT_ROWS = [
    "%s 0 0" % FOUNDER,
    "%s %s 0" % (HALF_SIRE, FOUNDER),
    "%s 0 %s" % (HALF_DAM, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (LAST, ORDINARY, HALF_DAM),
]

CUSTOM_SENTINEL = -999
CUSTOM_ROWS = [
    "%s %s %s" % (FOUNDER, CUSTOM_SENTINEL, CUSTOM_SENTINEL),
    "%s %s %s" % (HALF_SIRE, FOUNDER, CUSTOM_SENTINEL),
    "%s %s %s" % (HALF_DAM, CUSTOM_SENTINEL, FOUNDER),
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (LAST, ORDINARY, HALF_DAM),
]


def rows_to_ped(rows, pedformat="asd", **overrides):
    tmp = tempfile.mkdtemp(prefix="oldsave38_")
    path = os.path.join(tmp, "high1.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, **overrides)


def parse_oldsave(path):
    """
    Return the animal/sire/dam triples from an oldsave file.

    oldsave always writes animal, sire, dam as columns 0, 1, 2, then extra
    fields (birth year at least). Comment lines are ignored.
    """
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 3:
                raise AssertionError(
                    "oldsave row has fewer than 3 columns: %r" % line)
            rows.append((cols[0], cols[1], cols[2]))
    return rows


def written_as(rows):
    return {animal: (sire, dam) for animal, sire, dam in rows}


def snapshot(pedobj):
    return {
        "order": [a.originalID for a in pedobj.pedigree],
        "animalID": [a.animalID for a in pedobj.pedigree],
        "originalID": [a.originalID for a in pedobj.pedigree],
        "sireID": [a.sireID for a in pedobj.pedigree],
        "damID": [a.damID for a in pedobj.pedigree],
        "missing_parent": pedobj.kw["missing_parent"],
        "pedformat": pedobj.kw["pedformat"],
    }


def oldsave_triples(pedobj, **kwargs):
    with chdir_tmp() as tmp:
        out = os.path.join(tmp, "old.ped")
        ok = pedobj.oldsave(filename=out, **kwargs)
        rows = parse_oldsave(out) if os.path.exists(out) else []
    return ok, rows


# ===========================================================================
# Anti-vacuity. These must stay green on the broken tree AND the repaired
# tree: they prove the fixture can detect the wrap, not that the wrap is gone.
# ===========================================================================
class TestFixtureIsNotVacuous(unittest.TestCase):
    """The reproducer can actually see pedigree[-1] corruption."""

    def test_the_last_animal_is_not_the_missing_parent_sentinel(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        last = ped.pedigree[-1]
        missing = ped.kw["missing_parent"]
        self.assertEqual(int(LAST), int(last.originalID))
        self.assertNotEqual(int(last.originalID), int(missing))
        self.assertNotEqual(int(last.animalID), int(missing))

    def test_both_sire_and_dam_missing_slots_exist(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = ped.kw["missing_parent"]
        sire_missing = [a for a in ped.pedigree if a.sireID == missing]
        dam_missing = [a for a in ped.pedigree if a.damID == missing]
        self.assertGreaterEqual(len(sire_missing), 2)
        self.assertGreaterEqual(len(dam_missing), 2)
        # Founder (both), half-sire (dam only), half-dam (sire only).
        kinds = {
            (int(a.originalID), a.sireID == missing, a.damID == missing)
            for a in ped.pedigree
        }
        self.assertIn((int(FOUNDER), True, True), kinds)
        self.assertIn((int(HALF_SIRE), False, True), kinds)
        self.assertIn((int(HALF_DAM), True, False), kinds)
        self.assertIn((int(ORDINARY), False, False), kinds)

    def test_columns_1_and_2_are_the_parent_columns(self):
        """
        Animal 40 has two known parents. Whatever oldsave writes for that row
        in columns 1 and 2 must be those parents' original IDs -- otherwise
        we would be asserting about the wrong columns.
        """
        ped = rows_to_ped(DEFAULT_ROWS)
        _ok, rows = oldsave_triples(ped)
        by_id = written_as(rows)
        self.assertIn(ORDINARY, by_id)
        self.assertEqual((FOUNDER, HALF_SIRE), by_id[ORDINARY])

    def test_internal_ids_differ_from_original_ids_after_renumber(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        originals = [int(a.originalID) for a in ped.pedigree]
        internals = [int(a.animalID) for a in ped.pedigree]
        self.assertEqual([1, 2, 3, 4, 5], internals)
        self.assertEqual([10, 20, 30, 40, 999], originals)
        self.assertNotEqual(originals, internals)


# ===========================================================================
# THE REPRODUCER, INVERTED.
#
# On the phase baseline 2940f53 each of these asserted the BROKEN bytes and
# passed. They are kept rather than deleted -- a reproducer thrown away once
# it stops reproducing leaves nothing to stop the defect returning.
# ===========================================================================
class TestTheReproducerNoLongerReproduces(unittest.TestCase):
    """The Finding-38 reproducer, now inverted. Baseline: 2940f53."""

    def test_default_oldsave_no_longer_writes_the_last_animal_into_missing_slots(self):
        """Was: 10 999 999 / 20 10 999 / 30 999 10."""
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = str(ped.kw["missing_parent"])
        ok, rows = oldsave_triples(ped)
        self.assertTrue(ok)
        self.assertEqual(
            [
                (FOUNDER, missing, missing),
                (HALF_SIRE, FOUNDER, missing),
                (HALF_DAM, missing, FOUNDER),
                (ORDINARY, FOUNDER, HALF_SIRE),
                (LAST, ORDINARY, HALF_DAM),
            ],
            rows,
        )
        for animal, sire, dam in rows:
            if animal != LAST:
                self.assertNotEqual(LAST, sire)
                self.assertNotEqual(LAST, dam)

    def test_custom_sentinel_no_longer_makes_oldsave_return_false(self):
        """Was: return False and a header-only file (IndexError swallowed)."""
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        token = str(CUSTOM_SENTINEL)
        ok, rows = oldsave_triples(ped)
        self.assertTrue(ok)
        self.assertEqual((token, token), written_as(rows)[FOUNDER])


# ===========================================================================
# Desired contract. xfailed where current production fails it.
# ===========================================================================
class TestDesiredContract(unittest.TestCase):
    """OS-1 .. OS-7 and the unresolved-parent error."""

    def test_os1_both_parents_missing_stay_missing(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = str(ped.kw["missing_parent"])
        _ok, rows = oldsave_triples(ped)
        sire, dam = written_as(rows)[FOUNDER]
        self.assertEqual(missing, sire)
        self.assertEqual(missing, dam)

    def test_os2_known_sire_missing_dam(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = str(ped.kw["missing_parent"])
        _ok, rows = oldsave_triples(ped)
        sire, dam = written_as(rows)[HALF_SIRE]
        self.assertEqual(FOUNDER, sire)
        self.assertEqual(missing, dam)

    def test_os3_missing_sire_known_dam(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = str(ped.kw["missing_parent"])
        _ok, rows = oldsave_triples(ped)
        sire, dam = written_as(rows)[HALF_DAM]
        self.assertEqual(missing, sire)
        self.assertEqual(FOUNDER, dam)

    def test_os4_both_parents_known_are_preserved(self):
        """Green today -- this row never indexed pedigree[-1]."""
        ped = rows_to_ped(DEFAULT_ROWS)
        _ok, rows = oldsave_triples(ped)
        self.assertEqual((FOUNDER, HALF_SIRE), written_as(rows)[ORDINARY])
        self.assertEqual((ORDINARY, HALF_DAM), written_as(rows)[LAST])

    def test_os5_last_animal_never_appears_as_a_missing_parent(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        missing = str(ped.kw["missing_parent"])
        _ok, rows = oldsave_triples(ped)
        for animal, sire, dam in rows:
            if animal == LAST:
                continue
            self.assertNotEqual(
                LAST, sire,
                "%s wrote the last animal as sire" % animal)
            self.assertNotEqual(
                LAST, dam,
                "%s wrote the last animal as dam" % animal)
        founder_sire, founder_dam = written_as(rows)[FOUNDER]
        self.assertEqual(missing, founder_sire)
        self.assertEqual(missing, founder_dam)

    def test_os6_custom_missing_parent_is_honoured(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        token = str(CUSTOM_SENTINEL)
        ok, rows = oldsave_triples(ped)
        self.assertTrue(ok)
        by_id = written_as(rows)
        self.assertEqual((token, token), by_id[FOUNDER])
        self.assertEqual((FOUNDER, token), by_id[HALF_SIRE])
        self.assertEqual((token, FOUNDER), by_id[HALF_DAM])
        self.assertEqual((FOUNDER, HALF_SIRE), by_id[ORDINARY])
        self.assertEqual((ORDINARY, HALF_DAM), by_id[LAST])

    def test_os7_stable_identity_after_renumber(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        internals = [int(a.animalID) for a in ped.pedigree]
        originals = [int(a.originalID) for a in ped.pedigree]
        self.assertNotEqual(internals, originals)
        missing = str(ped.kw["missing_parent"])
        _ok, rows = oldsave_triples(ped)
        self.assertEqual(
            {
                FOUNDER: (missing, missing),
                HALF_SIRE: (FOUNDER, missing),
                HALF_DAM: (missing, FOUNDER),
                ORDINARY: (FOUNDER, HALF_SIRE),
                LAST: (ORDINARY, HALF_DAM),
            },
            written_as(rows),
        )
        # The file must not emit internal animalIDs in the original-ID path.
        emitted = {c for row in rows for c in row}
        self.assertTrue(emitted.isdisjoint({"1", "2", "3", "4", "5"}))

    def test_unresolved_known_parent_raises_structure_error(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        victim = ped.pedigree[3]  # original 40; both parents currently known
        victim.sireID = 99
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "old.ped")
            with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
                ped.oldsave(filename=out)


# ===========================================================================
# Already-correct neighbours. Not xfailed.
# ===========================================================================
class TestOldsaveDoesNotMutate(unittest.TestCase):
    """OS-8."""

    def test_os8_oldsave_mutates_neither_animals_nor_configuration(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        before = snapshot(ped)
        oldsave_triples(ped)
        self.assertEqual(before, snapshot(ped))

    def test_os8_custom_sentinel_path_also_mutates_nothing(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        before = snapshot(ped)
        oldsave_triples(ped)
        self.assertEqual(before, snapshot(ped))


class TestIdformatRAlreadyCorrect(unittest.TestCase):
    """The renumbered-ID path writes sireID/damID directly and never indexes."""

    def test_renumbered_ids_keep_the_in_memory_sentinel(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        ok, rows = oldsave_triples(ped, idformat="r")
        self.assertTrue(ok)
        by_id = written_as(rows)
        # animalID 1 is original 10, founders carry sentinel 0.
        self.assertEqual(("0", "0"), by_id["1"])
        self.assertEqual(("1", "0"), by_id["2"])
        self.assertEqual(("0", "1"), by_id["3"])
        self.assertEqual(("1", "2"), by_id["4"])
        self.assertEqual(("4", "3"), by_id["5"])

    def test_renumbered_ids_honour_a_custom_sentinel(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        ok, rows = oldsave_triples(ped, idformat="r")
        self.assertTrue(ok)
        token = str(CUSTOM_SENTINEL)
        by_id = written_as(rows)
        self.assertEqual((token, token), by_id["1"])
        self.assertEqual(("1", token), by_id["2"])
        self.assertEqual((token, "1"), by_id["3"])


class TestSaveAlreadyGuards(unittest.TestCase):
    """
    save(originalID=True) already checks missing_parent.  does not
    change it. Measured so a later edit of save() cannot be blamed on this
    phase, and so the two contracts are not silently conflated.

    save() writes hardcoded 0 for a missing parent, even when the configured
    sentinel is not 0. That is save()'s contract, not oldsave's.
    """

    def test_save_original_ids_write_zero_for_a_missing_parent(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "saved.ped")
            self.assertTrue(ped.save(filename=out, pedformat="asd",
                                     originalID=True))
            rows = parse_oldsave(out)
        self.assertEqual(
            [
                (FOUNDER, "0", "0"),
                (HALF_SIRE, FOUNDER, "0"),
                (HALF_DAM, "0", FOUNDER),
                (ORDINARY, FOUNDER, HALF_SIRE),
                (LAST, ORDINARY, HALF_DAM),
            ],
            rows,
        )

    def test_save_writes_hardcoded_zero_even_for_a_custom_sentinel(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "saved.ped")
            self.assertTrue(ped.save(filename=out, pedformat="asd",
                                     originalID=True))
            rows = parse_oldsave(out)
        self.assertEqual((FOUNDER, "0", "0"), rows[0])
        self.assertEqual((HALF_SIRE, FOUNDER, "0"), rows[1])
        self.assertEqual((HALF_DAM, "0", FOUNDER), rows[2])


def memory_graph(pedobj):
    """Stable-identity parentage: animals, sire edges, dam edges, missing slots."""
    missing = pedobj.kw["missing_parent"]
    by_id = {int(a.animalID): str(a.originalID) for a in pedobj.pedigree}
    animals = set(by_id.values())
    sire_edges, dam_edges, missing_slots = set(), set(), set()
    for animal in pedobj.pedigree:
        aid = str(animal.originalID)
        if animal.sireID == missing:
            missing_slots.add((aid, "sire"))
        else:
            sire_edges.add((aid, by_id[int(animal.sireID)]))
        if animal.damID == missing:
            missing_slots.add((aid, "dam"))
        else:
            dam_edges.add((aid, by_id[int(animal.damID)]))
    return animals, sire_edges, dam_edges, missing_slots


def file_graph(triples, missing_token):
    animals = {animal for animal, _s, _d in triples}
    sire_edges, dam_edges, missing_slots = set(), set(), set()
    token = str(missing_token)
    for animal, sire, dam in triples:
        if sire == token:
            missing_slots.add((animal, "sire"))
        else:
            sire_edges.add((animal, sire))
        if dam == token:
            missing_slots.add((animal, "dam"))
        else:
            dam_edges.add((animal, dam))
    return animals, sire_edges, dam_edges, missing_slots


class TestStructuralOracle(unittest.TestCase):
    """
    Input parentage and oldsave parentage, both in stable identity.

    Requires every known edge preserved, every missing slot still missing,
    no new edge, no sire/dam swap. Not a byte-identity check: oldsave
    always appends a birth year and that column is 's.
    """

    def test_default_sentinel_round_trips_structurally(self):
        ped = rows_to_ped(DEFAULT_ROWS)
        _ok, rows = oldsave_triples(ped)
        self.assertEqual(
            memory_graph(ped),
            file_graph(rows, ped.kw["missing_parent"]),
        )

    def test_custom_sentinel_round_trips_structurally(self):
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM_SENTINEL)
        _ok, rows = oldsave_triples(ped)
        self.assertEqual(
            memory_graph(ped),
            file_graph(rows, CUSTOM_SENTINEL),
        )

    def test_the_oracle_would_reject_the_pre_repair_wrap(self):
        """Anti-vacuity: the corrupt 2940f53 bytes fail the oracle."""
        ped = rows_to_ped(DEFAULT_ROWS)
        corrupt = [
            (FOUNDER, LAST, LAST),
            (HALF_SIRE, FOUNDER, LAST),
            (HALF_DAM, LAST, FOUNDER),
            (ORDINARY, FOUNDER, HALF_SIRE),
            (LAST, ORDINARY, HALF_DAM),
        ]
        self.assertNotEqual(
            memory_graph(ped),
            file_graph(corrupt, ped.kw["missing_parent"]),
        )


class TestFastReorderMissingparentClaim(unittest.TestCase):
    """
     pin, inverted where the wrappers were repaired.

    The public ``fast_reorder(..., missingparent=0)`` default is unchanged:
    omitting the argument with sentinel -999 still raises. The
    ``renumber()`` asdb wrapper now passes ``kw['missing_parent']``.
    """

    def test_the_signature_accepts_missingparent_defaulting_to_zero(self):
        from PyPedal import pyp_utils
        params = pyp_utils.fast_reorder.__defaults__
        self.assertEqual(("_new_reordered_", "no", False, 0), params)

    def test_omitting_missingparent_treats_a_custom_sentinel_as_a_real_id(self):
        from PyPedal import pyp_utils
        ped = rows_to_ped(
            ["1 -999 -999", "2 -999 -999", "3 1 2"],
            missing_parent=CUSTOM_SENTINEL,
            reorder=False, renumber=False, pedigree_is_renumbered=False,
        )
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            pyp_utils.fast_reorder(list(ped.pedigree))

    def test_passing_the_configured_sentinel_succeeds(self):
        from PyPedal import pyp_utils
        ped = rows_to_ped(
            ["1 -999 -999", "2 -999 -999", "3 1 2"],
            missing_parent=CUSTOM_SENTINEL,
            reorder=False, renumber=False, pedigree_is_renumbered=False,
        )
        out = pyp_utils.fast_reorder(
            list(ped.pedigree), missingparent=CUSTOM_SENTINEL)
        self.assertEqual(3, len(out))

    def test_asdb_fast_path_no_longer_raises_when_the_wrapper_omits_the_sentinel(self):
        """
        Was : NewPedigree.renumber() omitted missingparent.
        Closed on correctness/medium-3-fast-reorder-missingparent.
        """
        ped = rows_to_ped(
            ["1 -999 -999 1900", "2 -999 -999 1900", "3 1 2 1910"],
            pedformat="asdb",
            missing_parent=CUSTOM_SENTINEL,
            slow_reorder=False,
        )
        self.assertEqual(3, len(ped.pedigree))
        self.assertEqual(CUSTOM_SENTINEL, ped.kw["missing_parent"])

    def test_shipped_defaults_still_honour_a_custom_sentinel_on_asd(self):
        """
        Default slow_reorder=True takes reorder(), which does pass the
        sentinel. The claim is therefore live at the wrappers, not at the
        shipped asd load path.
        """
        ped = rows_to_ped(
            ["1 -999 -999", "2 -999 -999", "3 1 2"],
            missing_parent=CUSTOM_SENTINEL,
        )
        self.assertEqual(3, len(ped.pedigree))
        self.assertEqual(CUSTOM_SENTINEL, ped.kw["missing_parent"])
        for animal in ped.pedigree[:2]:
            self.assertEqual(CUSTOM_SENTINEL, animal.sireID)
            self.assertEqual(CUSTOM_SENTINEL, animal.damID)
