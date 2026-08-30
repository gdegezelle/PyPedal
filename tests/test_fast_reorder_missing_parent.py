"""
 -- internal ``fast_reorder`` callers must honour ``kw['missing_parent']``.

WHAT THIS FILE IS ABOUT
-----------------------
 put ``reorder()`` and ``fast_reorder()`` on one engine,
``_order_pedigree``, and added a trailing ``missingparent=0`` keyword to
``fast_reorder``. Three ``NewPedigree`` wrappers still call
``fast_reorder(self.pedigree)`` and therefore bind the engine's sentinel to
the public default ``0`` instead of the pedigree's configured value.

With the shipped default ``missing_parent = 0`` that is observationally
harmless. With a configured sentinel such as ``-999`` the omitted argument
makes ``_order_pedigree`` treat ``-999`` as a real parent ID and raise
``PyPedalPedigreeStructureError``.

The three wrappers, reconfirmed on ``e0917a0`` after :

    load()             reorder-without-renumber, ``slow_reorder=False``
    renumber()         birth-date/year fast path (``'b'`` or ``'y'`` in
                       pedformat and ``slow_reorder=False``)
    delete_animals()   reorder-without-renumber, ``slow_reorder=False``

``pyp_nrm`` and ``pyp_metrics`` already pass the configured sentinel. The
public ``fast_reorder(..., missingparent=0)`` default is kept for direct
callers that have no ``kw`` context.

WHAT IS *NOT* CLAIMED HERE
--------------------------
's ordering engine, tie-break and refusal contract are closed and
are not reopened.  /  (oldsave) is closed. 
(``delete_animals`` stale indices) is a separate defect: FR-5 proves only
that the configured sentinel reaches ``fast_reorder``, not that deletion
is otherwise correct.

The shipped ``slow_reorder=True`` default is not changed. An ordinary
``asd`` load with ``renumber=True`` takes ``reorder()``, which already
passes the sentinel, so it is not a  reproducer. Load-path tests
that claim the defect force the supported fast configuration and spy
``fast_reorder`` so they cannot pass by taking the slow wrapper.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with the repaired
contract marked ``xfail(strict=True)``. Those markers are removed here.
Characterisation of the omission is inverted, not deleted.

``self.subTest`` was kept out of the tests that were xfailed.
"""
import inspect
import os
import tempfile
import unittest
from unittest import mock

from PyPedal import pyp_errors, pyp_io, pyp_utils
from PyPedal.pyp_newclasses import load_pedigree

from _pedhelpers import chdir_tmp, load_corpus_from_path

CUSTOM = -999
DEFAULT = 0

FOUNDER = "10"
HALF_SIRE = "20"
HALF_DAM = "30"
ORDINARY = "40"

# Youngest first so the input is not a topological order. Distinctive
# original IDs so Finding-37 renumbering cannot be mistaken for identity.
CUSTOM_ROWS = [
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s" % (HALF_DAM, CUSTOM, FOUNDER),
    "%s %s %s" % (HALF_SIRE, FOUNDER, CUSTOM),
    "%s %s %s" % (FOUNDER, CUSTOM, CUSTOM),
]
CUSTOM_ASDB = [
    "%s %s %s 1910" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s %s %s 1905" % (HALF_DAM, CUSTOM, FOUNDER),
    "%s %s %s 1905" % (HALF_SIRE, FOUNDER, CUSTOM),
    "%s %s %s 1900" % (FOUNDER, CUSTOM, CUSTOM),
]
DEFAULT_ROWS = [
    "%s %s %s" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s 0 %s" % (HALF_DAM, FOUNDER),
    "%s %s 0" % (HALF_SIRE, FOUNDER),
    "%s 0 0" % FOUNDER,
]
DEFAULT_ASDB = [
    "%s %s %s 1910" % (ORDINARY, FOUNDER, HALF_SIRE),
    "%s 0 %s 1905" % (HALF_DAM, FOUNDER),
    "%s %s 0 1905" % (HALF_SIRE, FOUNDER),
    "%s 0 0 1900" % FOUNDER,
]

EXPECTED_STABLE_EDGES = {
    FOUNDER: ("MISSING", "MISSING"),
    HALF_SIRE: (FOUNDER, "MISSING"),
    HALF_DAM: ("MISSING", FOUNDER),
    ORDINARY: (FOUNDER, HALF_SIRE),
}


def rows_to_ped(rows, pedformat="asd", **overrides):
    tmp = tempfile.mkdtemp(prefix="m3_")
    path = os.path.join(tmp, "m3.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, **overrides)


def raw_custom(**overrides):
    opts = dict(
        missing_parent=CUSTOM,
        reorder=False,
        renumber=False,
        pedigree_is_renumbered=False,
    )
    opts.update(overrides)
    return rows_to_ped(CUSTOM_ROWS, **opts)


def raw_default(**overrides):
    opts = dict(
        missing_parent=DEFAULT,
        reorder=False,
        renumber=False,
        pedigree_is_renumbered=False,
    )
    opts.update(overrides)
    return rows_to_ped(DEFAULT_ROWS, **opts)


def stable(animal):
    return str(animal.originalID)


def parent_token(parent, by_id, missing):
    if str(parent) == str(missing):
        return "MISSING"
    return by_id.get(str(parent), "?%s" % parent)


def stable_graph(animals, missing):
    """
    Structural oracle. Keyed on original identity. A missing slot is the
    token MISSING, never a numeric stand-in, so default-0 and custom-999
    graphs compare equal after this normalisation.
    """
    by_id = {str(a.animalID): stable(a) for a in animals}
    edges = {}
    for animal in animals:
        edges[stable(animal)] = (
            parent_token(animal.sireID, by_id, missing),
            parent_token(animal.damID, by_id, missing),
        )
    return frozenset(stable(a) for a in animals), edges


def parent_violations(animals, missing):
    position = {str(a.animalID): i for i, a in enumerate(animals)}
    out = []
    for index, animal in enumerate(animals):
        for role, parent in (("sire", animal.sireID), ("dam", animal.damID)):
            if str(parent) == str(missing):
                continue
            if position.get(str(parent), index) >= index:
                out.append((stable(animal), role, str(parent)))
    return out


def founders_first(animals, missing):
    seen_non_founder = False
    for animal in animals:
        is_founder = (str(animal.sireID) == str(missing)
                      and str(animal.damID) == str(missing))
        if is_founder and seen_non_founder:
            return False
        if not is_founder:
            seen_non_founder = True
    return True


def assert_oracle(animals, missing, testcase):
    ids, edges = stable_graph(animals, missing)
    testcase.assertEqual(
        {FOUNDER, HALF_SIRE, HALF_DAM, ORDINARY}, ids)
    testcase.assertEqual(EXPECTED_STABLE_EDGES, edges)
    testcase.assertEqual([], parent_violations(animals, missing))
    testcase.assertTrue(founders_first(animals, missing))
    for animal in animals:
        if str(animal.sireID) == str(missing):
            testcase.assertEqual(missing, animal.sireID)
        if str(animal.damID) == str(missing):
            testcase.assertEqual(missing, animal.damID)


def called_missingparent(spy):
    args, kwargs = spy.call_args
    bound = inspect.signature(pyp_utils.fast_reorder).bind(*args, **kwargs)
    bound.apply_defaults()
    return bound.arguments["missingparent"]


def spy_fast_reorder():
    return mock.patch(
        "PyPedal.pyp_newclasses.pyp_utils.fast_reorder",
        wraps=pyp_utils.fast_reorder,
    )


def snapshot_animals(animals):
    return [
        (id(a), a.originalID, a.animalID, a.sireID, a.damID)
        for a in animals
    ]


# ===========================================================================
# Architecture and the public default. Green today; must stay green.
# ===========================================================================
class TestArchitectureAndPublicDefault(unittest.TestCase):
    """
    's surface: one engine, two wrappers, public default 0.
     must not change this.
    """

    def test_fast_reorder_signature_defaults_missingparent_to_zero(self):
        params = pyp_utils.fast_reorder.__defaults__
        self.assertEqual(("_new_reordered_", "no", False, 0), params)

    def test_reorder_signature_already_takes_missingparent(self):
        params = inspect.signature(pyp_utils.reorder).parameters
        self.assertIn("missingparent", params)
        self.assertEqual(0, params["missingparent"].default)

    def test_both_wrappers_delegate_to_the_shared_engine(self):
        self.assertIn(
            "_order_pedigree", pyp_utils.fast_reorder.__code__.co_names)
        self.assertIn(
            "_order_pedigree", pyp_utils.reorder.__code__.co_names)

    def test_canonical_key_is_missing_parent_and_defaults_to_zero(self):
        ped = raw_default()
        self.assertEqual(0, ped.kw["missing_parent"])
        self.assertIsInstance(ped.kw["missing_parent"], int)

    def test_omitting_the_public_argument_still_treats_minus_999_as_real(self):
        """
        The public default is 0. A direct caller that omits the argument
        is not a NewPedigree wrapper and is not this repair.
        """
        ped = raw_custom()
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError) as ctx:
            pyp_utils.fast_reorder(list(ped.pedigree))
        text = str(ctx.exception)
        self.assertIn("fast_reorder", text)
        self.assertIn("-999", text)


# ===========================================================================
# FR-1, FR-2, FR-6..FR-13 on the explicit / default-0 paths. Green today.
# ===========================================================================
class TestExplicitAndDefaultSentinelControls(unittest.TestCase):
    """Permanent controls.  is a propagation defect, not these."""

    def test_fr1_explicit_custom_sentinel_orders(self):
        ped = raw_custom()
        before = snapshot_animals(ped.pedigree)
        given = list(ped.pedigree)
        out = pyp_utils.fast_reorder(given, missingparent=CUSTOM)
        self.assertEqual(4, len(out))
        assert_oracle(out, CUSTOM, self)
        # FR-13: new list, same objects, input list and parent fields intact.
        self.assertIsNot(out, given)
        self.assertEqual(before, snapshot_animals(given))
        self.assertEqual({id(a) for a in given}, {id(a) for a in out})

    def test_fr2_default_zero_omit_and_explicit_agree(self):
        ped = raw_default()
        omitted = pyp_utils.fast_reorder(list(ped.pedigree))
        explicit = pyp_utils.fast_reorder(
            list(ped.pedigree), missingparent=DEFAULT)
        self.assertEqual(
            [a.originalID for a in omitted],
            [a.originalID for a in explicit])
        assert_oracle(omitted, DEFAULT, self)
        assert_oracle(explicit, DEFAULT, self)

    def test_fr6_founder_both_missing(self):
        out = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        founder = next(a for a in out if stable(a) == FOUNDER)
        self.assertEqual(CUSTOM, founder.sireID)
        self.assertEqual(CUSTOM, founder.damID)

    def test_fr7_known_sire_missing_dam(self):
        out = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        animal = next(a for a in out if stable(a) == HALF_SIRE)
        by_id = {str(a.animalID): stable(a) for a in out}
        self.assertEqual(FOUNDER, parent_token(animal.sireID, by_id, CUSTOM))
        self.assertEqual(CUSTOM, animal.damID)

    def test_fr8_missing_sire_known_dam(self):
        out = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        animal = next(a for a in out if stable(a) == HALF_DAM)
        by_id = {str(a.animalID): stable(a) for a in out}
        self.assertEqual(CUSTOM, animal.sireID)
        self.assertEqual(FOUNDER, parent_token(animal.damID, by_id, CUSTOM))

    def test_fr9_both_parents_known(self):
        out = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        animal = next(a for a in out if stable(a) == ORDINARY)
        by_id = {str(a.animalID): stable(a) for a in out}
        self.assertEqual(FOUNDER, parent_token(animal.sireID, by_id, CUSTOM))
        self.assertEqual(HALF_SIRE, parent_token(animal.damID, by_id, CUSTOM))

    def test_fr10_missing_slots_are_not_converted_into_edges(self):
        out = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        _ids, edges = stable_graph(out, CUSTOM)
        for sire, dam in edges.values():
            self.assertNotEqual(sire, "-999")
        self.assertEqual("MISSING", edges[FOUNDER][0])
        self.assertEqual("MISSING", edges[FOUNDER][1])
        self.assertEqual("MISSING", edges[HALF_SIRE][1])
        self.assertEqual("MISSING", edges[HALF_DAM][0])

    def test_fr11_unresolved_real_parent_still_refuses(self):
        # Drop HALF_SIRE after load so implicit-parent materialisation cannot
        # invent a record. ORDINARY still names that animal as dam.
        raw = raw_default()
        truncated = [a for a in raw.pedigree if stable(a) != HALF_SIRE]
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError) as ctx:
            pyp_utils.fast_reorder(truncated, missingparent=DEFAULT)
        self.assertIn("no animal with that ID", str(ctx.exception))

    def test_fr11_custom_sentinel_does_not_reclassify_a_real_parent(self):
        raw = raw_custom()
        truncated = [a for a in raw.pedigree if stable(a) != HALF_SIRE]
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError) as ctx:
            pyp_utils.fast_reorder(truncated, missingparent=CUSTOM)
        text = str(ctx.exception)
        self.assertIn("no animal with that ID", text)
        self.assertNotIn("-999 as its dam", text)

    def test_fr12_default_and_custom_agree_after_normalising_the_token(self):
        custom = pyp_utils.fast_reorder(
            list(raw_custom().pedigree), missingparent=CUSTOM)
        default = pyp_utils.fast_reorder(
            list(raw_default().pedigree), missingparent=DEFAULT)
        self.assertEqual(
            stable_graph(custom, CUSTOM),
            stable_graph(default, DEFAULT))
        self.assertEqual(
            [stable(a) for a in custom],
            [stable(a) for a in default])

    def test_shipped_asd_load_already_honours_a_custom_sentinel(self):
        """
        Default ``slow_reorder=True`` / ``renumber=True`` takes ``reorder()``,
        which already passes ``kw['missing_parent']``. Not a  path.
        """
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        assert_oracle(ped.pedigree, CUSTOM, self)

    def test_default_zero_fast_load_path_already_succeeds(self):
        ped = rows_to_ped(
            DEFAULT_ROWS,
            missing_parent=DEFAULT,
            reorder=True,
            renumber=False,
            pedigree_is_renumbered=False,
            slow_reorder=False,
        )
        assert_oracle(ped.pedigree, DEFAULT, self)

    def test_default_zero_fast_renumber_path_already_succeeds(self):
        ped = rows_to_ped(
            DEFAULT_ASDB,
            pedformat="asdb",
            missing_parent=DEFAULT,
            slow_reorder=False,
        )
        assert_oracle(ped.pedigree, DEFAULT, self)

    def test_options_dict_preserves_a_custom_sentinel(self):
        ped = raw_custom()
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])

    def test_ini_coerces_minus_999_to_int_and_it_survives_defaults(self):
        self.assertEqual(-999, pyp_io.coerce_ini_value("-999"))
        self.assertIsInstance(pyp_io.coerce_ini_value("-999"), int)
        with chdir_tmp() as tmp:
            pedfile = os.path.join(tmp, "m3.ped")
            ini = os.path.join(tmp, "m3.ini")
            with open(pedfile, "w", encoding="utf-8") as handle:
                handle.write("\n".join(CUSTOM_ROWS) + "\n")
            with open(ini, "w", encoding="utf-8") as handle:
                handle.write("\n".join([
                    "messages = quiet",
                    "pedfile = %s" % pedfile,
                    "pedformat = asd",
                    "sepchar = ' '",
                    "missing_parent = -999",
                    "renumber = 1",
                    "slow_reorder = 1",
                    "pedigree_summary = 0",
                ]) + "\n")
            ped = load_pedigree(optionsfile=ini)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        assert_oracle(ped.pedigree, CUSTOM, self)


# ===========================================================================
# The reproducer, inverted. Was: StructureError / bound missingparent=0.
# ===========================================================================
class TestTheReproducerNoLongerReproduces(unittest.TestCase):
    """
    Characterisation of the three wrappers, inverted. Baseline: e0917a0.
    """

    def test_load_fast_path_no_longer_raises_from_fast_reorder(self):
        """Was: PyPedalPedigreeStructureError naming -999 as a real sire."""
        ped = rows_to_ped(
            CUSTOM_ROWS,
            missing_parent=CUSTOM,
            reorder=True,
            renumber=False,
            pedigree_is_renumbered=False,
            slow_reorder=False,
        )
        self.assertEqual(4, len(ped.pedigree))
        assert_oracle(ped.pedigree, CUSTOM, self)

    def test_renumber_asdb_fast_path_no_longer_raises_from_fast_reorder(self):
        """Was: PyPedalPedigreeStructureError naming -999 as a real sire."""
        ped = rows_to_ped(
            CUSTOM_ASDB,
            pedformat="asdb",
            missing_parent=CUSTOM,
            slow_reorder=False,
        )
        self.assertEqual(4, len(ped.pedigree))
        assert_oracle(ped.pedigree, CUSTOM, self)

    def test_delete_animals_no_longer_binds_the_public_default(self):
        """
        Was: spy showed missingparent=0.  is still not under test.
        """
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        ped.kw["reorder"] = 1
        ped.kw["renumber"] = 0
        ped.kw["slow_reorder"] = False
        with spy_fast_reorder() as spy:
            ped.delete_animals([int(ORDINARY)])
        self.assertTrue(spy.called)
        self.assertEqual(CUSTOM, called_missingparent(spy))


class TestFR3LoadFastPath(unittest.TestCase):
    def test_load_fast_path_orders_with_custom_sentinel(self):
        with spy_fast_reorder() as spy:
            ped = rows_to_ped(
                CUSTOM_ROWS,
                missing_parent=CUSTOM,
                reorder=True,
                renumber=False,
                pedigree_is_renumbered=False,
                slow_reorder=False,
            )
        self.assertTrue(spy.called)
        self.assertEqual(CUSTOM, called_missingparent(spy))
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        assert_oracle(ped.pedigree, CUSTOM, self)
        before_ids = {stable(a) for a in ped.pedigree}
        after = snapshot_animals(ped.pedigree)
        again = pyp_utils.fast_reorder(
            list(ped.pedigree), missingparent=CUSTOM)
        self.assertEqual(before_ids, {stable(a) for a in again})
        self.assertEqual(after, snapshot_animals(again))


class TestFR4RenumberFastPath(unittest.TestCase):
    def test_renumber_asdb_fast_path_preserves_the_stable_graph(self):
        with spy_fast_reorder() as spy:
            ped = rows_to_ped(
                CUSTOM_ASDB,
                pedformat="asdb",
                missing_parent=CUSTOM,
                slow_reorder=False,
            )
        self.assertTrue(spy.called)
        self.assertEqual(CUSTOM, called_missingparent(spy))
        assert_oracle(ped.pedigree, CUSTOM, self)
        internals = [int(a.animalID) for a in ped.pedigree]
        originals = [int(a.originalID) for a in ped.pedigree]
        self.assertNotEqual(internals, originals)
        self.assertEqual(
            {FOUNDER, HALF_SIRE, HALF_DAM, ORDINARY},
            {str(a.originalID) for a in ped.pedigree})


class TestFR5DeleteAnimalsPropagation(unittest.TestCase):
    def test_delete_animals_passes_the_configured_sentinel(self):
        """
        Narrow seam test.  stale-index behaviour is not asserted.
        The load uses the shipped slow path so construction cannot vacuous
        the spy; flags are then flipped to the delete_animals fast wrapper.
        """
        ped = rows_to_ped(CUSTOM_ROWS, missing_parent=CUSTOM)
        ped.kw["reorder"] = 1
        ped.kw["renumber"] = 0
        ped.kw["slow_reorder"] = False
        with spy_fast_reorder() as spy:
            ped.delete_animals([int(ORDINARY)])
        self.assertTrue(spy.called)
        self.assertEqual(CUSTOM, called_missingparent(spy))


class TestFR12FastLoadEquivalence(unittest.TestCase):
    def test_fast_load_custom_matches_default_after_normalising(self):
        custom = rows_to_ped(
            CUSTOM_ROWS,
            missing_parent=CUSTOM,
            reorder=True,
            renumber=False,
            pedigree_is_renumbered=False,
            slow_reorder=False,
        )
        default = rows_to_ped(
            DEFAULT_ROWS,
            missing_parent=DEFAULT,
            reorder=True,
            renumber=False,
            pedigree_is_renumbered=False,
            slow_reorder=False,
        )
        self.assertEqual(
            stable_graph(custom.pedigree, CUSTOM),
            stable_graph(default.pedigree, DEFAULT))


class TestIniReachesFastReorder(unittest.TestCase):
    def test_ini_custom_sentinel_on_the_load_fast_path(self):
        with chdir_tmp() as tmp:
            pedfile = os.path.join(tmp, "m3.ped")
            ini = os.path.join(tmp, "m3.ini")
            with open(pedfile, "w", encoding="utf-8") as handle:
                handle.write("\n".join(CUSTOM_ROWS) + "\n")
            with open(ini, "w", encoding="utf-8") as handle:
                handle.write("\n".join([
                    "messages = quiet",
                    "pedfile = %s" % pedfile,
                    "pedformat = asd",
                    "sepchar = ' '",
                    "missing_parent = -999",
                    "reorder = 1",
                    "renumber = 0",
                    "pedigree_is_renumbered = 0",
                    "slow_reorder = 0",
                    "pedigree_summary = 0",
                ]) + "\n")
            with spy_fast_reorder() as spy:
                ped = load_pedigree(optionsfile=ini)
        self.assertEqual(CUSTOM, ped.kw["missing_parent"])
        self.assertTrue(spy.called)
        self.assertEqual(CUSTOM, called_missingparent(spy))
        assert_oracle(ped.pedigree, CUSTOM, self)


if __name__ == "__main__":
    unittest.main()
