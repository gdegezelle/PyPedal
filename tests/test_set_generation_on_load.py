"""
 -- ``load(set_generations=True)`` must invoke ``set_generation``
after the pedigree metadata it reads actually exists.

WHAT THIS FILE IS ABOUT
-----------------------
``pyp_utils.set_generation`` iterates ``range(pedobj.metadata.num_records)``.
``NewPedigree.__init__`` sets ``self.metadata = {}``. ``load()`` then calls
``set_generation`` when ``kw['set_generations']`` is true, *before* it
constructs ``PedigreeMetadata``. The routine's broad handler swallows the
``AttributeError``, logs it, and returns ``False``. ``load()`` ignores that
return value, so the caller gets a pedigree whose ``igen`` values are still
the initialisation sentinel.

The shipped default is ``set_generations = False``, so ordinary loads never
reach the defect.

The oracle for this phase is lifecycle equivalence, not generation
mathematics:

    load(set_generations=True)
        ==
    load(set_generations=False) + explicit pyp_utils.set_generation(...)

 owns ``gen_coeff`` / Pattie semantics.  established that
``igen`` and legacy ``gen`` are distinct; this file does not copy one into
the other.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with the repaired
contract marked ``xfail(strict=True)``. Those markers are removed here.
Characterisation of the lifecycle failure is inverted, not deleted.

WHAT IS *NOT* CLAIMED HERE
--------------------------
The numbers ``1, 2, 3`` on the chain fixture are what the *current*
``set_generation`` writes when it is allowed to run. They are not an
independent scientific convention.  is not adjudicated here.
``set_age``, birth-year sentinels, reorder, and 's reference
population are out of scope.
"""
import contextlib
import logging
import os
import tempfile
import unittest
from unittest import mock

from PyPedal import pyp_io, pyp_newclasses, pyp_utils
from PyPedal.pyp_newclasses import load_pedigree

from _pedhelpers import chdir_tmp, load_corpus, load_corpus_from_path

BASELINE = "1f90f39"

# Founder -> child -> grandchild. Input is already oldest-first so the
# values cannot be an artifact of reorder. Original IDs are 1..3 so the
# supported non-renumber path can still index animalID - 1.
CHAIN_ROWS = [
    "1 0 0",
    "2 1 0",
    "3 2 0",
]
# Explicit post-load set_generation on current master, measured 2026-08-23.
CHAIN_IGEN = {"1": 1, "2": 2, "3": 3}

# Mixed founder / half-founder / both-known parentage.
#   1, 2        founders
#   3           both parents known
#   4           known sire, missing dam
#   5           missing sire, known dam
#   6           both parents known, one of them a half-founder
MIXED_ROWS = [
    "1 0 0",
    "2 0 0",
    "3 1 2",
    "4 3 0",
    "5 0 3",
    "6 3 4",
]
MIXED_IGEN = {"1": 1, "2": 1, "3": 2, "4": 3, "5": 3, "6": 4}

CORPUS_NINE = (
    ("mrode.ped", "asd", " "),
    ("new_lacy.ped", "asd", " "),
    ("generations.ped", "asdbx", " "),
    ("hartlandclark.ped", "asdb", " "),
    ("boichard2a.ped", "asdg", " "),
    ("doug.ped", "ASDx", " "),
    ("new_ids.ped", "ASD", " "),
    ("horse.ped", "ASD", ","),
    ("userfield.ped", "asdu", " "),
)

GENERATION_FIELDS = ("igen", "gencoeff")


def rows_to_ped(rows, pedformat="asd", **overrides):
    tmp = tempfile.mkdtemp(prefix="m1_")
    path = os.path.join(tmp, "m1.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, **overrides)


def chain(**overrides):
    return rows_to_ped(CHAIN_ROWS, **overrides)


def mixed(**overrides):
    return rows_to_ped(MIXED_ROWS, **overrides)


def stable(animal):
    return str(animal.originalID)


def igen_by_original(ped):
    return {stable(a): a.igen for a in ped.pedigree}


def gen_by_original(ped):
    return {stable(a): a.gen for a in ped.pedigree}


def gencoeff_by_original(ped):
    return {stable(a): a.gencoeff for a in ped.pedigree}


def graph(ped):
    missing = ped.kw["missing_parent"]
    by_id = {str(a.animalID): stable(a) for a in ped.pedigree}
    edges = {}
    for animal in ped.pedigree:
        sire = ("MISSING" if str(animal.sireID) == str(missing)
                else by_id.get(str(animal.sireID), "?%s" % animal.sireID))
        dam = ("MISSING" if str(animal.damID) == str(missing)
               else by_id.get(str(animal.damID), "?%s" % animal.damID))
        edges[stable(animal)] = (sire, dam)
    return frozenset(by_id.values()), edges


def offspring_by_original(ped):
    out = {}
    for animal in ped.pedigree:
        out[stable(animal)] = (
            sorted(str(k) for k in animal.sons),
            sorted(str(k) for k in animal.daus),
            sorted(str(k) for k in animal.unks),
        )
    return out


def animal_fields(ped, skip=GENERATION_FIELDS):
    return [
        {k: repr(v) for k, v in sorted(a.__dict__.items()) if k not in skip}
        for a in ped.pedigree
    ]


def metadata_view(ped):
    md = ped.metadata
    return {
        "type": type(md).__name__,
        "num_records": md.num_records,
        "num_unique_gens": md.num_unique_gens,
        "unique_gens": sorted(repr(x) for x in md.unique_gen_list),
        "num_unique_founders": md.num_unique_founders,
        "num_implicit_parents": md.num_implicit_parents,
    }


@contextlib.contextmanager
def watch_load_set_generation():
    """
    Wrap the implementation ``load()`` actually calls.

    ``NewPedigree.load`` binds ``pyp_utils`` inside ``pyp_newclasses``, so
    the patch target is that name. The wrapper calls the live
    ``pyp_utils.set_generation``; a test that passed because a dummy was
    substituted would not prove the real routine ran.
    """
    real = pyp_utils.set_generation
    events = []

    def wrapper(pedobj):
        event = {
            "meta_type": type(pedobj.metadata).__name__,
            "has_num_records": hasattr(pedobj.metadata, "num_records"),
            "igen_before": [a.igen for a in pedobj.pedigree],
        }
        result = real(pedobj)
        event["result"] = result
        event["igen_after"] = [a.igen for a in pedobj.pedigree]
        events.append(event)
        return result

    with mock.patch(
        "PyPedal.pyp_newclasses.pyp_utils.set_generation",
        side_effect=wrapper,
    ) as spy:
        yield spy, events


def assert_chain_structure(ped, testcase):
    ids, edges = graph(ped)
    testcase.assertEqual({"1", "2", "3"}, ids)
    testcase.assertEqual(
        {"1": ("MISSING", "MISSING"),
         "2": ("1", "MISSING"),
         "3": ("2", "MISSING")},
        edges,
    )


def assert_path_equivalence(false_ped, true_ped, testcase):
    """
    load(False)+explicit vs load(True). Configuration may differ on the
    flag that selected the path. Everything the current set_generation
    is allowed to change must match; nothing else may move.
    """
    testcase.assertEqual(graph(false_ped), graph(true_ped))
    testcase.assertEqual(
        [stable(a) for a in false_ped.pedigree],
        [stable(a) for a in true_ped.pedigree],
    )
    testcase.assertEqual(
        [a.animalID for a in false_ped.pedigree],
        [a.animalID for a in true_ped.pedigree],
    )
    testcase.assertEqual(igen_by_original(false_ped), igen_by_original(true_ped))
    testcase.assertEqual(gen_by_original(false_ped), gen_by_original(true_ped))
    testcase.assertEqual(
        gencoeff_by_original(false_ped), gencoeff_by_original(true_ped))
    testcase.assertEqual(
        offspring_by_original(false_ped), offspring_by_original(true_ped))
    testcase.assertEqual(metadata_view(false_ped), metadata_view(true_ped))
    testcase.assertEqual(animal_fields(false_ped), animal_fields(true_ped))


# ===========================================================================
# Architecture that must stay true on both sides of the repair.
# ===========================================================================
class TestArchitectureAndOptionContract(unittest.TestCase):
    """SG-1 / option surface. Green on the baseline; must stay green."""

    def test_the_canonical_key_defaults_to_false(self):
        ped = chain()
        self.assertIn("set_generations", ped.kw)
        self.assertIs(False, ped.kw["set_generations"])

    def test_sg1_default_load_does_not_call_set_generation(self):
        with watch_load_set_generation() as (spy, events):
            ped = chain()
        self.assertFalse(spy.called)
        self.assertEqual([], events)
        self.assertIs(False, ped.kw["set_generations"])

    def test_default_load_leaves_igen_at_the_initialisation_sentinel(self):
        ped = chain()
        sentinel = ped.kw["missing_igen"]
        self.assertEqual([-999.0], [sentinel])
        self.assertEqual(
            {"1": sentinel, "2": sentinel, "3": sentinel},
            igen_by_original(ped),
        )

    def test_set_generation_writes_igen_and_does_not_write_gen(self):
        """
        : igen and gen are distinct. The existing routine writes
        igen. It does not write gen.  must not change that.
        """
        ped = chain()
        before_gen = gen_by_original(ped)
        self.assertTrue(pyp_utils.set_generation(ped))
        self.assertEqual(CHAIN_IGEN, igen_by_original(ped))
        self.assertEqual(before_gen, gen_by_original(ped))
        self.assertEqual(
            {"1": -999.0, "2": -999.0, "3": -999.0},
            gen_by_original(ped),
        )

    def test_an_options_dict_value_survives_setdefault(self):
        ped = chain(set_generations=True)
        self.assertTrue(ped.kw["set_generations"])

    def test_ini_coerces_one_to_int_and_it_is_truthy(self):
        self.assertEqual(1, pyp_io.coerce_ini_value("1"))
        self.assertTrue(pyp_io.coerce_ini_value("1"))
        self.assertIs(True, pyp_io.coerce_ini_value("true"))


# ===========================================================================
# THE REPRODUCER, INVERTED.
#
# On the phase baseline 1f90f39 each of these asserted the BROKEN behaviour
# and passed. They are kept rather than deleted.
# ===========================================================================
class TestTheReproducerNoLongerReproduces(unittest.TestCase):
    """The  reproducer, now inverted. Baseline: 1f90f39."""

    def test_sg2_set_generations_true_reaches_the_real_set_generation(self):
        with watch_load_set_generation() as (spy, events):
            ped = chain(set_generations=True)
        self.assertTrue(ped.kw["set_generations"])
        self.assertEqual(1, spy.call_count)
        self.assertEqual(1, len(events))

    def test_sg3_load_true_still_returns_a_pedigree(self):
        ped = chain(set_generations=True)
        self.assertIsInstance(ped, pyp_newclasses.NewPedigree)
        self.assertEqual(3, len(ped.pedigree))

    def test_the_call_happens_after_metadata_exists(self):
        """Was: metadata type was dict and num_records was missing."""
        with watch_load_set_generation() as (spy, events):
            chain(set_generations=True)
        self.assertEqual(1, len(events))
        self.assertEqual("PedigreeMetadata", events[0]["meta_type"])
        self.assertTrue(events[0]["has_num_records"])

    def test_set_generation_reports_success(self):
        """Was: assertIs(False, events[0]['result'])."""
        with watch_load_set_generation() as (spy, events):
            ped = chain(set_generations=True)
        self.assertIs(True, events[0]["result"])
        self.assertIsInstance(ped, pyp_newclasses.NewPedigree)

    def test_igen_is_assigned_from_the_sentinel_the_fixture_started_with(self):
        """
        Anti-vacuity: the animals begin at missing_igen. After the optional
        path they are not still there. A fixture that already carried 1/2/3
        would not have exercised the defect.
        """
        with watch_load_set_generation() as (spy, events):
            ped = chain(set_generations=True)
        sentinel = ped.kw["missing_igen"]
        self.assertEqual([sentinel, sentinel, sentinel], events[0]["igen_before"])
        self.assertEqual([1, 2, 3], events[0]["igen_after"])
        self.assertEqual(CHAIN_IGEN, igen_by_original(ped))

    def test_the_lifecycle_error_is_no_longer_logged(self):
        """Was: ERROR ... 'dict' object has no attribute 'num_records'."""
        logging.disable(logging.NOTSET)
        try:
            with self.assertLogs(level="ERROR") as captured:
                ped = chain(set_generations=True)
        finally:
            logging.disable(logging.CRITICAL)
        self.assertIsInstance(ped, pyp_newclasses.NewPedigree)
        text = "\n".join(captured.output)
        self.assertNotIn("unable to assign inferred generations", text)
        self.assertNotIn("num_records", text)

    def test_explicit_post_load_call_is_still_the_working_oracle(self):
        ped = chain(set_generations=False)
        self.assertIsInstance(ped.metadata, pyp_newclasses.PedigreeMetadata)
        self.assertEqual(3, ped.metadata.num_records)
        self.assertTrue(pyp_utils.set_generation(ped))
        self.assertEqual(CHAIN_IGEN, igen_by_original(ped))

    def test_sg7_metadata_exists_and_is_coherent_after_load_true(self):
        ped = chain(set_generations=True)
        self.assertIsInstance(ped.metadata, pyp_newclasses.PedigreeMetadata)
        self.assertEqual(3, ped.metadata.num_records)
        self.assertEqual(len(ped.pedigree), ped.metadata.num_records)
        # nug() reads gen, not igen. set_generation does not write gen.
        self.assertEqual(1, ped.metadata.num_unique_gens)
        self.assertEqual({-999.0}, ped.metadata.unique_gen_list)


# ===========================================================================
# Repaired load(True) contract. Was xfail(strict=True) on 1f90f39.
# ===========================================================================
class TestSG4LoadTrueEqualsExplicitPostLoad(unittest.TestCase):
    """
    Principal  oracle. load(True) must match load(False) followed
    by the current set_generation, on the same pedigree and configuration.
    """

    def test_default_renumber_chain(self):
        path_a = chain(set_generations=False)
        self.assertTrue(pyp_utils.set_generation(path_a))
        path_b = chain(set_generations=True)
        assert_path_equivalence(path_a, path_b, self)

    def test_supported_non_renumber_chain(self):
        opts = dict(renumber=False, pedigree_is_renumbered=False)
        path_a = chain(set_generations=False, **opts)
        self.assertTrue(pyp_utils.set_generation(path_a))
        path_b = chain(set_generations=True, **opts)
        assert_path_equivalence(path_a, path_b, self)


class TestSG5ChainValuesMatchThroughBothPaths(unittest.TestCase):
    """Lifecycle equivalence, not a Finding-34 convention claim."""

    def test_founder_child_grandchild_igen(self):
        path_a = chain(set_generations=False)
        self.assertTrue(pyp_utils.set_generation(path_a))
        path_b = chain(set_generations=True)
        self.assertEqual(CHAIN_IGEN, igen_by_original(path_a))
        self.assertEqual(CHAIN_IGEN, igen_by_original(path_b))
        assert_chain_structure(path_b, self)


class TestSG6MixedParentageEqualsExplicitCall(unittest.TestCase):

    def test_mixed_founder_half_founder_both_known(self):
        path_a = mixed(set_generations=False)
        self.assertTrue(pyp_utils.set_generation(path_a))
        path_b = mixed(set_generations=True)
        self.assertEqual(MIXED_IGEN, igen_by_original(path_a))
        self.assertEqual(MIXED_IGEN, igen_by_original(path_b))
        assert_path_equivalence(path_a, path_b, self)


class TestSG8OptionsDictTruePathAssignsIgen(unittest.TestCase):

    def test_options_dict_true_assigns_the_oracle_values(self):
        with watch_load_set_generation() as (spy, events):
            ped = chain(set_generations=True)
        self.assertEqual(1, spy.call_count)
        self.assertTrue(events[0]["has_num_records"])
        self.assertEqual("PedigreeMetadata", events[0]["meta_type"])
        self.assertIs(True, events[0]["result"])
        self.assertEqual(CHAIN_IGEN, igen_by_original(ped))


class TestSG9IniTruePathAssignsIgen(unittest.TestCase):

    def test_ini_set_generations_one_matches_explicit_post_load(self):
        path_a = chain(set_generations=False)
        self.assertTrue(pyp_utils.set_generation(path_a))
        with chdir_tmp() as tmp:
            pedfile = os.path.join(tmp, "m1.ped")
            ini = os.path.join(tmp, "m1.ini")
            with open(pedfile, "w", encoding="utf-8") as handle:
                handle.write("\n".join(CHAIN_ROWS) + "\n")
            with open(ini, "w", encoding="utf-8") as handle:
                handle.write("\n".join([
                    "messages = quiet",
                    "pedfile = %s" % pedfile,
                    "pedformat = asd",
                    "sepchar = ' '",
                    "set_generations = 1",
                    "renumber = 1",
                    "pedigree_summary = 0",
                ]) + "\n")
            path_b = load_pedigree(optionsfile=ini)
        self.assertTrue(path_b.kw["set_generations"])
        self.assertEqual(CHAIN_IGEN, igen_by_original(path_b))
        assert_path_equivalence(path_a, path_b, self)


# ===========================================================================
# Green on the baseline and required to stay green after the repair.
# ===========================================================================
class TestSG10Determinism(unittest.TestCase):

    def test_repeated_independent_default_loads_do_not_leak(self):
        first = igen_by_original(chain())
        second = igen_by_original(chain())
        self.assertEqual(first, second)
        self.assertEqual({"1": -999.0, "2": -999.0, "3": -999.0}, first)

    def test_repeated_independent_true_loads_are_deterministic(self):
        first = igen_by_original(chain(set_generations=True))
        second = igen_by_original(chain(set_generations=True))
        self.assertEqual(first, second)


class TestSG11NoUnrelatedFieldChanges(unittest.TestCase):
    """
    load(False) vs load(True) may differ only in the fields the current
    explicit set_generation writes. On the baseline that difference is
    empty because the optional path fails. After the repair it is igen.
    """

    def test_every_non_generation_field_matches(self):
        false_ped = chain(set_generations=False)
        true_ped = chain(set_generations=True)
        assert_chain_structure(false_ped, self)
        assert_chain_structure(true_ped, self)
        self.assertEqual(animal_fields(false_ped), animal_fields(true_ped))
        self.assertEqual(gen_by_original(false_ped), gen_by_original(true_ped))
        self.assertEqual(
            gencoeff_by_original(false_ped), gencoeff_by_original(true_ped))
        self.assertEqual(
            offspring_by_original(false_ped), offspring_by_original(true_ped))


class TestSG12DefaultCorpusUnchanged(unittest.TestCase):
    """
    The shipped default is set_generations=False. The normal 9-pedigree
    corpus must keep that baseline identity.
    """

    def test_each_corpus_pedigree_keeps_sentinel_igen_on_the_default_path(self):
        for name, pedformat, sepchar in CORPUS_NINE:
            ped = load_corpus(name, pedformat, sepchar=sepchar)
            self.assertIs(False, ped.kw["set_generations"], name)
            sentinel = ped.kw["missing_igen"]
            igens = [a.igen for a in ped.pedigree]
            self.assertTrue(len(igens) > 0, name)
            self.assertEqual([sentinel] * len(igens), igens, name)
            self.assertIsInstance(ped.metadata, pyp_newclasses.PedigreeMetadata)
            self.assertEqual(len(ped.pedigree), ped.metadata.num_records, name)


class TestMetadataDoesNotDependOnIgen(unittest.TestCase):
    """
    PedigreeMetadata.nug() reads animal.gen. set_generation writes igen.
    Moving the call across metadata construction cannot stale nug.
    """

    def test_nug_is_unchanged_by_an_explicit_set_generation(self):
        ped = chain()
        before = metadata_view(ped)
        self.assertTrue(pyp_utils.set_generation(ped))
        after = metadata_view(ped)
        self.assertEqual(before, after)
        self.assertEqual(CHAIN_IGEN, igen_by_original(ped))


class TestCallCountExactlyOnce(unittest.TestCase):
    """Performance / anti-duplication: the optional path is one call, not two."""

    def test_default_path_calls_zero_times(self):
        with watch_load_set_generation() as (spy, events):
            chain(set_generations=False)
        self.assertEqual(0, spy.call_count)
        self.assertEqual([], events)

    def test_true_path_calls_exactly_once(self):
        with watch_load_set_generation() as (spy, events):
            chain(set_generations=True)
        self.assertEqual(1, spy.call_count)
        self.assertEqual(1, len(events))


if __name__ == "__main__":
    unittest.main()
