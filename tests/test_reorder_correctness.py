"""
 -- pedigree reorder correctness.

WHAT THIS FILE IS ABOUT
-----------------------
``pyp_utils.reorder()`` establishes the invariant every single-forward-pass
algorithm in PyPedal depends on: **every known parent precedes its offspring**.

Until  it could not report failing to. It logged, broke out of its
loop, and returned the pedigree unchanged as a normal value; ``renumber()`` then
failed to resolve the offending parent in its incrementally-built ID map and
*silently set the parent to unknown*::

    # PyPedal/pyp_utils.py, before the repair
    try:
        animal.sireID = id_map[animal.sireID]
    except KeyError:
        animal.sireID = 0          # the edge is deleted; no log, no flag

So a pedigree that could not be ordered was silently stripped of exactly the
parent links that made it unorderable, and every downstream coefficient was then
computed on a graph the user never supplied. The original author diagnosed this
mechanism in ``CHANGES.txt`` on 06/26/2007, repaired the trigger, and left the
failure mode in place. See ``the algorithm notes``.

THE ORACLES
-----------
The suite's only pre-existing ordering test
(``test_correctness_invariants.py:485``) checks parent-before-offspring, and
that is precisely why this defect survived for so long: once ``renumber`` has
deleted the offending edge, **parent-before-offspring passes**. A pedigree that
"looks sorted" is not evidence that it is the pedigree that was loaded.

The load-bearing oracle here is therefore O1, graph preservation, expressed in
stable identity and compared for exact equality. O2 is kept, but as a
supporting check rather than the primary one.

WHAT IS *NOT* CLAIMED HERE
--------------------------
Nothing in this file adjudicates birth-date semantics.  is deferred
and untouched. ``TestReorderIgnoresChronology`` proves that reorder's output
does not depend on ``by``/``bd`` -- it says nothing whatever about what those
fields ought to mean, and it deliberately does not write a bare 1800 or 1900
into a pedigree file and call it "missing", because
``finding-36-birth-year-sentinel-provenance.md`` records that no such rule
exists in either codebase and that one must not be introduced.

Nor does this file reopen 's ``BOICHARD_TIE_BREAK = 'lowest_id'``.
That remains an adjudicated declared convention.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with every test
asserting the repaired behaviour marked ``xfail(strict=True)``, so the contract
was written down before any code satisfied it and nothing could pass by
accident. Those markers were removed in the commit that implemented the repair.

The characterisation tests that pinned the *broken* behaviour were kept rather
than deleted, restated as assertions that the defect no longer reproduces
(``TestTheReproducersNoLongerReproduce``). A reproducer thrown away once it
stops reproducing leaves nothing to stop the defect returning.

``self.subTest`` was deliberately kept out of the tests that were xfailed:
pytest 9 handles the subTest/xfail interaction differently from pytest 8, a
difference this project has already been bitten by
(``docs/-R2-HALF-FOUNDER-VERIFICATION.md`` section 2.4). Those tests still
use plain loops.
"""
import heapq
import os
import unittest

from PyPedal import pyp_errors, pyp_nrm, pyp_utils

from _pedhelpers import owned_temp_dir, chdir_tmp, load_corpus, load_corpus_from_path

# ---------------------------------------------------------------------------
# Fixtures and oracles
# ---------------------------------------------------------------------------

def rows_to_ped(rows, pedformat="asd", sepchar=" ", **overrides):
    """
    Load an inline pedigree, with every generated file confined to a tmpdir.

    Written through ``load_corpus_from_path`` rather than by hand so the
    repository-delta guard in ``conftest.py`` stays satisfied.
    """
    tmp = owned_temp_dir(prefix="pypedal_ro_")
    path = os.path.join(tmp, "ro.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, sepchar=sepchar, **overrides)


def raw_ped(rows, pedformat="asd", sepchar=" ", **overrides):
    """The same pedigree with reorder and renumber both off: input order."""
    opts = dict(renumber=False, reorder=False, pedigree_is_renumbered=False)
    opts.update(overrides)
    return rows_to_ped(rows, pedformat, sepchar, **opts)


def stable(animal, pedobj):
    """
    The identity that survives renumbering.

    For string-ID (``ASD``) pedigrees ``originalID`` holds the *hash*, not the
    source token (``pyp_newclasses.py:2585-2587``), so ``name`` is the only
    stable handle there.
    """
    if "A" in pedobj.kw["pedformat"]:
        return str(animal.name)
    return str(animal.originalID)


def graph_of(pedobj):
    """
    **Oracle A1.** The pedigree as a directed graph in stable identity:
    ``(animals, edges)`` with ``edges`` a set of ``(offspring, role, parent)``.

    Implementation-independent by construction -- it never calls the routine
    under test, and it is keyed on identity rather than on position, so a
    reordering is invisible to it and a *changed graph* is not.
    """
    missing = str(pedobj.kw["missing_parent"])
    by_id = {str(a.animalID): stable(a, pedobj) for a in pedobj.pedigree}
    animals, edges = set(), set()
    for animal in pedobj.pedigree:
        me = stable(animal, pedobj)
        animals.add(me)
        for role, parent in (("sire", animal.sireID), ("dam", animal.damID)):
            if str(parent) != missing:
                # An unresolvable parent is recorded verbatim rather than
                # dropped: silently forgetting it here would hide exactly the
                # defect this oracle exists to catch.
                edges.add((me, role, by_id.get(str(parent), "?%s" % parent)))
    return animals, edges


def order_of(pedobj):
    """The reordered sequence, in stable identity."""
    return [stable(a, pedobj) for a in pedobj.pedigree]


def parent_violations(pedobj):
    """**Oracle A2.** Known parent edges whose parent does not precede."""
    missing = str(pedobj.kw["missing_parent"])
    position = {str(a.animalID): i for i, a in enumerate(pedobj.pedigree)}
    out = []
    for index, animal in enumerate(pedobj.pedigree):
        for role, parent in (("sire", animal.sireID), ("dam", animal.damID)):
            if str(parent) == missing:
                continue
            if position.get(str(parent), index) >= index:
                out.append((stable(animal, pedobj), role, str(parent)))
    return out


def founders_first(pedobj):
    """**Oracle O3.** True when no founder follows a non-founder."""
    missing = str(pedobj.kw["missing_parent"])
    seen_non_founder = False
    for animal in pedobj.pedigree:
        is_founder = (str(animal.sireID) == missing
                      and str(animal.damID) == missing)
        if is_founder and seen_non_founder:
            return False
        if not is_founder:
            seen_non_founder = True
    return True


def convention_order(pedobj):
    """
    **Executable statement of the reorder convention.**

    Founders first in input order, then a stable topological order using
    original input position as the tie-break. ``pedobj`` must be loaded with
    reorder and renumber OFF, so ``pedobj.pedigree`` is still in input order.

    This encodes a *software convention*, not a scientific truth, which is
    exactly why writing it out as a specification is the right kind of oracle:
    there is no external authority to compare against, so the contract has to
    be stated somewhere it can be executed.
    """
    missing = str(pedobj.kw["missing_parent"])
    animals = list(pedobj.pedigree)
    position = {str(a.animalID): i for i, a in enumerate(animals)}
    by_id = {str(a.animalID): a for a in animals}

    indegree, children = {}, {}
    for animal in animals:
        key = str(animal.animalID)
        indegree[key] = 0
        children.setdefault(key, [])
    for animal in animals:
        key = str(animal.animalID)
        for parent in (animal.sireID, animal.damID):
            pkey = str(parent)
            if pkey != missing and pkey in by_id:
                indegree[key] += 1
                children[pkey].append(key)

    out, ready = [], []
    for animal in animals:
        key = str(animal.animalID)
        if (str(animal.sireID) == missing and str(animal.damID) == missing):
            out.append(stable(by_id[key], pedobj))
            for child in children[key]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    heapq.heappush(ready, position[child])
    # Anything that became eligible once the founders were emitted.
    for animal in animals:
        key = str(animal.animalID)
        if indegree[key] == 0 and not (str(animal.sireID) == missing
                                       and str(animal.damID) == missing):
            heapq.heappush(ready, position[key])

    emitted = set(out)
    while ready:
        key = str(animals[heapq.heappop(ready)].animalID)
        if stable(by_id[key], pedobj) in emitted:
            continue
        out.append(stable(by_id[key], pedobj))
        emitted.add(stable(by_id[key], pedobj))
        for child in children[key]:
            indegree[child] -= 1
            if indegree[child] == 0:
                heapq.heappush(ready, position[child])
    return out


# Fixture matrix RO-1 .. RO-9, the orderable cases.
ORDERABLE = {
    "RO-1 already ordered": ["1 0 0", "2 0 0", "3 1 2"],
    "RO-2 child before parents": ["3 1 2", "1 0 0", "2 0 0"],
    "RO-3 reversed chain": ["5 3 4", "4 1 2", "3 1 2", "2 0 0", "1 0 0"],
    "RO-4 branching": ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4", "6 3 4"],
    "RO-5 disconnected": ["3 1 2", "1 0 0", "2 0 0", "6 4 5", "4 0 0", "5 0 0"],
    "RO-6 half founder": ["3 1 0", "1 0 0", "4 0 2", "2 0 0"],
    "RO-7 founders only": ["1 0 0", "2 0 0", "3 0 0"],
    "RO-8 sparse ids": ["500 100 200", "100 0 0", "200 0 0", "900 500 200"],
}

# Structurally invalid input. RO-10 is covered separately because the duplicate
# has to be built as a genuine duplicate rather than a loader artefact.
UNORDERABLE = {
    "RO-11 self parent": ["1 1 0", "2 1 0"],
    "RO-12 two-node cycle": ["1 2 0", "2 1 0"],
    "RO-13 longer cycle": ["1 3 0", "2 1 0", "3 2 0"],
}


# ---------------------------------------------------------------------------
# 1. Characterisation -- what the code does TODAY. These pass on the baseline.
# ---------------------------------------------------------------------------

class TestTheReproducersNoLongerReproduce(unittest.TestCase):
    """
     section 2, inverted.

    Each test here was written as a *characterisation* of the defect, landed one
    commit ahead of the repair, and is kept -- restated as the repaired
    behaviour -- rather than deleted. A reproducer that is thrown away once it
    stops reproducing leaves nothing to stop the defect coming back, and these
    are cheap.
    """

    CYCLE = ["1 2 0", "2 1 0", "3 1 2"]

    def test_a_cycle_no_longer_has_an_edge_silently_deleted(self):
        """
        Both mutual edges are declared by the file. The load used to keep one
        and drop the other; it now refuses, so no partial graph is produced.
        """
        declared = graph_of(raw_ped(self.CYCLE))[1]
        self.assertIn(("1", "sire", "2"), declared)
        self.assertIn(("2", "sire", "1"), declared)
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            rows_to_ped(self.CYCLE)

    def test_the_public_entry_point_no_longer_answers_for_a_cycle(self):
        """
        What made this a release blocker rather than a nuisance: this call used
        to return ``{1: 0.0, 2: 0.0, 3: 0.25}`` -- a confident, well-formed
        answer for a pedigree in which animal 1 had been quietly orphaned.
        """
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            ped = rows_to_ped(self.CYCLE)
            pyp_nrm.inbreeding(ped, method="tabular", output=False)

    def test_a_self_parent_is_refused_rather_than_accepted(self):
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            rows_to_ped(["1 1 0", "2 1 0"])

    def test_reorder_is_idempotent(self):
        """The founder pre-pass used to reverse the founder block every call."""
        raw = raw_ped(["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"])
        once = pyp_utils.reorder(list(raw.pedigree), missingparent=0)
        once_ids = [a.originalID for a in once]
        twice = pyp_utils.reorder(list(once), missingparent=0)
        self.assertEqual(once_ids, [a.originalID for a in twice])

    def test_reorder_no_longer_replaces_the_animals_with_deep_copies(self):
        raw = raw_ped(["3 1 2", "1 0 0", "2 0 0"])
        given = list(raw.pedigree)
        before = sorted(id(a) for a in given)
        returned = pyp_utils.reorder(given, missingparent=0)
        self.assertIs(given, returned, "reorder is documented as in-place")
        self.assertEqual(before, sorted(id(a) for a in returned))

    def test_the_corpus_fixture_no_longer_duplicates_an_animal_id(self):
        """
         section 5b. ``'animal0'`` was materialised twice because the
        sire loop registered the hashed integer and the dam loop looked up the
        raw source token.
        """
        ped = load_corpus("new_ids.ped", renumber=False, reorder=False,
                          pedigree_is_renumbered=False)
        ids = [a.animalID for a in ped.pedigree]
        self.assertEqual(len(ids), len(set(ids)))

    def test_inbreeding_with_a_generation_limit_no_longer_returns_nothing(self):
        """
         section 5c. ``find_ancestors_g`` truncates, the truncation
        boundary dangles, and the blanket handler at ``pyp_nrm.py`` turned the
        resulting failure into ``{}`` -- a result shaped exactly like a real one.
        It refuses now. Whether a truncated pedigree *has* meaningful
        coefficients is a separate question this finding does not answer.
        """
        ped = load_corpus("hartlandclark.ped")
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_nrm.inbreeding(ped, method="tabular", gens=2, output=False)


class TestWhatAlreadyHoldsAndMustKeepHolding(unittest.TestCase):
    """
    The valid domain. These pass on the baseline and are the containment floor:
    the repair must not disturb any of them.
    """

    def test_graph_is_preserved_on_orderable_fixtures(self):
        for name, rows in ORDERABLE.items():
            with self.subTest(fixture=name):
                self.assertEqual(graph_of(raw_ped(rows)),
                                 graph_of(rows_to_ped(rows)))

    def test_parents_precede_offspring_on_orderable_fixtures(self):
        for name, rows in ORDERABLE.items():
            with self.subTest(fixture=name):
                self.assertEqual([], parent_violations(rows_to_ped(rows)))

    def test_founders_come_first_on_orderable_fixtures(self):
        for name, rows in ORDERABLE.items():
            with self.subTest(fixture=name):
                self.assertTrue(founders_first(rows_to_ped(rows)))

    def test_the_same_input_gives_the_same_order(self):
        for name, rows in ORDERABLE.items():
            with self.subTest(fixture=name):
                self.assertEqual(order_of(rows_to_ped(rows)),
                                 order_of(rows_to_ped(rows)))

    def test_fast_reorder_returns_a_new_list_of_the_same_objects(self):
        raw = raw_ped(["3 1 2", "1 0 0", "2 0 0"])
        given = list(raw.pedigree)
        snapshot = list(given)
        returned = pyp_utils.fast_reorder(given)
        self.assertIsNot(given, returned)
        self.assertEqual(snapshot, given, "input list must not be reordered")
        self.assertEqual(sorted(id(a) for a in snapshot),
                         sorted(id(a) for a in returned))


class TestReorderIgnoresChronology(unittest.TestCase):
    """
    RO-20. Reorder must depend on the parent graph and the tie-break, and on
    nothing else.

     is deferred and is NOT adjudicated here. In particular no
    variant writes a bare 1800 or 1900 into the pedigree file and calls it
    "missing": ``missing_byear``/``missing_bdate`` are *defaults applied when
    the birth column is absent* (``pyp_newclasses.py:2637-2646``), and
    ``finding-36-birth-year-sentinel-provenance.md`` records that no
    ``birth_year == 1900 -> unknown`` rule exists in either codebase and that
    one must not be introduced.

    Note that ``missing_byear`` alone does not reach ``animal.by`` -- ``by`` is
    derived from ``missing_bdate`` whenever a birth date is available -- so
    both are set together, as PyPedal 2.0.4 set them.
    """

    ROWS = ["5 3 4", "4 1 2", "3 1 2", "2 0 0", "1 0 0"]
    DATED = ["5 3 4 2010", "4 1 2 2001", "3 1 2 2000",
             "2 0 0 1991", "1 0 0 1990"]
    DATED_1800 = ["5 3 4 1804", "4 1 2 1803", "3 1 2 1802",
                  "2 0 0 1801", "1 0 0 1800"]
    DATED_1900 = ["5 3 4 1904", "4 1 2 1903", "3 1 2 1902",
                  "2 0 0 1901", "1 0 0 1900"]

    def variants(self):
        return [
            ("A known years", self.DATED, "asdb", {}),
            ("B no birth column", self.ROWS, "asd", {}),
            ("C real year 1800", self.DATED_1800, "asdb", {}),
            ("D real year 1900", self.DATED_1900, "asdb", {}),
        ]

    def test_the_variants_really_do_differ_in_birth_information(self):
        """Non-vacuity: if every variant carried the same ``by`` this would
        prove nothing."""
        seen = set()
        for _, rows, fmt, kw in self.variants():
            ped = rows_to_ped(rows, fmt, **kw)
            seen.add(tuple(a.by for a in ped.pedigree))
        self.assertGreater(len(seen), 1)

    def test_order_and_graph_do_not_depend_on_birth_information(self):
        expected_order, expected_graph = None, None
        for label, rows, fmt, kw in self.variants():
            with self.subTest(variant=label):
                ped = rows_to_ped(rows, fmt, **kw)
                if expected_order is None:
                    expected_order = order_of(ped)
                    expected_graph = graph_of(ped)
                self.assertEqual(expected_order, order_of(ped))
                self.assertEqual(expected_graph, graph_of(ped))
                self.assertEqual([], parent_violations(ped))

    def test_reorder_does_not_mutate_birth_fields(self):
        raw = raw_ped(self.DATED, "asdb")
        before = {a.animalID: (a.by, a.bd) for a in raw.pedigree}
        for animal in pyp_utils.reorder(list(raw.pedigree), missingparent=0):
            self.assertEqual(before[animal.animalID], (animal.by, animal.bd))


# ---------------------------------------------------------------------------
# 2. The contract under repair. Every test here is expected to FAIL until the
#    ordering engine lands, and the markers come off in that same commit.
# ---------------------------------------------------------------------------

class TestUnorderablePedigreesAreRefused(unittest.TestCase):
    """
    ``methodology.tex:10``: "PyPedal can reorder any pedigree unless there is
    an error in it that would prevent unambiguously placing parents before
    offspring." A cyclic pedigree is a documented error condition, not
    something reorder is licensed to resolve by deleting an edge.
    """

    def test_cycles_and_self_parents_refuse(self):
        for rows in UNORDERABLE.values():
            with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
                rows_to_ped(rows)

    def test_the_ordering_routines_refuse_directly(self):
        for rows in UNORDERABLE.values():
            raw = raw_ped(rows)
            with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
                pyp_utils.reorder(list(raw.pedigree), missingparent=0)
            with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
                pyp_utils.fast_reorder(list(raw.pedigree))

    def test_a_duplicate_animal_id_refuses(self):
        raw = raw_ped(["1 0 0", "2 0 0", "3 1 2"])
        animals = list(raw.pedigree)
        animals.append(animals[0])
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            pyp_utils.reorder(animals, missingparent=0)

    def test_a_dangling_parent_reference_refuses(self):
        raw = raw_ped(["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"])
        truncated = [a for a in raw.pedigree if str(a.animalID) != "3"]
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            pyp_utils.reorder(truncated, missingparent=0)

    def test_the_refusal_does_not_overclaim_which_animals_are_on_the_cycle(self):
        """
        The Kahn residue contains cyclic animals *and* their descendants. The
        message must not assert that every animal it lists lies on a cycle.
        """
        raw = raw_ped(["1 2 0", "2 1 0", "3 1 2"])
        try:
            pyp_utils.reorder(list(raw.pedigree), missingparent=0)
        except pyp_errors.PyPedalPedigreeStructureError as exc:
            text = str(exc).lower()
            # Animal 3 is a descendant of the 1<->2 cycle, not a member of it,
            # but it is genuinely unorderable and is genuinely listed. What the
            # message must not do is present the whole residue as the cycle.
            self.assertIn("unorderable", text)
            self.assertIn("descend", text)
            self.assertEqual(["1", "2", "3"],
                             sorted(str(a) for a in exc.animals))
        else:
            self.fail("expected a refusal")

    def test_a_generation_limited_request_refuses_instead_of_returning_nothing(self):
        ped = load_corpus("hartlandclark.ped")
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_nrm.inbreeding(ped, method="tabular", gens=2, output=False)


class TestTheOrderingIsADeterministicFixedPoint(unittest.TestCase):
    """D2 and D3. Plain loops, no ``subTest`` -- see the module docstring."""

    def test_reordering_twice_changes_nothing(self):
        for rows in ORDERABLE.values():
            raw = raw_ped(rows)
            once = pyp_utils.reorder(list(raw.pedigree), missingparent=0)
            once_ids = [a.animalID for a in once]
            twice = pyp_utils.reorder(list(once), missingparent=0)
            self.assertEqual(once_ids, [a.animalID for a in twice])

    def test_the_order_conforms_to_the_adjudicated_convention(self):
        for rows in ORDERABLE.values():
            self.assertEqual(convention_order(raw_ped(rows)),
                             order_of(rows_to_ped(rows)))

    def test_an_already_correct_pedigree_is_left_alone(self):
        """
        The founder block is currently reversed on every call, so a pedigree
        that is already in the adjudicated order does not survive one.
        """
        rows = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"]
        self.assertEqual(["1", "2", "3", "4", "5"], order_of(rows_to_ped(rows)))


class TestEveryInputPermutationIsIndependentlyCorrect(unittest.TestCase):
    """
    RO-16. Under an input-position tie-break the canonical sequence is a
    function of the graph **and** the input order, so different permutations
    may legitimately order incomparable animals differently. That is D4, and it
    is *not* permutation invariance -- this test does not assert that the
    permutations agree with each other.

    What every permutation must independently satisfy is the graph, the
    ordering invariants and its own tie-break.
    """

    BASE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2", "5 3 4"]

    def permutations(self):
        import itertools
        return list(itertools.permutations(self.BASE))[:12]

    def test_each_permutation_satisfies_every_ordering_property(self):
        for perm in self.permutations():
            rows = list(perm)
            raw = raw_ped(rows)
            ped = rows_to_ped(rows)
            self.assertEqual(graph_of(raw), graph_of(ped))
            self.assertEqual([], parent_violations(ped))
            self.assertTrue(founders_first(ped))
            self.assertEqual(convention_order(raw), order_of(ped))
            self.assertEqual(order_of(ped), order_of(rows_to_ped(rows)))

    def test_order_invariant_science_agrees_across_permutations(self):
        """
        A different claim from the one above, and deliberately separate: the
        reordered *sequence* may differ between permutations, but a coefficient
        keyed on stable identity may not.
        """
        expected = None
        for perm in self.permutations():
            with self.subTest(permutation=perm):
                ped = rows_to_ped(list(perm))
                back = {a.animalID: stable(a, ped) for a in ped.pedigree}
                got = {back[k]: round(v, 12) for k, v in
                       pyp_nrm.inbreeding(ped, method="tabular",
                                          output=False)["fx"].items()}
                if expected is None:
                    expected = got
                self.assertEqual(expected, got)


class TestAliasingIsPreserved(unittest.TestCase):
    """
    **Oracle A3.** The two public names have opposite aliasing contracts, and a
    shared engine is exactly the shape of change that could silently collapse
    them into one. Asserted with object identity, not equality.
    """

    def test_reorder_keeps_the_very_same_animal_objects(self):
        raw = raw_ped(["5 3 4", "4 1 2", "3 1 2", "2 0 0", "1 0 0"])
        given = list(raw.pedigree)
        before = sorted(id(a) for a in given)
        returned = pyp_utils.reorder(given, missingparent=0)
        self.assertIs(given, returned)
        self.assertEqual(before, sorted(id(a) for a in returned))


class TestTheImplicitParentIsMaterialisedOnce(unittest.TestCase):
    """
     section 5b, repaired narrowly: one implicit parent per unique
    source token. This does **not** make duplicate ``animalID``s acceptable --
    a genuine duplicate still refuses, which
    ``test_a_duplicate_animal_id_refuses`` above pins.
    """

    def test_new_ids_has_no_duplicate_animal_id(self):
        ped = load_corpus("new_ids.ped", renumber=False, reorder=False,
                          pedigree_is_renumbered=False)
        ids = [a.animalID for a in ped.pedigree]
        self.assertEqual(len(ids), len(set(ids)))

    def test_new_ids_keeps_every_declared_relationship(self):
        """
        Not marked xfail: the real parent links already survive today. The
        spurious duplicate costs a *record*, not a relationship. This is the
        containment floor the narrow de-duplication must not disturb, and it is
        stated here rather than left implicit because "the count changed"
        would otherwise be indistinguishable from "an animal was dropped".
        """
        ped = load_corpus("new_ids.ped")
        _, edges = graph_of(ped)
        self.assertIn(("'animal5'", "sire", "'animal2'"), edges)
        self.assertIn(("'animal5'", "dam", "'animal3'"), edges)
        self.assertIn(("'animal13'", "sire", "'animal7'"), edges)
        self.assertIn(("'animal13'", "dam", "'animal8'"), edges)


class TestTheCallersGetTheOrderingTheyRelyOn(unittest.TestCase):
    """
    Ordering is established by two different paths depending on ``kw``, and
    consumed by routines that assume it without checking. Both paths are
    exercised here, and the refusal is followed all the way out to the public
    API -- the defect was never that ``reorder`` alone was wrong, it was that
    nothing between ``reorder`` and the caller would say so.
    """

    ROWS = ["5 3 4", "4 1 2", "3 1 2", "2 0 0", "1 0 0"]

    def test_the_renumber_path_orders_and_the_index_contract_holds(self):
        """``renumber=True`` -- the shipped default, reorder nested inside."""
        ped = rows_to_ped(self.ROWS)
        self.assertEqual([], parent_violations(ped))
        self.assertTrue(founders_first(ped))
        for index, animal in enumerate(ped.pedigree):
            self.assertEqual(index + 1, int(animal.animalID))

    def test_the_standalone_reorder_path_orders_without_renumbering(self):
        """``reorder=True, renumber=False`` -- the other gate in ``load()``."""
        ped = rows_to_ped(self.ROWS, reorder=True, renumber=False,
                          pedigree_is_renumbered=False)
        self.assertEqual([], parent_violations(ped))
        self.assertTrue(founders_first(ped))

    def test_both_reorder_paths_agree(self):
        slow = rows_to_ped(self.ROWS, slow_reorder=True)
        fast = rows_to_ped(self.ROWS, slow_reorder=False)
        self.assertEqual(order_of(slow), order_of(fast))
        self.assertEqual(graph_of(slow), graph_of(fast))

    def test_a_non_zero_missing_parent_sentinel_is_honoured_on_both(self):
        """
        ``fast_reorder`` had no way to receive the configured sentinel and
        hardcoded 0, so flipping ``slow_reorder`` silently changed which
        animals counted as founders.
        """
        rows = ["5 3 4", "4 1 2", "3 1 2", "2 -999 -999", "1 -999 -999"]
        for slow in (True, False):
            ped = rows_to_ped(rows, missing_parent=-999, slow_reorder=slow)
            self.assertEqual([], parent_violations(ped))
            self.assertTrue(founders_first(ped))

    def test_the_refusal_reaches_the_public_api(self):
        """
        ``NewPedigree.renumber()`` used to collapse every failure into
        ``return False``, and all four of its call sites discard that.
        """
        with self.assertRaises(pyp_errors.PyPedalPedigreeStructureError):
            rows_to_ped(["1 2 0", "2 1 0"])

    def test_the_scientific_routines_still_agree_with_their_oracles(self):
        """
        The published values, keyed on stable identity. Reordering changes
        internal numbering; it must not change a coefficient.
        """
        mrode = rows_to_ped(["1 0 0", "2 0 0", "3 1 2", "4 1 0",
                             "5 4 3", "6 5 3"])
        back = {a.animalID: stable(a, mrode) for a in mrode.pedigree}
        fx = {back[k]: v for k, v in
              pyp_nrm.inbreeding(mrode, method="tabular",
                                 output=False)["fx"].items()}
        # Mrode's textbook answer for animal 5.
        self.assertAlmostEqual(0.125, fx["5"], places=12)

    def test_lacys_published_effective_founder_count_is_unchanged(self):
        from PyPedal import pyp_metrics
        with chdir_tmp():
            result = pyp_metrics.a_effective_founders_lacy(
                load_corpus("new_lacy.ped"))
        self.assertAlmostEqual(2.909090909090909,
                               float(result["fa_effective_founders"]),
                               places=9)


if __name__ == "__main__":
    unittest.main()
