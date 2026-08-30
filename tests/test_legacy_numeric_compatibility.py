"""
Class B -- MIGRATION-FIDELITY tests.

These are behavioural-compatibility tests. They do not establish mathematical
correctness.

Every value pinned here was measured to be identical under PyPedal 2.0.4 on
Python 2.7 and under this port (category A == category B in the audit's evidence
model), but has **no independent mathematical oracle**. Agreement between the
two implementations proves only that the migration preserved behaviour. Where
Python 2 was wrong -- and the audit found several such places -- it would prove
that the port preserved the error.

So: a failure here is a real regression signal and must be investigated. But a
value in this file must never be cited as evidence that PyPedal computes the
right answer, and must never be promoted into
``tests/test_correctness_invariants.py`` without independent evidence. Anything
that acquires an oracle moves to that file and stops being pinned here.

Deliberately **not** pinned: outputs the audit classified as wrong or
unadjudicated -- ``a_effective_founders_boichard`` (),
``a_effective_ancestors_definite`` / ``_indefinite`` (, ), and
``effective_founders_lacy`` on the half-founder pedigrees (). Pinning a
known-wrong output as regression truth would entrench it.
"""
import ast
import inspect
import textwrap
import unittest

from PyPedal import pyp_newclasses, pyp_nrm, pyp_metrics

from _pedhelpers import load_corpus, nrm_value


def _self_attribute_targets(target):
    if (
        isinstance(target, ast.Attribute)
        and isinstance(target.value, ast.Name)
        and target.value.id == "self"
    ):
        yield target.attr
    elif isinstance(target, (ast.Tuple, ast.List)):
        for elt in target.elts:
            yield from _self_attribute_targets(elt)


class TestNewAnimalAttributeParity(unittest.TestCase):
    """
    ``NewAnimal.__init__`` must assign everything PyPedal 2.0.4's does.

    This is the guard that would have caught **both**  (``userField``)
    and  (``ancestor``) at migration time, and it is what stops a third.
    Both were silent: the constructor simply omitted an assignment, leaving live
    consumers to raise ``AttributeError`` on a code path nothing exercised.

    The expected set is transcribed from 2.0.4 rather than computed, because that
    file is Python 2 and cannot be parsed by this interpreter. It was extracted
    from ``legacy/.../pyp_newclasses.py`` lines 3060-3244 (``class NewAnimal``
    through the end of its ``__init__``).

    A **superset** is asserted, not equality: the port legitimately adds
    attributes 2.0.4 never had. What must never happen is one going missing.
    """

    # PyPedal 2.0.4, class NewAnimal, __init__ only. 32 names.
    PY2_ATTRIBUTES = frozenset({
        "age", "alive", "alleles", "ancestor", "animalID", "bd", "breed", "by",
        "damID", "damName", "daus", "displayName", "fa", "fg", "founder", "gen",
        "gencoeff", "herd", "homozygosity", "igen", "name", "originalHerd",
        "originalID", "paddedID", "pedcomp", "renumberedID", "sex", "sireID",
        "sireName", "sons", "unks", "userField",
    })

    @staticmethod
    def _assigned_attributes(func):
        """Names bound as ``self.<name> = ...`` anywhere in ``func``."""
        # __init__ is indented inside its class; ast.parse needs it flush left.
        tree = ast.parse(textwrap.dedent(inspect.getsource(func)))
        found = set()
        for node in ast.walk(tree):
            targets = []
            if isinstance(node, ast.Assign):
                targets = node.targets
            elif isinstance(node, ast.AnnAssign):
                targets = [node.target]
            for target in targets:
                for attr_target in _self_attribute_targets(target):
                    found.add(attr_target)
        return found

    def test_port_assigns_every_attribute_python2_assigned(self):
        assigned = self._assigned_attributes(pyp_newclasses.NewAnimal.__init__)
        # Guard the guard: a parser that found nothing would pass vacuously.
        self.assertGreater(len(assigned), 20)

        missing = self.PY2_ATTRIBUTES - assigned
        self.assertEqual(
            set(), missing,
            "NewAnimal.__init__ no longer assigns %s, which PyPedal 2.0.4 "
            "assigns. Missing attributes used to fail only at first use."
            % sorted(missing))

    def test_the_two_repaired_attributes_specifically(self):
        """
        Named explicitly so a regression on either attribute points
        straight at the missing assignment rather than at a set difference.
        """
        assigned = self._assigned_attributes(pyp_newclasses.NewAnimal.__init__)
        self.assertIn("userField", assigned, "userField is no longer assigned")
        self.assertIn("ancestor", assigned, "ancestor is no longer assigned")


class TestPostLoadStructure(unittest.TestCase):
    """Record counts and founder counts after the full load pipeline."""

    EXPECTED = {
        # pedigree            n_animals
        "mrode.ped": 6,
        "new_lacy.ped": 7,
        "hartlandclark.ped": 15,
        "generations.ped": 13,
        "boichard2a.ped": 14,
    }

    def test_record_counts(self):
        for name, expected in self.EXPECTED.items():
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                self.assertEqual(expected, len(ped.pedigree))
                self.assertEqual(expected, ped.metadata.num_records)

    def test_founder_counts_are_stable(self):
        expected = {
            "mrode.ped": 2,
            "new_lacy.ped": 3,
            "hartlandclark.ped": 3,
            "generations.ped": 5,
            "boichard2a.ped": 4,
        }
        for name, want in expected.items():
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                founders = sum(1 for a in ped.pedigree
                               if getattr(a, "founder", "n") == "y")
                self.assertEqual(want, founders)
                self.assertEqual(want, ped.metadata.num_unique_founders)

    # String-ID (ASD) pedigrees. These counts are the PyPedal 2.0.4 values --
    # doug.ped loading 45 animals is the reference behaviour  broke.
    # The *invariant* that a string pedigree must never load empty is asserted
    # independently in test_correctness_invariants.py; this only pins the exact
    # counts.
    STRING_EXPECTED = (
        ("doug.ped", "ASDx", " ", 45),
        # 15 until . 'animal0' is this file's missing-parent
        # placeholder and is not PyPedal's, so it is materialised as a real
        # animal -- correctly -- but it used to be materialised twice, once by
        # the sire pass and once by the dam pass, because the passes registered
        # and probed `idmap` in different key domains. PyPedal 4 keeps a single
        # record; 2.0.4 duplicated it.
        ("new_ids.ped", "ASD", " ", 14),
        ("horse.ped", "ASD", ",", 16),
    )

    def test_string_id_record_counts(self):
        for name, pedformat, sepchar, expected in self.STRING_EXPECTED:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat, sepchar=sepchar)
                self.assertEqual(expected, len(ped.pedigree))

    def test_id_maps_round_trip(self):
        for name in self.EXPECTED:
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                self.assertEqual(len(ped.pedigree), len(ped.backmap))
                for original, renumbered in ped.idmap.items():
                    self.assertEqual(original, ped.backmap[renumbered])


class TestStorageBackendEquivalence(unittest.TestCase):
    """
    ``kw['matrix_type']`` selects a dense NumPy or a SciPy sparse NRM. The two
    are different code paths, and the audit found one comparison that had
    silently compared a dense matrix against a sparse one. They must agree
    elementwise.
    """

    def test_dense_and_sparse_nrm_agree(self):
        for name in ("mrode.ped", "new_lacy.ped", "hartlandclark.ped"):
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                dense = pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="dense")
                sparse = pyp_nrm.fast_a_matrix(ped.pedigree, ped.kw, method="sparse")
                n = len(ped.pedigree)
                for i in range(n):
                    for j in range(n):
                        self.assertAlmostEqual(
                            nrm_value(dense, i, j), nrm_value(sparse, i, j),
                            places=12, msg="A[%d,%d] differs in %s" % (i, j, name))


class TestPedigreeCompleteness(unittest.TestCase):
    """
    Completeness is a descriptive statistic with no closed-form oracle here; it
    is pinned purely as a behavioural regression guard.
    """

    # ``sum``/``nonfounder_sum`` are totals over animals and ``n``/
    # ``nonfounder_n`` are counts, so only these keys are per-animal
    # proportions bounded by 1.
    PROPORTION_KEYS = ("min", "max", "range", "average",
                       "nonfounder_min", "nonfounder_max",
                       "nonfounder_range", "nonfounder_average")

    def test_completeness_proportions_are_bounded(self):
        for name in ("mrode.ped", "hartlandclark.ped", "generations.ped"):
            with self.subTest(pedigree=name):
                ped = load_corpus(name)
                result = pyp_metrics.pedigree_completeness(ped, 3)
                self.assertIsInstance(result, dict)
                self.assertEqual(len(ped.pedigree), result["n"])
                for key in self.PROPORTION_KEYS:
                    self.assertGreaterEqual(float(result[key]), 0.0, key)
                    self.assertLessEqual(float(result[key]), 1.0, key)


if __name__ == "__main__":
    unittest.main()
