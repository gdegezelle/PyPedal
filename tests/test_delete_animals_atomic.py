"""RC4 an earlier revision: atomic delete_animals, typed refusal, NRM invalidation."""
import os
import time
import unittest
from unittest import mock

from PyPedal.pyp_errors import PyPedalPedigreeStructureError, PyPedalUsageError
from PyPedal.pyp_newclasses import NewAMatrix
from PyPedal.pyp_nrm import _matrix_value

from _pedhelpers import owned_temp_dir, load_corpus_from_path
from test_delete_animals import (
    ALL,
    FOUNDER,
    HALF_DAM,
    LEAF,
    ORDINARY,
    isolate,
    originals,
    rows_to_ped,
    snapshot_ids,
)


def form_nrm(ped):
    ped.nrm = NewAMatrix(dict(ped.kw))
    ped.nrm.kw["messages"] = "quiet"
    ped.nrm.form_a_matrix(ped.pedigree)
    return ped.nrm


def nrm_shape(ped):
    matrix = ped.nrm.nrm
    try:
        return matrix.shape
    except AttributeError:
        return matrix.get_shape()


class TestEmptyAndDeleteAll(unittest.TestCase):
    def test_empty_list_is_true_and_does_not_drop_an_existing_nrm(self):
        ped = isolate(rows_to_ped())
        form_nrm(ped)
        before = snapshot_ids(ped)
        marker = ped.nrm
        self.assertTrue(ped.delete_animals([]))
        self.assertIs(marker, ped.nrm)
        self.assertEqual(before, snapshot_ids(ped))

    def test_delete_all_leaves_an_empty_usable_pedigree(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.delete_animals(list(ALL)))
        self.assertEqual([], originals(ped))
        self.assertEqual(0, len(ped.pedigree))
        self.assertEqual(0, ped.metadata.num_records)
        self.assertEqual({}, ped.idmap)
        self.assertEqual({}, ped.backmap)
        self.assertIsNone(ped.nrm)

    def test_repeated_delete_of_removed_id_raises_without_further_mutation(self):
        ped = isolate(rows_to_ped())
        self.assertTrue(ped.delete_animals([LEAF]))
        after_first = snapshot_ids(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.delete_animals([LEAF])
        self.assertEqual(after_first, snapshot_ids(ped))


class TestNrmInvalidation(unittest.TestCase):
    def test_successful_delete_nulls_nrm_and_does_not_auto_reform(self):
        ped = rows_to_ped()
        form_nrm(ped)
        self.assertEqual((7, 7), nrm_shape(ped))
        with mock.patch.object(NewAMatrix, "form_a_matrix", wraps=NewAMatrix.form_a_matrix) as spy:
            ped.kw["form_nrm"] = True
            self.assertTrue(ped.delete_animals([LEAF]))
        self.assertIsNone(ped.nrm)
        self.assertFalse(spy.called)

    def test_explicit_reform_matches_the_new_pedigree(self):
        ped = rows_to_ped()
        form_nrm(ped)
        self.assertTrue(ped.delete_animals([LEAF]))
        self.assertIsNone(ped.nrm)
        form_nrm(ped)
        self.assertEqual((6, 6), nrm_shape(ped))
        # Mrode-style founder pair in this fixture: animals 10 and 20 are
        # unrelated after load+renumber; remaining diagonal is 1.0.
        self.assertAlmostEqual(float(_matrix_value(ped.nrm.nrm, 0, 0)), 1.0, places=6)

    def test_refused_parent_delete_does_not_drop_nrm(self):
        ped = rows_to_ped()
        form_nrm(ped)
        marker = ped.nrm
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.delete_animals([ORDINARY])
        self.assertIs(marker, ped.nrm)

    def test_addanimal_invalidates_existing_nrm(self):
        ped = rows_to_ped()
        form_nrm(ped)
        new_id = max(int(k) for k in ped.idmap) + 1
        self.assertTrue(ped.addanimal(new_id, FOUNDER, HALF_DAM))
        self.assertIsNone(ped.nrm)

    def test_delanimal_invalidates_existing_nrm(self):
        ped = isolate(rows_to_ped())
        form_nrm(ped)
        self.assertTrue(ped.delanimal(LEAF))
        self.assertIsNone(ped.nrm)


class TestFormNrmFlagDoesNotAllocate(unittest.TestCase):
    def test_form_nrm_true_does_not_call_form_a_matrix_on_delete(self):
        ped = isolate(rows_to_ped())
        ped.kw["form_nrm"] = True
        with mock.patch.object(NewAMatrix, "form_a_matrix") as spy:
            self.assertTrue(ped.delete_animals([LEAF]))
        self.assertFalse(spy.called)
        self.assertIsNone(ped.nrm)


class TestPerformanceSanity(unittest.TestCase):
    def test_structural_delete_is_linear_and_does_not_form_nrm(self):
        n = 1200
        k = 40
        rows = ["%s 0 0" % i for i in range(1, n + 1)]
        tmp = owned_temp_dir(prefix="del_perf_")
        path = os.path.join(tmp, "wide.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(rows) + "\n")
        ped = load_corpus_from_path(path, "asd")
        isolate(ped)
        ped.kw["form_nrm"] = True
        requested = list(range(n - k + 1, n + 1))
        with mock.patch.object(NewAMatrix, "form_a_matrix") as spy:
            started = time.perf_counter()
            self.assertTrue(ped.delete_animals(requested))
            elapsed = time.perf_counter() - started
        self.assertFalse(spy.called)
        self.assertIsNone(ped.nrm)
        self.assertEqual(n - k, len(ped.pedigree))
        self.assertLess(elapsed, 5.0)
