"""RC4 an earlier revision: list_likely_same_animals and merge_animals."""
import copy
import os
import tempfile
import unittest
from unittest import mock

from PyPedal.pyp_errors import PyPedalPedigreeStructureError, PyPedalUsageError
from PyPedal.pyp_newclasses import NewAMatrix
from PyPedal.pyp_utils import list_duplicates, list_likely_same_animals

from _pedhelpers import load_corpus_from_path


POJKA, SUMO = 1, 2
FREDDY_KEEP, FREDDY_DROP = 4, 5
NINA, QUINTEN, QUIRINE = 10, 11, 12

FREDDY_ROWS = [
    "1 0 0 Pojka Unknown 01011800",
    "2 0 0 Sumo Unknown 01011800",
    "4 1 2 Freddy S-2015/23 04252015",
    "5 1 2 Freddy S-2015/23 04252015",
    "10 4 0 Nina Unknown 01012018",
    "11 4 0 Quinten Unknown 01012019",
    "12 5 0 Quirine Unknown 01022019",
]


def rows_to_ped(rows, pedformat="asdnub", **overrides):
    tmp = tempfile.mkdtemp(prefix="merge_")
    path = os.path.join(tmp, "merge.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    options = {"messages": "quiet"}
    options.update(overrides)
    return load_corpus_from_path(path, pedformat, **options)


def by_original(ped, oid):
    return next(a for a in ped.pedigree if a.originalID == oid)


def originals(ped):
    return [a.originalID for a in ped.pedigree]


def snapshot(ped):
    return {
        "originals": originals(ped),
        "objects": [id(a) for a in ped.pedigree],
        "slots": {
            a.originalID: (a.sireID, a.damID, a.animalID, a.name, a.bd, a.userField)
            for a in ped.pedigree
        },
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "n": len(ped.pedigree),
        "meta": ped.metadata.num_records,
        "nrm": getattr(ped, "nrm", None),
        "kw": copy.deepcopy(dict(ped.kw)),
    }


def form_nrm(ped):
    ped.nrm = NewAMatrix(dict(ped.kw))
    ped.nrm.kw["messages"] = "quiet"
    ped.nrm.form_a_matrix(ped.pedigree)
    return ped.nrm


class TestFreddyControl(unittest.TestCase):
    def test_both_freddy_records_load_as_distinct_identities(self):
        ped = rows_to_ped(FREDDY_ROWS)
        self.assertEqual(
            [POJKA, SUMO, FREDDY_KEEP, FREDDY_DROP, NINA, QUINTEN, QUIRINE],
            originals(ped),
        )
        self.assertEqual("Freddy", by_original(ped, FREDDY_KEEP).name)
        self.assertEqual("Freddy", by_original(ped, FREDDY_DROP).name)
        self.assertNotEqual(
            by_original(ped, FREDDY_KEEP).animalID,
            by_original(ped, FREDDY_DROP).animalID,
        )
        self.assertEqual([], list_duplicates(ped))

    def test_likely_same_is_strong_when_userfield_is_declared_unique(self):
        ped = rows_to_ped(FREDDY_ROWS, unique_external_field="userField")
        groups = list_likely_same_animals(ped)
        strong = [g for g in groups if g["strength"] == "strong"]
        self.assertEqual(1, len(strong))
        self.assertEqual(
            {FREDDY_KEEP, FREDDY_DROP}, set(strong[0]["animals"]))
        self.assertIn("userField", strong[0]["evidence"])

    def test_delete_of_drop_freddy_refuses_while_quirine_survives(self):
        ped = rows_to_ped(FREDDY_ROWS)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.delete_animals([FREDDY_DROP])
        self.assertEqual(before, snapshot(ped))

    def test_merge_redirects_all_offspring_and_drops_id_5(self):
        ped = rows_to_ped(FREDDY_ROWS)
        form_nrm(ped)
        self.assertTrue(ped.merge_animals(FREDDY_KEEP, FREDDY_DROP))
        self.assertNotIn(FREDDY_DROP, originals(ped))
        self.assertNotIn(FREDDY_DROP, ped.idmap)
        kept = by_original(ped, FREDDY_KEEP)
        self.assertEqual("Freddy", kept.name)
        keep_aid = kept.animalID
        self.assertEqual(keep_aid, by_original(ped, NINA).sireID)
        self.assertEqual(keep_aid, by_original(ped, QUINTEN).sireID)
        self.assertEqual(keep_aid, by_original(ped, QUIRINE).sireID)
        live = {a.animalID for a in ped.pedigree}
        for animal in ped.pedigree:
            if animal.sireID != ped.kw["missing_parent"]:
                self.assertIn(animal.sireID, live)
            if animal.damID != ped.kw["missing_parent"]:
                self.assertIn(animal.damID, live)
        leftover = [
            kid
            for animal in ped.pedigree
            for kid in list(animal.sons) + list(animal.daus) + list(animal.unks)
        ]
        self.assertTrue(all(kid in live for kid in leftover))
        self.assertEqual(6, ped.metadata.num_records)
        self.assertIsNone(ped.nrm)


class TestCandidateStrength(unittest.TestCase):
    def test_same_name_alone_is_not_a_candidate(self):
        rows = [
            "1 0 0 Freddy Unknown 01011800",
            "2 0 0 Other Unknown 01011800",
            "3 2 0 Freddy Unknown 05011990",
        ]
        ped = rows_to_ped(rows)
        self.assertEqual([], list_likely_same_animals(ped))

    def test_same_parents_and_date_without_strong_id_is_heuristic(self):
        rows = [
            "1 0 0 Pojka Unknown 01011800",
            "2 0 0 Sumo Unknown 01011800",
            "4 1 2 Freddy Unknown 04252015",
            "5 1 2 Freddy Unknown 04252015",
        ]
        ped = rows_to_ped(rows)
        groups = list_likely_same_animals(ped)
        self.assertEqual(1, len(groups))
        self.assertEqual("heuristic", groups[0]["strength"])
        self.assertEqual({4, 5}, set(groups[0]["animals"]))
        self.assertIn("name", groups[0]["evidence"])
        self.assertIn("parents", groups[0]["evidence"])

    def test_different_date_is_not_an_exact_strong_match(self):
        rows = [
            "1 0 0 Pojka Unknown 01011800",
            "2 0 0 Sumo Unknown 01011800",
            "4 1 2 Freddy S-2015/23 04252015",
            "5 1 2 Freddy S-2015/99 05012016",
        ]
        ped = rows_to_ped(rows, unique_external_field="userField")
        groups = list_likely_same_animals(ped)
        self.assertFalse(any(g["strength"] == "strong" for g in groups))
        self.assertFalse(
            any(
                set(g["animals"]) == {4, 5} and g["strength"] == "heuristic"
                for g in groups
            )
        )

    def test_different_external_id_is_not_strong(self):
        rows = [
            "1 0 0 Pojka Unknown 01011800",
            "2 0 0 Sumo Unknown 01011800",
            "4 1 2 Freddy S-2015/23 04252015",
            "5 1 2 Freddy S-2015/99 04252015",
        ]
        ped = rows_to_ped(rows, unique_external_field="userField")
        groups = list_likely_same_animals(ped)
        strong = [g for g in groups if g["strength"] == "strong"]
        self.assertEqual([], strong)
        heuristic = [g for g in groups if g["strength"] == "heuristic"]
        self.assertEqual(1, len(heuristic))
        self.assertEqual({4, 5}, set(heuristic[0]["animals"]))

    def test_undeclared_userfield_is_not_treated_as_registration(self):
        ped = rows_to_ped(FREDDY_ROWS)
        groups = list_likely_same_animals(ped)
        self.assertFalse(any(g["strength"] == "strong" for g in groups))
        heuristic = [g for g in groups if g["strength"] == "heuristic"]
        self.assertEqual(1, len(heuristic))

    def test_detection_does_not_mutate(self):
        ped = rows_to_ped(FREDDY_ROWS, unique_external_field="userField")
        before = snapshot(ped)
        list_likely_same_animals(ped)
        self.assertEqual(before, snapshot(ped))


class TestMergeConflicts(unittest.TestCase):
    def test_conflicting_sire_refuses_without_mutation(self):
        rows = [
            "1 0 0 A Unknown 01011800",
            "2 0 0 B Unknown 01011800",
            "3 0 0 C Unknown 01011800",
            "4 1 2 Keep Unknown 04252015",
            "5 3 2 Drop Unknown 04252015",
        ]
        ped = rows_to_ped(rows)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.merge_animals(4, 5)
        self.assertEqual(before, snapshot(ped))

    def test_conflicting_dam_refuses_without_mutation(self):
        rows = [
            "1 0 0 A Unknown 01011800",
            "2 0 0 B Unknown 01011800",
            "3 0 0 C Unknown 01011800",
            "4 1 2 Keep Unknown 04252015",
            "5 1 3 Drop Unknown 04252015",
        ]
        ped = rows_to_ped(rows)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.merge_animals(4, 5)
        self.assertEqual(before, snapshot(ped))

    def test_conflicting_birth_date_refuses_without_mutation(self):
        rows = [
            "1 0 0 Pojka Unknown 01011800",
            "2 0 0 Sumo Unknown 01011800",
            "4 1 2 Freddy Unknown 04252015",
            "5 1 2 Freddy Unknown 05012016",
        ]
        ped = rows_to_ped(rows)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.merge_animals(4, 5)
        self.assertEqual(before, snapshot(ped))

    def test_ancestor_relationship_refuses_without_mutation(self):
        ped = rows_to_ped(FREDDY_ROWS)
        before = snapshot(ped)
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.merge_animals(FREDDY_KEEP, NINA)
        self.assertEqual(before, snapshot(ped))
        with self.assertRaises(PyPedalPedigreeStructureError):
            ped.merge_animals(NINA, FREDDY_KEEP)
        self.assertEqual(before, snapshot(ped))

    def test_keep_equals_drop_is_usage_error(self):
        ped = rows_to_ped(FREDDY_ROWS)
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            ped.merge_animals(FREDDY_KEEP, FREDDY_KEEP)
        self.assertEqual(before, snapshot(ped))

    def test_missing_keep_field_copies_known_drop_name(self):
        rows = [
            "1 0 0 Pojka Unknown 01011800",
            "2 0 0 Sumo Unknown 01011800",
            "4 1 2 Unknown_Name Unknown 04252015",
            "5 1 2 Freddy Unknown 04252015",
        ]
        ped = rows_to_ped(rows)
        self.assertTrue(ped.merge_animals(4, 5))
        self.assertEqual("Freddy", by_original(ped, 4).name)
        self.assertNotIn(5, originals(ped))

    def test_fa_is_not_copied_without_per_animal_provenance(self):
        ped = rows_to_ped(FREDDY_ROWS)
        keep = by_original(ped, FREDDY_KEEP)
        drop = by_original(ped, FREDDY_DROP)
        keep.fa = 0.0
        drop.fa = 0.42
        self.assertTrue(ped.merge_animals(FREDDY_KEEP, FREDDY_DROP))
        self.assertEqual(0.0, by_original(ped, FREDDY_KEEP).fa)


class TestMergeNrmInvalidation(unittest.TestCase):
    def test_merge_nulls_nrm_and_does_not_auto_reform(self):
        ped = rows_to_ped(FREDDY_ROWS)
        form_nrm(ped)
        with mock.patch.object(NewAMatrix, "form_a_matrix") as spy:
            ped.kw["form_nrm"] = True
            self.assertTrue(ped.merge_animals(FREDDY_KEEP, FREDDY_DROP))
        self.assertIsNone(ped.nrm)
        self.assertFalse(spy.called)
        form_nrm(ped)
        matrix = ped.nrm.nrm
        try:
            shape = matrix.shape
        except AttributeError:
            shape = matrix.get_shape()
        self.assertEqual((6, 6), shape)
