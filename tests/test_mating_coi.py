"""RC4 an earlier revision -- read-only ``mating_coi``.

Prospective offspring inbreeding is A_ij / 2 for distinct parents and
(1 + F_i) / 2 for self-mating. Production calculation never inserts a
phantom child. Invalid IDs and unsupported gens raise PyPedalUsageError
rather than returning -999.9.
"""
import copy
import os
import unittest

from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_newclasses import NewAMatrix

from _pedhelpers import owned_temp_dir, chdir_tmp, load_corpus, load_corpus_from_path, load_example, load_griffon_1871_1890

STUD_ROWS = [
    "100 0 0",
    "200 0 0",
    "300 100 200",
]


def studbook(**overrides):
    tmp = owned_temp_dir(prefix="mating_stud_")
    path = os.path.join(tmp, "stud.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(STUD_ROWS) + "\n")
    return load_corpus_from_path(path, "asd", **overrides)


def write_pedigree(text, pedformat="asd", **overrides):
    tmp = owned_temp_dir(prefix="mating_fix_")
    path = os.path.join(tmp, "fixture.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return load_corpus_from_path(path, pedformat, **overrides)


def form_nrm(ped):
    ped.nrm = NewAMatrix(dict(ped.kw))
    ped.nrm.kw["messages"] = "quiet"
    ped.nrm.form_a_matrix(ped.pedigree)
    return ped.nrm


def nrm_payload(ped):
    nrm = getattr(ped, "nrm", None)
    if nrm is None:
        return None
    matrix = getattr(nrm, "nrm", None)
    if matrix is None:
        return ("attached", id(nrm), None)
    if hasattr(matrix, "toarray"):
        data = matrix.toarray().tolist()
    else:
        data = [list(row) for row in matrix]
    return ("attached", id(nrm), data)


def snapshot(ped):
    return {
        "n": len(ped.pedigree),
        "objects": [id(animal) for animal in ped.pedigree],
        "animal_ids": [int(animal.animalID) for animal in ped.pedigree],
        "original_ids": [animal.originalID for animal in ped.pedigree],
        "sire_ids": [animal.sireID for animal in ped.pedigree],
        "dam_ids": [animal.damID for animal in ped.pedigree],
        "renumbered_ids": [animal.renumberedID for animal in ped.pedigree],
        "fa": [animal.fa for animal in ped.pedigree],
        "founder": [animal.founder for animal in ped.pedigree],
        "gen": [animal.gen for animal in ped.pedigree],
        "igen": [animal.igen for animal in ped.pedigree],
        "idmap": dict(ped.idmap),
        "backmap": dict(ped.backmap),
        "namemap": dict(ped.namemap),
        "namebackmap": dict(ped.namebackmap),
        "offspring": [
            (
                animal.originalID,
                dict(animal.sons),
                dict(animal.daus),
                dict(animal.unks),
            )
            for animal in ped.pedigree
        ],
        "meta": ped.metadata.num_records,
        "kw": copy.deepcopy(dict(ped.kw)),
        "nrm": nrm_payload(ped),
        "f_computed": ped.kw.get("f_computed"),
        "pedigree_is_renumbered": ped.kw.get("pedigree_is_renumbered"),
    }


class TestMrodeScientificControls(unittest.TestCase):
    def test_unrelated_founders_are_zero(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(1, 2, ped)
            got = pyp_metrics.mating_coi(1, 2, ped)
        self.assertEqual(0.0, a_ij)
        self.assertEqual(0.0, got)
        self.assertIsInstance(got, float)

    def test_mrode_1_by_6_is_one_eighth(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(1, 6, ped)
            got = pyp_metrics.mating_coi(1, 6, ped)
        self.assertAlmostEqual(0.25, a_ij, places=12)
        self.assertAlmostEqual(0.125, got, places=12)

    def test_half_siblings_3_by_4(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(3, 4, ped)
            got = pyp_metrics.mating_coi(3, 4, ped)
        self.assertAlmostEqual(0.25, a_ij, places=12)
        self.assertAlmostEqual(0.125, got, places=12)

    def test_parent_offspring_1_by_3(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(1, 3, ped)
            got = pyp_metrics.mating_coi(1, 3, ped)
        self.assertAlmostEqual(0.5, a_ij, places=12)
        self.assertAlmostEqual(0.25, got, places=12)

    def test_inbred_parent_pair_is_not_a_simple_degree(self):
        """3 × 5 is parent × inbred offspring: A_35 = 0.625, not 0.5."""
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(3, 5, ped)
            got = pyp_metrics.mating_coi(3, 5, ped)
        self.assertAlmostEqual(0.625, a_ij, places=12)
        self.assertAlmostEqual(0.3125, got, places=12)
        self.assertNotAlmostEqual(0.25, got, places=12)

    def test_half_founder_parent_offspring_matches_a_matrix(self):
        ped = load_corpus("mrode.ped")
        animal_4 = ped.pedigree[3]
        self.assertEqual(0, int(animal_4.damID))
        self.assertNotEqual("y", animal_4.founder)
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(1, 4, ped)
            got = pyp_metrics.mating_coi(1, 4, ped)
        self.assertAlmostEqual(0.5, a_ij, places=12)
        self.assertAlmostEqual(0.25, got, places=12)

    def test_half_founder_unrelated_to_the_unknown_side(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(2, 4, ped)
            got = pyp_metrics.mating_coi(2, 4, ped)
        self.assertEqual(0.0, a_ij)
        self.assertEqual(0.0, got)


class TestFullSiblingControl(unittest.TestCase):
    def test_full_sibs_are_one_quarter(self):
        ped = write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 1 2\n")
        form_nrm(ped)
        with chdir_tmp():
            a_ij = pyp_metrics.relationship(3, 4, ped)
            nrm_ij = pyp_nrm._matrix_value(ped.nrm.nrm, 2, 3)
            got = pyp_metrics.mating_coi(3, 4, ped)
        self.assertAlmostEqual(0.5, a_ij, places=12)
        self.assertAlmostEqual(0.5, nrm_ij, places=12)
        self.assertAlmostEqual(0.25, got, places=12)


class TestSelfMating(unittest.TestCase):
    def test_founder_self_is_one_half(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 1, ped)
        self.assertAlmostEqual(0.5, got, places=12)

    def test_inbred_mrode_5_self_is_not_relationship_shortcut(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            shortcut = pyp_metrics.relationship(5, 5, ped)
            got = pyp_metrics.mating_coi(5, 5, ped)
            fx = pyp_nrm.inbreeding(ped, method="tabular", output=False)["fx"][5]
        self.assertEqual(1.0, shortcut)
        self.assertAlmostEqual(0.125, fx, places=12)
        self.assertAlmostEqual(0.5625, got, places=12)
        self.assertNotAlmostEqual(shortcut / 2.0, got, places=12)

    def test_self_matches_attached_nrm_diagonal(self):
        ped = load_corpus("mrode.ped")
        form_nrm(ped)
        a_55 = pyp_nrm._matrix_value(ped.nrm.nrm, 4, 4)
        self.assertAlmostEqual(1.125, a_55, places=12)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(5, 5, ped)
        self.assertAlmostEqual(a_55 / 2.0, got, places=12)
        self.assertAlmostEqual(0.5625, got, places=12)

    def test_self_without_nrm_matches_one_plus_f_over_two(self):
        ped = load_corpus("mrode.ped")
        self.assertIsNone(getattr(ped, "nrm", None))
        with chdir_tmp():
            got = pyp_metrics.mating_coi(5, 5, ped)
        self.assertIsNone(getattr(ped, "nrm", None))
        self.assertAlmostEqual(0.5625, got, places=12)


class TestGensContract(unittest.TestCase):
    def test_gens_zero_and_minus_one_agree(self):
        ped = load_corpus("mrode.ped")
        with chdir_tmp():
            a = pyp_metrics.mating_coi(1, 6, ped, 0)
            b = pyp_metrics.mating_coi(1, 6, ped, -1)
        self.assertAlmostEqual(a, b, places=12)
        self.assertAlmostEqual(0.125, a, places=12)

    def test_gens_one_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError) as raised:
            pyp_metrics.mating_coi(1, 6, ped, 1)
        self.assertNotIn("-999.9", str(raised.exception))
        self.assertEqual(before, snapshot(ped))

    def test_gens_three_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(1, 6, ped, 3)
        self.assertEqual(before, snapshot(ped))


class TestIdValidation(unittest.TestCase):
    def test_unknown_current_id_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(99999, 1, ped)
        self.assertEqual(before, snapshot(ped))

    def test_zero_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(0, 1, ped)

    def test_negative_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(-1, 1, ped)

    def test_out_of_range_raises(self):
        ped = load_corpus("mrode.ped")
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(7, 1, ped)

    def test_non_integral_raises(self):
        ped = load_corpus("mrode.ped")
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(1.5, 2, ped)
        self.assertEqual(before, snapshot(ped))

    def test_original_id_is_not_translated(self):
        ped = studbook()
        self.assertEqual({100: 1, 200: 2, 300: 3}, ped.idmap)
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(100, 200, ped)
        self.assertEqual(before, snapshot(ped))

    def test_current_ids_on_the_studbook_compute(self):
        ped = studbook()
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped)
        self.assertEqual(0.0, got)

    def test_unrenumbered_pedigree_raises(self):
        ped = load_corpus(
            "mrode.ped",
            renumber=False,
            pedigree_is_renumbered=False,
        )
        self.assertFalse(ped.kw.get("pedigree_is_renumbered"))
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(1, 2, ped)
        self.assertEqual(before, snapshot(ped))


class TestSexIsNotRequired(unittest.TestCase):
    def test_same_sex_pair_is_a_valid_calculation(self):
        ped = write_pedigree("1 0 0 M\n2 0 0 M\n3 1 2 F\n", pedformat="asdx")
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 2, ped)
        self.assertEqual(0.0, got)


class TestZeroMutation(unittest.TestCase):
    def test_distinct_mating_leaves_state(self):
        ped = load_corpus("mrode.ped")
        form_nrm(ped)
        before = snapshot(ped)
        with chdir_tmp():
            pyp_metrics.mating_coi(1, 6, ped)
        self.assertEqual(before, snapshot(ped))

    def test_self_mating_leaves_state_and_does_not_form_nrm(self):
        ped = load_corpus("mrode.ped")
        self.assertIsNone(getattr(ped, "nrm", None))
        before = snapshot(ped)
        with chdir_tmp():
            pyp_metrics.mating_coi(5, 5, ped)
        self.assertEqual(before, snapshot(ped))
        self.assertIsNone(getattr(ped, "nrm", None))

    def test_error_paths_leave_state(self):
        ped = load_corpus("mrode.ped")
        form_nrm(ped)
        before = snapshot(ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(99999, 1, ped)
        with self.assertRaises(PyPedalUsageError):
            pyp_metrics.mating_coi(1, 2, ped, gens=2)
        self.assertEqual(before, snapshot(ped))

    def test_one_hundred_repeated_calls_leave_state(self):
        ped = load_corpus("mrode.ped")
        form_nrm(ped)
        before = snapshot(ped)
        with chdir_tmp():
            for _ in range(120):
                pyp_metrics.mating_coi(1, 6, ped)
                pyp_metrics.mating_coi(5, 5, ped)
        self.assertEqual(before, snapshot(ped))

    def test_no_phantom_animal_on_studbook(self):
        ped = studbook()
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(1, 3, ped)
        self.assertAlmostEqual(0.25, got, places=12)
        self.assertEqual(before, snapshot(ped))
        self.assertEqual(3, len(ped.pedigree))


class TestHypotheticalChildCrossCheck(unittest.TestCase):
    def test_isolated_child_inbreeding_matches_mating_coi(self):
        original = load_corpus("mrode.ped")
        with chdir_tmp():
            expected = pyp_metrics.mating_coi(1, 6, original)
        child = write_pedigree(
            "# independent copy with a hypothetical child\n"
            "1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n7 1 6\n"
        )
        with chdir_tmp():
            fx = pyp_nrm.inbreeding(child, method="tabular", output=False)["fx"]
        self.assertAlmostEqual(expected, float(fx[7]), places=12)
        self.assertEqual(6, len(original.pedigree))

    def test_full_sib_child_cross_check(self):
        parents = write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 1 2\n")
        with chdir_tmp():
            expected = pyp_metrics.mating_coi(3, 4, parents)
        child = write_pedigree("1 0 0\n2 0 0\n3 1 2\n4 1 2\n5 3 4\n")
        with chdir_tmp():
            fx = pyp_nrm.inbreeding(child, method="tabular", output=False)["fx"]
        self.assertAlmostEqual(0.25, expected, places=12)
        self.assertAlmostEqual(expected, float(fx[5]), places=12)


class TestNrmPolicy(unittest.TestCase):
    def test_uses_attached_nrm_without_forming_another(self):
        ped = load_corpus("mrode.ped")
        form_nrm(ped)
        nrm_id = id(ped.nrm)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(3, 4, ped)
        self.assertIs(ped.nrm, ped.nrm)
        self.assertEqual(nrm_id, id(ped.nrm))
        self.assertAlmostEqual(0.125, got, places=12)

    def test_none_nrm_does_not_attach_one(self):
        ped = load_corpus("mrode.ped")
        self.assertIsNone(getattr(ped, "nrm", None))
        with chdir_tmp():
            pyp_metrics.mating_coi(1, 2, ped)
        self.assertIsNone(getattr(ped, "nrm", None))


class TestPerformanceSanity(unittest.TestCase):
    def test_repeated_pair_calls_do_not_form_dense_nrm(self):
        rows = ["%s 0 0" % i for i in range(1, 81)]
        rows.extend("%s %s 0" % (80 + i, i) for i in range(1, 41))
        ped = write_pedigree("\n".join(rows) + "\n")
        self.assertIsNone(getattr(ped, "nrm", None))
        with chdir_tmp():
            for i in range(1, 21):
                pyp_metrics.mating_coi(i, i + 1, ped)
        self.assertIsNone(getattr(ped, "nrm", None))
        self.assertEqual(120, len(ped.pedigree))


class TestGriffonDoesNotGrow(unittest.TestCase):
    def test_griffon_pair_does_not_append(self):
        ped = load_griffon_1871_1890()
        ids = [int(a.animalID) for a in ped.pedigree]
        before = snapshot(ped)
        with chdir_tmp():
            got = pyp_metrics.mating_coi(ids[0], ids[1], ped)
        self.assertIsInstance(got, float)
        self.assertNotEqual(-999.9, got)
        self.assertEqual(before, snapshot(ped))
