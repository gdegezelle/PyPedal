"""Dict-compatible analysis result types (4.1-A2).

Returned objects remain dicts. Existing key access is unchanged; named
properties are convenience accessors only.
"""
import copy
import json
import pickle
import unittest

from _pedhelpers import chdir_tmp, load_corpus

from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_app import _format_inbreeding
from PyPedal.pyp_results import (
    EffectiveFoundersResult,
    InbreedingResult,
    MatingCoIGroupResult,
)

LACY_APPENDIX_A_FE = 2.909090909090909


def _assert_mapping_compat(test, result, expected, cls):
    test.assertIsInstance(result, dict)
    test.assertIsInstance(result, cls)
    test.assertEqual(result, expected)
    test.assertEqual(dict(result), expected)
    test.assertEqual(len(result), len(expected))
    test.assertEqual(list(result), list(expected))
    for key, value in expected.items():
        test.assertIn(key, result)
        test.assertEqual(result[key], value)
        test.assertEqual(result.get(key), value)
    copied = copy.copy(result)
    test.assertIsInstance(copied, cls)
    test.assertEqual(copied, expected)
    deep = copy.deepcopy(result)
    test.assertIsInstance(deep, cls)
    test.assertEqual(deep, expected)
    pickled = pickle.loads(pickle.dumps(result))
    test.assertIsInstance(pickled, cls)
    test.assertEqual(pickled, expected)


class TestInbreedingResult(unittest.TestCase):
    def test_mapping_and_property_contract(self):
        with chdir_tmp():
            result = pyp_nrm.inbreeding(
                load_corpus("mrode.ped"), method="tabular", output=False)
        self.assertEqual(0.125, result["fx"][5])
        self.assertEqual(0.125, result.fx[5])
        self.assertIs(result.fx, result["fx"])
        self.assertIs(result.metadata, result["metadata"])
        self.assertIn("all", result.metadata)
        self.assertNotIn("rel_dict", result)
        self.assertIsNone(result.rel_dict)
        expected = dict(result)
        _assert_mapping_compat(self, result, expected, InbreedingResult)
        self.assertEqual(
            json.dumps(result, sort_keys=True),
            json.dumps(expected, sort_keys=True),
        )
        deep = copy.deepcopy(result)
        deep.fx[5] = 1.0
        self.assertEqual(0.125, result.fx[5])

    def test_rel_dict_is_present_only_when_requested(self):
        with chdir_tmp():
            with_rels = pyp_nrm.inbreeding(
                load_corpus("mrode.ped"), method="tabular", rels=1,
                output=False)
            without = pyp_nrm.inbreeding(
                load_corpus("mrode.ped"), method="tabular", output=False)
        self.assertIn("rel_dict", with_rels)
        self.assertIs(with_rels.rel_dict, with_rels["rel_dict"])
        self.assertGreater(with_rels.rel_dict["r_count"], 0)
        self.assertNotIn("rel_dict", without)
        self.assertIsNone(without.rel_dict)
        self.assertIsInstance(with_rels, InbreedingResult)
        self.assertIsInstance(without, InbreedingResult)

    def test_gui_formatter_still_accepts_the_dict_subclass(self):
        with chdir_tmp():
            result = pyp_nrm.inbreeding(
                load_corpus("mrode.ped"), method="tabular", output=False)
        text = _format_inbreeding(result)
        self.assertIn("0.125000", text)

    def test_native_methods_all_return_inbreeding_result(self):
        for method in ("tabular", "vanraden", "meu_luo", "mod_meu_luo"):
            with self.subTest(method=method):
                with chdir_tmp():
                    result = pyp_nrm.inbreeding(
                        load_corpus("mrode.ped"), method=method, output=False)
                self.assertIsInstance(result, InbreedingResult)
                self.assertIsInstance(result, dict)
                self.assertEqual(0.125, result["fx"][5])


class TestEffectiveFoundersResult(unittest.TestCase):
    def test_both_entry_points_share_the_class_and_keys(self):
        with chdir_tmp():
            scalable = pyp_metrics.effective_founders_lacy(
                load_corpus("new_lacy.ped"), output=False)
            dense = pyp_metrics.a_effective_founders_lacy(
                load_corpus("new_lacy.ped"), output=False)
        for result in (scalable, dense):
            self.assertIsInstance(result, EffectiveFoundersResult)
            self.assertEqual(
                ["fa_animal_count", "fa_founder_count",
                 "fa_descendant_count", "fa_effective_founders"],
                list(result),
            )
            self.assertAlmostEqual(
                LACY_APPENDIX_A_FE, result["fa_effective_founders"], places=12)
            self.assertAlmostEqual(
                LACY_APPENDIX_A_FE, result.fa_effective_founders, places=12)
            self.assertEqual(
                result.fa_effective_founders, result["fa_effective_founders"])
            self.assertEqual(result.fa_founder_count, result["fa_founder_count"])
            self.assertEqual(result.fa_animal_count, result["fa_animal_count"])
            self.assertEqual(
                result.fa_descendant_count, result["fa_descendant_count"])
            _assert_mapping_compat(self, result, dict(result),
                                   EffectiveFoundersResult)
            self.assertEqual(
                json.dumps(result, sort_keys=True),
                json.dumps(dict(result), sort_keys=True),
            )

    def test_output_true_and_false_return_the_same_class(self):
        with chdir_tmp():
            on = pyp_metrics.effective_founders_lacy(load_corpus("new_lacy.ped"))
            off = pyp_metrics.effective_founders_lacy(
                load_corpus("new_lacy.ped"), output=False)
        self.assertIsInstance(on, EffectiveFoundersResult)
        self.assertIsInstance(off, EffectiveFoundersResult)
        self.assertEqual(on["fa_effective_founders"],
                         off["fa_effective_founders"])
        self.assertTrue(on)
        self.assertTrue(off)


class TestMatingCoIGroupResult(unittest.TestCase):
    def test_mapping_and_property_contract(self):
        with chdir_tmp():
            result = pyp_metrics.mating_coi_group(
                [(4, 3), (1, 2)], load_corpus("mrode.ped"))
        self.assertIsInstance(result, MatingCoIGroupResult)
        self.assertAlmostEqual(0.125, result["matings"][(4, 3)], places=12)
        self.assertAlmostEqual(0.125, result.matings[(4, 3)], places=12)
        self.assertIs(result.matings, result["matings"])
        self.assertIs(result.metadata, result["metadata"])
        self.assertEqual(0.0, result["matings"][(1, 2)])
        _assert_mapping_compat(self, result, dict(result), MatingCoIGroupResult)
        with self.assertRaises(TypeError):
            json.dumps(result)
        with self.assertRaises(TypeError):
            json.dumps(dict(result))


if __name__ == "__main__":
    unittest.main()
