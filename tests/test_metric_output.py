"""Optional analysis ``.dat`` output (4.1-A1).

Default ``output=True`` keeps historical files. ``output=False`` still
computes and returns the same scientific result, but writes none of the
routine's analysis files. Goldens were captured from v4.0.1 before the
writer extraction.
"""
import inspect
import os
import unittest

from _pedhelpers import chdir_tmp, load_corpus

from PyPedal import pyp_metrics, pyp_nrm

GOLDEN = os.path.join(os.path.dirname(__file__), "fixtures", "metric_dat")


def _dat_names(ped):
    parent = os.path.dirname(ped.kw["filetag"])
    base = os.path.basename(ped.kw["filetag"])
    return sorted(
        name for name in os.listdir(parent)
        if name.startswith(base) and name.endswith(".dat")
    )


def _bytes(ped, suffix):
    with open(ped.kw["filetag"] + suffix, "rb") as handle:
        return handle.read()


def _golden(name):
    with open(os.path.join(GOLDEN, name), "rb") as handle:
        return handle.read()


def _strip_timestamp(payload):
    """``effective_founders_lacy`` stamps ``pyp_nice_time()`` into the header."""
    lines = payload.splitlines(keepends=True)
    return b"".join(
        line for line in lines if b"Created by " not in line
    )


class TestOutputFalseSuppressesFiles(unittest.TestCase):
    def test_a_effective_founders_lacy(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            result = pyp_metrics.a_effective_founders_lacy(ped, output=False)
            self.assertIn("fa_effective_founders", result)
            self.assertGreater(result["fa_effective_founders"], 0.0)
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(os.path.exists(ped.kw["filetag"] + "_fe_lacy_.dat"))

    def test_effective_founders_lacy(self):
        with chdir_tmp():
            ped = load_corpus("new_lacy.ped")
            before = _dat_names(ped)
            result = pyp_metrics.effective_founders_lacy(ped, output=False)
            self.assertAlmostEqual(2.91, result["fa_effective_founders"], places=2)
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(os.path.exists(ped.kw["filetag"] + "_fe_lacy.dat"))

    def test_a_effective_founders_boichard(self):
        with chdir_tmp():
            ped = load_corpus("boichard2a.ped", "asdg")
            before = _dat_names(ped)
            value = pyp_metrics.a_effective_founders_boichard(ped, output=False)
            self.assertEqual(4.0, value)
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(os.path.exists(ped.kw["filetag"] + "_fe_boichard_.dat"))

    def test_a_effective_ancestors_definite(self):
        with chdir_tmp():
            ped = load_corpus("boichard2a.ped", "asdg")
            before = _dat_names(ped)
            value = pyp_metrics.a_effective_ancestors_definite(ped, output=False)
            self.assertEqual(2.0, value)
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(
                os.path.exists(ped.kw["filetag"] + "_fa_boichard_definite_.dat"))

    def test_a_effective_ancestors_indefinite(self):
        with chdir_tmp():
            ped = load_corpus("boichard2a.ped", "asdg")
            before = _dat_names(ped)
            bounds = pyp_metrics.a_effective_ancestors_indefinite(
                ped, n=25, output=False)
            self.assertEqual(2, len(bounds))
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(
                os.path.exists(ped.kw["filetag"] + "_fa_boichard_indefinite_.dat"))

    def test_a_coefficients(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            result = pyp_metrics.a_coefficients(ped, output=False)
            self.assertIsInstance(result, dict)
            self.assertEqual(before, _dat_names(ped))
            for suffix in (
                "_rel_to_pop_.dat",
                "_population_coefficients_.dat",
                "_individual_coefficients_.dat",
            ):
                self.assertFalse(os.path.exists(ped.kw["filetag"] + suffix))

    def test_theoretical_ne_from_metadata(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            ok = pyp_metrics.theoretical_ne_from_metadata(ped, output=False)
            self.assertAlmostEqual(ok, 4.8)
            self.assertIsInstance(ok, float)
            self.assertEqual(before, _dat_names(ped))
            self.assertFalse(
                os.path.exists(ped.kw["filetag"] + "_ne_from_metadata_.dat"))

    def test_a_decompose(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            D, T = pyp_nrm.a_decompose(ped, output=False)
            self.assertEqual((6, 6), D.shape)
            self.assertEqual((6, 6), T.shape)
            self.assertEqual(before, _dat_names(ped))

    def test_a_inverse_dnf(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            inverse = pyp_nrm.a_inverse_dnf(ped, output=False)
            self.assertEqual((6, 6), inverse.shape)
            self.assertEqual(before, _dat_names(ped))

    def test_a_inverse_df(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            before = _dat_names(ped)
            inverse = pyp_nrm.a_inverse_df(ped, output=False)
            self.assertEqual((6, 6), inverse.shape)
            self.assertEqual(before, _dat_names(ped))


class TestOutputTrueMatchesDefaultAndGoldens(unittest.TestCase):
    def test_lacy_dense_default_matches_explicit_and_golden(self):
        with chdir_tmp():
            default_ped = load_corpus("mrode.ped")
            default = pyp_metrics.a_effective_founders_lacy(default_ped)
            explicit_ped = load_corpus("mrode.ped")
            explicit = pyp_metrics.a_effective_founders_lacy(
                explicit_ped, output=True)
            self.assertEqual(default, explicit)
            payload = _bytes(explicit_ped, "_fe_lacy_.dat")
            self.assertEqual(payload, _bytes(default_ped, "_fe_lacy_.dat"))
            self.assertEqual(_golden("mrode_fe_lacy_.dat"), payload)

    def test_lacy_scalable_default_matches_explicit(self):
        with chdir_tmp():
            default_ped = load_corpus("new_lacy.ped")
            default = pyp_metrics.effective_founders_lacy(default_ped)
            explicit_ped = load_corpus("new_lacy.ped")
            explicit = pyp_metrics.effective_founders_lacy(
                explicit_ped, output=True)
            self.assertEqual(default, explicit)
            self.assertEqual(
                _strip_timestamp(_bytes(default_ped, "_fe_lacy.dat")),
                _strip_timestamp(_bytes(explicit_ped, "_fe_lacy.dat")),
            )
            self.assertAlmostEqual(
                2.91, default["fa_effective_founders"], places=2)

    def test_boichard_founders_matches_golden(self):
        with chdir_tmp():
            default_ped = load_corpus("boichard2a.ped", "asdg")
            default = pyp_metrics.a_effective_founders_boichard(default_ped)
            explicit_ped = load_corpus("boichard2a.ped", "asdg")
            explicit = pyp_metrics.a_effective_founders_boichard(
                explicit_ped, output=True)
            self.assertEqual(default, explicit)
            payload = _bytes(explicit_ped, "_fe_boichard_.dat")
            self.assertEqual(payload, _bytes(default_ped, "_fe_boichard_.dat"))
            self.assertEqual(_golden("boichard2a_fe_boichard_.dat"), payload)

    def test_boichard_definite_matches_golden(self):
        with chdir_tmp():
            default_ped = load_corpus("boichard2a.ped", "asdg")
            default = pyp_metrics.a_effective_ancestors_definite(default_ped)
            explicit_ped = load_corpus("boichard2a.ped", "asdg")
            explicit = pyp_metrics.a_effective_ancestors_definite(
                explicit_ped, output=True)
            self.assertEqual(default, explicit)
            payload = _bytes(explicit_ped, "_fa_boichard_definite_.dat")
            self.assertEqual(
                payload, _bytes(default_ped, "_fa_boichard_definite_.dat"))
            self.assertEqual(
                _golden("boichard2a_fa_boichard_definite_.dat"), payload)

    def test_boichard_indefinite_matches_golden(self):
        with chdir_tmp():
            default_ped = load_corpus("boichard2a.ped", "asdg")
            default = pyp_metrics.a_effective_ancestors_indefinite(
                default_ped, n=25)
            explicit_ped = load_corpus("boichard2a.ped", "asdg")
            explicit = pyp_metrics.a_effective_ancestors_indefinite(
                explicit_ped, n=25, output=True)
            self.assertEqual(default, explicit)
            payload = _bytes(explicit_ped, "_fa_boichard_indefinite_.dat")
            self.assertEqual(
                payload, _bytes(default_ped, "_fa_boichard_indefinite_.dat"))
            self.assertEqual(
                _golden("boichard2a_fa_boichard_indefinite_.dat"), payload)

    def test_a_coefficients_matches_goldens(self):
        with chdir_tmp():
            default_ped = load_corpus("mrode.ped")
            default = pyp_metrics.a_coefficients(default_ped)
            explicit_ped = load_corpus("mrode.ped")
            explicit = pyp_metrics.a_coefficients(explicit_ped, output=True)
            self.assertEqual(default, explicit)
            for suffix, golden in (
                ("_rel_to_pop_.dat", "mrode_rel_to_pop_.dat"),
                ("_population_coefficients_.dat",
                 "mrode_population_coefficients_.dat"),
                ("_individual_coefficients_.dat",
                 "mrode_individual_coefficients_.dat"),
            ):
                payload = _bytes(explicit_ped, suffix)
                self.assertEqual(payload, _bytes(default_ped, suffix))
                self.assertEqual(_golden(golden), payload)

    def test_theoretical_ne_matches_golden_and_returns_float(self):
        with chdir_tmp():
            default_ped = load_corpus("mrode.ped")
            default = pyp_metrics.theoretical_ne_from_metadata(default_ped)
            explicit_ped = load_corpus("mrode.ped")
            explicit = pyp_metrics.theoretical_ne_from_metadata(
                explicit_ped, output=True)
            self.assertAlmostEqual(default, 4.8)
            self.assertAlmostEqual(explicit, 4.8)
            self.assertEqual(default, explicit)
            payload = _bytes(explicit_ped, "_ne_from_metadata_.dat")
            self.assertEqual(
                payload, _bytes(default_ped, "_ne_from_metadata_.dat"))
            self.assertEqual(_golden("mrode_ne_from_metadata_.dat"), payload)

    def test_nrm_writers_match_goldens(self):
        with chdir_tmp():
            default_ped = load_corpus("mrode.ped")
            pyp_nrm.a_decompose(default_ped)
            pyp_nrm.a_inverse_dnf(default_ped)
            pyp_nrm.a_inverse_df(default_ped)
            explicit_ped = load_corpus("mrode.ped")
            pyp_nrm.a_decompose(explicit_ped, output=True)
            pyp_nrm.a_inverse_dnf(explicit_ped, output=True)
            pyp_nrm.a_inverse_df(explicit_ped, output=True)
            for suffix, golden in (
                ("_a_decompose_d_.dat", "mrode_a_decompose_d_.dat"),
                ("_a_decompose_t_.dat", "mrode_a_decompose_t_.dat"),
                ("_a_inverse_dnf_a_inv.dat", "mrode_a_inverse_dnf_a_inv.dat"),
                ("_a_inverse_dnf_d_inv.dat", "mrode_a_inverse_dnf_d_inv.dat"),
                ("_a_inverse_df_a_inv.dat", "mrode_a_inverse_df_a_inv.dat"),
                ("_a_inverse_df_l.dat", "mrode_a_inverse_df_l.dat"),
                ("_a_inverse_df_d_inv.dat", "mrode_a_inverse_df_d_inv.dat"),
            ):
                payload = _bytes(explicit_ped, suffix)
                self.assertEqual(payload, _bytes(default_ped, suffix), suffix)
                self.assertEqual(_golden(golden), payload, suffix)


class TestOutputFalseKeepsScience(unittest.TestCase):
    def test_lacy_and_boichard_results_match_with_and_without_files(self):
        with chdir_tmp():
            on_lacy = pyp_metrics.effective_founders_lacy(
                load_corpus("new_lacy.ped"))
            off_lacy = pyp_metrics.effective_founders_lacy(
                load_corpus("new_lacy.ped"), output=False)
            self.assertEqual(on_lacy, off_lacy)

            on_boichard = pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"))
            off_boichard = pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"), output=False)
            self.assertEqual(on_boichard, off_boichard)

            on_fa = pyp_metrics.a_effective_ancestors_definite(
                load_corpus("boichard2a.ped", "asdg"))
            off_fa = pyp_metrics.a_effective_ancestors_definite(
                load_corpus("boichard2a.ped", "asdg"), output=False)
            self.assertEqual(on_fa, off_fa)


class TestFastACoefficientsFileIo(unittest.TestCase):
    def test_default_file_io_writes(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            pyp_metrics.fast_a_coefficients(ped)
            self.assertTrue(
                os.path.exists(ped.kw["filetag"] + "_population_coefficients_.dat"))

    def test_file_io_false_still_suppresses_even_when_output_true(self):
        with chdir_tmp():
            ped = load_corpus("mrode.ped")
            ped.kw["file_io"] = False
            before = _dat_names(ped)
            pyp_metrics.fast_a_coefficients(ped, output=True)
            self.assertEqual(before, _dat_names(ped))

    def test_file_io_true_writes_unless_output_false(self):
        with chdir_tmp():
            on_ped = load_corpus("mrode.ped")
            on_ped.kw["file_io"] = True
            on = pyp_metrics.fast_a_coefficients(on_ped)
            self.assertTrue(
                os.path.exists(on_ped.kw["filetag"] + "_population_coefficients_.dat"))

            off_ped = load_corpus("mrode.ped")
            off_ped.kw["file_io"] = True
            off = pyp_metrics.fast_a_coefficients(off_ped, output=False)
            self.assertEqual(on, off)
            self.assertFalse(
                os.path.exists(off_ped.kw["filetag"] + "_population_coefficients_.dat"))


class TestBallouDoesNotWriteInbreedingDat(unittest.TestCase):
    def test_internal_inbreeding_call_suppresses_file(self):
        with chdir_tmp():
            ped = load_corpus("ballou_fig1.ped")
            self.assertFalse(ped.kw.get("f_computed", False))
            fa = pyp_metrics.ballou_ancestral_inbreeding(ped)
            self.assertTrue(fa)
            self.assertTrue(ped.kw.get("f_computed", False))
            self.assertFalse(os.path.exists(ped.kw["filetag"] + "_inbreeding.dat"))
            self.assertNotIn(
                os.path.basename(ped.kw["filetag"]) + "_inbreeding.dat",
                _dat_names(ped),
            )


class TestBoichardOutputIsKeywordOnly(unittest.TestCase):
    def test_output_does_not_become_a_fourth_positional(self):
        with chdir_tmp():
            ped = load_corpus("boichard2a.ped", "asdg")
            with self.assertRaises(TypeError):
                pyp_metrics.a_effective_founders_boichard(ped, None, None, [7, 8])

    def test_output_is_keyword_only_after_reference(self):
        for routine in (
            pyp_metrics.a_effective_founders_boichard,
            pyp_metrics.a_effective_ancestors_definite,
            pyp_metrics.a_effective_ancestors_indefinite,
        ):
            parameters = inspect.signature(routine).parameters
            positional = [
                name for name, param in parameters.items()
                if param.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
            ]
            self.assertEqual(["pedobj", "a", "gen"], positional[:3], routine.__name__)
            self.assertEqual(
                inspect.Parameter.KEYWORD_ONLY, parameters["output"].kind,
                routine.__name__)
            self.assertIs(True, parameters["output"].default, routine.__name__)
            self.assertEqual(
                inspect.Parameter.KEYWORD_ONLY, parameters["reference"].kind,
                routine.__name__)


if __name__ == "__main__":
    unittest.main()
