import os
import unittest

from PyPedal import pyp_metrics, pyp_nrm
from PyPedal.pyp_app import _format_inbreeding

from _pedhelpers import load_example, load_griffon_test_small


class TestCorePedigreeWorkflow(unittest.TestCase):
    def test_load_lacy_pedigree(self):
        ped = load_example(
            "new_lacy.ped",
            {
                "pedformat": "asd",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
                "pedname": "Lacy test",
            },
        )
        self.assertEqual(len(ped.pedigree), 7)
        self.assertGreater(ped.metadata.num_unique_founders, 0)
        self.assertTrue(ped.metadata.stringme())

    def test_load_pedigree_alias(self):
        ped = load_example(
            "new_graphics.ped",
            {
                "pedformat": "asdgy",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
            },
        )
        self.assertTrue(ped)
        self.assertGreater(len(ped.pedigree), 0)

    def test_effective_founders_lacy(self):
        # load_example, not load_pedigree: effective_founders_lacy writes
        # <filetag>_fe_lacy.dat, and filetag follows the pedfile path.
        ped = load_example(
            "new_lacy.ped",
            {
                "pedformat": "asd",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
                "pedigree_is_renumbered": True,
            },
        )
        result = pyp_metrics.effective_founders_lacy(ped)
        fe = result.get("fa_effective_founders", -999.9)
        self.assertGreater(fe, 0)
        self.assertEqual(2.91, round(fe, 2))

    def test_effective_founders_griffon_small(self):
        ped = load_griffon_test_small()
        result = pyp_metrics.effective_founders_lacy(ped)
        fe = result.get("fa_effective_founders", -999.9)
        self.assertNotEqual(fe, -999.9)
        self.assertGreater(fe, 0)

    def test_inbreeding_returns_coefficients(self):
        ped = load_example(
            "new_lacy.ped",
            {
                "pedformat": "asd",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
            },
        )
        result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
        if isinstance(result, tuple):
            result = result[0]
        self.assertIsInstance(result, dict)
        self.assertTrue("fx" in result or "metadata" in result)
        summary = _format_inbreeding(result)
        self.assertIn("Inbreeding", summary)
        self.assertIn("f_count", summary)
        self.assertIn("Coefficients by animal", summary)
        self.assertTrue(all(float(v) == 0.0 for v in result["fx"].values()))
        self.assertIn("No inbreeding in this pedigree", summary)
        self.assertNotIn("new_lacy.ped", summary)
        self.assertNotIn("mrode.ped", summary)
        self.assertNotIn("hartlandclark.ped", summary)

    def test_inbreeding_mrode_has_nonzero_coi(self):
        ped = load_example(
            "mrode.ped",
            {
                "pedformat": "asd",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
                "pedigree_is_renumbered": True,
            },
        )
        result = pyp_nrm.inbreeding(ped, method="meu_luo", output=False)
        fx = {int(k): float(v) for k, v in result["fx"].items()}
        self.assertAlmostEqual(fx[5], 0.125, places=3)
        self.assertGreater(fx[6], 0.0)

    def _load_mrode(self):
        # load_example, not load_pedigree: fast_a_coefficients and
        # a_coefficients write three <filetag>_*.dat files apiece.
        return load_example(
            "mrode.ped",
            {
                "pedformat": "asd",
                "sepchar": " ",
                "messages": "quiet",
                "renumber": True,
                "pedigree_is_renumbered": True,
            },
        )

    def test_inbreeding_vanraden_mrode(self):
        ped = self._load_mrode()
        result = pyp_nrm.inbreeding(ped, method="vanraden", output=False)
        fx = {int(k): float(v) for k, v in result["fx"].items()}
        self.assertAlmostEqual(fx[5], 0.125, places=3)
        self.assertGreater(fx[6], 0.0)

    def test_inbreeding_tabular_mrode(self):
        ped = self._load_mrode()
        result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
        fx = {int(k): float(v) for k, v in result["fx"].items()}
        self.assertAlmostEqual(fx[5], 0.125, places=3)
        self.assertGreater(fx[6], 0.0)

    def test_inbreeding_rels_does_not_crash(self):
        ped = self._load_mrode()
        result = pyp_nrm.inbreeding(ped, method="tabular", rels=1, output=False)
        self.assertIn("rel_dict", result)
        self.assertGreater(result["rel_dict"]["r_count"], 0)

    def test_find_ancestors_accepts_legacy_accumulator(self):
        from PyPedal import pyp_network

        ped = self._load_mrode()
        graph = pyp_network.ped_to_graph(ped)
        two_arg = pyp_network.find_ancestors(graph, 5)
        three_arg = pyp_network.find_ancestors(graph, 5, [])
        self.assertTrue(two_arg)
        self.assertEqual(sorted(two_arg), sorted(three_arg))
        limited = pyp_network.find_ancestors_g(graph, 5, {}, 1)
        self.assertIsInstance(limited, dict)
        self.assertTrue(limited)

    def test_relationship_and_coefficients_mrode(self):
        ped = self._load_mrode()
        self.assertAlmostEqual(pyp_metrics.relationship(4, 3, ped), 0.25, places=3)
        self.assertAlmostEqual(pyp_metrics.relationship(5, 5, ped), 1.0, places=3)
        coeffs = pyp_metrics.fast_a_coefficients(ped, storage="sparse")
        self.assertAlmostEqual(float(coeffs[5]), 0.125, places=3)

    def test_package_import_and_db_aliases(self):
        import PyPedal
        from PyPedal import pyp_db

        self.assertTrue(hasattr(PyPedal, "pyp_nrm"))
        self.assertTrue(callable(pyp_db.connectToDatabase))
        self.assertTrue(callable(pyp_db.doesTableExist))
        self.assertTrue(callable(pyp_db.deleteTable))

    def test_a_coefficients_sparse_mrode(self):
        ped = self._load_mrode()
        ped.kw["filetag"] = os.path.join(os.path.dirname(__file__), "_tmp_mrode")
        try:
            coeffs = pyp_metrics.a_coefficients(ped)
            self.assertAlmostEqual(float(coeffs[5]), 0.125, places=3)
        finally:
            for suffix in (
                "_rel_to_pop_.dat",
                "_population_coefficients_.dat",
                "_individual_coefficients_.dat",
            ):
                path = ped.kw["filetag"] + suffix
                if os.path.exists(path):
                    os.remove(path)

    def test_founders_by_year_griffon(self):
        from PyPedal import pyp_demog

        ped = load_griffon_test_small()
        by_year = pyp_demog.founders_by_year(ped)
        self.assertTrue(by_year)
        self.assertGreater(sum(by_year.values()), 0)

    def test_griffon_birth_year_from_mmddyyyy(self):
        ped = load_griffon_test_small()
        years = {animal.by for animal in ped.pedigree if animal.by is not None}
        self.assertTrue(any(year >= 1900 for year in years))


if __name__ == "__main__":
    unittest.main()
