#!/usr/bin/env python3

###############################################################################
# NAME: pyp_tests.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
#   Provide unit tests for PyPedal.
###############################################################################

"""
Unit tests for the PyPedal library.
"""

import os
import unittest

from . import pyp_metrics, pyp_newclasses


class PyPedalMetricsTestCases(unittest.TestCase):
    def setUp(self):
        """Set up the base path for test files."""
        self.base_path = os.path.join(os.path.dirname(__file__), "examples")

    def testMetricsMinMaxF(self):
        options = {
            "renumber": False,
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "hartlandclark.ped"),
            "pedformat": "asdb",
            "sepchar": " ",
            "pedigree_is_renumbered": True,
            "form_nrm": "1",
        }
        example = pyp_newclasses.NewPedigree(options)
        example.load()
        high_coi, low_coi = pyp_metrics.min_max_f(example, n=5)
        print(high_coi)
        print(low_coi)

    def testMetricsEffectiveFoundersLacy(self):
        options = {
            "renumber": False,
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "new_lacy.ped"),
            "pedformat": "asd",
            "sepchar": " ",
            "pedigree_is_renumbered": True,
        }
        example = pyp_newclasses.NewPedigree(options)
        example.load()

        # Call the function and extract the effective founder number
        result = pyp_metrics.effective_founders_lacy(example)
        fe = result.get("fa_effective_founders", -999.9)  # Default to -999.9 if key not found

        # Assert the result
        self.assertEqual(2.91, round(fe, 2))

    def testMetricsEffectiveFoundersBoichardA(self):
        options = {
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "boichard2a.ped"),
            "pedname": "Boichard Pedigree (Family 1 only)",
            "pedformat": "asdg",
            "pedigree_is_renumbered": True,
            "filetag": "example2a",
        }
        example2a = pyp_newclasses.NewPedigree(options)
        example2a.load()
        fe = pyp_metrics.a_effective_founders_boichard(example2a)
        self.assertEqual(round(4.0, 1), round(fe, 1))

    @unittest.skip("Example file boichard2b.ped is not in this repository")
    def testMetricsEffectiveFoundersBoichardB(self):
        options = {
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "boichard/boichard2b.ped"),
            "pedname": "Boichard Pedigree (Family 2 only)",
            "pedformat": "asdg",
            "pedigree_is_renumbered": True,
            "filetag": "example2b",
        }
        example2b = pyp_newclasses.NewPedigree(options)
        example2b.load()
        fe = pyp_metrics.a_effective_founders_boichard(example2b)
        self.assertEqual(round(2.0, 1), round(fe, 1))

    def testMetricsEffectiveAncestorsDefiniteBoichardA(self):
        options = {
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "boichard2a.ped"),
            "pedname": "Boichard Pedigree (Family 1 only)",
            "pedformat": "asdg",
            "pedigree_is_renumbered": True,
            "filetag": "example2a",
        }
        example2a = pyp_newclasses.NewPedigree(options)
        example2a.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2a)
        self.assertEqual(round(2.0, 2), round(fa, 2))

    @unittest.skip("Example file boichard2b.ped is not in this repository")
    def testMetricsEffectiveAncestorsDefiniteBoichardB(self):
        options = {
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "boichard/boichard2b.ped"),
            "pedname": "Boichard Pedigree (Family 2 only)",
            "pedformat": "asdg",
            "pedigree_is_renumbered": True,
            "filetag": "example2a",
        }
        example2b = pyp_newclasses.NewPedigree(options)
        example2b.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2b)
        self.assertEqual(round(2.0, 2), round(fa, 2))

    def testMetricsEffectiveAncestorsDefiniteBoichardC(self):
        options = {
            "messages": "quiet",
            "pedfile": os.path.join(self.base_path, "boichard2.ped"),
            "pedname": "Boichard Pedigree (Family 1 and 2)",
            "pedformat": "asdg",
            "pedigree_is_renumbered": True,
            "filetag": "example2",
        }
        example2 = pyp_newclasses.NewPedigree(options)
        example2.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2)
        self.assertEqual(round(2.94, 2), round(fa, 2))


class PyPedalNrmTestCases(unittest.TestCase):
    pass


class PyPedalUtilsTestCases(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
