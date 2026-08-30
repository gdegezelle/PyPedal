"""
The Lacy (1989) oracle against the paper's published values.

Validates the ORACLE, not PyPedal -- production is adjudicated separately, and
only after the instrument is calibrated.

Two independent sets of published numbers are available, which is unusual and
worth using in full:

* **Appendix A** (pp.121-123), one worked pedigree with f_e = 2.91 and
  f_g = 2.18. That pedigree is this repository's ``new_lacy.ped``.
* **Table 1** (p.116), seven rows of (f_e, f_g, r_1, r_2) over two founder
  pairs with varying numbers of offspring. H is also published but is a
  10,000-replicate simulation figure and is not an exact target.

Together they pin the effective founder number, the founder genome equivalent,
and the allele-retention closed form at 28 published values.
"""
import os
import tempfile
import unittest

from _pedhelpers import corpus
from oracles import LACY_MODES, lacy_f_e, lacy_n_half_founders

import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "oracles"))
from oracle_lacy import (                     # noqa: E402
    DEFAULT_MODE,
    f_e_lacy,
    f_g_lacy,
    retention,
)
from oracle_meuwissen_luo import read_pedigree, renumber   # noqa: E402

# Published to two decimals; tolerance is the rounding half-width.
TWO_DP = 5e-3 + 1e-9
THREE_DP = 5e-4 + 1e-9


def _two_founder_pairs(x1, x2):
    """
    Table 1's design: two founder pairs, the first with ``x1`` first-generation
    descendants and the second with ``x2``. f = 4 in every row.
    """
    rows = [(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]
    nid = 5
    for _ in range(x1):
        rows.append((nid, 1, 2))
        nid += 1
    for _ in range(x2):
        rows.append((nid, 3, 4))
        nid += 1
    return rows


class TestThePaperIsTheDefault(unittest.TestCase):

    def test_default_mode_is_phantom(self):
        """
        Lacy p.113 defines half-founder handling, and it is the phantom
        treatment. The oracle's default must be the paper's rule, not the
        historical PyPedal one.
        """
        self.assertEqual("phantom", DEFAULT_MODE)

    def test_all_three_modes_are_still_available(self):
        """
        The wrong treatments are kept as comparators so the change is
        explicable. Removing them would leave nothing to compare against.
        """
        self.assertEqual(("lacy", "half", "phantom"), LACY_MODES)


class TestAppendixA(unittest.TestCase):
    """
    Lacy Appendix A, article pp.121-123, Figure 3 with Tables 2 and 3.
    """

    def test_new_lacy_is_the_appendix_a_pedigree(self):
        """
        Three founders and four descendants, with no half-founder -- so all
        three modes agree here and the fixture cannot discriminate between
        them. That is worth asserting: it is why this pedigree never settled
        .
        """
        self.assertEqual(0, lacy_n_half_founders(corpus("new_lacy.ped")))
        values = set()
        for mode in LACY_MODES:
            f_e, gate = lacy_f_e(corpus("new_lacy.ped"), mode=mode)
            self.assertTrue(gate)
            values.add(round(f_e, 12))
        self.assertEqual(1, len(values))

    def test_published_f_e(self):
        f_e, gate = lacy_f_e(corpus("new_lacy.ped"))
        self.assertAlmostEqual(2.91, f_e, delta=TWO_DP)
        self.assertTrue(gate, "q must sum to one without being normalised")

    def test_published_f_g(self):
        rows = read_pedigree(corpus("new_lacy.ped"))
        ped, _back = renumber(rows)
        f_g, _r = f_g_lacy(ped)
        self.assertAlmostEqual(2.18, f_g, delta=TWO_DP)

    def test_published_retention(self):
        """Appendix A: "the fraction of alleles retained ... (r_i) is 0.75"."""
        rows = read_pedigree(corpus("new_lacy.ped"))
        ped, _back = renumber(rows)
        for founder, r_i in retention(ped).items():
            with self.subTest(founder=founder):
                self.assertAlmostEqual(0.75, r_i, places=9)

    def test_f_g_never_exceeds_f_e(self):
        """
        p.116: "The number of founder genome equivalents is always less than
        the number of founder equivalents", because it discounts alleles lost
        to drift as well as unequal contributions.
        """
        rows = read_pedigree(corpus("new_lacy.ped"))
        ped, _back = renumber(rows)
        f_e, _q, _c, _d = f_e_lacy(ped)
        f_g, _r = f_g_lacy(ped)
        self.assertLessEqual(f_g, f_e + 1e-9)


class TestTable1(unittest.TestCase):
    """
    Lacy Table 1, article p.116 -- "Effects of founder contributions on founder
    and founder genome equivalents". Two founder pairs, f = 4 throughout.

    H is published too (0.812 ... 0.851) and is deliberately not tested: the
    footnote says it "was determined from 10,000 computer simulations", so it
    is a Monte-Carlo figure, not an exact target.
    """

    # (offspring of pair 1, offspring of pair 2, f_e, f_g, r_1, r_2)
    PUBLISHED = (
        (2, 2, 4.00, 3.00, 0.750, 0.750),
        (3, 3, 4.00, 3.50, 0.875, 0.875),
        (4, 4, 4.00, 3.75, 0.938, 0.938),
        (8, 8, 4.00, 3.98, 0.996, 0.996),
        (2, 4, 3.60, 3.21, 0.750, 0.938),
        (4, 8, 3.60, 3.54, 0.938, 0.996),
        (8, 16, 3.60, 3.60, 0.996, 1.00),
    )

    def test_published_f_e(self):
        for x1, x2, want, _fg, _r1, _r2 in self.PUBLISHED:
            with self.subTest(offspring=(x1, x2)):
                ped = _two_founder_pairs(x1, x2)
                got, _q, _c, _d = f_e_lacy(ped)
                self.assertAlmostEqual(want, got, delta=TWO_DP)

    def test_published_f_g(self):
        for x1, x2, _fe, want, _r1, _r2 in self.PUBLISHED:
            with self.subTest(offspring=(x1, x2)):
                ped = _two_founder_pairs(x1, x2)
                got, _r = f_g_lacy(ped)
                self.assertAlmostEqual(want, got, delta=TWO_DP)

    def test_published_retention(self):
        """r_i = 1 - .5^x, the Table 1 footnote's closed form."""
        for x1, x2, _fe, _fg, want1, want2 in self.PUBLISHED:
            with self.subTest(offspring=(x1, x2)):
                ped = _two_founder_pairs(x1, x2)
                r = retention(ped)
                self.assertAlmostEqual(want1, r[1], delta=THREE_DP)
                self.assertAlmostEqual(want1, r[2], delta=THREE_DP)
                self.assertAlmostEqual(want2, r[3], delta=THREE_DP)
                self.assertAlmostEqual(want2, r[4], delta=THREE_DP)

    def test_equal_contributions_give_f_e_equal_to_the_founder_count(self):
        """
        p.114: "If all founders contribute equally to the descendant
        population, the founder equivalent number is equal to the actual number
        of founders." The first four rows are exactly that case.
        """
        for x1, x2, want, _fg, _r1, _r2 in self.PUBLISHED[:4]:
            with self.subTest(offspring=(x1, x2)):
                got, _q, contributors, _d = f_e_lacy(_two_founder_pairs(x1, x2))
                self.assertAlmostEqual(4.0, got, places=9)
                self.assertEqual(4, len(contributors))

    def test_f_e_never_exceeds_the_founder_count(self):
        """
        The bound withdrawn in d1129c4 and reinstated by p.114. It holds here
        with k counted as the paper counts it.
        """
        for x1, x2, _fe, _fg, _r1, _r2 in self.PUBLISHED:
            with self.subTest(offspring=(x1, x2)):
                got, _q, contributors, _d = f_e_lacy(_two_founder_pairs(x1, x2))
                self.assertLessEqual(got, len(contributors) + 1e-9)


class TestHalfFounderModesDisagreeAsDocumented(unittest.TestCase):
    """
    On a pedigree WITH a half-founder the three modes must differ, and only the
    paper's must produce a probability vector. If they ever agreed, the
    comparators would have stopped comparing anything.
    """

    PED = "1 0 0\n2 0 0\n3 1 2\n4 3 0\n5 3 4\n6 3 4\n"

    def _pedfile(self):
        tmp = tempfile.mkdtemp(prefix="pypedal_lacy_")
        path = os.path.join(tmp, "half.ped")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.PED)
        return path

    def test_only_phantom_yields_a_probability_vector(self):
        path = self._pedfile()
        self.assertEqual(1, lacy_n_half_founders(path))
        gates = {mode: lacy_f_e(path, mode=mode)[1] for mode in LACY_MODES}
        self.assertTrue(gates["phantom"])
        self.assertFalse(gates["lacy"])
        self.assertFalse(gates["half"])

    def test_the_modes_give_different_answers(self):
        path = self._pedfile()
        values = {round(lacy_f_e(path, mode=mode)[0], 9) for mode in LACY_MODES}
        self.assertEqual(3, len(values),
                         "the comparators must still differ from the paper's rule")


if __name__ == "__main__":
    unittest.main()
