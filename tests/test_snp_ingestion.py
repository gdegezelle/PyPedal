"""
SNP file ingestion.

These tests cover reading, validating and attaching genotypes.
They do not compute genomic inbreeding coefficients. VanRaden (2008)
Method 1 fixes the centring convention, the scaling denominator and the
reference allele frequencies; callers must know those conventions.
"""
import os
import tempfile
import unittest

from PyPedal import pyp_errors, pyp_newclasses, pyp_snp


def write(lines, name):
    tmp = tempfile.mkdtemp(prefix="pypedal_snp_")
    path = os.path.join(tmp, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


# Four animals, ten loci. Hand-checkable and small enough to reason about.
PEDIGREE = ["1 0 0", "2 0 0", "3 1 2", "4 1 2"]
GENOTYPES = [
    "1 chip1 10 0120120120",
    "2 chip1 10 1201201201",
    "3 chip1 10 2012012012",
    "4 chip1 10 0000000000",
]


def load(pedigree_rows=PEDIGREE, genotype_rows=GENOTYPES, **overrides):
    pedfile = write(pedigree_rows, "ped.ped")
    options = {
        "pedfile": pedfile,
        "pedformat": "asd",
        "sepchar": " ",
        "messages": "quiet",
        "pedigree_summary": 0,
        "renumber": True,
    }
    if genotype_rows is not None:
        options["snpfile"] = write(genotype_rows, "geno.txt")
    options.update(overrides)
    cwd = os.getcwd()
    try:
        os.chdir(os.path.dirname(pedfile))
        return pyp_newclasses.loadPedigree(options=options)
    finally:
        os.chdir(cwd)


class TestSnpFileIsRead(unittest.TestCase):

    def test_genotypes_are_attached_to_the_pedigree(self):
        ped = load()
        self.assertIsNot(ped.snp, False,
                         "kw['snpfile'] was set but pedobj.snp was never populated")
        self.assertEqual(4, len(ped.snp))
        self.assertEqual(10, len(ped.snp["genotype"].iloc[0]))

    def test_without_snpfile_the_pedigree_still_loads(self):
        """Ingestion is opt-in; nothing about a plain pedigree load changes."""
        ped = load(genotype_rows=None)
        self.assertIs(ped.snp, False)
        self.assertEqual(4, len(ped.pedigree))

    def test_grm_construction_becomes_reachable(self):
        """
        The point of the stage. Not an assertion about the VALUES -- only that
        the guard no longer short-circuits, which it always did before.
        """
        ped = load()
        grm = pyp_snp.form_grm_from_snp(ped)
        self.assertIsNot(grm, False,
                         "form_grm_from_snp still returns False, so SNP data is "
                         "not reaching it")
        self.assertEqual((4, 4), grm.shape)

    def test_comments_are_ignored(self):
        rows = ["# a comment line"] + GENOTYPES
        ped = load(genotype_rows=rows)
        self.assertEqual(4, len(ped.snp))


class TestGenotypeValidation(unittest.TestCase):
    """
    Each check is about the data, not about any estimator, so it holds whatever
    is later computed from the genotypes.
    """

    def test_inconsistent_genotype_lengths_raise(self):
        rows = list(GENOTYPES)
        rows[2] = "3 chip1 8 20120120"
        with self.assertRaises(pyp_errors.PyPedalInputError) as caught:
            load(genotype_rows=rows)
        self.assertIn("inconsistent lengths", str(caught.exception))

    def test_non_dosage_characters_raise(self):
        rows = list(GENOTYPES)
        rows[1] = "2 chip1 10 12012012A1"
        with self.assertRaises(pyp_errors.PyPedalInputError) as caught:
            load(genotype_rows=rows)
        self.assertIn("other than 0, 1 and 2", str(caught.exception))

    def test_duplicate_animals_raise(self):
        rows = list(GENOTYPES) + ["1 chip1 10 2222222222"]
        with self.assertRaises(pyp_errors.PyPedalInputError) as caught:
            load(genotype_rows=rows)
        self.assertIn("more than one genotype", str(caught.exception))

    def test_an_empty_snp_file_raises(self):
        path = write(["# nothing but a comment"], "geno.txt")
        ped = load(genotype_rows=None)
        with self.assertRaises(pyp_errors.PyPedalInputError):
            pyp_snp.load_snp_file(ped, snpfile=path)

    def test_a_missing_snp_file_raises_input_error(self):
        ped = load(genotype_rows=None)
        with self.assertRaises(pyp_errors.PyPedalInputError):
            pyp_snp.load_snp_file(ped, snpfile="/nonexistent/path/geno.txt")

    def test_no_snpfile_configured_raises_configuration_error(self):
        ped = load(genotype_rows=None)
        with self.assertRaises(pyp_errors.PyPedalConfigurationError):
            pyp_snp.load_snp_file(ped)

    def test_a_partially_genotyped_pedigree_is_allowed(self):
        """
        Genotyping a subset is normal practice. It warns rather than failing --
        what is unacceptable is doing it silently.
        """
        ped = load(genotype_rows=GENOTYPES[:2])
        self.assertEqual(2, len(ped.snp))
        self.assertEqual(4, len(ped.pedigree))


class TestPandasCompatibility(unittest.TestCase):

    def test_readers_do_not_use_the_removed_delim_whitespace_keyword(self):
        """
        ``delim_whitespace`` is deprecated in pandas 2.2 and removed in 3.0.
        All three AGIL readers used it, so they were a live breakage on the
        dependency floor the project is moving to, independent of any genomic
        decision.
        """
        source = os.path.join(os.path.dirname(pyp_snp.__file__), "pyp_snp.py")
        with open(source, encoding="utf-8") as handle:
            self.assertNotIn("delim_whitespace", handle.read())

    def test_reading_emits_no_future_warnings(self):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            ped = load()
        self.assertEqual(4, len(ped.snp))


if __name__ == "__main__":
    unittest.main()
