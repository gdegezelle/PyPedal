"""filetag is the pedigree path with only the final extension stripped."""
import os
import tempfile
import unittest

from PyPedal.pyp_newclasses import NewPedigree


def _filetag(pedfile):
    ped = NewPedigree(
        {
            "pedfile": pedfile,
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )
    return ped.kw["filetag"]


class TestFiletagDerivation(unittest.TestCase):
    def test_plain_name_keeps_the_stem(self):
        self.assertEqual("mrode", _filetag("mrode.ped"))

    def test_dot_slash_keeps_the_relative_prefix(self):
        self.assertEqual("./mrode", _filetag("./mrode.ped"))

    def test_parent_directory_keeps_the_relative_prefix(self):
        self.assertEqual("../dir/mrode", _filetag("../dir/mrode.ped"))

    def test_dotted_directory_is_not_truncated(self):
        self.assertEqual("/tmp/data.v1/mrode", _filetag("/tmp/data.v1/mrode.ped"))

    def test_two_relative_pedigrees_do_not_share_untitled_filetag(self):
        mrode = _filetag("./mrode.ped")
        lacy = _filetag("./lacy.ped")
        self.assertEqual("./mrode", mrode)
        self.assertEqual("./lacy", lacy)
        self.assertNotEqual(mrode, lacy)
        self.assertNotEqual("untitled_pedigree", mrode)
        self.assertNotEqual("untitled_pedigree", lacy)

    def test_two_relative_pedigrees_cannot_collide_on_output_names(self):
        tmp = tempfile.mkdtemp(prefix="pypedal_filetag_")
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            open("mrode.ped", "w", encoding="utf-8").write("1 0 0\n")
            open("lacy.ped", "w", encoding="utf-8").write("1 0 0\n")
            mrode = _filetag("./mrode.ped")
            lacy = _filetag("./lacy.ped")
            mrode_out = mrode + "_inbreeding.dat"
            lacy_out = lacy + "_inbreeding.dat"
            self.assertNotEqual(mrode_out, lacy_out)
            open(mrode_out, "w", encoding="utf-8").write("mrode\n")
            open(lacy_out, "w", encoding="utf-8").write("lacy\n")
            with open(mrode_out, encoding="utf-8") as handle:
                self.assertEqual("mrode\n", handle.read())
            with open(lacy_out, encoding="utf-8") as handle:
                self.assertEqual("lacy\n", handle.read())
        finally:
            os.chdir(cwd)
