"""filetag is the pedigree path with only the final extension stripped."""
import unittest

from PyPedal.pyp_newclasses import filetag_from_pedfile


class TestFiletagDerivation(unittest.TestCase):
    def test_plain_name_keeps_the_stem(self):
        self.assertEqual("mrode", filetag_from_pedfile("mrode.ped"))

    def test_dot_slash_keeps_the_relative_prefix(self):
        self.assertEqual("./mrode", filetag_from_pedfile("./mrode.ped"))

    def test_parent_directory_keeps_the_relative_prefix(self):
        self.assertEqual("../dir/mrode", filetag_from_pedfile("../dir/mrode.ped"))

    def test_dotted_directory_is_not_truncated(self):
        self.assertEqual(
            "/tmp/data.v1/mrode", filetag_from_pedfile("/tmp/data.v1/mrode.ped")
        )

    def test_two_relative_pedigrees_do_not_share_untitled_filetag(self):
        mrode = filetag_from_pedfile("./mrode.ped")
        lacy = filetag_from_pedfile("./lacy.ped")
        self.assertEqual("./mrode", mrode)
        self.assertEqual("./lacy", lacy)
        self.assertNotEqual(mrode, lacy)
        self.assertNotEqual("untitled_pedigree", mrode)
        self.assertNotEqual("untitled_pedigree", lacy)

    def test_two_relative_pedigrees_cannot_collide_on_output_names(self):
        mrode_out = filetag_from_pedfile("./mrode.ped") + "_inbreeding.dat"
        lacy_out = filetag_from_pedfile("./lacy.ped") + "_inbreeding.dat"
        self.assertNotEqual(mrode_out, lacy_out)
        self.assertNotEqual(
            filetag_from_pedfile("./mrode.ped"),
            filetag_from_pedfile("./lacy.ped"),
        )
