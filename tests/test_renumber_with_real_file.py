import os
import shutil
import unittest
from PyPedal.pyp_newclasses import NewPedigree
from PyPedal.pyp_utils import renumber
from _pedhelpers import owned_temp_dir


class TestRenumberWithRealPedFile(unittest.TestCase):
    def setUp(self):
        """Set up paths and load the test pedigree file."""
        self.base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../PyPedal/examples"))
        self.pedfile_path = os.path.join(self.base_dir, "new_graphics.ped")

        # renumber(io='yes') with no filetag writes _renumbered__renum.ped and
        # _renumbered__id_map.map relative to the *working directory*. Those are
        # matched by the `_renumbered_*` line in .gitignore, so for a long time
        # they were rewritten into the repository root on every test run and
        # `git status` never said a word about it. Run from a temporary
        # directory, and write the explicit-output-dir case there too, so this
        # suite leaves no trace in the source tree.
        self._tmp = owned_temp_dir(prefix="pypedal_test_renumber_")
        self._cwd = os.getcwd()
        os.chdir(self._tmp)
        self.output_dir = self._tmp
        local_ped = os.path.join(self._tmp, "new_graphics.ped")
        shutil.copy(self.pedfile_path, local_ped)

        self.options = {
            "pedfile": local_ped,
            "pedformat": "asdxy",
            "sepchar": " ",
            "messages": "verbose",
            "renumber": 0,
            "pedigree_is_renumbered": False,
            "form_nrm": True,
        }

        # Load the test pedigree
        self.mock_pedigree = NewPedigree(self.options)
        self.mock_pedigree.load()

    def tearDown(self):
        os.chdir(self._cwd)
        shutil.rmtree(self._tmp, ignore_errors=True)

    # @unittest.skip("Test is disabled temporarily.")
    def test_renumber_with_real_file(self):
        """Test the renumber function using a real .ped file."""
        result = renumber(self.mock_pedigree, io="yes", debug=True, cleanmap=False)
        self.assertIsInstance(result, NewPedigree)
        self.assertTrue(self.mock_pedigree.kw.get("pedigree_is_renumbered", False))

    # @unittest.skip("Test is disabled temporarily.")
    def test_renumber_writes_output(self):
        """Test that renumbered files are written to the correct location."""
        file_prefix = os.path.splitext(os.path.basename(self.options["pedfile"]))[0]
        filetag = f"{file_prefix}_renumbered_"
        output_dir = self.output_dir  # a temporary directory; see setUp

        renum_file = os.path.join(output_dir, f"{filetag}_renum.ped")
        id_map_file = os.path.join(output_dir, f"{filetag}_id_map.map")


        print(f"[TEST DEBUG] Expected output directory: {output_dir}")
        print(f"[TEST DEBUG] Expected renumbered file: {renum_file}")
        print(f"[TEST DEBUG] Expected ID map file: {id_map_file}")

        # Run renumber with the correct filetag and output_dir
        result = renumber(self.mock_pedigree, io="yes", debug=True, cleanmap=False, output_dir=output_dir, filetag=filetag)

        self.assertIsInstance(result, NewPedigree)

        # Verify files were created in the correct directory
        self.assertTrue(os.path.exists(renum_file), f"Expected file {renum_file} not found.")
        self.assertTrue(os.path.exists(id_map_file), f"Expected file {id_map_file} not found.")

        # Verify renumbered file contents
        with open(renum_file, "r", encoding="utf-8") as file:
            content = file.readlines()
            self.assertGreater(len(content), 2, "Renumbered file content seems too small.")

        # Cleanup
        if os.path.exists(renum_file):
            os.remove(renum_file)
        if os.path.exists(id_map_file):
            os.remove(id_map_file)



if __name__ == "__main__":
    unittest.main()
