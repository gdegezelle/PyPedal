import os
import tempfile
import unittest

import pytest

from PyPedal.pyp_utils import fast_reorder
from PyPedal.pyp_newclasses import NewPedigree
from PyPedal.pyp_newclasses import NewAnimal

from _pedhelpers import (
    GRIFFON_TEST_SMALL_IDS,
    chdir_tmp,
    load_canonical_griffon,
    write_canonical_griffon_subset,
)


def validate_order(pedigree):
    """
    Validate that all parents appear before their offspring in the reordered pedigree.
    :param pedigree: List of NewAnimal objects representing the reordered pedigree.
    :return: True if the order is valid, raises ValueError otherwise.
    """
    id_to_index = {animal.animalID: idx for idx, animal in enumerate(pedigree)}

    for idx, animal in enumerate(pedigree):
        if animal.sireID in id_to_index and id_to_index[animal.sireID] >= idx:
            raise ValueError(f"Sire {animal.sireID} of animal {animal.animalID} appears after its offspring.")
        if animal.damID in id_to_index and id_to_index[animal.damID] >= idx:
            raise ValueError(f"Dam {animal.damID} of animal {animal.animalID} appears after its offspring.")
    return True


def calculate_statistics(pedigree):
    """
    Calculate basic statistics of the pedigree for integrity checks.
    :param pedigree: List of NewAnimal objects.
    :return: A dictionary of statistics.
    """
    stats = {
        "total_animals": len(pedigree),
        "unique_sires": len(set(a.sireID for a in pedigree if a.sireID > 0)),
        "unique_dams": len(set(a.damID for a in pedigree if a.damID > 0)),
        "founders": sum(1 for a in pedigree if a.sireID == -999 and a.damID == -999),
    }
    return stats


def write_pedigree_to_file(pedigree, filename):
    """
    Write the pedigree to a file for manual inspection.
    :param pedigree: List of NewAnimal objects representing the pedigree.
    :param filename: Path to the output file.
    """
    with open(filename, 'w', encoding='utf-8') as file:
        for animal in pedigree:
            file.write(f"{animal.animalID},{animal.sireID},{animal.damID}\n")


class TestFastReorder(unittest.TestCase):
    def test_small_pedigree(self):
        """
        Test fast_reorder with a small, simple pedigree.

        an earlier revision stage D2: this was skipped as "disabled temporarily". Run
        unskipped it passes its assertions -- including the birth-year check,
        which was suspected of being ill-founded because ``new_graphics.ped`` is
        loaded as ``asdgy``. That suspicion was wrong; the assertion holds.

        What it did do was write ``_new_reordered__reordered.ped`` and
        ``_new_reordered__id_map.map`` into the working directory and then
        assert on those relative paths, so it both depended on where pytest was
        invoked from and littered the repository. Wrapped in ``chdir_tmp()``
        like its ``_asdxb`` sibling, which is the same fix D1 applied there.
        """
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../PyPedal/examples"))
        pedfile_path = os.path.join(base_dir, "new_graphics.ped")
        options = {
            "messages": "verbose",
            "renumber": 0,
            "pedfile": pedfile_path,
            "pedformat": "asdgy",
            "sepchar": " ",
        }
        example = NewPedigree(options)
        example.load()

        with chdir_tmp():
            # Reorder the pedigree
            reordered_pedigree = fast_reorder(
                example.pedigree,
                debug=True,
                io="yes",
            )

            # Validate the order
            try:
                validate_order(reordered_pedigree)
            except ValueError as e:
                self.fail(f"Reordering validation failed: {e}")

            # Verify integrity
            self.assertEqual(len(example.pedigree), len(reordered_pedigree))

            # Expected filenames to match fast_reorder() output
            expected_output_file = f"_new_reordered__reordered.ped"
            expected_id_map_file = f"_new_reordered__id_map.map"

            # Check if output files exist
            self.assertTrue(os.path.exists(expected_output_file), f"Expected output file {expected_output_file} not found!")
            self.assertTrue(os.path.exists(expected_id_map_file), f"Expected ID map file {expected_id_map_file} not found!")

        # Birth years must survive reordering; compare by ID, not list position
        original_by = {animal.animalID: animal.by for animal in example.pedigree}
        for reordered in reordered_pedigree:
            if reordered.animalID in original_by:
                self.assertEqual(
                    original_by[reordered.animalID],
                    reordered.by,
                    f"Birthyear mismatch for animal {reordered.animalID}",
                )

    # @unittest.skip("Test is disabled temporarily.")
    def test_small_pedigree_asdxb(self):
        """
        Test fast_reorder with a small, simple pedigree.
        """
        tmp = tempfile.mkdtemp(prefix="pypedal_test_")
        pedfile_path = os.path.join(tmp, "griffon_test_small.ped")
        write_canonical_griffon_subset(pedfile_path, GRIFFON_TEST_SMALL_IDS)
        options = {
            "messages": "quiet",
            "renumber": 0,
            "pedfile": pedfile_path,
            "pedformat": "asdxb",
            "sepchar": ",",
        }
        example = NewPedigree(options)
        example.load()

        # fast_reorder(io='yes') writes _new_reordered__reordered.ped and
        # _new_reordered__id_map.map relative to the *working directory*, and
        # the assertions below are on those relative paths. Run the block inside
        # a temporary directory so the test neither depends on where pytest was
        # invoked from nor litters the repository.
        with chdir_tmp():
            # Reorder the pedigree
            reordered_pedigree = fast_reorder(
                example.pedigree,
                debug=True,
                io="yes",
            )

            # Validate the order
            try:
                validate_order(reordered_pedigree)
                print("Reordering is valid for small pedigree!")
            except ValueError as e:
                self.fail(f"Reordering validation failed: {e}")

            # Verify integrity
            self.assertEqual(len(example.pedigree), len(reordered_pedigree))

            # Expected filenames to match fast_reorder() output
            expected_output_file = f"_new_reordered__reordered.ped"
            expected_id_map_file = f"_new_reordered__id_map.map"

            # Check if output files exist
            self.assertTrue(os.path.exists(expected_output_file), f"Expected output file {expected_output_file} not found!")
            self.assertTrue(os.path.exists(expected_id_map_file), f"Expected ID map file {expected_id_map_file} not found!")

        # Birth years must survive reordering; compare by ID, not list position
        original_by = {animal.animalID: animal.by for animal in example.pedigree}
        for reordered in reordered_pedigree:
            if reordered.animalID in original_by:
                self.assertEqual(
                    original_by[reordered.animalID],
                    reordered.by,
                    f"Birthyear mismatch for animal {reordered.animalID}",
                )

    @pytest.mark.integration
    def test_large_pedigree_reorder(self):
        """
        Test reordering with a large pedigree file.

        THE ONLY SCALE TEST IN THE SUITE. Everything else runs on pedigrees of
        6 to 45 animals, so nothing else would catch a reorder regression that
        only appears at ~100k records.

        Skipped from an earlier revision stage D2 until , with an accurate reason:
        it did not complete the *load* of this 98,016-record file in twelve
        minutes. That was the quadratic ordering, and its own skip message asked
        for the test back behind the 'integration' marker once the cost was
        understood. Measured after the repair, the load takes about 1.5 seconds.

        Loaded through ``load_example`` and reordered inside ``chdir_tmp`` so
        that ``io='yes'`` writes land in a temporary directory: this test wrote
        ``_new_reordered__reordered.ped`` and ``_new_reordered__id_map.map``
        into whatever directory pytest was invoked from, which the
        repository-delta guard now catches.

        Identity columns and true chronology (``asdxb``). The corrected
        canonical sample loads with recorded dates; this test is still
        reorder-at-scale, not a dated-load oracle.
        """
        example = load_canonical_griffon({
            "messages": "quiet",
            "renumber": True,
            "pedformat": "asdxb",
            "sepchar": ",",
            "pedigree_summary": 0,
        })
        self.assertEqual(len(example.pedigree), 98001,
                           "the scale test must actually be at scale")

        # The graph, in stable identity, before anything is reordered again.
        before = sorted(
            "%s|%s|%s" % (a.originalID, role, p)
            for a in example.pedigree
            for role, p in (("s", a.sireID), ("d", a.damID)))

        with chdir_tmp():
            reordered_pedigree = fast_reorder(example.pedigree, io="yes")

        try:
            validate_order(reordered_pedigree)
        except ValueError as e:
            self.fail(f"Reordering validation failed: {e}")

        # Not a single animal gained, lost or re-parented.
        self.assertEqual(before, sorted(
            "%s|%s|%s" % (a.originalID, role, p)
            for a in reordered_pedigree
            for role, p in (("s", a.sireID), ("d", a.damID))))

        original_stats = calculate_statistics(example.pedigree)
        reordered_stats = calculate_statistics(reordered_pedigree)
        self.assertEqual(original_stats, reordered_stats, "Statistics mismatch after reordering!")

        # Reordering an already-ordered pedigree of this size is a fixed point.
        again = fast_reorder(reordered_pedigree)
        self.assertEqual([a.animalID for a in reordered_pedigree],
                         [a.animalID for a in again])


    def test_empty_pedigree(self):
        """
        Test fast_reorder with an empty pedigree.

        an earlier revision stage D2: this was skipped as "disabled temporarily" with no
        further reason. Run unskipped, it passes in under half a second. There
        was nothing wrong with it -- it looks like collateral damage from a
        blanket disable of the whole class, and it is the degenerate-input
        contract for fast_reorder, which is worth having.
        """
        empty_pedigree = []
        reordered = fast_reorder(empty_pedigree, debug=True)
        self.assertEqual(reordered, [])

    def test_circular_reference(self):
        """
        Test fast_reorder with a circular reference in parentage.
        """
        from PyPedal.pyp_newclasses import NewAnimal

        # Define mock data for circular reference
        locations = {"animal": 0, "sire": 1, "dam": 2, "birthyear": 3}
        data = [
            ["1", "2", "3", "2020"],  # Animal 1 with sire 2 and dam 3
            ["2", "3", "1", "2019"],  # Animal 2 with sire 3 and dam 1 (circular reference)
            ["3", "1", "2", "2018"],  # Animal 3 with sire 1 and dam 2 (circular reference)
        ]
        mykw = {"missing_parent": -999, "missing_byear": 0}

        # Create mock pedigree
        pedigree = [NewAnimal(locations, row, mykw) for row in data]

        # This should raise a ValueError due to circular reference
        with self.assertRaises(ValueError):
            fast_reorder(pedigree, debug=True)

    def test_pad_id(self):
        """
        Test the pad_id method of NewAnimal to ensure correct padded ID generation.
        """
        locations = {
            'animal': 0,
            'sire': 1,
            'dam': 2,
            'birthyear': 3,
        }
        data = [123, 456, 789, 2005]  # Example data: [animalID, sireID, damID, birthyear]
        mykw = {
            'missing_parent': -1,
            'missing_name': 'unknown',
            'missing_byear': 1900,
        }
        
        # Initialize the NewAnimal object
        animal = NewAnimal(locations, data, mykw)
        
        # Generate the padded ID
        padded_id = animal.pad_id()
        
        # Expected padded ID
        expected_padded_id = "200500000000123"  # Birth year + padded animal ID
        
        self.assertEqual(padded_id, expected_padded_id)



if __name__ == "__main__":
    unittest.main()
