"""
``NewAnimal.ancestor`` must exist on every animal ().

PyPedal 2.0.4 initialises the attribute unconditionally in ``NewAnimal.__init__``
(``pyp_newclasses.py:3206``)::

    self.paddedID = self.pad_id()
    self.ancestor = 0
    self.sons = {}

The Python 3 port dropped that line. Nothing else initialises it:
``pyp_utils.set_ancestor_flag`` only ever writes ``1``, and only onto animals that
*are* a sire or dam of someone else (``pyp_utils.py:91`` and ``:98``). Every animal
that is not a parent therefore had no ``ancestor`` attribute at all, which is why
``pyp_db.populate_table`` -- which reads ``record.ancestor`` at ``pyp_db.py:214``
-- could not run.

This is the same defect class as  (``userField``): an assignment present
in 2.0.4 and absent from the port, with live consumers left in place. Those two are
the *only* two omissions in ``NewAnimal.__init__``; the guard that asserts so lives
in ``test_legacy_numeric_compatibility.py``.
"""
import unittest

from PyPedal import pyp_utils

from _pedhelpers import load_corpus


class TestAncestorIsInitialised(unittest.TestCase):

    def test_every_animal_has_ancestor_zero_on_a_fresh_load(self):
        """
        The whole of : the attribute exists, and it starts at 0.

        ``set_ancestor_flag`` has not been called, so no animal should be
        flagged -- and, crucially, none should be *missing* the attribute.
        """
        ped = load_corpus("new_lacy.ped")
        for animal in ped.pedigree:
            self.assertTrue(
                hasattr(animal, "ancestor"),
                f"animal {animal.originalID} has no 'ancestor' attribute; "
                "PyPedal 2.0.4 sets it at pyp_newclasses.py:3206")
            self.assertEqual(
                0, animal.ancestor,
                f"animal {animal.originalID} starts with ancestor="
                f"{animal.ancestor!r}, expected 0")

    def test_set_ancestor_flag_marks_parents_and_leaves_the_rest_at_zero(self):
        """
        The half that only holds once the initialiser exists.

        Before the repair the non-parents had no attribute to read, so the
        "and the rest are still 0" assertion could not even be expressed.
        """
        ped = load_corpus("new_lacy.ped")

        missing = ped.kw["missing_parent"]
        expected_parents = set()
        for animal in ped.pedigree:
            if animal.sireID != missing:
                expected_parents.add(int(animal.sireID))
            if animal.damID != missing:
                expected_parents.add(int(animal.damID))
        # A vacuous parent set would make the rest of this test meaningless.
        self.assertTrue(expected_parents)

        pyp_utils.set_ancestor_flag(ped)

        flagged = {a.animalID for a in ped.pedigree if a.ancestor == 1}
        self.assertEqual(expected_parents, flagged)

        for animal in ped.pedigree:
            if animal.animalID not in expected_parents:
                self.assertEqual(
                    0, animal.ancestor,
                    f"animal {animal.animalID} is not a parent but was left "
                    f"with ancestor={animal.ancestor!r}")


if __name__ == "__main__":
    unittest.main()
