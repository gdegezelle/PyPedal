#!/usr/bin/env python3

#
# PyPedal - a software package for pedigree analysis
#
# Copyright (C) 2001-2024  John B. Cole (john.b.cole@gmail.com)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from PyPedal import pyp_newclasses
from PyPedal import pyp_metrics
from PyPedal import pyp_nrm
from PyPedal import pyp_utils
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == "__main__":

    print(f"Starting new_methods.py at {pyp_nice_time()}")

    example = pyp_newclasses.load_pedigree(options_file="new_format.ini")
    pyp_utils.assign_offspring(example)

    inbr = pyp_nrm.inbreeding(example, method="vanraden", rels=1)
    print(f"inbr:\n{inbr}")
    print(f"fx:\n{inbr.get('fx')}")
    print(f"rel_dict:\n{inbr.get('rel_dict')}")

    # Default NRM is SciPy sparse. NewAMatrix.save() uses ndarray.tofile(),
    # which that matrix type does not provide. Print the in-memory matrix.
    example.nrm = pyp_newclasses.NewAMatrix(example.kw)
    example.nrm.form_a_matrix(example.pedigree)
    example.nrm.printme()
    nrm = example.nrm.nrm
    if hasattr(nrm, "toarray"):
        nrm = nrm.toarray()
    print(nrm[1, 4])
    print(nrm[4, 1])

    print(f"Calling related_animals() at {pyp_nice_time()}")
    list_a = pyp_metrics.related_animals(example.pedigree[4].animalID, example)
    print(list_a)

    print(f"Calling related_animals() at {pyp_nice_time()}")
    list_b = pyp_metrics.related_animals(example.pedigree[13].animalID, example)
    print(list_b)

    print(f"Calling common_ancestors() at {pyp_nice_time()}")
    list_r = pyp_metrics.common_ancestors(
        example.pedigree[4].animalID, example.pedigree[13].animalID, example
    )
    print(list_r)

    print('\tTesting NewPedigree::addanimal() at %s' % ( pyp_nice_time() ))
    if not example.kw.get('pedigree_is_renumbered'):
        example.renumber()
    _added = example.addanimal(15,10,11)
    if _added:
       print('\t\tAdded animal 15 to %s' % ( example.kw['pedname'] ))

    print('\tTesting NewPedigree::delete_animals() at %s' % ( pyp_nice_time() ))
    if _added:
        _deleted = example.delete_animals([15])
        if _deleted:
           print('\t\tDeleted animal 15 from %s' % ( example.kw['pedname'] ))

    print('\tTesting pyp_metrics.mating_coi() at %s' % ( pyp_nice_time() ))
    _a = int(example.pedigree[0].animalID)
    _b = int(example.pedigree[1].animalID)
    print('\t\tmating_coi(%s, %s) = %s' % (
        _a, _b, pyp_metrics.mating_coi(_a, _b, example)
    ))

    desc1 = pyp_metrics.descendants(5, example, {})
    print(f"desc1: {desc1}")

    descf = pyp_metrics.founder_descendants(example)
    print(f"descf: {descf}")

    print(f"Stopping new_methods.py at {pyp_nice_time()}")
