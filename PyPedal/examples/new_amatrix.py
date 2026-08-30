#!/usr/bin/env python3

# Not an active successful example. new_amatrix.ini / new_amatrix.ped is
# asdgb historical evidence: dam 2047 is year 1997 while named offspring
# dated load. See docs/manual/recipes.md.

from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

import os
import numpy
import time

if __name__ == "__main__":

    print("Starting test new_amatrix.py at %s" % (pyp_nice_time()))

    options = os.path.join(os.path.dirname(__file__), 'new_amatrix.ini')

    example = pyp_newclasses.load_pedigree(options_file=options, debug_load=True)

    # amatrix = pyp_newclasses.NewAMatrix(example.kw)
    # amatrix.form_a_matrix(example.pedigree)

    # Here's how to save a matrix to a binary file.
    #    amatrix.save('boichard2_pedigree.bin')

    # Here's how to load a matrix from a binary file.
    #    amatrix2 = pyp_newclasses.NewAMatrix(example.kw)
    #    amatrix2.load('boichard2_pedigree.bin')

    # Calculate coefficients of inbreeding on this pedigree.
    #print(f"\tEntering pyp_nrm.inbreeding() at {pyp_nice_time()}")
    #example_inbreeding = pyp_nrm.inbreeding(example, method="vanraden")
    #print(f"\tReturning from pyp_nrm.inbreeding() at {pyp_nice_time()}")
    #print(example_inbreeding["metadata"])
    
    
    #print("Stopping test new_amatrix.py at %s" % (pyp_nice_time()))
