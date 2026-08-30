#!/usr/bin/env python3

#######################################################################
# NAME: example_a_matrix.py
# VERSION: 2.0.0 (29 August 2013)
# AUTHOR: John B. Cole (john.cole@ars.usda.gov)
# LICENSE: LGPL
#######################################################################

# Import (load) the pieces of PyPedal that we're going to need
# for this example.
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_utils

import os

ini_file = os.path.join(
    os.path.dirname(__file__), 'new_graphics.ini')

# Create a dictionary to hold options, such as the name of the
# input file, and its format. new_amatrix.ped/asdgb is historical
# evidence with an impossible dam year; this example loads the valid
# asdgb graphics pedigree instead.
options = {}
options["pedfile"] = "new_graphics.ped"
options["pedname"] = "new_graphics"
options["pedformat"] = "asdgb"
options["messages"] = "verbose"
options["renumber"] = True

if __name__ == "__main__":

    print(f"Starting pypedal.py at {pyp_utils.pyp_nice_time()}")

    # example is a PyPedal NewPedigree object that stores the
    # individual animal records, as well as other information about the
    # pedigree.
    # Most things are done by calling methods on example, or by passing
    # it to functions.
    example = pyp_newclasses.load_pedigree(options_file=ini_file)

    # We're going to create an instance of the NewAMatrix class, which we
    # can use to store relationship matrices, and then create the A matrix
    # from the pedigree stored in example. Note that we can only store
    # matrices that we can allocate in memory!
    amatrix = pyp_newclasses.NewAMatrix(example.kw)
    amatrix.form_a_matrix(example.pedigree)

    # Here's how to save the A matrix to a binary file.
    # amatrix.save(os.path.join(os.path.dirname(
    #     __file__), "boichard2_pedigree.bin"))

    # Here's how to load the A matrix from a binary file.
    # amatrix2 = pyp_newclasses.NewAMatrix(example.kw)
    # amatrix2.load(os.path.join(os.path.dirname(
    #     __file__), "boichard2_pedigree.bin"))

    # Now we're going to calculate coefficients of inbreeding for the individuals
    # in this pedigree. There are a few different algorithms for doing this in
    # PyPedal.
    # VanRaden's method is slow, but it will work (eventually) on just about any
    # pedigree you give it. It also is the only method that will allow you to
    # compute summary statistics about coefficients of relationship.
    print(f"\tEntering pyp_nrm.inbreeding() at {pyp_utils.pyp_nice_time()}")
    example_inbreeding = pyp_nrm.inbreeding(example, method="vanraden")
    print(
        f"\tReturning from pyp_nrm.inbreeding() at {pyp_utils.pyp_nice_time()}")
    print(example_inbreeding["metadata"])

    # If all you want are fast calculations of inbreeding coefficients, then you
    # should use the Meuwissen & Luo or modified Meuwissen & Luo methods.
    print(
        f"\tStarted computing inbreeding using Meuwissen and Luo's method  at {pyp_utils.pyp_nice_time()}")
    test_inbreeding_ml = pyp_nrm.inbreeding(example, method="meu_luo")
    print(
        f"\tFinished computing inbreeding using Meuwissen and Luo's method  at {pyp_utils.pyp_nice_time()}")
    print(test_inbreeding_ml["metadata"])

    # All methods should produce the same coefficients of inbreeding.
    print(
        f"\n\tStarted computing inbreeding using Meuwissen and Luo's modified method  at {pyp_utils.pyp_nice_time()}")
    test_inbreeding_qu = pyp_nrm.inbreeding(example, method="mod_meu_luo")
    print(
        f"\tFinished computing inbreeding using Meuwissen and Luo's modified method  at {pyp_utils.pyp_nice_time()}")
    print(test_inbreeding_qu["metadata"])
