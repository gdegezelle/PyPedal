#!/usr/bin/env python3

###############################################################################
# NAME: new_graphics.py
# VERSION: 2.0.0b5 (19DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm

import os
test_file = test_file = os.path.join(
    os.path.dirname(__file__), 'new_renumbering.ini')

if __name__ == "__main__":
    example = pyp_newclasses.load_pedigree(options_file=test_file)

    print("-" * 80)
    print(f"INFO: Pedigree ID map: {example.idmap}")
    print("-" * 80)
    print(f"INFO: Pedigree ID reverse map: {example.backmap}")
    print("-" * 80)
    my_inbreeding = pyp_nrm.inbreeding(example)
    print(f"INFO: Coefficients of inbreeding: {my_inbreeding}")

    # pyp_graphics.draw_pedigree(
    #     example, gfilename="new_renumbering", gtitle="My  Pedigree"
    # )
