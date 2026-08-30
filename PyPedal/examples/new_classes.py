#!/usr/bin/env python3

###############################################################################
# NAME: new_classes.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

import os
ini_file = os.path.join(
    os.path.dirname(__file__), 'new_classes.ini')

if __name__ == "__main__":

    print(f"Starting pypedal.py at {pyp_nice_time()}")

    example = pyp_newclasses.load_pedigree(
        options_file=ini_file
    )

    if example.kw["messages"] == "verbose":
        print(
            f"[INFO]: Forming numerator relationship matrix at {pyp_nice_time()}")

    my_a = pyp_nrm.fast_a_matrix_r(example.pedigree, example.kw)

    if example.kw["messages"] == "verbose":
        print(f"[INFO]: Visualizing NRM sparsity at {pyp_nice_time()}")

    pyp_graphics.rmuller_spy_matrix_pil(my_a, fname="boichard2_spy.png")

    if example.kw["messages"] == "verbose":
        print(f"[INFO]: Visualizing NRM in pseudocolor at {pyp_nice_time()}")

    pyp_graphics.rmuller_pcolor_matrix_pil(my_a, fname="boichard2_pcolor.png")

    pyp_graphics.draw_pedigree(
        example,
        gfilename="boichard2_pedigree",
        gtitle="",
        gorient="p",
        gdirec="RL",
        gfontsize=12,
    )

    print(f"Stopping pypedal.py at {pyp_nice_time()}")
