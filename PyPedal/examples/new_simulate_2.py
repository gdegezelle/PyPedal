#!/usr/bin/env python3

###############################################################################
# NAME: new_simulate.py
# VERSION: 2.0.0b15 (08JUNE2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_graphics
from PyPedal import pyp_network
import networkx

if __name__ == "__main__":

    for i in range(1):
        print("Simulating and drawing pedigree %d" % (i))
        example = pyp_newclasses.load_pedigree(options_file="new_simulate.ini")
        print(example)
