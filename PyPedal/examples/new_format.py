#!/usr/bin/env python3

###############################################################################
# NAME: new_format.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses

import os
ini_file = os.path.join(
    os.path.dirname(__file__), 'new_format.ini')

if __name__ == "__main__":

    # Example taken from Boichard et al. (1997), Figure 2 / Table II.
    example2a = pyp_newclasses.load_pedigree(options_file=ini_file)
    example2a.metadata.printme()
