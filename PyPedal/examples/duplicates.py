#!/usr/bin/env python3

#######################################################################
# NAME: duplicates.py
# VERSION: 2.0.0b15 (18SEPTEMBER2006); updated for PyPedal 4.0
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
#######################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_utils

options = {
    "pedname": "Test Duplicate Handling",
    "pedformat": "asd",
    "pedfile": "duplicates.ped",
    "sepchar": ",",
    "messages": "verbose",
    "debug": True,
    "assign_sexes": True,
    "renumber": False,
    "reorder": False,
}

if __name__ == "__main__":

    example = pyp_newclasses.load_pedigree(options=options, debug_load=True)

    print("Duplicate detection (not DUPLICATE_REDIRECT / merge):")
    for d in pyp_utils.list_duplicates(example):
        print(d)
    print(
        "PyPedal 4.0 refuses to reorder/renumber a pedigree that contains "
        "duplicate animal IDs. There is no DUPLICATE_REDIRECT."
    )
