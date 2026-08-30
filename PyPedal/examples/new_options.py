#!/usr/bin/env python3

###############################################################################
# NAME: new_methods.py
# VERSION: 2.0.0a5 (12DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == "__main__":

    print("Starting pypedal.py at %s" % (pyp_nice_time()))

    print("=" * 80)

    # These two loads are expected to fail: no pedfile and no option file.
    print("This load should fail")
    try:
        myped1 = pyp_newclasses.load_pedigree()
        print(myped1)
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}")

    print("-" * 80)

    options = {}
    print("This load should fail")
    try:
        myped2 = pyp_newclasses.load_pedigree(options)
        print(myped2)
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}")

    print("-" * 80)

    # This should work because we are providing a dictionary of
    # options and no configuration file name.
    options = {}
    options["pedfile"] = "new_lacy.ped"
    options["pedformat"] = "asd"
    options["pedname"] = "Lacy Pedigree"
    myped3 = pyp_newclasses.load_pedigree(options)
    print(myped3)

    print("-" * 80)

    # This should work because we are providing an empty dictionary
    # and the name of a valid configuration file.
    options = {}
    myped4 = pyp_newclasses.load_pedigree(options_file="new_options.ini")
    print(myped4)

    print("=" * 80)

    print("Stopping pypedal.py at %s" % (pyp_nice_time()))
