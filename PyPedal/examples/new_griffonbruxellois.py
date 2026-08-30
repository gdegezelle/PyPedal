#!/usr/bin/env python3

"""Load the curated Griffon Bruxellois 2026 export with true asdxb chronology."""

import os

from PyPedal import pyp_newclasses, pyp_utils

if __name__ == "__main__":
    pedfile = os.path.join(
        os.path.dirname(__file__), "griffonbruxellois_2026_pyp.ped"
    )
    print("Starting at %s" % pyp_utils.pyp_nice_time())
    ped = pyp_newclasses.load_pedigree(
        options={
            "pedfile": pedfile,
            "pedformat": "asdxb",
            "sepchar": ",",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )
    print("n = %d" % len(ped.pedigree))
    print("implicit parents = %d" % ped.metadata.num_implicit_parents)
    print("unknown birth years = %d" % ped.metadata.num_unknown_birth_years)
    print("Finished at %s" % pyp_utils.pyp_nice_time())
