#!/usr/bin/env python3

###############################################################################
# NAME: new_mating_coi.py
# VERSION: 4.0.0-rc4
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################
# Read-only test mating on Mrode (2005) Table 2.1. No phantom offspring,
# no output files, no mate-selection engine.
###############################################################################

from PyPedal import pyp_metrics
from PyPedal.pyp_newclasses import load_pedigree
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == "__main__":
    print("Starting new_mating_coi.py at %s" % (pyp_nice_time(),))
    ped = load_pedigree(
        options={
            "pedfile": "mrode.ped",
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
        }
    )

    print("mating_coi(1, 2) =", pyp_metrics.mating_coi(1, 2, ped))
    print("mating_coi(1, 6) =", pyp_metrics.mating_coi(1, 6, ped))
    print("mating_coi(3, 4) =", pyp_metrics.mating_coi(3, 4, ped))
    print("mating_coi(5, 5) =", pyp_metrics.mating_coi(5, 5, ped))

    group = pyp_metrics.mating_coi_group([(1, 2), (1, 6), (3, 4)], ped)
    print("group matings:", group["matings"])
    print("group metadata:", group["metadata"])
    print("Finished new_mating_coi.py at %s" % (pyp_nice_time(),))
