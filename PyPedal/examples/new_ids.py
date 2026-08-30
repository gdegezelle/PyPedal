#!/usr/bin/env python3

###############################################################################
# NAME: new_ids.py
# VERSION: 2.0.0b5 (14DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses

if __name__ == "__main__":

    example = pyp_newclasses.load_pedigree(options_file="new_ids.ini")

    for _p in example.pedigree:
        print(
            _p.animalID,
            _p.originalID,
            _p.sireID,
            _p.damID,
            _p.name,
            _p.sex,
            _p.gen,
            getattr(_p, "userField", ""),
        )
