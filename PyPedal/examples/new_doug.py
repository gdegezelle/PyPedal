#!/usr/bin/env python3

###############################################################################
# NAME: new_doug.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal.pyp_utils import pyp_nice_time

import os
ini_file = os.path.join(
    os.path.dirname(__file__), 'new_doug.ini')

if __name__ == "__main__":

    example = pyp_newclasses.load_pedigree(options_file=ini_file)

    if example.kw["messages"] == "verbose":
        print(f"[INFO]: Calling pyp_graphics.draw_pedigree() at {pyp_nice_time()}")

    pyp_graphics.draw_pedigree(
        example,
        gfilename="doug_below",
        gtitle="Doug the German Shepherd (B)",
        gorient="p",
        gname=1,
        gdirec="",
        gfontsize=12,
        garrow=0,
        gtitloc="b",
    )

    pyp_graphics.draw_pedigree(
        example,
        gfilename="doug_above",
        gtitle="Doug the German Shepherd (A)",
        gorient="p",
        gname=1,
        gdirec="",
        gfontsize=12,
        garrow=0,
        gtitloc="t",
    )

    pyp_graphics.draw_pedigree(
        example,
        gfilename="doug_p_rl_notitle",
        gtitle="",
        gorient="p",
        gname=1,
        gdirec="RL",
        gfontsize=12,
    )
