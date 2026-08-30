#!/usr/bin/env python3

###############################################################################
# NAME: new_reporting.py
# VERSION: 4.0.0-rc4
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################
# Supported 4.0 workflow: load a small pedigree and write headless PDF
# reports with ReportLab (the reports extra). No GUI is required.
###############################################################################

import tempfile
from pathlib import Path

from PyPedal import pyp_newclasses
from PyPedal import pyp_reports

if __name__ == "__main__":

    example = pyp_newclasses.load_pedigree(options_file="new_reporting.ini")
    out_dir = Path(tempfile.mkdtemp(prefix="pypedal_pdf_reports_"))

    metadata_path = pyp_reports.pdf_pedigree_metadata(
        example,
        titlepage=1,
        reporttitle="Mrode pedigree metadata",
        reportfile=str(out_dir / "mrode_metadata.pdf"),
    )
    three_path = pyp_reports.pdf_three_gen_ped(
        5,
        example,
        reportfile=str(out_dir / "mrode_three_generation.pdf"),
    )
    print(f"metadata PDF: {metadata_path}")
    print(f"three-generation PDF: {three_path}")
    print(f"reports written under {out_dir}")
