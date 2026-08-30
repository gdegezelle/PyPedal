# Notices

## Authorship and maintenance

PyPedal was originally written by **John B. Cole**.

The current maintainer of this Python 3 line is **Geert Degezelle**.

This manual describes PyPedal 4.0. It is newly written documentation of
the current library. It is not a relicensed copy of the 2012 PyPedal 2
manual.

## Library license

The PyPedal library is distributed under the **GNU Lesser General Public
License version 2.1 or later**. The license text is the `LICENSE` file
at the repository root (`LGPL-2.1-or-later`).

File-level modification notices for original PyPedal source that was
changed for the Python 3 release are in the affected modules. This page
does not replace them and does not assert a new copyright ownership
claim.

## Contributors

People known to have contributed code, files, testing, or suggestions
are listed in `PyPedal/AUTHORS.txt` and `PyPedal/CREDITS.txt` in the
repository. Those files are product credit, not a legal copyright
assignment. Historical testers and correspondents named there include
(among others) Kathy Hanford, Paul VanRaden, Gregor Gorjanc, Bradley J.
Heins, and others recorded by the original project.

## Third-party material that remains in the product

**Graphics helpers.** Some matrix visualisation routines in
`PyPedal/pyp_graphics.py` were taken from Rick Muller’s recipes in the
ASPN Python Cookbook and are used under the terms noted in that source
file. They were not written by John B. Cole. The original comments in
`pyp_graphics.py` must be kept.

**Newfoundland example.** The Newfoundland pedigree example distributed
with PyPedal was taken from the Newfoundland Dog database
(http://www.newfoundlanddog-database.net/) and was used with permission
in the original project documentation. The example file remains in
`PyPedal/examples/` in a checkout or sdist.

**ADOdb.** PyPedal 2.x bundled ADOdb for Python. PyPedal 4 uses only the
standard library `sqlite3` module. No ADOdb-derived source remains in
the package.

## This manual

Small names, formulas, API identifiers, scientific terminology, and
bibliographic citations are used as facts. Long passages from the 2.0
LaTeX manual are not reproduced here.
