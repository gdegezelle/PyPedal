# Third-party material in PyPedal

The PyPedal library is licensed under LGPL-2.1-or-later (`LICENSE`).
This file records attribution for material that is not original PyPedal
source, and historical dependencies that are **not** distributed.

## Rick Muller / ASPN Python Cookbook (graphics helpers)

Some matrix visualisation routines in `PyPedal/pyp_graphics.py`
(`rmuller_spy_matrix_pil`, `rmuller_pcolor_matrix_pil`, `rmuller_get_color`)
were taken from Rick Muller’s recipes in the ASPN Python Cookbook
(http://aspn.activestate.com/ASPN/Cookbook/Python/) and are used under
the Cookbook terms noted in that source file. They were not written by
John B. Cole. The comments next to those functions must be kept.

## Newfoundland Dog database example

The Newfoundland pedigree example distributed with PyPedal was taken from
the Newfoundland Dog database (http://www.newfoundlanddog-database.net/)
and was used with permission in the original project documentation.
The example remains in `PyPedal/examples/` in a checkout or sdist. Wheels
do not ship examples.

## ADOdb (historical, not distributed)

PyPedal 2.x bundled ADOdb for Python. PyPedal 4 uses only the standard
library `sqlite3` module. No ADOdb-derived source is distributed, and
ADOdb’s license is not a PyPedal project license.

## Ordinary Python dependencies

NumPy, pandas, SciPy, NetworkX, and optional extras (matplotlib, Graphviz
bindings, ReportLab, CustomTkinter) are external packages. They are not
bundled as source in this repository.
