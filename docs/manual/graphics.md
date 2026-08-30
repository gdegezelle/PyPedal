# Graphics

PyPedal can draw a pedigree as a Graphviz graph and can plot simple
summaries with matplotlib. Drawing is optional. Core inbreeding does not
require it.

PDF pedigree **reports** are a separate ReportLab extra. See
[PDF reports](pdf-reports.md).

## Extras

| Extra | Packages | Used for |
|---|---|---|
| `graphics` | matplotlib, pydot, graphviz, Pillow | matplotlib helpers; `draw_pedigree` (pydot) |
| `graphviz-extra` | pygraphviz | `new_draw_pedigree` |

You also need the Graphviz **`dot`** program on `PATH`. Python bindings
without the Graphviz binaries cannot write the image.

NetworkX is a **core** dependency. `pyp_network.ped_to_graph` builds an
in-memory directed graph for ancestor queries. It is not a drawing API.

```bash
python -m pip install -e ".[graphics]"
```

## Pedigree drawings

`pyp_graphics.draw_pedigree` uses pydot. The pedigree should be
renumbered (the default load).

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_graphics

work = Path(tempfile.mkdtemp())
pedfile = work / "mrode.ped"
pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n4 1 0\n5 4 3\n6 5 2\n")

ped = load_pedigree(
    options={
        "pedfile": str(pedfile),
        "pedformat": "asd",
        "messages": "quiet",
        "pedigree_summary": 0,
    }
)

stem = work / "mrode_graph"
ok = pyp_graphics.draw_pedigree(
    ped, gfilename=str(stem), gformat="png", gdot=0, gtitle="mrode"
)
png = Path(str(stem) + ".png")
print(ok, png.exists())
```

When pydot and Graphviz are installed that prints `1 True`. `gfilename`
is a **stem**; the function appends `.png` (or `.jpg`, `.ps`). Return
value is `1` on success and `0` on failure.

`pyp_graphics.new_draw_pedigree` is the pygraphviz variant (`graphviz-extra`).

Some matrix-spy and colour-mesh helpers in `pyp_graphics` come from Rick
Muller’s ASPN Python Cookbook recipes. See [Notices](notices.md).

The 2.0 wxPython drawing windows are historical. They are not part of
PyPedal 4.0.
