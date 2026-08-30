# PDF reports

PyPedal can write **headless PDF pedigree reports** with ReportLab. This is
not Graphviz drawing and not the desktop app. No GUI is required.

Install the `reports` extra from a local checkout:

```bash
python -m pip install -e ".[reports]"
```

Calling a PDF function without ReportLab raises `PyPedalDependencyError`.
`import PyPedal` does not require ReportLab.

## Metadata report

```python
import tempfile
from pathlib import Path

from PyPedal.pyp_newclasses import load_pedigree
from PyPedal import pyp_reports

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

path = pyp_reports.pdf_pedigree_metadata(
    ped,
    titlepage=1,
    reporttitle="Mrode pedigree metadata",
    reportauthor="PyPedal 4.0",
    reportfile=str(work / "mrode_metadata.pdf"),
)
print(path)
```

`pdf_pedigree_metadata` writes current pedigree metadata as paginated
text and returns the output path. If `reportfile` is omitted, the default
name is `{default_report}_metadata.pdf`. An existing target is
overwritten. The parent directory must already exist.

Unknown birth chronology is counted as unknown. Years such as 1800 or
1900 are real years when they appear.

## Three-generation pedigree

```python
path = pyp_reports.pdf_three_gen_ped(
    5,
    ped,
    reportfile=str(work / "mrode_three_generation.pdf"),
)
print(path)
```

This is the historical **15-slot** layout: the subject, two parents, four
grandparents, and eight great-grandparents. It is not a seven-box
pedigree.

The subject is a **current `animalID`**. Names are displayed when present;
they are not used to select the subject. A single ID or an iterable of
IDs is accepted. Every ID is validated before the PDF is opened. An
invalid ID raises `PyPedalUsageError` and does not create a partial file.
Each subject is one page. Missing parents render as `(Unknown Parent)`.

The default filename, when `reportfile` is omitted, is
`three_generation_pedigrees.pdf`.

## Compatibility names

These names call the same functions:

| Compatibility name | Modern name |
|---|---|
| `pdfPedigreeMetadata` | `pdf_pedigree_metadata` |
| `pdf3GenPed` | `pdf_three_gen_ped` |
| `pdf_3_gen_ped` | `pdf_three_gen_ped` |

Prefer snake_case in new scripts.

`pdf_mean_metric_by` writes a simple metric-by-group PDF. It is not a
pedigree drawing.

The CustomTkinter app has **no PDF menu**. Use these library functions.
