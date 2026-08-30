# Desktop application

PyPedal 4.0 ships a **CustomTkinter** desktop app. It is a thin front end
over the same `NewPedigree` load path and analysis functions documented
in this manual. It does not reimplement the science in the GUI.

The 2.0.4 **wxPython** application is historical. This chapter does not
claim feature parity with it.

## Install and launch

```bash
python -m pip install -e ".[gui]"
python -m PyPedal
# or
pypedal
```

`pypedal-gui` is registered as the same entry point.

CustomTkinter is imported lazily. The library remains importable without
the GUI extra. Starting the app without CustomTkinter prints a message to
install `PyPedal[gui]`.

The window needs a graphical display. Headless environments do not
instantiate the app window.

## What it does

1. **Open a pedigree** — file picker, `pedformat` (default `asdxb`),
   separator (empty means space; type `,` for CSV — surrounding spaces
   are ignored), optional renumber (default on).
2. **Metadata** — pedigree summary text.
3. **List animals** — current `animalID`, sire, dam, birth year (blank
   when unknown), sex, name.
4. **Inbreeding** — `pyp_nrm.inbreeding(..., method="meu_luo")` (not
   tabular).
5. **Effective founders** — Lacy’s *f<sub>e</sub>*.
6. **Inbreeding by year** — mean `meu_luo` *F* grouped by recorded year,
   skipping missing years.
7. **About** — version string.

Analyses that can take a while run on a background thread. Typed errors
are shown as messages; they do not shut the process down with a zero
exit status.

## What it does not do

There is no GUI for genomic matrices, gene dropping, Boichard *f<sub>a</sub>*,
drawings, or PDF reports. Use the Python API for those. Headless PDFs
are [library functions](pdf-reports.md), not app menus.

Inbreeding in the app is **`meu_luo` only**, which is the appropriate
default for large files.
