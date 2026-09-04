# Desktop application

PyPedal 4.2 ships a **PySide6** desktop application. It is a native front
end over the same load path and analysis functions documented in this
manual. It does not reimplement the science in the GUI.

The 2.0.4 **wxPython** GUI is historical. The 4.0–4.1 **CustomTkinter**
application is also historical.

## Install and launch

```bash
python -m pip install -e ".[gui]"
python -m PyPedal
# or
pypedal
```

`pypedal` and `pypedal-gui` are the same desktop. `gui` is a
compatibility alias for the `desktop` extra; both install PySide6.
PyPedal 4.2 does not provide a command-oriented analysis CLI.

PySide6 is imported only when the window starts. The library remains
importable without the GUI extra. Starting the app without PySide6
prints a message to install `PyPedal[gui]`.

`--version` prints the package version and does not open a window.

The window needs a graphical display. Headless environments do not
instantiate the app window.

## Opening a pedigree

Use **File → Open…** (the platform Open shortcut: Command-O on macOS,
Ctrl-O elsewhere). There is no Open toolbar and no in-window Open
button.

The native file chooser is followed by format options:

- **Format** — `pedformat` (default `asd`)
- **Separator** — empty means a space; type `,` for CSV. Surrounding
  spaces are ignored, so `, ` is still a comma.
- **Renumber** — default on

**Open Recent** lists successful loads. A missing recent file is
reported and removed from the list. A failed load keeps the previous
pedigree.

## What the window shows

1. **Metadata** — compact inspector of the loaded pedigree (counts,
   file, format). Not a dump of `stringme()`.
2. **Animals** — virtual table of every animal: current animal ID, sire,
   dam, birth year, sex, name, and inbreeding *F* after an inbreeding
   run. Search and column sort do not copy the pedigree into a second
   list.
3. **Inbreeding** — Meuwissen–Luo (`meu_luo`) coefficients for the
   whole pedigree. Results appear in a table; they are not written as
   `.dat` files on Run.
4. **Inbreeding by Year** — mean *F* by recorded birth year, using the
   cached inbreeding result. It does not rerun Meuwissen–Luo when
   coefficients are already available.
5. **Effective Founders** — Lacy’s *f<sub>e</sub>*.
6. **Relationship** — coefficient of relationship for two animals. Search
   by **display name**, **original ID**, or **current animal ID**, then
   choose a result. Duplicate names are not unique identities; the match
   list includes each matching animal so you can pick the right one. This is not a
   pedigree-wide matrix.
7. **Mating** — prospective offspring inbreeding for an explicit pair,
   or for an explicit small group of pairs you add. It does not mate
   every animal with every other animal. The pair selectors are the same
   name/ID search used on Relationship. Sex is shown in the match list;
   animals are not swapped automatically.
8. **Population** — theoretical Ne from pedigree metadata.

**File → Save Pedigree As…** writes the loaded pedigree to a path you
choose. Analysis pages that have a result can **export** that result as
UTF-8 CSV or text to a path you choose. Run never writes those files
automatically.

**File → Export** can write a metadata PDF and a three-generation
pedigree PDF when the `reports` extra (ReportLab) is installed. You
choose the destination. Missing ReportLab is reported as a dependency
error. PDFs are not generated automatically.

**About PyPedal** (application menu on macOS) shows the version, original
author, maintainer, and GNU LGPL-2.1-or-later.

## Progress and errors

Loads and long analyses run on a background thread. The window stays
usable enough to show progress; Open and analysis actions are disabled
until the job finishes. Typed errors appear in a dialog; they do not
shut the process down with a zero exit status. There is no cancel
button. Single-pair Relationship and Mating use an exact selected
calculation and are interactive on the curated Griffon sample; they do
not construct a dense or ancestor sub-NRM for one pair.

## Animal identifiers in Relationship and Mating

Relationship and Mating let you **search** for an animal by display/call
name, original ID, or current (renumbered) animal ID. Choosing a result
commits that animal. Typing in the search box without choosing a result
does not select an animal, and Compute stays disabled until both sides
have an explicit choice.

Names are **not** unique identities. Two animals can share a call name;
the match list shows original ID, sex, birth year, and current ID so you
can tell them apart. The scientific functions still receive current
`animalID` values.

Pedigrees need an `n` column only for name search. Original ID and
current ID selection work on files without names (including the
scientific Griffon `asdxb` sample).

A pedigree reload clears selected animals because current IDs may
change.

## What it does not do

There is no GUI for genomic matrices, gene dropping, Boichard
*f<sub>a</sub>*, or drawings. Use the Python API for those.
