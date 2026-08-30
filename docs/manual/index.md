# Introduction

PyPedal is a Python package for **pedigree analysis**. It loads a recorded
animal pedigree, checks common data problems, and computes measures that
animal breeders and population geneticists use every day: inbreeding,
additive relationships, founder and ancestor contributions, and related
summaries. Optional extras add pedigree drawings, headless PDF reports, and
a small desktop application.

The original program was written by **John B. Cole**. PyPedal 4 is a
Python 3 reimplementation of PyPedal 2.0.4, maintained by Geert Degezelle.
This manual describes **PyPedal 4.0**. The current line is a local
**4.0.0-rc8** publication candidate. It is not tagged, not published to
PyPI, and not a final 4.0.0 release.

PyPedal is a research tool. It works from pedigree structure, not from
individual genotypes, except for an optional genomic relationship matrix.
It is not a general genealogy product.

## Who this manual is for

You can use PyPedal if you can run a short Python script. You do not need
to know how the library is implemented. The first chapters walk through
loading a pedigree and calculating inbreeding. Later chapters cover data
preparation, common analyses, reports, and large files. Reference pages at
the end list options, format codes, and scientific citations.

## What you can do with PyPedal 4.0

- Load a pedigree from a text file (integer IDs or unique string identities)
- Check missing parents, duplicate IDs, and inconsistent chronology
- Calculate inbreeding coefficients (*F*)
- Calculate additive (numerator) relationships
- Evaluate the inbreeding of a prospective offspring (`mating_coi`)
  without adding that animal to the pedigree
- Estimate effective founders, effective ancestors, and founder genomes
- Simulate gene dropping
- Write PDF reports and Graphviz drawings (optional extras)
- Open a pedigree in a CustomTkinter desktop app (optional extra)

PyPedal 4.0 requires **Python 3.12 or newer**. See
[Installation](installation.md).

## A first look

The six-animal pedigree from Mrode (2005) Table 2.1 is used throughout
this manual. Animal 5 is inbred. Its inbreeding coefficient is **0.125**.
That example is written inline so it works from an installed wheel; you
do not need the `PyPedal/examples` directory. Full walkthroughs are in
[Your first pedigree](first-pedigree.md) and
[Your first analysis](first-analysis.md).

## How to read this manual

1. **Getting started** — install PyPedal, load a pedigree, run inbreeding.
2. **Preparing data** — files, IDs, dates, checks, and export.
3. **Common analyses** — inbreeding, relationships, matings, founders.
4. **Output** — PDFs, drawings, SQLite, the desktop app.
5. **Large pedigrees** — which algorithm to use, and a ~98,000-animal
   sample.
6. **Recipes** — copy-and-paste workflows.
7. **Reference** — options, format codes, objects, glossary, citations.

Scientific methods are attributed to their published sources. See
[References](references.md) and [Notices](notices.md).
