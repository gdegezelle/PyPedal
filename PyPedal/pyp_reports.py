###############################################################################
# NAME: pyp_reports.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
"""
Headless PDF pedigree reports (ReportLab).

Primary modern APIs:

    pdf_pedigree_metadata
    pdf_three_gen_ped

Historical names ``pdfPedigreeMetadata``, ``pdf3GenPed``, and the pypedal3
spelling ``pdf_3_gen_ped`` are compatibility aliases. They are supported, not
deprecated for 4.0.

ReportLab is optional (the ``reports`` extra). Importing this module does not
require it. Calling a PDF API without ReportLab raises
``PyPedalDependencyError``.

Layout helpers ``_pdfInitialize``, ``_pdfDrawPageFrame``, and
``_pdfCreateTitlePage`` follow the historical 2.0.4 ReportLab geometry
(page frame, title wrapping, 15-slot three-generation pedigree). They are
ported to Python 3 and the current object model; missing parents are never
resolved through ``pedigree[-1]``.
"""

from __future__ import annotations

import logging
import os
import textwrap
from collections.abc import Iterable

from . import pyp_db, pyp_errors, pyp_utils

# Global mappings used by meanMetricBy (surviving metric-report helper).
metric_to_column = {'fa': 'coi'}
byvar_to_column = {
    'by': 'birthyear',
    'gen': 'generation',
    'sex': 'sex'
}

_UNKNOWN_PARENT_LABEL = '(Unknown Parent)'
_DEFAULT_THREE_GEN_FILENAME = 'three_generation_pedigrees.pdf'
_TITLE_WRAP = 26
_AUTHOR_WRAP = 50

# 15 ancestral slots: proband + 2 parents + 4 grandparents + 8 great-grandparents.
_THREE_GEN_SLOTS = (
    'a',
    's', 'd',
    'ss', 'sd', 'ds', 'dd',
    'sss', 'ssd', 'sds', 'sdd', 'dss', 'dsd', 'dds', 'ddd',
)


def _require_reportlab():
    """Import ReportLab or raise PyPedalDependencyError.

    Always imports inside the function so a missing extra cannot leak
    ``ImportError`` / ``NameError`` from module import, and so tests can
    block ReportLab in a subprocess without uninstalling it.
    """
    try:
        from reportlab.lib.pagesizes import A4, letter
        from reportlab.lib.units import cm, inch
        from reportlab.pdfgen import canvas
    except ImportError as exc:
        raise pyp_errors.PyPedalDependencyError(
            "PDF reports require ReportLab. Install the reports extra: "
            "python -m pip install -e '.[reports]'"
        ) from exc
    return canvas, letter, A4, inch, cm


def _pdf_text(value):
    """Coerce a value to a ReportLab Times-Roman (WinAnsi) string."""
    if value is None:
        return ''
    text = str(value)
    try:
        text.encode('latin-1')
        return text
    except UnicodeEncodeError:
        return text.encode('latin-1', 'replace').decode('latin-1')


def _resolve_output_path(reportfile, default_name):
    """Return the output path string. Empty/None uses *default_name*."""
    if reportfile is None or reportfile == '':
        return default_name
    return os.fspath(reportfile)


def _is_missing_parent(parent_id, missing):
    return parent_id is None or parent_id == missing


def _animals_by_id(pedobj):
    return {animal.animalID: animal for animal in pedobj.pedigree}


def _as_current_animal_id(value):
    """Coerce one subject to a current animalID. Names are not IDs."""
    if isinstance(value, bool) or value is None:
        raise pyp_errors.PyPedalUsageError(
            'pdf_three_gen_ped subject IDs must be current animalID integers, '
            f'not {value!r}.'
        )
    if isinstance(value, str):
        raise pyp_errors.PyPedalUsageError(
            'pdf_three_gen_ped subject IDs must be current animalID values; '
            'names are not IDs and are not searched.'
        )
    if isinstance(value, int):
        return value
    shape = getattr(value, 'shape', None)
    if shape == ():
        try:
            return int(value)
        except (TypeError, ValueError) as exc:
            raise pyp_errors.PyPedalUsageError(
                f'pdf_three_gen_ped cannot use {value!r} as a current animalID.'
            ) from exc
    raise pyp_errors.PyPedalUsageError(
        f'pdf_three_gen_ped cannot use {value!r} as a current animalID.'
    )


def _subject_id_list(animal_id):
    """Normalise a single current animalID or an iterable of them."""
    if isinstance(animal_id, (str, bytes)):
        raise pyp_errors.PyPedalUsageError(
            'pdf_three_gen_ped subject IDs must be current animalID values; '
            'names are not IDs and are not searched.'
        )
    if isinstance(animal_id, bool) or animal_id is None:
        raise pyp_errors.PyPedalUsageError(
            'pdf_three_gen_ped requires at least one current animalID.'
        )
    if isinstance(animal_id, int) or getattr(animal_id, 'shape', None) == ():
        return [_as_current_animal_id(animal_id)]
    if isinstance(animal_id, Iterable):
        ids = [_as_current_animal_id(item) for item in animal_id]
        if not ids:
            raise pyp_errors.PyPedalUsageError(
                'pdf_three_gen_ped requires at least one current animalID.'
            )
        return ids
    return [_as_current_animal_id(animal_id)]


def _validate_subjects(pedobj, subject_ids):
    """Validate every current animalID before any PDF is opened."""
    by_id = _animals_by_id(pedobj)
    missing = pedobj.kw.get('missing_parent', 0)
    invalid = []
    for animal_id in subject_ids:
        if _is_missing_parent(animal_id, missing) or animal_id not in by_id:
            invalid.append(animal_id)
    if invalid:
        raise pyp_errors.PyPedalUsageError(
            'pdf_three_gen_ped subject IDs must be current animalID values '
            f'present in the pedigree; invalid: {invalid!r}.'
        )
    return by_id


def _parents_of(animal_id, by_id, missing):
    """Return (sireID, damID) without indexing pedigree[-1]."""
    if _is_missing_parent(animal_id, missing) or animal_id not in by_id:
        return missing, missing
    animal = by_id[animal_id]
    sire = animal.sireID
    dam = animal.damID
    if _is_missing_parent(sire, missing):
        sire = missing
    if _is_missing_parent(dam, missing):
        dam = missing
    return sire, dam


def _three_gen_places(subject_id, by_id, missing):
    """Populate the historical 15-slot ancestral dictionary."""
    places = {'a': subject_id}
    places['s'], places['d'] = _parents_of(places['a'], by_id, missing)
    places['ss'], places['sd'] = _parents_of(places['s'], by_id, missing)
    places['ds'], places['dd'] = _parents_of(places['d'], by_id, missing)
    places['sss'], places['ssd'] = _parents_of(places['ss'], by_id, missing)
    places['sds'], places['sdd'] = _parents_of(places['sd'], by_id, missing)
    places['dss'], places['dsd'] = _parents_of(places['ds'], by_id, missing)
    places['dds'], places['ddd'] = _parents_of(places['dd'], by_id, missing)
    return places


def _slot_caption(place_id, by_id, missing):
    if _is_missing_parent(place_id, missing) or place_id not in by_id:
        return _UNKNOWN_PARENT_LABEL, ''
    animal = by_id[place_id]
    name = animal.name if animal.name not in (None, '') else str(animal.originalID)
    return _pdf_text(name), _pdf_text(f'({animal.originalID})')


def _metadata_report_lines(pedobj):
    """Human-readable metadata lines for the PDF (current PedigreeMetadata)."""
    metadata = getattr(pedobj, 'metadata', None)
    if metadata is None:
        raise pyp_errors.PyPedalUsageError(
            'pdf_pedigree_metadata requires a loaded pedigree with metadata.'
        )
    lines = metadata.stringme().split('\n')
    unknown_years = getattr(metadata, 'num_unknown_birth_years', None)
    if unknown_years is not None:
        lines.append(f'\tUnknown birth years:\t{unknown_years}')
    year_list = getattr(metadata, 'unique_year_list', None)
    if year_list:
        recorded = sorted(year for year in year_list if year is not None)
        if recorded:
            lines.append(
                '\tRecorded birth years:\t' + ', '.join(str(year) for year in recorded)
            )
    return lines


def _new_canvas(outfile, pdf_settings):
    canvas_mod, _letter, _a4, _inch, _cm = _require_reportlab()
    canv = canvas_mod.Canvas(
        outfile,
        pagesize=pdf_settings['_pdfCalcs']['_page'],
        invariant=1,
    )
    canv.setPageCompression(1)
    return canv


##
# meanMetricBy() returns a dictionary of means keyed by levels of the 'byvar'
# for summary statistics.
def meanMetricBy(pedobj, metric='fa', byvar='by', createpdf=False, conn=None):
    conn_created = False
    result_dict = {}
    if conn is None:
        conn = pyp_db.connectToDatabase(pedobj)
        conn_created = True

    try:
        metric = metric_to_column.get(metric, 'fa')
        byvar = byvar_to_column.get(byvar, 'birthyear')

        if pyp_db.doesTableExist(pedobj, conn=conn):
            sql = (
                f"SELECT {byvar}, AVG({metric}) FROM {pedobj.kw['database_table']} "
                f"WHERE {byvar} IS NOT NULL GROUP BY {byvar} ORDER BY {byvar}"
            )
            cursor = conn.cursor()
            cursor.execute(sql)
            for row in cursor:
                result_dict[row[0]] = row[1]
            conn.commit()
        else:
            logging.error('Table does not exist for pyp_reports/meanMetricBy()')

        if conn_created:
            conn.close()

        if createpdf:
            mmbPdfTitle = f"{pedobj.kw['default_report']}_mean_metric_{metric}_{byvar}"
            success = pdfMeanMetricBy(pedobj, result_dict, 1, mmbPdfTitle)
            if success:
                logging.info('PDF report generated successfully.')
            else:
                logging.error('Failed to generate PDF report.')

    except pyp_errors.PyPedalError:
        raise
    except Exception as e:
        logging.error(f"Error in meanMetricBy: {e}")

    return result_dict


##
# pdfMeanMetricBy() generates a PDF report summarizing metrics.
def pdfMeanMetricBy(pedobj, results, titlepage=0, reporttitle='', reportauthor='', reportfile=''):
    _require_reportlab()
    try:
        _pdfOutfile = _resolve_output_path(
            reportfile,
            f"{pedobj.kw['default_report']}_mean_metric_by.pdf",
        )
        _pdfSettings = _pdfInitialize(pedobj)
        canv = _new_canvas(_pdfOutfile, _pdfSettings)

        if titlepage:
            reporttitle = reporttitle or f"Mean Metric By Report\n{pedobj.kw['pedname']}"
            _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)
        _pdfDrawPageFrame(canv, _pdfSettings)

        canv.setFont("Times-Bold", 12)
        tx = canv.beginText(
            _pdfSettings['_pdfCalcs']['_left_margin'],
            _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * _pdfSettings['_pdfCalcs']['_unit']
        )
        for _k, _v in results.items():
            _line = f"\t{_k}:\t{_v}"
            tx.textLine(_pdf_text(_line))
            if tx.getY() < _pdfSettings['_pdfCalcs']['_bottom_margin'] + \
                    0.5 * _pdfSettings['_pdfCalcs']['_unit']:
                canv.drawText(tx)
                canv.showPage()
                _pdfDrawPageFrame(canv, _pdfSettings)
                canv.setFont('Times-Roman', 12)
                tx = canv.beginText(
                    _pdfSettings['_pdfCalcs']['_left_margin'],
                    _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * _pdfSettings['_pdfCalcs']['_unit']
                )
        if tx:
            canv.drawText(tx)
            canv.showPage()
        canv.save()
        return 1
    except pyp_errors.PyPedalError:
        raise
    except Exception as e:
        logging.error(f"Error in pdfMeanMetricBy: {e}")
        return 0


def pdf_pedigree_metadata(
    pedobj,
    titlepage=0,
    reporttitle='',
    reportauthor='',
    reportfile='',
):
    """Write a PDF of current PedigreeMetadata and return the output path.

    Parameters
    ----------
    pedobj
        A loaded ``NewPedigree``.
    titlepage
        If true, a title page is written first.
    reporttitle, reportauthor
        Optional title-page text. An empty title uses the pedigree name.
    reportfile
        Output path (``str`` or ``os.PathLike``). Empty uses
        ``{default_report}_metadata.pdf``. An existing file is overwritten.
        The parent directory must already exist.
    """
    _require_reportlab()
    if getattr(pedobj, 'metadata', None) is None:
        raise pyp_errors.PyPedalUsageError(
            'pdf_pedigree_metadata requires a loaded pedigree with metadata.'
        )
    _pdfOutfile = _resolve_output_path(
        reportfile,
        f"{pedobj.kw['default_report']}_metadata.pdf",
    )
    if pedobj.kw.get('messages') == 'verbose':
        print(f'Writing metadata report to {_pdfOutfile}')
    logging.info('Writing metadata report to %s', _pdfOutfile)

    _pdfSettings = _pdfInitialize(pedobj)
    canv = _new_canvas(_pdfOutfile, _pdfSettings)

    if titlepage:
        if reporttitle == '':
            reporttitle = f"Metadata for Pedigree\n{pedobj.kw['pedname']}"
        _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)

    _pdfDrawPageFrame(canv, _pdfSettings)
    canv.setFont("Times-Bold", 12)
    tx = canv.beginText(
        _pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * _pdfSettings['_pdfCalcs']['_unit'],
    )
    for line in _metadata_report_lines(pedobj):
        tx.textLine(_pdf_text(line))
        if tx.getY() < _pdfSettings['_pdfCalcs']['_bottom_margin'] + \
                0.5 * _pdfSettings['_pdfCalcs']['_unit']:
            canv.drawText(tx)
            canv.showPage()
            _pdfDrawPageFrame(canv, _pdfSettings)
            canv.setFont('Times-Roman', 12)
            tx = canv.beginText(
                _pdfSettings['_pdfCalcs']['_left_margin'],
                _pdfSettings['_pdfCalcs']['_top_margin'] -
                0.5 * _pdfSettings['_pdfCalcs']['_unit'],
            )
    if tx:
        canv.drawText(tx)
        canv.showPage()
    canv.save()
    return _pdfOutfile


def pdf_three_gen_ped(
    animal_id,
    pedobj,
    titlepage=0,
    reporttitle='',
    reportauthor='',
    reportfile='',
):
    """Draw a 15-slot three-generation pedigree PDF and return the output path.

    Subject IDs are current ``animalID`` values (after the default renumber).
    A single ID or an iterable of IDs is accepted; every ID is validated
    before the output file is opened. Invalid IDs raise
    ``PyPedalUsageError`` and leave no new PDF.

    Each subject occupies one page: the proband, two parents, four
    grandparents, and eight great-grandparents (15 slots). Missing parents
    are labelled ``(Unknown Parent)``.
    """
    _require_reportlab()
    subject_ids = _subject_id_list(animal_id)
    by_id = _validate_subjects(pedobj, subject_ids)
    missing = pedobj.kw.get('missing_parent', 0)

    _pdfOutfile = _resolve_output_path(reportfile, _DEFAULT_THREE_GEN_FILENAME)
    if pedobj.kw.get('messages') == 'verbose':
        print(f'Writing 3GenPed to {_pdfOutfile}')
    logging.info('Writing 3GenPed to %s', _pdfOutfile)

    _pdfSettings = _pdfInitialize(pedobj)
    canv = _new_canvas(_pdfOutfile, _pdfSettings)

    if titlepage:
        if reporttitle == '':
            reporttitle = 'Three-generation Pedigrees'
        _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)

    for subject_id in subject_ids:
        _draw_three_gen_page(canv, _pdfSettings, pedobj, by_id, subject_id, missing)
    canv.save()
    return _pdfOutfile


def _draw_three_gen_page(canv, pdf_settings, pedobj, by_id, subject_id, missing):
    """Render one 15-slot pedigree page. Geometry follows PyPedal 2.0.4."""
    places = _three_gen_places(subject_id, by_id, missing)
    _pdfDrawPageFrame(canv, pdf_settings)

    calcs = pdf_settings['_pdfCalcs']
    sill_width = calcs['_frame_width'] * 0.25
    sixteenth = calcs['_frame_height'] / 16.0
    sixtyfourth = calcs['_frame_height'] / 64.0
    x = calcs['_left_margin']
    y = calcs['_bottom_margin']

    os_xy = {key: {} for key in _THREE_GEN_SLOTS}
    os_xy['a']['x'] = x
    os_xy['a']['y'] = y + (8 * sixteenth)
    os_xy['d']['x'] = x + sill_width
    os_xy['d']['y'] = y + (4 * sixteenth)
    os_xy['s']['x'] = x + sill_width
    os_xy['s']['y'] = y + (12 * sixteenth)
    os_xy['dd']['x'] = x + (2 * sill_width)
    os_xy['dd']['y'] = y + (2 * sixteenth)
    os_xy['ds']['x'] = x + (2 * sill_width)
    os_xy['ds']['y'] = y + (6 * sixteenth)
    os_xy['sd']['x'] = x + (2 * sill_width)
    os_xy['sd']['y'] = y + (10 * sixteenth)
    os_xy['ss']['x'] = x + (2 * sill_width)
    os_xy['ss']['y'] = y + (14 * sixteenth)
    os_xy['ddd']['x'] = x + (3 * sill_width)
    os_xy['ddd']['y'] = y + (1 * sixteenth)
    os_xy['dds']['x'] = x + (3 * sill_width)
    os_xy['dds']['y'] = y + (3 * sixteenth)
    os_xy['dsd']['x'] = x + (3 * sill_width)
    os_xy['dsd']['y'] = y + (5 * sixteenth)
    os_xy['dss']['x'] = x + (3 * sill_width)
    os_xy['dss']['y'] = y + (7 * sixteenth)
    os_xy['sdd']['x'] = x + (3 * sill_width)
    os_xy['sdd']['y'] = y + (9 * sixteenth)
    os_xy['sds']['x'] = x + (3 * sill_width)
    os_xy['sds']['y'] = y + (11 * sixteenth)
    os_xy['ssd']['x'] = x + (3 * sill_width)
    os_xy['ssd']['y'] = y + (13 * sixteenth)
    os_xy['sss']['x'] = x + (3 * sill_width)
    os_xy['sss']['y'] = y + (15 * sixteenth)

    subject = by_id[subject_id]
    header = _pdf_text(f'Pedigree for {subject.name} ({subject.originalID})')
    canv.setFont('Times-Bold', 12)
    canv.drawString(x, calcs['_top_margin'] - 0.25 * sixteenth, header)
    canv.setLineWidth(1)

    for key in _THREE_GEN_SLOTS:
        canv.line(
            os_xy[key]['x'], os_xy[key]['y'],
            os_xy[key]['x'] + sill_width, os_xy[key]['y'],
        )
        line1, line2 = _slot_caption(places[key], by_id, missing)
        canv.setFont('Times-Bold', 12)
        canv.drawString(os_xy[key]['x'], os_xy[key]['y'] + 2, line1)
        canv.setFont('Times-Roman', 12)
        if line2:
            canv.drawString(os_xy[key]['x'], os_xy[key]['y'] - sixtyfourth, line2)
        if key == 'a':
            _draw_subject_details(canv, subject, pedobj, x, y, sill_width, sixteenth)

    canv.line(os_xy['s']['x'], os_xy['d']['y'], os_xy['s']['x'], os_xy['s']['y'])
    canv.line(os_xy['ds']['x'], os_xy['dd']['y'], os_xy['ds']['x'], os_xy['ds']['y'])
    canv.line(os_xy['ss']['x'], os_xy['sd']['y'], os_xy['ss']['x'], os_xy['ss']['y'])
    canv.line(os_xy['dds']['x'], os_xy['ddd']['y'], os_xy['dds']['x'], os_xy['dds']['y'])
    canv.line(os_xy['dss']['x'], os_xy['dsd']['y'], os_xy['dss']['x'], os_xy['dss']['y'])
    canv.line(os_xy['sds']['x'], os_xy['sdd']['y'], os_xy['sds']['x'], os_xy['sds']['y'])
    canv.line(os_xy['sss']['x'], os_xy['ssd']['y'], os_xy['sss']['x'], os_xy['sss']['y'])
    canv.showPage()


def _draw_subject_details(canv, subject, pedobj, x, y, sill_width, sixteenth):
    """Footer fields that remain current on NewAnimal."""
    canv.setFont('Times-Roman', 12)
    missing_herd = pedobj.kw.get('missing_herd', 'Unknown_Herd')
    herd = getattr(subject, 'originalHerd', missing_herd)
    if herd in ('u', None, ''):
        herd = missing_herd
    breed = getattr(subject, 'breed', pedobj.kw.get('missing_breed', 'unknown'))
    inbreed = getattr(subject, 'fa', pedobj.kw.get('missing_inbreeding', 0.0))
    pedcomp = getattr(subject, 'pedcomp', pedobj.kw.get('missing_pedcomp', -999.0))
    missing_pedcomp = pedobj.kw.get('missing_pedcomp', -999.0)

    canv.drawString(x, y + 0.25 * sixteenth, 'Pedigree completeness:')
    if pedcomp == missing_pedcomp:
        canv.drawString(x + sill_width, y + 0.25 * sixteenth, 'unknown')
    else:
        try:
            canv.drawString(x + sill_width, y + 0.25 * sixteenth, f'{float(pedcomp):5.3f}')
        except (TypeError, ValueError):
            canv.drawString(x + sill_width, y + 0.25 * sixteenth, _pdf_text(pedcomp))
    canv.drawString(x, y + 0.5 * sixteenth, 'Inbreeding:')
    canv.drawString(x + sill_width, y + 0.5 * sixteenth, _pdf_text(inbreed))
    canv.drawString(x, y + 0.75 * sixteenth, 'Breed:')
    canv.drawString(x + sill_width, y + 0.75 * sixteenth, _pdf_text(breed))
    canv.drawString(x, y + sixteenth, 'Herd:')
    canv.drawString(x + sill_width, y + sixteenth, _pdf_text(herd))


###############################################################################
# _pdfDrawPageFrame() was taken from the procedure drawPageFrame() included in
# the demo program odyssey.py in the ReportLab distribution.  _pdfInitialize()
# is rolled together using some of the code in odyssey.py as an example, as
# well as some of my own work.
###############################################################################

##
# _pdfInitialize() returns a dictionary of metadata that is used for report
# generation.
# @param pedobj A PyPedal pedigree object.
# @retval A dictionary of metadata that is used for report generation.
def _pdfInitialize(pedobj):
    """
    _pdfInitialize() returns a dictionary of metadata that is used for report
    generation.
    """
    _canvas_mod, letter, A4, inch, cm = _require_reportlab()
    _pdfSettings = {
        '_pdfTitle': pedobj.kw.get('pedname', 'Unknown Title'),
        '_pdfPageinfo': pedobj.kw.get('filetag', 'Unknown Info'),
    }
    unit = inch if pedobj.kw.get('default_unit', 'inch') == 'inch' else cm
    if pedobj.kw.get('paper_size', 'letter') == 'letter':
        page = letter
    else:
        page = A4
    page_width, page_height = page
    _pdfCalcs = {
        '_unit': unit,
        '_page': page,
        '_page_width': page_width,
        '_page_height': page_height,
        '_top_margin': page_height - inch,
        '_bottom_margin': inch,
        '_left_margin': inch,
        '_right_margin': page_width - inch,
    }
    _pdfCalcs['_frame_width'] = _pdfCalcs['_right_margin'] - _pdfCalcs['_left_margin']
    _pdfCalcs['_frame_height'] = _pdfCalcs['_top_margin'] - _pdfCalcs['_bottom_margin']
    _pdfSettings['_pdfCalcs'] = _pdfCalcs
    return _pdfSettings


##
# _pdfDrawPageFrame() nicely frames page contents and includes the
# document title in a header and the page number in a footer.
# @param canv An instance of a ReportLab Canvas object.
# @param _pdfSettings An options dictionary created by _pdfInitialize().
# @retval None
def _pdfDrawPageFrame(canv, _pdfSettings):
    """
    _pdfDrawPageFrame() nicely frames page contents and includes the
    document title in a header and the page number in a footer.
    """
    calcs = _pdfSettings['_pdfCalcs']
    canv.line(
        calcs['_left_margin'], calcs['_top_margin'],
        calcs['_right_margin'], calcs['_top_margin'],
    )
    canv.line(
        calcs['_left_margin'], calcs['_bottom_margin'],
        calcs['_right_margin'], calcs['_bottom_margin'],
    )
    canv.setFont('Times-Italic', 12)
    canv.drawString(
        calcs['_left_margin'],
        calcs['_top_margin'] + 2,
        _pdf_text(_pdfSettings.get('_pdfTitle', '')),
    )
    timestamp = pyp_utils.pyp_nice_time()
    if timestamp:
        canv.drawString(
            calcs['_right_margin'] - 1.85 * calcs['_unit'],
            calcs['_top_margin'] + 2,
            _pdf_text(timestamp),
        )
    canv.drawCentredString(
        0.5 * calcs['_page_width'],
        0.5 * calcs['_unit'],
        f"Page {canv.getPageNumber()}",
    )


def _centred_wrapped_lines(canv, page_width, start_y, unit, text, wrap_width, font, size):
    """Draw centred wrapped lines; return the next y position."""
    canv.setFont(font, size)
    y = start_y
    for chunk in str(text).split('\n'):
        wrapped = textwrap.wrap(chunk, wrap_width, break_long_words=True) or ['']
        for line in wrapped:
            canv.drawCentredString(0.5 * page_width, y, _pdf_text(line))
            y -= unit
    return y


##
# _pdfCreateTitlePage() adds a title page to a ReportLab canvas object.
# @param canv An instance of a ReportLab Canvas object.
# @param _pdfSettings An options dictionary created by _pdfInitialize().
# @param reporttitle Title of report; if '', _pdfTitle is used.
# @param reportauthor Author/preparer of report.
# @retval None
def _pdfCreateTitlePage(canv, _pdfSettings, reporttitle='', reportauthor=''):
    """
    _pdfCreateTitlePage() adds a title page to a ReportLab canvas object.
    """
    _pdfDrawPageFrame(canv, _pdfSettings)
    calcs = _pdfSettings['_pdfCalcs']
    title_y = 7 * calcs['_unit']
    title = reporttitle if reporttitle != '' else _pdfSettings['_pdfTitle']
    title_y = _centred_wrapped_lines(
        canv, calcs['_page_width'], title_y, calcs['_unit'],
        title, _TITLE_WRAP, 'Times-Bold', 36,
    )
    if reportauthor != '':
        _centred_wrapped_lines(
            canv, calcs['_page_width'], title_y, calcs['_unit'],
            reportauthor, _AUTHOR_WRAP, 'Times-Bold', 18,
        )
    canv.showPage()


# Historical / pypedal3 compatibility names. Supported in 4.0; no warning.
pdfPedigreeMetadata = pdf_pedigree_metadata
pdf3GenPed = pdf_three_gen_ped
pdf_3_gen_ped = pdf_three_gen_ped
pdf_mean_metric_by = pdfMeanMetricBy
mean_metric_by = meanMetricBy
