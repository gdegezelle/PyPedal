#!/usr/bin/env python3

###############################################################################
# NAME: pyp_io.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################

from configparser import ConfigParser, Error
import logging, pickle, time
from typing import Optional, Dict
import numpy as np

from . import pyp_chronology, pyp_utils, pyp_nrm, pyp_errors
from pathlib import Path

global LINE1, LINE2
LINE1 = '=' * 80
LINE2 = '-' * 80

logging.basicConfig(level=logging.INFO)

def a_inverse_to_file(pedobj, ainv=''):
    """
    Write the inverse of a relationship matrix to a file using pickle.
    """
    try:
        logging.info('Entered a_inverse_to_file()')
        if not ainv:
            ainv = pyp_nrm.a_inverse_df(pedobj)
        filetag = pedobj.kw['filetag']
        output_file = f"{filetag}_a_inverse_pickled.pkl"
        with open(output_file, 'wb') as f:
            pickle.dump(ainv, f)
        logging.info(f"A-inverse saved to {output_file}")
        return True
    except Exception as e:
        logging.error(f"Failed to write A-inverse to file: {e}")
        return False


def a_inverse_from_file(inputfile):
    """
    Read the inverse of a relationship matrix from a file using pickle.
    """
    try:
        logging.info('Entered a_inverse_from_file()')
        with open(inputfile, 'rb') as f:
            a_inv = pickle.load(f)
        logging.info('Successfully loaded A-inverse')
        return a_inv
    except Exception as e:
        logging.error(f"Failed to load A-inverse: {e}")
        return np.zeros((1, 1), dtype=float)


def dissertation_pedigree_to_file(pedobj):
    """
    Write a pedigree in 'asdxfg' format to a file.
    """
    try:
        logging.info('Entered dissertation_pedigree_to_file()')
        length = len(pedobj.pedigree)
        outputfile = f"{pedobj.kw['filetag']}_diss.ped"
        with open(outputfile, 'w') as f:
            f.write('# DISSERTATION pedigree produced by PyPedal.\n')
            f.write('% asdbxfg\n')
            for animal in pedobj.pedigree:
                f.write(f"{animal.animalID},{animal.sireID},{animal.damID},"
                        f"{pyp_chronology.format_year_token(animal.by)},"
                        f"{animal.sex},{animal.fa},{animal.gen}\n")
        logging.info(f"Dissertation pedigree written to {outputfile}")
        return True
    except Exception as e:
        logging.error(f"Failed to write dissertation pedigree: {e}")
        return False


def dissertation_pedigree_to_pedig_format(pedobj):
    """
    Format a pedigree for Didier Boichard's 'pedig' suite of programs and write to a file.
    """
    try:
        logging.info('Entered dissertation_pedigree_to_pedig_format()')
        outputfile = f"{pedobj.kw['filetag']}_pedig.ped"
        with open(outputfile, 'w') as f:
            for animal in pedobj.pedigree:
                sex = 1 if animal.sex.lower() == 'm' else 2
                f.write(f"{animal.animalID} {animal.sireID} {animal.damID} "
                        f"{pyp_chronology.format_year_token(animal.by)} {sex} 1 1\n")
        logging.info(f"Pedigree written to {outputfile}")
        return True
    except Exception as e:
        logging.error(f"Failed to write pedigree format: {e}")
        return False


def dissertation_pedigree_to_pedig_interest_format(pedobj):
    """
    Format a pedigree for the 'parente' program's studied individuals file.
    """
    try:
        logging.info('Entered dissertation_pedigree_to_pedig_interest_format()')
        outputfile = f"{pedobj.kw['filetag']}_parente.ped"
        with open(outputfile, 'w') as f:
            for animal in pedobj.pedigree:
                f.write(f"{animal.animalID} 1\n")
        logging.info(f"Pedigree interest format written to {outputfile}")
        return True
    except Exception as e:
        logging.error(f"Failed to write pedigree interest format: {e}")
        return False


def dissertation_pedigree_to_pedig_format_mask(pedobj: object) -> bool:
    """
    Takes a pedigree in 'asdbxfg' format, formats it into the form used by Didier Boichard's 'pedig' suite of programs,
    and writes it to a file. This function masks the generation ID with a fake birth year and writes the fake birth
    year to the file instead of the true birth year. This is an attempt to fool PEDIG to get f_e, f_a et al. by generation.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.

    Returns
    -------
    bool
        True on success, False on failure.
    """
    logging.info("Entered dissertation_pedigree_to_pedig_format_mask()")
    
    try:
        length = len(pedobj.pedigree)
        outputfile = f"{pedobj.kw['filetag']}_pedig_mask.ped"
        
        with open(outputfile, 'w') as aout:
            for l in range(length):
                # Mask generations
                mygen = float(pedobj.pedigree[l].gen)
                if 0 < mygen <= 1.25:
                    _gen = 10
                elif 1.25 < mygen <= 1.75:
                    _gen = 15
                elif 1.75 < mygen <= 2.25:
                    _gen = 20
                elif 2.25 < mygen <= 2.75:
                    _gen = 25
                elif 2.75 < mygen <= 3.25:
                    _gen = 30
                elif 3.25 < mygen <= 3.75:
                    _gen = 35
                elif 3.75 < mygen <= 4.25:
                    _gen = 40
                elif 4.25 < mygen <= 4.75:
                    _gen = 45
                elif 4.75 < mygen <= 5.25:
                    _gen = 50
                elif 5.25 < mygen <= 5.75:
                    _gen = 55
                elif 5.75 < mygen <= 6.25:
                    _gen = 60
                elif 6.25 < mygen <= 6.75:
                    _gen = 65
                elif 6.75 < mygen <= 7.25:
                    _gen = 70
                elif 7.25 < mygen <= 7.75:
                    _gen = 75
                else:
                    _gen = 0
                
                _maskgen = 1950 + _gen

                # Convert sexes
                sex = 1 if pedobj.pedigree[l].sex.lower() == 'm' else 2
                
                # Write to file
                aout.write(f"{pedobj.pedigree[l].animalID} {pedobj.pedigree[l].sireID} "
                           f"{pedobj.pedigree[l].damID} {_maskgen} {sex} 1 1\n")
        
        logging.info("Exited dissertation_pedigree_to_pedig_format_mask()")
        return True
    
    except Exception as e:
        logging.error(f"Error in dissertation_pedigree_to_pedig_format_mask: {e}")
        return False


#: Section synthesised for legacy option files that have no header at all.
DEFAULT_INI_SECTION = 'default'


def coerce_ini_value(value: str):
    """
    Convert an INI value string to the type PyPedal's options actually need.

    ConfigParser returns every value as a string, and PyPedal's `kw` dict is
    read with plain truth tests. ``renumber = 0`` therefore arrived as the
    string ``'0'``, which is TRUE in Python, so turning an option off in an
    option file turned it on.

    Coercion order matters: booleans before integers, because ``True``/``False``
    must not fall through to the int branch, and integers before floats so that
    ``rounds = 100`` stays an int.
    """
    if not isinstance(value, str):
        return value
    text = value.strip()

    # Strip one layer of matching quotes, so messages = 'verbose' compares
    # equal to messages = verbose rather than to "'verbose'".
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        return text[1:-1]

    lowered = text.lower()
    if lowered in ('true', 'yes', 'on'):
        return True
    if lowered in ('false', 'no', 'off'):
        return False
    if lowered in ('none', 'null'):
        return None

    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def read_ini_file(kwfile: str, section: Optional[str] = None) -> Dict[str, object]:
    """
    Read a PyPedal option file into a FLAT, type-coerced options dictionary.

    Three things the previous implementation got wrong:

      * it rejected section-less files, which is the format 17 of the 26
        example option files actually use;
      * it returned ``{section: {key: value}}``, a shape nothing downstream
        consumes -- every PyPedal option is read as ``kw['renumber']``, never
        ``kw['analysis']['renumber']``;
      * it left every value a string, so falsy options were true.

    Sections are still accepted, and are flattened, because in most option
    files they are presentational grouping rather than namespacing -- PyPedal
    option names are global and every consumer reads kw['renumber'].

    A few files use sections differently, to hold several alternative
    configurations in one place: examples/new_inbreeding2multiple.ini defines a
    [noinbreeding] pedigree and a [horse] pedigree. Flattening those would
    silently merge two pedigrees into one. So a key defined in more than one
    section with different values raises unless `section` selects which one is
    wanted.

    Parameters
    ----------
    kwfile : str
        Path to the option file.
    section : str, optional
        Read only this section. Use it for files that hold several alternative
        configurations. When omitted, all sections are flattened together.
    """
    parser = ConfigParser()
    # Preserve the case of option names: PyPedal has camelCase keys such as
    # debugLoad, and ConfigParser lowercases them by default.
    parser.optionxform = str

    with open(kwfile, 'r', encoding='utf-8') as handle:
        text = handle.read()

    if not _has_section_header(text):
        # Legacy section-less file: synthesise a header so ConfigParser will
        # accept it, rather than rejecting a supported input format.
        text = '[%s]\n%s' % (DEFAULT_INI_SECTION, text)

    parser.read_string(text)

    if section is not None:
        if not parser.has_section(section):
            raise ValueError(
                "Option file %r has no section [%s]. Available sections: %s"
                % (kwfile, section, ", ".join(parser.sections()) or "(none)"))
        wanted = [section]
    else:
        wanted = parser.sections()

    flat: Dict[str, object] = {}
    origin: Dict[str, str] = {}
    for name in wanted:
        for key, value in parser.items(name):
            coerced = coerce_ini_value(value)
            if key in flat and flat[key] != coerced:
                raise ValueError(
                    "Option %r is defined in both [%s] and [%s] of %r with "
                    "different values (%r vs %r). This file holds several "
                    "alternative configurations rather than one; pass "
                    "section='%s' or section='%s' to choose. Available "
                    "sections: %s"
                    % (key, origin[key], name, kwfile, flat[key], coerced,
                       origin[key], name, ", ".join(parser.sections())))
            flat[key] = coerced
            origin[key] = name
    return flat


def _has_section_header(text: str) -> bool:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith('[') and stripped.endswith(']'):
            return True
    return False


def read_ini_to_dict(kwfile: str) -> Optional[Dict[str, dict]]:
    """
    Return a configuration dictionary from an INI file.

    Parameters
    ----------
    kwfile : str
        The path to the INI file.

    Returns
    -------
    dict or None
        A dictionary containing the configuration, structured by sections.
        Returns None if an error occurs.
    """
    config = ConfigParser()

    try:
        # Read the INI file
        config.read(kwfile)

        if not config.sections():
            raise ValueError(f"INI file '{kwfile}' has no sections or is empty.")

        # Build a dictionary structured by sections
        kw = {section: dict(config.items(section)) for section in config.sections()}
        return kw
    
    except (Error, ValueError) as e:
        logging.error(f"Error reading INI file '{kwfile}': {e}")
        return None
    

def pyp_file_header(ofhandle, caller="Unknown PyPedal routine"):
    """
    Writes a header to a page of PyPedal output.

    :param ofhandle: A file handle opened for writing.
    :type ofhandle: file-like object
    :param caller: A string indicating the name of the calling routine, defaults to "Unknown PyPedal routine".
    :type caller: str, optional
    :returns: None
    :rtype: None
    """
    try:
        ofhandle.write(f"{'-' * 80}\n")
        ofhandle.write(f"Created by {caller} at {pyp_utils.pyp_nice_time()}\n")
        ofhandle.write(f"{'-' * 80}\n")
    except Exception as e:  # Explicitly catch exceptions and log them
        print(f"Error writing header: {e}")


def pyp_file_footer(ofhandle, caller: str = "Unknown PyPedal routine") -> None:
    """
    Writes a footer to a page of PyPedal output.

    Parameters
    ----------
    ofhandle : file-like object
        A writable file-like object (e.g., a file handle).
    caller : str, optional
        A string indicating the name of the calling routine, by default "Unknown PyPedal routine".
    """
    try:
        ofhandle.write(f"{'-' * 80}\n")
    except Exception as e:
        print(f"Error writing footer: {e}")


def render_title(title_string: str, title_level: int = 1, output_type: str = "html") -> str:
    """
    Renders page titles (produces HTML output by default).

    Parameters
    ----------
    title_string : str
        The string to be enclosed in HTML "H" tags or underlined.
    title_level : int, optional
        Size to be attached to "H" tags, e.g., "H1". Defaults to 1.
    output_type : str, optional
        The output type ("html" or "text"). Defaults to "html".

    Returns
    -------
    str
        Rendered title string in the specified format.
    """
    # Validate title_level
    if not (1 <= title_level <= 3):
        title_level = 1

    # Generate the title based on the output type
    if output_type.lower() in ["html", "h"]:
        rendered_title = f"<H{title_level}>{title_string}</H{title_level}>\n"
    else:
        underline = "=" * len(title_string)
        rendered_title = f"{title_string}\n{underline}"

    return rendered_title


def render_body_text(text_string: str, output_type: str = "html") -> str:
    """
    Renders page contents (produces HTML output by default).

    Parameters
    ----------
    text_string : str
        The string to be rendered with either a trailing newline or enclosed in HTML "P" tags.
    output_type : str, optional
        The output type ("html" or "text"). Defaults to "html".

    Returns
    -------
    str
        Rendered body text in the specified format.
    """
    # Generate the body text based on the output type
    if output_type.lower() in ["html", "h"]:
        rendered_body_text = f"<p>{text_string}</p>"
    else:
        rendered_body_text = f"{text_string}\n"

    return rendered_body_text


def pickle_pedigree(pedobj, filename: str = "") -> int:
    """
    Pickles a pedigree object.

    Parameters
    ----------
    pedobj : object
        An instance of a PyPedal pedigree object.
    filename : str, optional
        The name of the file to which the pedigree object should be pickled. 
        If not provided, the filename is generated using the object's filetag.

    Returns
    -------
    int
        1 on success, 0 otherwise.
    """
    logging.info("Entered pickle_pedigree()")
    try:
        # Determine the pickle filename
        if not filename:
            pfn = f"{pedobj.kw['filetag']}.pkl"
        else:
            pfn = f"{filename}.pkl"

        # Pickle the object
        with open(pfn, "wb") as outfile:
            pickle.dump(pedobj, outfile)
        
        logging.info("Pickled pedigree %s to file %s", pedobj.kw.get("pedname", "Unknown"), pfn)
        if pedobj.kw.get("messages") == "verbose":
            print(f"Pickled pedigree {pedobj.kw.get('pedname', 'Unknown')} to file {pfn}")
        return 1
    except Exception as e:
        logging.error("Unable to pickle pedigree %s to file %s: %s", pedobj.kw.get("pedname", "Unknown"), pfn, str(e))
        return 0
    finally:
        logging.info("Exited pickle_pedigree()")


def unpickle_pedigree(filename: str = ""):
    """
    Reads a pickled pedigree from a file and returns the unpacked pedigree object.

    Parameters
    ----------
    filename : str
        The name of the pickle file.

    Returns
    -------
    object
        An instance of a NewPedigree object on success, or 0 if the operation fails.
    """
    logging.info("Entered unpickle_pedigree()")
    try:
        if not filename:
            logging.error("No filename provided for pedigree unpickling!")
            return 0

        # Ensure the file has a `.pkl` extension
        if not filename.endswith(".pkl"):
            filename = f"{filename}.pkl"
            logging.info("No file extension provided for %s. An extension (.pkl) was added.", filename)

        # Load the pickled object
        with open(filename, "rb") as infile:
            my_pedigree = pickle.load(infile)

        logging.info("Unpickled pedigree %s from file %s", my_pedigree.kw.get("pedname", "Unknown"), filename)
        return my_pedigree

    except Exception as e:
        logging.error("Unable to unpickle pedigree from file %s: %s", filename, str(e))
        return 0
    finally:
        logging.info("Exited unpickle_pedigree()")


def summary_inbreeding(f_metadata: dict) -> str:
    """
    Returns a string representation of the data contained in the 'metadata' dictionary 
    from the output of pyp_nrm/pyp_inbreeding().

    Parameters
    ----------
    f_metadata : dict
        Dictionary of inbreeding metadata.

    Returns
    -------
    str
        A string representation of the inbreeding statistics, or '0' if an error occurs.
    """
    try:
        summary = []
        line1 = "=" * 80
        line2 = "-" * 80

        summary.append(line1)
        summary.append("Inbreeding Statistics")
        summary.append(line1)
        summary.append("All animals:")
        summary.append(line2)

        # Process all animals
        for k, v in f_metadata.get("all", {}).items():
            summary.append(f"\t{k}\t{v}")

        summary.append(line1)
        summary.append("Animals with non-zero CoI:")
        summary.append(line2)

        # Process animals with non-zero CoI
        for k, v in f_metadata.get("nonzero", {}).items():
            summary.append(f"\t{k}\t{v}")

        summary.append(line1)

        return "\n".join(summary)

    except Exception as e:
        print(f"An error occurred while generating the summary: {e}")
        return "0"


def save_ijk(pedobj, nrm_filename: str) -> int:
    """
    Saves an NRM to a file in the form "animal A" "animal B" "rAB".

    Parameters
    ----------
    pedobj : object
        The pedigree object to which the NRM is attached.
    nrm_filename : str
        The file to which the matrix should be written.

    Returns
    -------
    int
        A save status indicator (0: failed, 1: success).
    """
    try:
        if pedobj.kw.get('messages') == 'verbose':
            print(f"[INFO]: Saving A-matrix to file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")

        logging.info("Saving A-matrix to file %s", nrm_filename)

        with open(nrm_filename, "w") as file:
            for i in range(pedobj.metadata.num_records):
                for j in range(i, pedobj.metadata.num_records):
                    line = (
                        f"{pedobj.backmap[i + 1]} {pedobj.backmap[j + 1]} "
                        f"{pyp_nrm._matrix_value(pedobj.nrm.nrm, i, j)}\n"
                    )
                    file.write(line)

        if pedobj.kw.get('messages') == 'verbose':
            print(f"[INFO]: A-matrix successfully saved to file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")
        logging.info("A-matrix successfully saved to file %s", nrm_filename)
        return 1

    except Exception as e:
        logging.error("Failed to save A-matrix to file %s. Error: %s", nrm_filename, e)
        return 0


def load_from_gedcom(
    infilename: str,
    messages: str = 'verbose',
    standalone: int = 1,
    missing_sex: str = 'u',
    missing_parent: int = 0,
    missing_name: str = 'Unknown Name',
    missing_byear: str = None,
    debug: bool = False
) -> str:
    """
    Reads and parses pedigree data conforming to a subset of the GEDCOM 5.5 specification.

    Parameters
    ----------
    infilename : str
        The file to which the matrix should be written.
    messages : str, optional
        Controls output to the screen, by default 'verbose'.
    standalone : int, optional
        Uses logging if called by a NewPedigree method, by default 1.
    missing_sex : str, optional
        Value assigned to an animal with unknown sex, by default 'u'.
    missing_parent : int, optional
        Value assigned to unknown parents, by default 0.
    missing_name : str, optional
        Name assigned by default, by default 'Unknown Name'.
    missing_byear : str, optional
        Ignored. Unknown GEDCOM dates are omitted / stored as None.
        Retained only so older callers do not break.
    debug : bool, optional
        Flag turning debugging messages on or off, by default False.

    Returns
    -------
    str
        A string containing the pedigree format code. 'xxxx' if there was a problem.
    """
    known_tags = {'BIRT', 'CHIL', 'DATE', 'FAM', 'FAMC', 'FAMS', 'HUSB', 'INDI', 'NAME', 'SEX', 'WIFE'}
    indi, fam = {}, {}
    fam2husb, fam2wife, fam2chil, indi2name, indi2sex, indi2famc, indi2fams, indi2birth = {}, {}, {}, {}, {}, {}, {}, {}

    try:
        if standalone == 0:
            logging.info(f"Opening GEDCOM pedigree file {infilename}.")

        with open(infilename, 'r') as infile:
            inlines = infile.readlines()

        all_done = False
        zero_mark = 0

        while not all_done:
            next_zero = False
            current_record = []

            for l in range(zero_mark, len(inlines)):
                try:
                    line = inlines[l].strip()
                    if line.startswith("0") and l > zero_mark:
                        next_zero = True
                        zero_mark = l
                        break
                    else:
                        current_record.append(line)
                except IndexError:
                    all_done = True

            if debug:
                print("-" * 80)
                print(current_record)
                print("-" * 80)

            first_line = current_record[0].split()
            if first_line[-1] == 'INDI':
                _indi, _name, _sex, _famc, _fams, _byr = "", "", "", "", [], ""
                any_sexes, get_birth = False, False

                for line in current_record[1:]:
                    linelist = line.split()
                    if not linelist:
                        continue

                    tag = linelist[1].upper() if len(linelist) > 1 else ""
                    if tag in known_tags:
                        if debug:
                            print(f"Processing tag {tag}")

                        if tag == "NAME":
                            _name = " ".join(linelist[2:])
                        elif tag == "SEX":
                            _sex = linelist[2]
                            any_sexes = True
                        elif tag == "FAMC":
                            _famc = linelist[2][1:-1]
                        elif tag == "FAMS":
                            _fams.append(linelist[2][1:-1])
                        elif tag == "BIRT":
                            get_birth = True
                        elif tag == "DATE" and get_birth:
                            _byr = linelist[-1]
                            get_birth = False

                if not any_sexes:
                    _sex = missing_sex.upper()
                indi2sex[_indi] = _sex
                indi2birth[_indi] = _byr or None
                indi2famc[_indi] = _famc
                indi2fams[_indi] = _fams
                indi2name[_indi] = _name or missing_name

            elif first_line[-1] == 'FAM':
                _fam, _husb, _wife, _child = "", "", "", []

                for line in current_record[1:]:
                    linelist = line.split()
                    if not linelist:
                        continue

                    tag = linelist[1].upper() if len(linelist) > 1 else ""
                    if tag in known_tags:
                        if debug:
                            print(f"Processing tag {tag}")

                        if tag == "HUSB":
                            _husb = linelist[2][1:-1]
                        elif tag == "WIFE":
                            _wife = linelist[2][1:-1]
                        elif tag == "CHIL":
                            _child.append(linelist[2][1:-1])

                fam2husb[_fam] = _husb
                fam2wife[_fam] = _wife
                fam2chil[_fam] = {ch: ch for ch in _child}

        assembled = {}
        for i in indi2sex.keys():
            assembled[i] = {
                'indi': i,
                'sex': indi2sex[i],
                'birth': indi2birth[i],
                'name': indi2name[i],
                'sire': fam2husb.get(indi2famc.get(i, ""), missing_parent),
                'dam': fam2wife.get(indi2famc.get(i, ""), missing_parent),
            }

        if messages == 'verbose':
            print(f"[INFO]: Successfully imported pedigree from the GEDCOM file {infilename}.")

        logging.info(f"Successfully imported pedigree from the GEDCOM file {infilename}.")

    except Exception as e:
        if messages == 'verbose':
            print(f"[ERROR]: Unable to import pedigree from the GEDCOM file {infilename}. Error: {e}")
        logging.error(f"Unable to import pedigree from the GEDCOM file {infilename}. Error: {e}")
        return "xxxx"

    try:
        outfilename = f"{infilename}.tmp"
        pedformat = save_from_gedcom(outfilename, assembled)
    except Exception as e:
        logging.error(f"Error saving GEDCOM file: {e}")
        return "xxxx"

    return pedformat


def save_from_gedcom(outfilename: str, assembled: dict) -> str:
    """
    Takes pedigree data parsed by load_from_gedcom() and writes it to a text file 
    in an ASD format that PyPedal can easily read.

    Parameters
    ----------
    outfilename : str
        The file to which the records should be written.
    assembled : dict
        A dictionary of records read from a GEDCOM input file.

    Returns
    -------
    str
        A string containing the pedigree format code. 'xxxx' if there was a problem.
    """
    pedformat = 'xxxx'
    try:
        with open(outfilename, 'w') as ofh:
            for _i, record in assembled.items():
                pedformat = 'ASDxbu'
                birth = record['birth']
                if birth is None or str(birth).strip() == '':
                    birth_out = pyp_chronology.TEXT_MISSING_TOKEN
                else:
                    text = str(birth).strip()
                    if text.isdigit() and len(text) in (7, 8):
                        birth_out = text
                    elif text.isdigit():
                        birth_out = pyp_chronology.format_year_token(text)
                    else:
                        birth_out = text
                outstring = '{},{},{},{},{},{}\n'.format(
                    record['indi'],
                    record['sire'],
                    record['dam'],
                    record['sex'],
                    birth_out,
                    record['name'],
                )
                ofh.write(outstring)
        logging.info(f"Saved GEDCOM pedigree to the file {outfilename}!")
    except Exception as e:
        logging.error(f"Unable to save GEDCOM pedigree to the file {outfilename}! Error: {e}")
    return pedformat


def save_to_gedcom(pedobj, outfilename: str) -> int:
    """
    Writes a PyPedal NewPedigree object to a file in GEDCOM 5.5 format.

    Parameters
    ----------
    pedobj : object
        An instance of a PyPedal NewPedigree object.
    outfilename : str
        The file to which the data should be written.

    Returns
    -------
    int
        A save status indicator (1 for success, 0 for failure).
    """
    try:
        with open(outfilename, 'w') as ofh:
            # Write file header
            ofh.write("0 HEAD\n")
            ofh.write("1 SOUR PYPEDAL\n")
            ofh.write("2 VERS V2.0\n")
            ofh.write("2 CORP USDA-ARS-BA-ANRI-AIPL\n")
            ofh.write("1 DEST PYPEDAL\n")
            ofh.write(f"1 DATE {time.strftime('%d %b %Y', time.localtime())}\n")
            ofh.write(f"1 FILE {pedobj.kw['pedfile']}\n")
            ofh.write("1 GEDC\n")
            ofh.write("2 VERS 5.5\n")
            ofh.write("2 FORM Lineage-Linked\n")
            ofh.write("1 CHAR ASCII\n")

            indi = {}
            fam = {}
            par2spouses = {}

            # Sweep pedigree for family and individual mappings
            for p in pedobj.pedigree:
                if p.sireID != pedobj.kw['missing_parent'] or p.damID != pedobj.kw['missing_parent']:
                    par2spouses.setdefault(p.sireID, [])
                    par2spouses.setdefault(p.damID, [])
                    if p.sireName == pedobj.kw['missing_name'] and p.damName == pedobj.kw['missing_name']:
                        _spouses = 'F0_0'
                    elif p.sireName == pedobj.kw['missing_name']:
                        _spouses = f"F{p.damName}"
                    elif p.damName == pedobj.kw['missing_name']:
                        _spouses = f"F{p.sireName}"
                    else:
                        _spouses = f"F{p.sireName}_{p.damName}"

                    if p.sireID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.sireID]:
                            par2spouses[p.sireID].append(_spouses)
                    if p.damID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.damID]:
                            par2spouses[p.damID].append(_spouses)

            # Create INDI and FAM records
            for p in pedobj.pedigree:
                if p.sireName == pedobj.kw['missing_name'] and p.damName == pedobj.kw['missing_name']:
                    _fam = 'F0_0'
                elif p.sireName == pedobj.kw['missing_name']:
                    _fam = f"F{p.damName}"
                elif p.damName == pedobj.kw['missing_name']:
                    _fam = f"F{p.sireName}"
                else:
                    _fam = f"F{p.sireName}_{p.damName}"

                if p.animalID not in indi:
                    indi[p.animalID] = f"0 @{pedobj.namebackmap[pedobj.backmap[p.animalID]]}@ INDI\n"
                    indi[p.animalID] += f"1 SEX {p.sex.upper()}\n"

                    if 'n' in pedobj.kw['pedformat']:
                        indi[p.animalID] += f"1 NAME {p.name}\n"
                    elif 'u' in pedobj.kw['pedformat']:
                        indi[p.animalID] += f"1 NAME {p.userField}\n"

                    if 'y' in pedobj.kw['pedformat'] and p.by is not None:
                        indi[p.animalID] += "1 BIRT\n"
                        indi[p.animalID] += f"2 DATE {p.by}\n"
                    elif 'b' in pedobj.kw['pedformat'] and p.bd is not None:
                        indi[p.animalID] += "1 BIRT\n"
                        indi[p.animalID] += (
                            f"2 DATE {pyp_chronology.format_date_token(p.bd)}\n"
                        )
                    elif 'b' in pedobj.kw['pedformat'] and p.by is not None:
                        indi[p.animalID] += "1 BIRT\n"
                        indi[p.animalID] += f"2 DATE {p.by}\n"
                    if _fam != 'F0_0':
                        indi[p.animalID] += f"1 FAMC @{_fam}@\n"

                    if p.animalID in par2spouses:
                        for _p2s in par2spouses[p.animalID]:
                            if _p2s != 'F0_0':
                                indi[p.animalID] += f"1 FAMS @{_p2s}@\n"

                if _fam not in fam and _fam[3:] != '0_0':
                    fam[_fam] = f"0 @{_fam}@ FAM\n"
                if 'HUSB' not in fam.get(_fam, '') and p.sireName != pedobj.kw['missing_name']:
                    fam[_fam] += f"1 HUSB @{p.sireName}@\n"
                if 'WIFE' not in fam.get(_fam, '') and p.damName != pedobj.kw['missing_name']:
                    fam[_fam] += f"1 WIFE @{p.damName}@\n"
                fam[_fam] += f"1 CHIL @{pedobj.namebackmap[pedobj.backmap[p.animalID]]}@\n"

            # Write INDI and FAM records
            for record in indi.values():
                ofh.write(record)
            for record in fam.values():
                ofh.write(record)

            # Write footer
            ofh.write("0 TRLR\n")

        if pedobj.kw['messages'] == 'verbose':
            print(f"[INFO]: Successfully exported pedigree to the GEDCOM file {outfilename}!")
        logging.info(f"Successfully exported pedigree to the GEDCOM file {outfilename}!")
        return 1
    except Exception as e:
        if pedobj.kw['messages'] == 'verbose':
            print(f"[ERROR]: Unable to save pedigree to the GEDCOM file {outfilename}! Error: {e}")
        logging.error(f"Unable to export pedigree to the GEDCOM file {outfilename}! Error: {e}")
        return 0


def save_newanimals_to_file(animal_list, filename: str, kw: dict, new_animal_attr: dict) -> bool:
    """
    save_newanimals_to_file() takes a list of PyPedal NewAnimal objects as input and writes them to a pedigree file.

    Parameters
    ----------
    animal_list : list
        A list of PyPedal NewAnimal objects.
    filename : str
        The name of the file to which the animals should be written.
    kw : dict
        A dictionary of pedigree options.
    new_animal_attr : dict
        A mapping of pedigree format codes to attributes of NewAnimal objects.

    Returns
    -------
    bool
        True if the file was saved successfully, False otherwise.
    """
    if not animal_list:
        if kw.get('messages') == 'verbose':
            print('[WARNING]: There were no animals in the list passed to save_newanimals_to_file()!')
        logging.warning('There were no animals in the list passed to save_newanimals_to_file()!')
        return False

    try:
        with open(filename, 'w') as ofh:
            for animal in animal_list:
                output_line = []

                for pf in kw.get('pedformat', []):
                    if pf == 'a':
                        value = animal.animalID
                    elif pf == 'A':
                        value = animal.name
                    elif pf == 's':
                        value = animal.sireID if animal.sireID != kw.get('missing_parent') else kw.get('missing_parent')
                    elif pf == 'S':
                        value = animal.sireName if animal.sireID != kw.get('missing_parent') else kw.get('missing_parent')
                    elif pf == 'd':
                        value = animal.damID if animal.damID != kw.get('missing_parent') else kw.get('missing_parent')
                    elif pf == 'D':
                        value = animal.damName if animal.damID != kw.get('missing_parent') else kw.get('missing_parent')
                    else:
                        formatted = pyp_chronology.format_pedigree_field(
                            pf, animal, kw.get('pedformat', '')
                        )
                        if formatted is not None:
                            value = formatted
                        else:
                            value = getattr(animal, new_animal_attr.get(pf, ''), '')

                    output_line.append(str(value))

                sepchar = kw.get('sepchar', ',')
                ofh.write(f"{sepchar.join(output_line)}\n")
        return True
    except Exception as e:
        logging.error('Error writing to file %s: %s', filename, e)
        return False
