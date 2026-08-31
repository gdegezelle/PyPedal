#!/usr/bin/env python3

###############################################################################
# NAME: pyp_snp.py
# VERSION: originally 2.0.0; current package version is PyPedal.__version__
# AUTHOR: John B. Cole, PhD
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
# FUNCTIONS:
#   read_agil_chromosome_data()
#   read_agil_genotypes_txt()
#   read_agil_true_frequency()
#   form_p_matrix_from_snp()
#   form_m_matrix_from_snp()
#   form_grm_from_snp()
#   compute_genomic_inbreeding_from_grm()
#   compute_genomic_homozygosity_from_snp()
#   renumber_snp_ids()
#   generate_random_genotype()
###############################################################################

## @package pyp_snp
# pyp_snp contains several procedures for working with single nucleotide
# polymorphism (SNP) genotype data.

import logging

import numpy as np
import pandas as pd
import random

from . import pyp_errors, pyp_validate
from .pyp_errors import PyPedalUsageError


logger = logging.getLogger(__name__)

def read_agil_chromosome_data(filename='chromosome.data'):
    """
    Load SNP marker information from the chromosome.data file used by AGIL and CDCB.
    Reads the first 5 columns: SNP name, chromosome number, within-chromosome marker
    number, overall marker number, and location in base pairs.
    """
    try:
        df = pd.read_csv(
            filename,
            sep=r'\s+',
            usecols=[0, 1, 2, 3, 4],
            names=['snp_name', 'chrome', 'within', 'overall', 'location']
        )
        print(f'[INFO]: Successfully read {len(df)} SNP records from {filename}.')
        logger.info(f'Successfully read {len(df)} SNP records from {filename}.')
        return df
    except Exception as e:
        logger.error(f'Could not access the file {filename}: {e}')
        print(f'[ERROR]: Could not access the file {filename}.')
        return False


def read_agil_genotypes_txt(filename='genotypes.txt'):
    """
    Load SNP genotypes from the genotypes.txt file used by AGIL and CDCB.
    """
    try:
        df = pd.read_csv(
            filename,
            sep=r'\s+',
            dtype={'genotype': str},
            names=['animalID', 'chip_type', 'n_snps', 'genotype']
        )
        print(f'[INFO]: Successfully read {len(df)} genotype records from {filename}.')
        logger.info(f'Successfully read {len(df)} genotype records from {filename}.')
        return df
    except Exception as e:
        logger.error(f'Could not access the file {filename}: {e}')
        print(f'[ERROR]: Could not access the file {filename}.')
        return False


def read_agil_true_frequency(filename='true.frequency'):
    """
    Load SNP frequency data from the true.frequency file used by AGIL and CDCB.
    """
    try:
        df = pd.read_csv(
            filename,
            sep=r'\s+',
            names=['snp_name', 'overall', 'frequency']
        )
        print(f'[INFO]: Successfully read {len(df)} SNP frequencies from {filename}.')
        logger.info(f'Successfully read {len(df)} SNP frequencies from {filename}.')
        return df
    except Exception as e:
        logger.error(f'Could not access the file {filename}: {e}')
        print(f'[ERROR]: Could not access the file {filename}.')
        return False


def validate_base_frequencies(base_frequencies, n_loci, caller):
    """
    Check and align caller-supplied base-population allele frequencies.

    VanRaden (2008) p.4416 is explicit that ``P`` "should be from the unselected
    base population rather than from those that appear after selection or
    inbreeding", and p.4419 that "frequency estimation can bias the genomic
    inbreeding coefficients". So supplying them is the source-faithful path --
    but source-faithful is not the same as blindly trusted, and a frequency
    vector silently one element short, or accidentally in the wrong order, would
    produce a plausible and wrong G.

    Accepted forms:

    * a sequence of length ``n_loci``, positional;
    * a mapping from locus index (0-based) to frequency, **aligned by key**;
    * a pandas Series or DataFrame from :func:`read_agil_true_frequency`, whose
      ``frequency`` column is used and whose index is treated as locus identity.

    Returns a plain NumPy array in locus order. Raises ``PyPedalUsageError`` on
    anything that would leave the alignment or the values in doubt; after
    validation the values are used verbatim -- no renormalisation, no clipping.
    """
    if isinstance(base_frequencies, pd.DataFrame):
        if 'frequency' not in base_frequencies.columns:
            raise PyPedalUsageError(
                '%s: a DataFrame of base frequencies must have a "frequency" '
                'column; got %r.' % (caller, list(base_frequencies.columns)))
        base_frequencies = base_frequencies['frequency']

    if isinstance(base_frequencies, pd.Series):
        base_frequencies = {int(k): float(v)
                            for k, v in base_frequencies.items()}

    if isinstance(base_frequencies, dict):
        missing = [i for i in range(n_loci) if i not in base_frequencies]
        if missing:
            raise PyPedalUsageError(
                '%s: base_frequencies is keyed by locus, but %d of %d loci have '
                'no entry (first missing: %r). Alignment must be explicit; '
                'PyPedal will not fill gaps.'
                % (caller, len(missing), n_loci, missing[:5]))
        extra = [k for k in base_frequencies if k not in range(n_loci)]
        if extra:
            raise PyPedalUsageError(
                '%s: base_frequencies contains %d key(s) that are not loci of '
                'this genotype set (first: %r). That means the frequencies '
                'describe a different marker panel.'
                % (caller, len(extra), extra[:5]))
        values = np.array([base_frequencies[i] for i in range(n_loci)],
                          dtype=float)
    else:
        values = np.asarray(base_frequencies, dtype=float).ravel()
        if values.shape[0] != n_loci:
            raise PyPedalUsageError(
                '%s: expected exactly one base frequency per locus -- %d loci, '
                '%d frequencies. A positional vector of the wrong length cannot '
                'be aligned, and guessing would silently mis-centre every '
                'marker.' % (caller, n_loci, values.shape[0]))

    if not np.all(np.isfinite(values)):
        raise PyPedalUsageError(
            '%s: base frequencies must all be finite; %d are NaN or infinite.'
            % (caller, int((~np.isfinite(values)).sum())))
    if np.any(values < 0.0) or np.any(values > 1.0):
        raise PyPedalUsageError(
            '%s: base frequencies are allele frequencies and must lie in '
            '[0, 1]; the range supplied is [%r, %r].'
            % (caller, float(values.min()), float(values.max())))
    return values


def form_p_matrix_from_snp(pedobj, base_frequencies=None, debug=False):
    """
    Allele frequencies for the marker panel, one per locus.

    ``base_frequencies`` supplies them directly, which is what VanRaden (2008)
    p.4416 calls for: frequencies "from the unselected base population rather
    than from those that appear after selection or inbreeding". They are
    validated (see :func:`validate_base_frequencies`) and then used verbatim.

    Omitting it falls back to ESTIMATING them from the genotyped sample. That is
    the paper's own "simple estimate" and is retained for convenience and
    backward compatibility, but it is an **estimate of** the quantity the
    equation names, not that quantity: p.4419 states that "frequency estimation
    can bias the genomic inbreeding coefficients". The fallback logs that it was
    taken, so the choice appears in the run record instead of being implicit.
    """
    if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
        logger.error('No SNP data available in the pedigree.')
        print('[ERROR]: No SNP data available in the pedigree.')
        return False

    n_loci = len(pedobj.snp.iloc[0, 3])
    if base_frequencies is not None:
        P = validate_base_frequencies(base_frequencies, n_loci,
                                      'form_p_matrix_from_snp')
    else:
        logger.info(
            'form_p_matrix_from_snp(): no base_frequencies supplied, so allele '
            'frequencies are being ESTIMATED from the genotyped sample. '
            'VanRaden (2008) p.4416 calls for frequencies from the unselected '
            'base population, and p.4419 notes that estimation can bias genomic '
            'inbreeding coefficients.')
        P = np.zeros(n_loci)
        for a in range(len(pedobj.snp)):
            for s in range(n_loci):
                P[s] += float(pedobj.snp.iloc[a, 3][s])
        P /= (2.0 * len(pedobj.snp))

    if debug:
        print('P:', P)
    return P


def form_m_matrix_from_snp(pedobj, scale_m=True, base_frequencies=None, debug=False):
    """
    Form the centred marker matrix used to build the genomic relationship matrix.

    Genotypes are stored as allele counts (0, 1, 2). VanRaden (2008) p.4416
    codes ``M`` as -1/0/1 and centres it with ``P``, column *i* being
    ``2(p_i - 0.5)``, so that ``Z = M - P = counts - 2p``.

    That subtraction is split across two functions here: this one subtracts
    ``p`` and :func:`form_grm_from_snp` subtracts it again, which together give
    the paper's ``Z``. **The value returned here is therefore an intermediate,
    not the paper's M and not its Z.** It is documented rather than
    restructured, because the composed result is correct and changing the
    decomposition would change this function's public return value.

    Parameters
    ----------
    scale_m : bool
        Must be True. ``scale_m=False`` skips the centring, leaving
        ``Z = counts - p``, which is **not** VanRaden's Z; under it ``G_jj - 1``
        is not a genomic inbreeding coefficient. It was an option in name only
        and is now refused -- a trap rather than a choice.
    """
    if not scale_m:
        raise PyPedalUsageError(
            'form_m_matrix_from_snp: scale_m=False is not a supported option. '
            'VanRaden (2008) p.4416 centres the marker matrix by P, column i = '
            '2(p_i - 0.5); skipping it leaves Z = counts - p instead of '
            'counts - 2p, so the resulting matrix is not a genomic relationship '
            'matrix and its diagonal minus one is not a genomic inbreeding '
            'coefficient. The centring is not optional.')
    if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
        logger.error('No SNP data available in the pedigree.')
        print('[ERROR]: No SNP data available in the pedigree.')
        return False

    P = form_p_matrix_from_snp(pedobj, base_frequencies=base_frequencies, debug=debug)
    M = np.zeros((len(pedobj.snp), len(pedobj.snp.iloc[0, 3])))

    for a in range(len(pedobj.snp)):
        for s in range(len(pedobj.snp.iloc[0, 3])):
            M[a, s] = pedobj.snp.iloc[a, 3][s]
            if scale_m:
                M[a, s] -= P[s]

    if debug:
        print('M:', M)
    return M


def form_grm_from_snp(pedobj, scale_m=True, method=1, base_frequencies=None, debug=False):
    """
    Genomic relationship matrix by VanRaden (2008) Method 1, article p.4416::

        G = ZZ' / (2 * sum(p_i * (1 - p_i)))

    with ``Z = M - P``. "Division by 2*sum(p_i(1-p_i)) scales G to be analogous
    to the numerator relationship matrix A."

    The genomic inbreeding coefficient of individual *j* is ``G_jj - 1``. Note
    that it is bounded below by -1 but has **no generic finite upper bound**:
    an individual homozygous for rare alleles can have an arbitrarily large
    value, so the pedigree-inbreeding validators must not be applied to it.

    Parameters
    ----------
    scale_m : bool
        Must be True; see :func:`form_m_matrix_from_snp`.
    method : int
        Only Method 1 is implemented. Methods 2 and 3 of the paper are not.

    Raises
    ------
    PyPedalUsageError
        If ``scale_m`` is False, or ``method`` is not 1.
    PyPedalValidationError
        If every locus is monomorphic, so that ``2*sum(p(1-p))`` is zero. That
        is not a convention: Method 1 is mathematically undefined there, and the
        paper is silent because there is nothing to state. The guard is derived
        algebraically, not quoted.
    """
    if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
        logger.error('No SNP data available in the pedigree.')
        print('[ERROR]: No SNP data available in the pedigree.')
        return False

    P = form_p_matrix_from_snp(pedobj, base_frequencies=base_frequencies, debug=debug)
    M = form_m_matrix_from_snp(pedobj, scale_m=scale_m, base_frequencies=base_frequencies, debug=debug)

    if method != 1:
        raise PyPedalUsageError(
            'form_grm_from_snp: only VanRaden (2008) Method 1 is implemented, '
            'not method=%r. Methods 2 and 3 of that paper are different '
            'estimators, not variations, and returning a Method 1 matrix for a '
            'Method 2 request would be worse than refusing.' % (method,))

    sum_freq = float(np.sum(P * (1 - P)))
    if sum_freq <= 0.0:
        raise pyp_validate.PyPedalValidationError(
            'form_grm_from_snp: 2*sum(p(1-p)) is zero because every locus is '
            'monomorphic, so VanRaden Method 1 is undefined for this input. '
            'This is algebra, not a convention -- there is no genomic '
            'relationship to measure when no marker varies.')

    Z = M - P
    G = Z @ Z.T / (2 * sum_freq)
    if debug:
        print('G:', G)
    return G


def compute_genomic_inbreeding_from_grm(pedobj, grm=None, base_frequencies=None,
                                        store=True):
    """
    Genomic inbreeding coefficients from a genomic relationship matrix.

    VanRaden (2008) article p.4416, verbatim: "The genomic inbreeding
    coefficient for individual j is simply G_jj - 1."

    Restored so PyPedal can both parse the ``G`` pedformat column and compute it.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree carrying SNP genotypes.
    grm : optional
        A precomputed genomic relationship matrix. Built with
        :func:`form_grm_from_snp` when omitted.
    base_frequencies : optional
        Base-population allele frequencies, passed through when ``grm`` is
        built here. See :func:`form_p_matrix_from_snp`.
    store : bool
        Write each value onto ``NewAnimal.genomic_inbreeding`` and set
        ``kw['g_computed']``.

    Returns
    -------
    dict
        ``{animalID: F_g}``.

    Notes
    -----
    ``F_g >= -1``, attainably, and has **no generic finite upper bound**. The
    pedigree-inbreeding validators must not be applied to it; see
    :func:`pyp_validate.check_genomic_inbreeding`.
    """
    if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
        raise pyp_errors.PyPedalInputError(
            'compute_genomic_inbreeding_from_grm: the pedigree carries no SNP '
            'genotypes, so there is no genomic relationship matrix to take a '
            'diagonal from.')

    if grm is None:
        grm = form_grm_from_snp(pedobj, base_frequencies=base_frequencies)

    diagonal = np.asarray(grm).diagonal() if hasattr(grm, 'diagonal') \
        else np.diagonal(np.asarray(grm))
    if len(diagonal) != len(pedobj.snp):
        raise PyPedalUsageError(
            'compute_genomic_inbreeding_from_grm: the matrix is %d x %d but %d '
            'animals are genotyped. A relationship matrix built from a '
            'different set of animals cannot be indexed by this one.'
            % (len(diagonal), len(diagonal), len(pedobj.snp)))

    result = {}
    for row in range(len(pedobj.snp)):
        animal_id = pedobj.snp.iloc[row, 0]
        result[animal_id] = float(diagonal[row]) - 1.0

    pyp_validate.check_genomic_inbreeding(
        pedobj, result, 'compute_genomic_inbreeding_from_grm')

    if store:
        by_id = {str(k): v for k, v in result.items()}
        for animal in pedobj.pedigree:
            for key in (str(animal.animalID), str(animal.originalID)):
                if key in by_id:
                    animal.genomic_inbreeding = by_id[key]
                    break
        pedobj.kw['g_computed'] = True

    return result


def compute_genomic_homozygosity_from_snp(pedobj, missing_code=None, store=True):
    """
    Per-individual genomic homozygosity: the proportion of typed loci at which
    an individual carries two copies of the same allele.

    Independent of VanRaden -- this is a direct count, not an estimator, and its
    contract is genuinely ``[0, 1]``. That is worth stating beside
    :func:`compute_genomic_inbreeding_from_grm`, whose value shares the word
    "genomic" and has no upper bound at all.

    Genotypes are allele dosages, so 0 and 2 are the homozygotes and 1 the
    heterozygote.

    Parameters
    ----------
    missing_code : str, optional
        A character marking an untyped locus, excluded from both numerator and
        denominator. The loader currently admits only ``0``, ``1`` and ``2``, so
        this matters for genotypes assembled by other means.
    store : bool
        Write each value onto ``NewAnimal.homozygosity``.

    Returns
    -------
    dict
        ``{animalID: homozygosity}``. An individual with **no typed loci** gets
        ``kw['missing_homozygosity']``, not 0.0: zero typed loci is a missing
        result, and reporting it as complete homozygosity-of-zero would be a
        measurement that was never made.
    """
    if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
        raise pyp_errors.PyPedalInputError(
            'compute_genomic_homozygosity_from_snp: the pedigree carries no SNP '
            'genotypes.')

    missing_value = pedobj.kw.get('missing_homozygosity', -999.0)
    result = {}
    for row in range(len(pedobj.snp)):
        animal_id = pedobj.snp.iloc[row, 0]
        genotype = pedobj.snp.iloc[row, 3]
        typed = [c for c in genotype if missing_code is None or c != missing_code]
        if not typed:
            logger.warning(
                'compute_genomic_homozygosity_from_snp(): animal %s has no '
                'typed loci, so its homozygosity is undefined and is reported '
                'as the missing value rather than as 0.0.', animal_id)
            result[animal_id] = missing_value
            continue
        homozygous = sum(1 for c in typed if c in ('0', '2'))
        result[animal_id] = homozygous / float(len(typed))

    for animal_id, value in result.items():
        if value == missing_value:
            continue
        if not 0.0 <= value <= 1.0:
            raise pyp_validate.PyPedalValidationError(
                'compute_genomic_homozygosity_from_snp: homozygosity for animal '
                '%s is %r, outside [0, 1]. It is a proportion of typed loci.'
                % (animal_id, value))

    if store:
        by_id = {str(k): v for k, v in result.items()}
        for animal in pedobj.pedigree:
            for key in (str(animal.animalID), str(animal.originalID)):
                if key in by_id:
                    animal.homozygosity = by_id[key]
                    break

    return result


# def renumber_snp_ids(pedobj):
#     """
#     Renumber SNP animal IDs in the SNP dataframe to match pedigree IDs.
#     """
#     if pedobj.snp is False or pedobj.snp is None or getattr(pedobj.snp, "empty", True):
#         logger.error('No SNP data available to renumber IDs.')
#         print('[ERROR]: No SNP data available to renumber IDs.')
#         return

#     for p in pedobj.pedigree:
#         pedobj.snp.loc[:, 'animalID'] = pedobj.snp['animalID'].replace(
#             to_replace=p.originalID, value=p.animalID)


def renumber_snp_ids(pedobj):
    """
    Renumbers animal IDs in the SNP dataframe attached to a pedigree object.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object with SNP data and a pedigree.

    Returns
    -------
    None
    """
    # Check if SNP data is available and valid
    if pedobj.snp is False or pedobj.snp is None:
        logger.error(
            "pyp_snp/renumber_snp_ids(): SNP data is not available. No IDs need to be renumbered."
        )
        if pedobj.kw.get('debug_messages', False):
            print(
                "[ERROR]: pyp_snp/renumber_snp_ids(): SNP data is not available. No IDs need to be renumbered."
            )
        return

    if isinstance(pedobj.snp, pd.DataFrame) and not pedobj.snp.empty:
        if pedobj.kw.get('debug_messages', False):
            print(
                "[INFO]: pyp_snp/renumber_snp_ids(): Renumbering animal IDs in the SNP dataframe."
            )
        logger.info("pyp_snp/renumber_snp_ids(): Renumbering animal IDs in the SNP dataframe")

        # Replace `originalID` with `animalID` in the SNP dataframe.
        #
        # One mapped replace on the column, not a loop of chained-assignment
        # `inplace=True` calls. Chained assignment through a column selection
        # raises a FutureWarning in pandas 2.2 and "will never work" in 3.0 --
        # the intermediate object is always a copy -- so the old form was a
        # silent no-op waiting to happen on the dependency floor this project is
        # moving to. Keys are stringified because the genotype file's animalID
        # column is read as text while originalID may be either.
        mapping = {str(animal.originalID): animal.animalID
                   for animal in pedobj.pedigree}
        pedobj.snp['animalID'] = (
            pedobj.snp['animalID'].astype(str).map(mapping).fillna(
                pedobj.snp['animalID'])
        )
    else:
        logger.error(
            "pyp_snp/renumber_snp_ids(): Invalid or empty SNP dataframe. No renumbering performed."
        )
        if pedobj.kw.get('debug_messages', False):
            print(
                "[ERROR]: pyp_snp/renumber_snp_ids(): Invalid or empty SNP dataframe. No renumbering performed."
            )




def load_snp_file(pedobj, snpfile=None, sepchar=None):
    """
    Read a SNP genotype file and attach it to ``pedobj.snp``.

    Until this loader ran, ``kw['snpfile']`` was stored but never read.
    ``self.snp = False`` was the only assignment, so every entry point that
    opens with the ``pedobj.snp is False`` guard always took the failure path and
    ``form_grm_from_snp`` was unreachable in production.

    File format -- whitespace-separated, one line per animal::

        animalID  chip_type  n_snps  genotype

    ``genotype`` is a string of per-locus dosages, one character each, counting
    copies of the reference allele: ``0``, ``1`` or ``2``. This is the layout
    ``read_agil_genotypes_txt`` already expects.

    Parameters
    ----------
    pedobj : NewPedigree
        Loaded pedigree. Animals are matched on ``originalID``.
    snpfile : str, optional
        Path to read. Defaults to ``pedobj.kw['snpfile']``.
    sepchar : str, optional
        Ignored; the format is whitespace-separated. Accepted so callers can
        pass ``kw['snp_sepchar']`` without special-casing.

    Returns
    -------
    pandas.DataFrame
        The genotypes, also attached to ``pedobj.snp``.

    Raises
    ------
    PyPedalConfigurationError
        No SNP file was configured.
    PyPedalInputError
        The file is unreadable, empty, or its genotypes are inconsistent.
    """
    if snpfile is None:
        snpfile = pedobj.kw.get('snpfile', False)
    if not snpfile:
        raise pyp_errors.PyPedalConfigurationError(
            "No SNP file configured. Set kw['snpfile'] or pass snpfile=.")

    try:
        snp = pd.read_csv(
            snpfile,
            sep=r'\s+',
            dtype={'animalID': str, 'genotype': str},
            names=['animalID', 'chip_type', 'n_snps', 'genotype'],
            comment='#',
        )
    except OSError as exc:
        raise pyp_errors.PyPedalInputError(
            'Could not read the SNP file %r: %s' % (snpfile, exc)) from exc
    except Exception as exc:
        raise pyp_errors.PyPedalInputError(
            'Could not parse the SNP file %r: %s. Expected whitespace-separated '
            'columns: animalID chip_type n_snps genotype.'
            % (snpfile, exc)) from exc

    if snp.empty:
        raise pyp_errors.PyPedalInputError(
            'The SNP file %r contains no genotype records.' % snpfile)

    _validate_genotypes(snp, snpfile)
    _warn_about_unmatched_animals(pedobj, snp, snpfile)

    pedobj.snp = snp
    logger.info('Read %d genotype records of %d loci from %s',
                 len(snp), len(snp['genotype'].iloc[0]), snpfile)
    return snp


def _validate_genotypes(snp, snpfile):
    """
    Fail loudly on genotypes that would silently corrupt a GRM.

    Every check here is about the *data*, not about any estimator, so it applies
    regardless of what is later computed from the genotypes.
    """
    missing = snp['genotype'].isna()
    if bool(missing.any()):
        raise pyp_errors.PyPedalInputError(
            'The SNP file %r has %d record(s) with no genotype column.'
            % (snpfile, int(missing.sum())))

    lengths = snp['genotype'].str.len()
    if lengths.nunique() != 1:
        counts = lengths.value_counts().to_dict()
        raise pyp_errors.PyPedalInputError(
            'Genotypes in %r have inconsistent lengths %r. Every animal must be '
            'typed at the same loci, in the same order, or the columns of the '
            'genotype matrix do not correspond to the same markers.'
            % (snpfile, counts))

    bad = snp['genotype'].str.contains(r'[^012]', regex=True, na=False)
    if bool(bad.any()):
        offenders = snp.loc[bad, 'animalID'].head(5).tolist()
        raise pyp_errors.PyPedalInputError(
            'Genotypes in %r contain characters other than 0, 1 and 2 -- first '
            'offending animal(s): %r. Dosages count copies of the reference '
            'allele and must be 0, 1 or 2.' % (snpfile, offenders))

    duplicated = snp['animalID'].duplicated()
    if bool(duplicated.any()):
        offenders = snp.loc[duplicated, 'animalID'].head(5).tolist()
        raise pyp_errors.PyPedalInputError(
            'The SNP file %r has more than one genotype for animal(s) %r. '
            'Which one to use is undefined.' % (snpfile, offenders))


def _warn_about_unmatched_animals(pedobj, snp, snpfile):
    """
    Report animals present in one source and not the other.

    Deliberately a warning, not an error: genotyping a subset of a pedigree is
    normal practice. What is not acceptable is doing it silently, which is how
    a GRM ends up describing different animals than the caller believes.
    """
    pedigree_ids = {str(a.originalID) for a in pedobj.pedigree}
    snp_ids = {str(i) for i in snp['animalID']}
    only_snp = snp_ids - pedigree_ids
    only_ped = pedigree_ids - snp_ids
    if only_snp:
        logger.warning(
            '%d animal(s) in %s are not in the pedigree, e.g. %r',
            len(only_snp), snpfile, sorted(only_snp)[:5])
    if only_ped:
        logger.warning(
            '%d animal(s) in the pedigree have no genotype in %s, e.g. %r',
            len(only_ped), snpfile, sorted(only_ped)[:5])


def generate_random_genotype(string_length):
    """
    Generate a random string of genotypes (0, 1, 2).
    """
    return ''.join(random.choices(['0', '1', '2'], k=string_length))
