#!/usr/bin/env python3

###############################################################################
# NAME: pyp_metrics.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################

"""
pyp_metrics.py
Metrics calculations for PyPedal.
"""

import copy
import logging

import numbers
import operator
import random
import warnings

from typing import Iterable, Optional, Dict, Tuple, List, Union

import numpy
import numpy as np

from . import pyp_chronology, pyp_io, pyp_network, pyp_nrm, pyp_utils, pyp_validate
from .pyp_errors import (
    PyPedalError, PyPedalInternalError,
    PyPedalPedigreeStructureError, PyPedalUsageError, PyPedalValidationError)


logger = logging.getLogger(__name__)

def _current_animal_id(pedobj, animal_id):
    """Return the live (possibly renumbered) ID for an original or current ID."""
    mapped = pedobj.idmap.get(animal_id)
    if mapped is not None:
        return mapped
    return animal_id


def _most_recent_generation(gens, routine):
    """
    Return the generation ID with the largest NUMERIC value.

    ``NewAnimal.gen`` holds whatever the pedigree file's ``g`` column contained,
    as a string, or the ``missing_gen`` default. The Boichard routines used to
    pick their default reference population with ``sorted(gens, reverse=True)[0]``,
    which orders strings lexicographically: ``'9' > '12'``. On any pedigree with
    ten or more generations that silently selected generation 9.

    That contradicts the documented contract. ``pyp-metrics.tex:14-38`` states,
    for all three Boichard routines, that "the most recent generation -- the
    generation with the largest generation ID -- will be used as the reference
    population".

    The original value is returned, not the parsed number, so downstream
    comparisons against ``individual.gen`` keep working unchanged.

    Non-numeric generation IDs raise rather than falling back to some other
    ordering. "Largest generation ID" has no meaning for labels that are not
    ordered quantities, and quietly picking a lexicographic or insertion order
    would be inventing a contract PyPedal has never documented. No pedigree
    shipped with PyPedal uses non-numeric generations.
    """
    parsed = []
    unparseable = []
    for value in gens:
        try:
            parsed.append((float(value), value))
        except (TypeError, ValueError):
            unparseable.append(value)
    if unparseable:
        raise PyPedalError(
            '%s: the pedigree contains generation IDs that are not numeric '
            '(%r), so "the generation with the largest generation ID" is not '
            'defined. PyPedal will not guess an ordering. Supply numeric '
            'generations, or select one explicitly with gen=.'
            % (routine, sorted(str(v) for v in unparseable)))
    return max(parsed, key=lambda pair: pair[0])[1]


def _current_animal_index(pedobj, animal_id):
    """Return the 0-based pedigree index for an original or current ID."""
    return int(_current_animal_id(pedobj, animal_id)) - 1


def min_max_f(
    pedobj,
    a: Optional[np.ndarray] = None,
    n: int = 10,
    forma: str = 'dense'
) -> Union[Tuple[List[Tuple[str, float]], List[Tuple[str, float]]], bool]:
    """
    Given a pedigree or relationship matrix, return a list of individuals
    with the n largest and n smallest coefficients of inbreeding.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object with attributes `pedigree` and `kw`.
    a : np.ndarray, optional
        A numerator relationship matrix. If None, it will be generated.
    n : int, optional
        Number of coefficients to return (e.g., 10 smallest/largest). Default is 10.
    forma : str, optional
        Specify whether dense or sparse matrices should be used. Default is 'dense'.

    Returns
    -------
    tuple
        Lists of the n largest and the n smallest CoI in the pedigree, or False on failure.
    """
    logger.info('Entered min_max_f()')

    # Validate 'forma'
    if forma not in ['dense', 'sparse']:
        logger.warning(f"Invalid 'forma' value: {forma}. Defaulting to 'dense'.")
        forma = 'dense'

    try:
        # Generate or validate relationship matrix `a`
        if not pedobj.kw['form_nrm'] and a is None:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw, method=forma)
            if pedobj.kw.get('debug_messages'):
                logger.debug("Matrix `a` created in min_max_f():")
                logger.debug(a)

        # Calculate individual coefficients of inbreeding
        individual_coi = fast_a_coefficients(pedobj, a=a)
        if not individual_coi:
            logger.error("Failed to compute individual CoI in min_max_f()")
            return False
        
        # Convert dictionary to sorted list
        # Cast coefficients to native Python floats
        sorted_coi = sorted(
            [(str(k), float(v)) for k, v in individual_coi.items()],
            key=lambda x: x[1]
        )

        # Adjust `n` if there are not enough individuals
        total_individuals = len(sorted_coi)
        if n > total_individuals:
            old_n = n
            n = max(1, total_individuals // 2)  # Ensure at least one result
            logger.info(
                f"Requested {old_n} high/low CoI values, but only {total_individuals} available. "
                f"Adjusted `n` to {n}."
            )

        # Extract the n smallest and n largest CoI values
        low_coi = sorted_coi[:n]
        high_coi = sorted_coi[-n:]

        if pedobj.kw.get('debug_messages'):
            logger.debug(f"Low CoI: {low_coi}")
            logger.debug(f"High CoI: {high_coi}")

        logger.info('Exited min_max_f()')
        return high_coi, low_coi
    
    except Exception as e:
        logger.error(f"Error in min_max_f: {e}")
        logger.info('Exited min_max_f() with failure')
        return False


# ---------------------------------------------------------------------------
# Lacy (1989) half-founder handling
#
# Lacy RC (1989) "Analysis of Founder Representation in Pedigrees", Zoo Biol
# 8:111-123, article p.113, defines what to do with an animal that has one
# known parent: "the unknown parent is considered a founder of the captive
# population". One distinct phantom founder per unknown parent, confirmed by
# the okapi (p.117) and Goeldi's monkey (p.118) analyses, which show them as
# separate founder sources rather than one pooled bucket.
#
# PyPedal 4 supports only that phantom treatment. Historical strict/absorb
# (half=) partitions do not yield a founder probability vector when
# half-founders are present, so public APIs refuse them.
# ---------------------------------------------------------------------------

#: Sentinel for "argument not supplied". Distinguishing an omitted ``half``
#: from an explicit ``half=False`` keeps conflicting ``mode`` + ``half``
#: calls detectable rather than guessed.
_UNSET = object()

#: Supported public Lacy mode. Lacy (1989) p.113: one phantom founder per
#: unknown parent; the half-founder itself remains a descendant.
LACY_MODES = ('phantom',)
LACY_DEFAULT_MODE = 'phantom'

# Founder-row contribution batches. Internal only; not a public option.
_LACY_SOURCE_BATCH = 256

_LACY_SUPPORTED_MESSAGE = (
    "PyPedal 4 supports Lacy effective-founder calculations only with "
    "mode='phantom'. The historical strict/absorb (half=...) variants do not "
    "form valid Lacy founder probability partitions. Use the default "
    "mode='phantom'."
)


def _resolve_lacy_mode(mode, half, routine):
    """
    Resolve ``mode`` / ``half`` to the supported phantom treatment, or refuse.

    Public callers may omit both arguments or pass ``mode='phantom'``. Any
    historical ``half=`` value, ``mode='strict'``, ``mode='absorb'``, an
    unknown mode, or supplying both arguments raises ``PyPedalUsageError``.
    """
    if mode is not _UNSET and half is not _UNSET:
        raise PyPedalUsageError(
            '%s: pass mode= or half=, not both. %s'
            % (routine, _LACY_SUPPORTED_MESSAGE))

    if half is not _UNSET:
        raise PyPedalUsageError(
            '%s: half=%r is not supported. %s'
            % (routine, half, _LACY_SUPPORTED_MESSAGE))

    if mode is _UNSET:
        return LACY_DEFAULT_MODE

    if mode == LACY_DEFAULT_MODE:
        return LACY_DEFAULT_MODE

    raise PyPedalUsageError(
        '%s: mode=%r is not supported. %s'
        % (routine, mode, _LACY_SUPPORTED_MESSAGE))


def lacy_phantom_slots(pedobj):
    """
    THE SHARED PHANTOM ENUMERATION. One entry per unknown parent of an animal
    that has at least one known parent -- i.e. per phantom founder Lacy p.113
    calls for.

    Returns ``[(animal_id, 'sire'|'dam'), ...]``, deterministic in pedigree
    order. Both public entry points consume this, so they cannot disagree about
    what the founder set is; they previously disagreed about far more than that.
    """
    missing = str(pedobj.kw['missing_parent'])
    slots = []
    for animal in pedobj.pedigree:
        sire_missing = str(animal.sireID) == missing
        dam_missing = str(animal.damID) == missing
        if sire_missing and dam_missing:
            continue                       # a real founder, not a half-founder
        if sire_missing:
            slots.append((int(animal.animalID), 'sire'))
        if dam_missing:
            slots.append((int(animal.animalID), 'dam'))
    return slots


class _LacyRecord:
    """
    Minimal animal record for the phantom-extended pedigree.

    ``fast_a_matrix`` reads only ``animalID``, ``sireID`` and ``damID`` and
    indexes by position, so a full ``NewAnimal`` is unnecessary here -- and
    building one would mean inventing values for fields the NRM never touches.
    """

    __slots__ = ('animalID', 'sireID', 'damID')

    def __init__(self, animal_id, sire_id, dam_id):
        self.animalID = animal_id
        self.sireID = sire_id
        self.damID = dam_id


def _lacy_extended_pedigree(pedobj):
    """
    Build the pedigree Lacy p.113 describes: every unknown parent of a
    non-founder becomes a real founder with a record of its own.

    Phantoms take IDs ``1..k`` and the original animals shift up by ``k``, which
    preserves topological order because the originals were already sorted.

    Returns ``(records, k, slots)``.

    WHY EXTEND RATHER THAN DERIVE. A phantom's relationship to a descendant is
    NOT half its half-founder's: ``A(half_founder, d)`` also counts paths that
    reach ``d`` from the half-founder's KNOWN parent, which the phantom has no
    part in. On mrode.ped the true ratio is 0.4, not 0.5. The additive
    relationships therefore have to be solved on the completed pedigree.
    """
    missing = int(pedobj.kw['missing_parent'])
    slots = lacy_phantom_slots(pedobj)
    k = len(slots)
    fill = {slot: index + 1 for index, slot in enumerate(slots)}

    records = [_LacyRecord(index + 1, missing, missing) for index in range(k)]
    for animal in pedobj.pedigree:
        animal_id = int(animal.animalID)
        sire_id = int(animal.sireID)
        dam_id = int(animal.damID)
        sire_id = sire_id + k if sire_id != missing else fill.get((animal_id, 'sire'), missing)
        dam_id = dam_id + k if dam_id != missing else fill.get((animal_id, 'dam'), missing)
        records.append(_LacyRecord(animal_id + k, sire_id, dam_id))
    return records, k, slots


def _lacy_partition(pedobj, mode):
    """
    Split the pedigree into (founder ids, phantom slots, descendant ids).

    Public Lacy APIs only pass ``mode='phantom'``. The ``strict`` / ``absorb``
    branches remain as unused internal layout of those historical partitions;
    they are not a supported product domain and are not claimed to match
    each other or Lacy (1989).
    """
    missing = str(pedobj.kw['missing_parent'])
    full, halves, rest = [], [], []
    for animal in pedobj.pedigree:
        animal_id = int(animal.animalID)
        sire_missing = str(animal.sireID) == missing
        dam_missing = str(animal.damID) == missing
        if sire_missing and dam_missing:
            full.append(animal_id)
        elif sire_missing or dam_missing:
            halves.append(animal_id)
        else:
            rest.append(animal_id)

    if mode == 'absorb':
        return sorted(full + halves), [], sorted(rest)
    descendants = sorted(halves + rest)
    if mode == 'strict':
        return full, [], descendants
    return full, lacy_phantom_slots(pedobj), descendants


def _lacy_source_contributions(sire_i, dam_i, source_idx, desc_idx):
    """
    Sum of A(source, descendant) over descendants for each founder source.

    This is the founder row of the tabular NRM: c[source] = 1 and
    c[i] = 0.5 c[sire] + 0.5 c[dam]. For a founder that equals the gene
    contribution Lacy uses for p_k. Sources are processed in batches so
    memory stays O(n * batch), not O(n * n_sources).
    """
    n_src = len(source_idx)
    if n_src == 0:
        return []
    n = len(sire_i)
    sire_a = np.asarray(sire_i, dtype=np.intp)
    dam_a = np.asarray(dam_i, dtype=np.intp)
    desc_a = np.asarray(desc_idx, dtype=np.intp)
    out = np.zeros(n_src, dtype=np.float64)
    batch = _LACY_SOURCE_BATCH
    for start in range(0, n_src, batch):
        cols = np.asarray(source_idx[start:start + batch], dtype=np.intp)
        width = int(cols.size)
        contrib = np.zeros((n, width), dtype=np.float64)
        contrib[cols, np.arange(width)] = 1.0
        for i in range(n):
            sire = int(sire_a[i])
            dam = int(dam_a[i])
            if sire >= 0 and dam >= 0:
                contrib[i] = 0.5 * contrib[sire] + 0.5 * contrib[dam]
            elif sire >= 0:
                contrib[i] = 0.5 * contrib[sire]
            elif dam >= 0:
                contrib[i] = 0.5 * contrib[dam]
            self_cols = np.flatnonzero(cols == i)
            if self_cols.size:
                contrib[i, self_cols] = 1.0
        if desc_a.size:
            out[start:start + width] = contrib[desc_a].sum(axis=0)
    return [float(value) for value in out]


def _lacy_family_contributions(pedobj, mode):
    """
    Founder-source contributions over the pedigree ``mode`` implies.

    This is the memory-efficient route ``effective_founders_lacy`` exists to
    provide: it never forms a relationship matrix over the whole pedigree.
    For the supported phantom partition it agrees with the population-wise
    route because both consume the same founder set and phantom enumeration.

    Returns ``(contributions, n_sources, n_descendants)`` with ``contributions``
    a list aligned to the founder sources, real founders first and phantoms in
    slot order.
    """
    missing = int(pedobj.kw['missing_parent'])
    fs, phantom_slots, ds = _lacy_partition(pedobj, mode)

    if phantom_slots:
        records, shift, _slots = _lacy_extended_pedigree(pedobj)
        sources = [f + shift for f in fs] + list(range(1, len(phantom_slots) + 1))
        descendants = [d + shift for d in ds]
    else:
        records = [_LacyRecord(int(x.animalID), int(x.sireID), int(x.damID))
                   for x in pedobj.pedigree]
        sources = list(fs)
        descendants = list(ds)

    id_to_index = {record.animalID: index for index, record in enumerate(records)}
    sire_i = [-1] * len(records)
    dam_i = [-1] * len(records)
    for index, record in enumerate(records):
        if record.sireID != missing:
            sire_i[index] = id_to_index.get(record.sireID, -1)
        if record.damID != missing:
            dam_i[index] = id_to_index.get(record.damID, -1)

    source_idx = [id_to_index[source] for source in sources]
    desc_idx = [id_to_index[descendant] for descendant in descendants]
    contributions = _lacy_source_contributions(sire_i, dam_i, source_idx, desc_idx)
    return contributions, len(sources), len(descendants)


def _write_a_effective_founders_lacy(
    pedobj, lenped, n_f, fs, phantom_slots, n_d, ds, f_e
):
    """Write ``{filetag}_fe_lacy_.dat`` for :func:`a_effective_founders_lacy`."""
    outputfile = f"{pedobj.kw['filetag']}_fe_lacy_.dat"
    with open(outputfile, 'w', encoding = 'utf-8') as aout:
        line = "=" * 60
        aout.write(f"{line}\n")
        aout.write(f"{lenped} animals\n")
        aout.write(f"{n_f} founder sources: {fs} plus phantoms {phantom_slots}\n")
        aout.write(f"{n_d} descendants: {ds}\n")
        aout.write(f"effective number of founders: {f_e}\n")
        aout.write(f"{line}\n")


def _write_effective_founders_lacy(
    pedobj, caller, mode, n_sources, f_contribs, f_contribs_weighted,
    f_contribs_sum, f_contribs_weighted_sum_sq, f_e,
):
    """Write ``{filetag}_fe_lacy.dat`` for :func:`effective_founders_lacy`."""
    outputfile = f"{pedobj.kw['filetag']}_fe_lacy.dat"
    with open(outputfile, "w", encoding = "utf-8") as aout:
        pyp_io.pyp_file_header(aout, caller)
        aout.write(f"{len(pedobj.pedigree)} animals in pedigree\n")
        aout.write(f"{n_sources} founder sources (mode={mode!r})\n")
        aout.write(f"Founder contributions: {f_contribs}\n")
        aout.write(f"Proportional founder contributions: {f_contribs_weighted}\n")
        aout.write(f"Sum of founder contributions: {f_contribs_sum}\n")
        aout.write(f"Sum of squared proportional contributions: {f_contribs_weighted_sum_sq}\n")
        aout.write(f"Effective founder number: {f_e}\n")
        pyp_io.pyp_file_footer(aout, caller)


def _write_boichard_founders_output(
    pedobj, q, n_f, founders, n_d, explicit, gens, gen, ngen, population, f_e
):
    """Write ``{filetag}_fe_boichard_.dat``."""
    outputfile = f"{pedobj.kw['filetag']}_fe_boichard_.dat"
    with open(outputfile, 'w', encoding='utf-8') as aout:
        line = '=' * 60 + '\n'
        aout.write(line)
        aout.write(f'q: {[q[i] for i in sorted(q)]}\n')
        aout.write(f'{n_f} founders: {sorted(founders)}\n')
        aout.write(f'{n_d} descendants\n')
        # On the explicit path there IS no generation, so none is
        # fabricated to keep the legacy wording usable. Reporting a
        # generation label here would be the same lie the `gen` coupling
        # already tells: it would name a selector that played no part.
        if explicit is None:
            aout.write(f'generations: {gens}\n')
            aout.write(f'{ngen} animals in generation {gen}\n')
        else:
            aout.write(f'{ngen} animals in the explicit reference '
                       f'population: {sorted(population)}\n')
        aout.write(f'effective number of founders: {f_e}\n')
        aout.write(line)


def _write_boichard_definite_output(
    pedobj, n_f, n_d, population, explicit, gens, gen, f_a, ancestors,
    phantoms, order,
):
    """Write ``{filetag}_fa_boichard_definite_.dat``."""
    outputfile = f"{pedobj.kw['filetag']}_fa_boichard_definite_.dat"
    with open(outputfile, 'w', encoding='utf-8') as aout:
        line = '=' * 60 + '\n'
        aout.write(line)
        aout.write(f'{n_f} candidate ancestors\n')
        aout.write(f'{n_d} animals in the population under study: {sorted(population)}\n')
        # The line above is already truthful on both paths and is kept. Only
        # the generation provenance below it changes: on the explicit path
        # there is no generation to report, and inventing one would misdescribe
        # how the population was chosen.
        if explicit is None:
            aout.write(f'generations: {gens}\n')
            aout.write(f'{n_d} animals in generation {gen}\n')
        else:
            aout.write('reference population supplied explicitly '
                       '(no generation selector)\n')
        aout.write(f'effective number of ancestors: {f_a}\n')
        aout.write(f'ancestors: {_boichard_label_ancestors(ancestors, phantoms)}\n')
        aout.write('ancestor contributions: '
                   f'{_boichard_label_contributions(order, phantoms)}\n')
        if phantoms:
            aout.write(_boichard_phantom_note(phantoms))
        aout.write(line)


def _write_boichard_indefinite_output(
    pedobj, n_taken, order, n_founders, c, taken, phantoms, f_l, f_u
):
    """Write ``{filetag}_fa_boichard_indefinite_.dat``."""
    outputfile = f"{pedobj.kw['filetag']}_fa_boichard_indefinite_.dat"
    with open(outputfile, 'w', encoding='utf-8') as aout:
        line = '=' * 60 + '\n'
        aout.write(line)
        aout.write(f'{n_taken} ancestors taken of {len(order)} with positive contribution\n')
        aout.write(f'{n_founders} founders\n')
        aout.write(f'cumulated contribution c: {c}\n')
        aout.write('ancestors taken: '
                   f'{_boichard_label_contributions(taken, phantoms)}\n')
        aout.write(f'lower bound f_l: {f_l}\n')
        aout.write(f'upper bound f_u: {f_u}\n')
        if phantoms:
            aout.write(_boichard_phantom_note(phantoms))
        aout.write(line)


def _write_a_coefficients_output(
    pedobj, a, lenped, f_n, f_sum, f_avg, fnz_n, fnz_sum, fnz_avg,
    r_n, r_sum, r_avg, rnz_n, rnz_sum, rnz_avg,
):
    """Write the three ``a_coefficients`` analysis files."""
    outputfile2 = f"{pedobj.kw['filetag']}_rel_to_pop_.dat"
    with open(outputfile2, 'w', encoding = 'utf-8') as aout2:
        aout2.write("# Average relationship to population (renumbered ID, r)\n")
        for row in range(lenped):
            r_pop_avg = sum(
                pyp_nrm._matrix_value(a, row, col) for col in range(lenped) if row != col
            ) / lenped
            aout2.write(f"{pedobj.pedigree[row].animalID} {r_pop_avg:.6f}\n")

    outputfile = f"{pedobj.kw['filetag']}_population_coefficients_.dat"
    with open(outputfile, 'w') as aout:
        line = "=" * 60
        aout.write(f"{line}\n")
        aout.write("# Population average coefficients of inbreeding and relationship\n")
        aout.write("f_n: {}\nf_sum: {:.6f}\nf_avg: {:.6f}\n".format(f_n, f_sum, f_avg))
        aout.write("fnz_n: {}\nfnz_sum: {:.6f}\nfnz_avg: {:.6f}\n".format(fnz_n, fnz_sum, fnz_avg))
        aout.write("r_n: {}\nr_sum: {:.6f}\nr_avg: {:.6f}\n".format(r_n, r_sum, r_avg))
        aout.write("rnz_n: {}\nrnz_sum: {:.6f}\nrnz_avg: {:.6f}\n".format(rnz_n, rnz_sum, rnz_avg))
        aout.write(f"{line}\n")

    outputfile = f"{pedobj.kw['filetag']}_individual_coefficients_.dat"
    with open(outputfile, 'w', encoding = 'utf-8') as aout:
        aout.write("# Individual coefficients of inbreeding\n")
        aout.write("# animalID f_a\n")
        for row in range(lenped):
            aout.write(
                f"{pedobj.pedigree[row].animalID}\t{pyp_nrm._coi_from_matrix(a, row):.4f}\n"
            )


def _write_fast_a_coefficients_output(
    pedobj, a, lenped, f_n, f_sum, f_avg, fnz_n, fnz_sum, fnz_avg,
    r_n, r_sum, r_avg, rnz_n, rnz_sum, rnz_avg,
):
    """Write the three ``fast_a_coefficients`` analysis files when ``file_io`` is on."""
    outputfile2 = f"{pedobj.kw['filetag']}_rel_to_pop_.dat"
    with open(outputfile2, 'w') as aout2:
        aout2.write("# Average relationship to population (renumbered ID, r)\n")
        for row in range(lenped):
            r_pop_avg = sum(
                pyp_nrm._matrix_value(a, row, col) for col in range(lenped) if col != row
            ) / lenped
            aout2.write(f"{pedobj.pedigree[row].animalID} {r_pop_avg:.6f}\n")

    outputfile = f"{pedobj.kw['filetag']}_population_coefficients_.dat"
    with open(outputfile, 'w') as aout:
        line = "=" * 60
        aout.write(f"{line}\n")
        aout.write("# Population average coefficients of inbreeding and relationship (fast_a_coefficients)\n")
        aout.write(f"f_n: {f_n}\nf_sum: {f_sum:.6f}\nf_avg: {f_avg:.6f}\n")
        aout.write(f"fnz_n: {fnz_n}\nfnz_sum: {fnz_sum:.6f}\nfnz_avg: {fnz_avg:.6f}\n")
        aout.write(f"r_n: {r_n}\nr_sum: {r_sum:.6f}\nr_avg: {r_avg:.6f}\n")
        aout.write(f"rnz_n: {rnz_n}\nrnz_sum: {rnz_sum:.6f}\nrnz_avg: {rnz_avg:.6f}\n")
        aout.write(f"{line}\n")

    outputfile = f"{pedobj.kw['filetag']}_individual_coefficients_.dat"
    with open(outputfile, 'w') as aout:
        aout.write("# Individual coefficients of inbreeding\n")
        aout.write("# animalID f_a\n")
        for row in range(lenped):
            aout.write(
                f"{pedobj.pedigree[row].animalID}\t{pyp_nrm._coi_from_matrix(a, row):.4f}\n"
            )


def _write_theoretical_ne_output(pedobj, ns, nd, ne):
    """Write ``{filetag}_ne_from_metadata_.dat``."""
    outputfile = f"{pedobj.kw['filetag']}_ne_from_metadata_.dat"
    with open(outputfile, 'w') as aout:
        line = "=" * 60 + '\n'
        aout.write(line)
        aout.write("# Theoretical effective population size (N_e)\n")
        aout.write("#   n_sires = number of sires\n")
        aout.write("#   n_dams = number of dams\n")
        aout.write("#   n_e = effective population size\n")
        aout.write(f"n_sires: {ns}\n")
        aout.write(f"n_dams: {nd}\n")
        aout.write(f"n_e: {ne}\n")
        aout.write(line)


def a_effective_founders_lacy(
    pedobj,
    a=None,
    mode=_UNSET,
    half=_UNSET,
    output=True,
) -> Dict[str, float]:
    """
    Calculate the number of effective founders in a pedigree using the exact method of Lacy.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : np.ndarray, optional
        A numerator relationship matrix. If None, it will be generated.
    mode : str, optional
        Half-founder treatment. The only supported value is ``'phantom'``
        (the default): each unknown parent is itself a founder, and the
        half-founder remains a descendant (Lacy 1989 p.113). Historical
        ``'strict'`` and ``'absorb'`` values raise ``PyPedalUsageError``.
    half : bool, optional
        Not supported. Historical ``half=True`` / ``half=False`` raise
        ``PyPedalUsageError``. Use the default ``mode='phantom'``.
    output : bool, optional
        If True (the default), write ``{filetag}_fe_lacy_.dat``. If False,
        perform the calculation and return the same dictionary without
        writing that analysis file.

    Returns
    -------
    dict
        A dictionary of results, including the effective founder number.

    Raises
    ------
    PyPedalUsageError
        If ``mode`` is not ``'phantom'``, if ``half`` is supplied, or if both
        ``mode`` and ``half`` are supplied.
    PyPedalError
        On an unexpected internal failure. Typed PyPedal errors are
        re-raised. ``{"fa_effective_founders": -999.9}`` is not a
        supported failure return.
    """
    mode = _resolve_lacy_mode(mode, half, 'a_effective_founders_lacy')
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered a_effective_founders_lacy()')

    try:
        # Generate numerator relationship matrix if not provided
        if a is None:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
            # try:
                
            # except Exception:
            #     return -999.9

        lenped = len(pedobj.pedigree)
        fs, phantom_slots, ds = _lacy_partition(pedobj, mode)
        n_f, n_d = len(fs) + len(phantom_slots), len(ds)
        if n_d == 0:
            raise PyPedalError(
                'a_effective_founders_lacy: the pedigree has no descendants, so '
                'there is no population whose founder representation could be '
                'measured.')

        # Contribution of each founder source to each descendant.
        p = numpy.zeros((n_d, n_f), dtype=float)
        if phantom_slots:
            # The relationships must be solved on the COMPLETED pedigree; a
            # phantom's contribution cannot be derived from the original NRM.
            # Any `a` the caller supplied is for the original pedigree and does
            # not describe the phantoms, so it is deliberately not used here.
            records, shift, _slots = _lacy_extended_pedigree(pedobj)
            a_ext = pyp_nrm.fast_a_matrix(records, pedobj.kw)
            for row, descendant in enumerate(ds):
                for col, founder in enumerate(fs):
                    p[row, col] = pyp_nrm._matrix_value(
                        a_ext, founder + shift - 1, descendant + shift - 1)
                for offset in range(len(phantom_slots)):
                    p[row, len(fs) + offset] = pyp_nrm._matrix_value(
                        a_ext, offset, descendant + shift - 1)
        else:
            for row, descendant in enumerate(ds):
                for col, founder in enumerate(fs):
                    p[row, col] = pyp_nrm._matrix_value(a, founder - 1, descendant - 1)

        # Proportional contributions: Lacy p.114 divides each founder's total
        # contribution by the number of living descendants.
        p_sums = [np.sum(p[:, col]) for col in range(n_f)]
        rel_p = [sum_val / n_d for sum_val in p_sums]
        rel_p_sq = [val ** 2 for val in rel_p]

        # With the phantom founder set the paper specifies, the contributions
        # form a genuine probability vector. Nothing normalises them here.
        if mode == 'phantom':
            pyp_validate.check_contribution_vector(
                pedobj, rel_p, 'a_effective_founders_lacy: founder contributions')

        # Calculate effective number of founders
        sum_rel_p_sq = sum(rel_p_sq)
        f_e = 1.0 / sum_rel_p_sq if sum_rel_p_sq > 0 else 0.0

        # Verbose output if enabled
        if pedobj.kw.get('messages') == 'verbose':
            print("=" * 60)
            print(f"animals:\t{len(fs) + len(ds)}")
            print(f"founders:\t{n_f}")
            print(f"descendants:\t{n_d}")
            print(f"f_e:\t\t{f_e:.3f}")
            print("=" * 60)

        # Prepare output dictionary
        out_dict = {
            "fa_animal_count": len(fs) + len(ds),
            "fa_founder_count": n_f,
            "fa_descendant_count": n_d,
            "fa_effective_founders": f_e,
        }

        if output:
            _write_a_effective_founders_lacy(
                pedobj, lenped, n_f, fs, phantom_slots, n_d, ds, f_e)

        if pedobj.kw.get('debug_messages'):
            logger.info('Exited a_effective_founders_lacy()')

        return out_dict

    except PyPedalError:
        raise
    except Exception as e:
        logger.error("Error in a_effective_founders_lacy: %s", e)
        if pedobj.kw.get('debug_messages'):
            logger.info('Exited a_effective_founders_lacy() with failure')
        raise PyPedalError(
            'a_effective_founders_lacy: unexpected failure: %s' % e
        ) from e


def effective_founders_lacy(pedobj, mode=_UNSET, half=_UNSET, output=True) -> Dict[str, Union[int, float]]:
    """
    Calculate the number of effective founders in a pedigree using the exact method of Lacy.

    This version of the routine a_effective_founders_lacy() is designed to work
    with larger pedigrees: it walks founder-row contributions instead of
    forming a population-wise relationship matrix.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    mode : str, optional
        Half-founder treatment; see :func:`a_effective_founders_lacy`. The
        only supported value is ``'phantom'`` (the default).
    half : bool, optional
        Not supported. Historical ``half=True`` / ``half=False`` raise
        ``PyPedalUsageError``. Use the default ``mode='phantom'``.
    output : bool, optional
        If True (the default), write ``{filetag}_fe_lacy.dat``. If False,
        perform the calculation and return the same dictionary without
        writing that analysis file.

    Returns
    -------
    dict
        A dictionary of results, including the effective founder number.

    Raises
    ------
    PyPedalUsageError
        If ``mode`` is not ``'phantom'``, if ``half`` is supplied, or if both
        ``mode`` and ``half`` are supplied.
    PyPedalError
        On an unexpected internal failure. Typed PyPedal errors
        (including ``PyPedalValidationError``) are re-raised.
        ``{"fa_effective_founders": -999.9}`` is not a supported
        failure return.
    """
    mode = _resolve_lacy_mode(mode, half, 'effective_founders_lacy')
    logger.info("Entered effective_founders_lacy()")

    caller = "pyp_metrics.effective_founders_lacy"
    out_dict = {}

    try:
        # Ensure the pedigree is renumbered
        if not pedobj.kw.get("pedigree_is_renumbered", False):
            logger.info("[NOTE]: The pedigree is not renumbered. Renumbering...")
            pedobj.kw["renumber"] = True
            pedobj.renumber()

        # Founder-source contributions via the founder-row recurrence. This
        # shares _lacy_partition and the phantom enumeration with
        # a_effective_founders_lacy, so the two entry points cannot disagree
        # about what the founder set is -- which, before this, they did.
        _contribs, _n_sources, _n_desc = _lacy_family_contributions(pedobj, mode)
        _f_contribs = {index: value for index, value in enumerate(_contribs)}
        _f_contribs_sum = sum(_contribs)

        # Lacy p.114 divides each founder's contribution by the number of living
        # descendants. With the founder set the paper specifies the two
        # normalisations coincide, because the contributions then sum to exactly
        # the descendant count -- which is the same fact as q summing to one.
        _f_contribs_weighted = {
            index: (value / _n_desc if _n_desc else 0.0)
            for index, value in _f_contribs.items()
        }
        _f_contribs_weighted_sum_sq = sum(
            value ** 2 for value in _f_contribs_weighted.values())

        if mode == 'phantom':
            pyp_validate.check_contribution_vector(
                pedobj, list(_f_contribs_weighted.values()),
                'effective_founders_lacy: founder contributions')

        # Calculate the effective number of founders
        _f_e = 1.0 / _f_contribs_weighted_sum_sq if _f_contribs_weighted_sum_sq > 0 else 0.0

        # Verbose output if enabled
        if pedobj.kw.get("messages") == "verbose":
            print(f"\tFounder contributions: {_f_contribs}")
            print(f"\tProportional founder contributions: {_f_contribs_weighted}")
            print(f"\tSum of founder contributions: {_f_contribs_sum}")
            print(f"\tSum of squared proportional contributions: {_f_contribs_weighted_sum_sq}")
            print(f"Effective founder number (f_e): {_f_e}")

        if output:
            _write_effective_founders_lacy(
                pedobj, caller, mode, _n_sources, _f_contribs, _f_contribs_weighted,
                _f_contribs_sum, _f_contribs_weighted_sum_sq, _f_e)

        # Populate the output dictionary
        out_dict = {
            "fa_animal_count": len(pedobj.pedigree),
            "fa_founder_count": _n_sources,
            "fa_descendant_count": _n_desc,
            "fa_effective_founders": _f_e,
        }

        # Postcondition. f_e is the reciprocal of a sum of
        # squared founder contributions and cannot be below 1.
        #
        # No contribution-vector check is applied here on purpose:
        # _f_contribs_weighted is divided by its own total, so "sums to 1" is
        # true by construction and would be a vacuous assertion. The meaningful
        # unforced-sum check lives in the independent oracle
        # independent Lacy oracle, where the vector is not normalised.
        # `_f_contribs_weighted` is the contribution vector f_e was computed
        # from, so its length is exactly the k in f_e <= k.
        pyp_validate.check_effective_number(
            pedobj, float(_f_e), 'effective_founders_lacy: f_e',
            n_contributors=len(_f_contribs_weighted))

        logger.info("Exited effective_founders_lacy()")
        return out_dict

    except PyPedalError:
        raise
    except Exception as e:
        logger.error("Error in effective_founders_lacy: %s", e)
        logger.info("Exited effective_founders_lacy() with failure")
        raise PyPedalError(
            'effective_founders_lacy: unexpected failure: %s' % e
        ) from e


def boichard_probabilities_of_gene_origin(pedobj, reference):
    """
    Appendix A of Boichard, Maignel & Verrier (1997), article p.22, verbatim.

    Returns ``(q, founders)``: the probability of gene origin for every animal,
    keyed by animalID, and the list of founder IDs the paper's own definition
    selects.

    The algorithm::

        2. initialise q with 1 for animals in the population under study
        3. process the pedigree YOUNGEST to OLDEST:
               q(sire(i)) += 0.5*q(i);  q(dam(i)) += 0.5*q(i)
        4. if an animal is a 'half founder' -- one known parent and one unknown
           -- multiply ITS contribution by 0.5. Divide the vector q by N.

    TWO THINGS THE ORDER OF STEPS 3 AND 4 DECIDES, and both were wrong before:

    * The known parent of a half-founder receives half of the animal's FULL
      contribution. The halving in step 4 happens *after* the whole propagation
      pass, not during it. Halving first passes 0.25*q(i) up the pedigree
      instead of 0.5*q(i), which is what the port did.
    * A half-founder IS a founder. Article p.7: "A founder is defined as an
      ancestor with unknown parents. Note that when an animal has only one
      known parent, the animal is considered as a founder." The port required
      BOTH parents to be unknown, so half-founders were excluded from the
      founder set entirely and their retained contribution vanished from
      sum(q^2) -- leaving q summing to less than 1.

    The paper adds that halving the half-founder's own contribution "is
    equivalent to considering the unknown parent as a founder". The two
    descriptions give the same multiset of contributions, so f_e is identical
    either way; they differ only in which label the 0.5 is filed under. Lacy
    (1989) p.113 takes the phantom-parent view of the same fact.
    """
    missing = int(pedobj.kw['missing_parent'])
    animals = list(pedobj.pedigree)
    ids = [int(x.animalID) for x in animals]
    index_of = {aid: i for i, aid in enumerate(ids)}
    sire = [int(x.sireID) for x in animals]
    dam = [int(x.damID) for x in animals]

    q = np.zeros(len(animals), dtype=float)
    for animal_id in reference:
        q[index_of[int(animal_id)]] = 1.0

    # Step 3 -- YOUNGEST to OLDEST, with the FULL q(i) propagated.
    for i in range(len(animals) - 1, -1, -1):
        if sire[i] != missing:
            q[index_of[sire[i]]] += 0.5 * q[i]
        if dam[i] != missing:
            q[index_of[dam[i]]] += 0.5 * q[i]

    # Step 4 -- and only now is the half-founder's own contribution halved.
    for i in range(len(animals)):
        if (sire[i] == missing) != (dam[i] == missing):
            q[i] *= 0.5

    q /= float(len(reference))

    # Founders per the paper's definition (p.7): unknown parents, and animals
    # with exactly one known parent.
    founders = [ids[i] for i in range(len(animals))
                if sire[i] == missing or dam[i] == missing]
    return {ids[i]: float(q[i]) for i in range(len(animals))}, founders


def a_effective_founders_boichard(pedobj, a: Optional[np.ndarray] = None, gen: Optional[int] = None,
                                  *, reference: Optional[Iterable[int]] = None, output: bool = True) -> float:
    """
    Effective number of founders, ``f_e = 1 / sum(q_k^2)`` over founders.

    Implements Appendix A and equation 1 of Boichard, Maignel & Verrier (1997),
    *Genet Sel Evol* 29:5-23, article pp.7 and 22, via
    :func:`boichard_probabilities_of_gene_origin`.

    REPAIRED AGAINST THE SOURCE. Two independent defects, both in the
    half-founder path and both shared with PyPedal 2.0.4:

    * the halving was applied *before* propagation, so a half-founder's known
      parent received a quarter of the contribution instead of a half;
    * half-founders were excluded from the founder set, because
      ``NewAnimal.founder`` is ``'y'`` only when both parents are unknown, so
      their retained contribution was dropped from the sum entirely.

    Together these made ``q`` sum to less than 1 on any pedigree containing a
    half-founder, which is not a probability vector, so ``1 / sum(q^2)`` was not
    an effective number of founders. Pedigrees without half-founders are
    unaffected and their values are unchanged.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : numpy.ndarray, optional
        Accepted for API compatibility and **ignored**. Appendix A propagates
        contributions directly and needs no relationship matrix; the previous
        implementation built one on every call and never read it.
    gen : int, optional
        Generation to use as the population under study. Defaults to the most
        recent generation, compared numerically rather than lexicographically. Mutually
        exclusive with ``reference``.
    reference : iterable of int, keyword-only, optional
        The population under study stated directly: the *N* animals carrying
        the gene pool of interest (Appendix A step 1, article p.22). Elements
        are real, current ``animalID`` integers from this pedigree; original
        and string IDs are not translated, so map them through
        ``pedobj.idmap`` first. Order is irrelevant and the set is
        canonicalised. An empty set, a duplicate, an unknown ID, the
        missing-parent sentinel, a Reading-C phantom ID, a bool, and any
        non-integral value are each refused with a
        :class:`~PyPedal.pyp_errors.PyPedalUsageError`.

        When supplied, ``gen`` is not read and ``require_generations`` is not
        consulted -- there is no generation being used for selection, so the
        guard would have nothing scientifically relevant to validate. R is
        never inferred: not from ``igen``, pedigree depth, terminal status,
        birth year, sex or parent completeness.
    output : bool, optional, keyword-only
        If True (the default), write ``{filetag}_fe_boichard_.dat``. If False,
        perform the calculation and return the same ``f_e`` without writing
        that analysis file.

    Returns
    -------
    float
        The effective number of founders.
    """
    routine = 'a_effective_founders_boichard'
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered a_effective_founders_boichard()')

    # Argument-domain checks first, then the pedigree-domain guards in their
    # existing order. Only the SELECTION of R changes here; this routine's
    # validation domain is unchanged, and in particular it still has no
    # antichain guard -- it has never had one, because Appendix A propagates
    # contributions and never selects ancestors.
    explicit = _boichard_select_reference(pedobj, gen, reference, routine)
    if explicit is None:
        pyp_validate.require_generations(pedobj, routine)
    _boichard_require_topological(pedobj, routine)

    try:
        # `reference` is what the CALLER passed and may be None; `population`
        # is the R actually analysed. Keeping them separate stops one name
        # meaning two things, and keeps the parameter's declared type honest.
        if explicit is None:
            population, gen, most_recent, gens = _boichard_reference_population(
                pedobj, gen, routine)
        else:
            population = explicit

        q, founders = boichard_probabilities_of_gene_origin(pedobj, population)
        contributions = [q[animal_id] for animal_id in founders]

        # Appendix A step 4: dividing by N is what makes founder contributions
        # sum to one. Nothing normalises them, so this is a real check on the
        # propagation and on the half-founder handling -- it is exactly what
        # failed before the repair.
        pyp_validate.check_contribution_vector(
            pedobj, contributions, '%s: founder contributions' % routine)

        sum_sq = sum(value ** 2 for value in contributions)
        f_e = 0.0 if sum_sq == 0.0 else 1.0 / sum_sq

        n_f = len(founders)
        n_d = len(pedobj.pedigree) - n_f
        ngen = len(population)

        if pedobj.kw.get('messages') == 'verbose':
            print('=' * 60)
            print(f'animals:\t{len(pedobj.pedigree)}')
            print(f'founders:\t{n_f}')
            print(f'descendants:\t{n_d}')
            print(f'f_e:\t\t{f_e:.3f}')
            print('=' * 60)

        if output:
            _write_boichard_founders_output(
                pedobj, q, n_f, founders, n_d, explicit, gens if explicit is None else None,
                gen if explicit is None else None, ngen, population, f_e)

        if pedobj.kw.get('debug_messages'):
            logger.info('Exited a_effective_founders_boichard()')

        # f_e <= k, with k the number of founders the routine actually summed
        # over. That count now includes half-founders, so the bound is no
        # longer tighter than the mathematics allows.
        pyp_validate.check_effective_number(
            pedobj, float(f_e), '%s: f_e' % routine, n_contributors=n_f)

        return float(f_e)

    except pyp_validate.PyPedalValidationError:
        raise
    except PyPedalError:
        raise
    except Exception as e:
        # This returned -999.9, which is not an effective
        # founder number and which a caller storing the result in a report or a
        # database has no way to distinguish from a computed value.
        logger.error(f"Error in a_effective_founders_boichard: {e}")
        if pedobj.kw.get('debug_messages'):
            logger.info('Exited a_effective_founders_boichard() with failure')
        raise PyPedalError(
            "a_effective_founders_boichard() failed: %s: %s. No effective "
            "founder number was computed." % (type(e).__name__, e)) from e


# ---------------------------------------------------------------------------
# Boichard Appendix B -- the shared marginal-contribution engine
#
# Boichard, Maignel & Verrier (1997), Genet Sel Evol 29:5-23, Appendix B
# (article pp.22-23). The effective number of ancestors and the pp.9-10 bounds
# are the SAME sequence, exactly and truncated, so both are built on this one
# generator. Keeping them separate is how they drifted apart before: the
# bounded routine ended up computing one expression twice and calling the
# results f_l and f_u.
#
# The three residual semantics (tie-breaking, half-founders, antichain
# reference populations) are NOT in the paper; they are handled by a declared
# convention or a refusal here, never by silent inference.
# ---------------------------------------------------------------------------

#: Lowest-ID tie-break. Boichard p.10 states that where two ancestors hold the same
#: marginal contribution "the final result may depend on the chosen one", and
#: gives no rule. This is a PyPedal CONVENTION chosen for determinism -- lowest
#: animalID, i.e. the oldest animal after renumbering. It is not the paper's.
#:
#: It is deliberately an explicit code path rather than whatever numpy.argmax
#: happens to do with ties, so that it is testable and cannot change under us.
#:
#: Measured on the paper's own Table I: the tie-break changes only WHICH animal
#: is credited, never f_e, f_a or p_n, so no statistic depends on this choice.
#: The paper's own table resolves its three ties inconsistently by ID, so no
#: ID rule reproduces it.
BOICHARD_TIE_BREAK = 'lowest_id'


def _boichard_require_topological(pedobj, routine):
    """
    Parents must precede offspring in ``pedobj.pedigree``.

    Appendix B steps 5 and 6 are single passes in opposite directions, and both
    are only correct on a topologically ordered pedigree. PyPedal reorders by
    default, so this holds in normal use; it is checked rather than assumed
    because a violation would not raise, it would quietly compute a different
    number.
    """
    missing = str(pedobj.kw['missing_parent'])
    seen = set()
    for animal in pedobj.pedigree:
        for parent in (animal.sireID, animal.damID):
            if str(parent) != missing and str(parent) not in seen:
                raise PyPedalError(
                    '%s: animal %s has parent %s, which does not appear earlier '
                    'in the pedigree. Appendix B propagates in single passes and '
                    'requires parents to precede offspring. Load the pedigree '
                    'with reorder enabled.' % (routine, animal.animalID, parent))
        seen.add(str(animal.animalID))


# ---------------------------------------------------------------------------
# Half-founders in Appendix B, by pedigree completion
#
# This routine used to refuse. It no longer does. The source text did not
# change; the supported domain did.
#
# SOURCE-EXPLICIT, Appendix A step 4, article p.22:
#
#     "(4) if an animal is a 'half founder' (ie, with one known parent and one
#     unknown parent), multiply its contribution by 0.5. This is equivalent to
#     considering the unknown parent as a founder. Divide the vector q by N, so
#     that founder contributions sum to 1."
#
# SOURCE-EXPLICIT, p.7: "A founder is defined as an ancestor with unknown
# parents. Note that when an animal has only one known parent, the animal is
# considered as a founder."
#
# APPENDIX B IS SILENT. It neither restates the halving nor excludes it. Step 8
# restates only the normalisation -- Appendix A step 4 is [halve] + [divide by
# N], Appendix B step 8 is [divide by N] alone. Silence is not an exclusion,
# and writing "Appendix B says to complete the pedigree" would be false.
#
# WHAT WAS ADJUDICATED, and at what level:
#
#     Reading C -- pedigree-level phantom completion before an unchanged
#     Appendix B -- is MATHEMATICALLY IMPLIED / INDEPENDENTLY SUPPORTED.
#     It is NOT source-explicit Appendix-B text.
#
# It was selected against invariants the paper states about ITSELF, not against
# invariants of ours: p.7 founder contributions sum to one; p.8 "the marginal
# contributions over all ancestors sum to one" and "the number of ancestors
# with a positive contribution is less than or equal to the total number of
# founders"; p.17 "Consequently f_a is always lower than or equal to f_e". The
# two rival readings each violate some of these and Reading C violates none --
# measured on four half-founder
# fixtures plus a half-founder-free control that proves the scorer is not
# biased towards completion.
#
# This is NOT the gene-drop half-founder rule and NOT Lacy's phantom-founder
# rule, even though all three mechanically give an unknown parental side its
# own founder. They rest on different sources of different strengths and none
# justifies another.
# ---------------------------------------------------------------------------

def _boichard_label(animal_id, phantoms):
    """
    Render one selected ancestor for a human reader.

    A synthetic ID must never be presented as though it were a real pedigree
    animal. The integer is kept, because the report structure is a list of
    IDs and dropping it would make the report untraceable, but it is never
    printed bare: ``501`` becomes ``501=phantom (dam of 4)``.
    """
    slot = phantoms.get(animal_id)
    if slot is None:
        return str(animal_id)
    return '%d=phantom (%s of %s)' % (animal_id, slot[1], slot[0])


def _boichard_label_ancestors(ancestors, phantoms):
    return '[%s]' % ', '.join(_boichard_label(a, phantoms) for a in ancestors)


def _boichard_label_contributions(order, phantoms):
    return '{%s}' % ', '.join('%s: %r' % (_boichard_label(a, phantoms), v)
                              for a, v in order)


def _boichard_phantom_note(phantoms):
    """
    Say in the report itself what the synthetic IDs are and where they came
    from, so a reader who has only the .dat file is not left guessing why it
    names animals the pedigree does not contain.
    """
    listed = ', '.join(_boichard_label(pid, phantoms)
                       for pid in sorted(phantoms))
    return (
        '\n'
        'NOTE: %d phantom founder(s) were created for this analysis, one per\n'
        'unknown parental slot: %s.\n'
        'They are NOT in the pedigree, were not added to it, and have no\n'
        'animal record. They exist only for the duration of this analysis.\n'
        'Boichard Appendix A step 4 (p.22) states that halving a half-founder\n'
        'contribution "is equivalent to considering the unknown parent as a\n'
        'founder"; applying that as a completion of the pedigree ahead of an\n'
        'unchanged Appendix B (half-founder pedigree completion) is mathematically implied and\n'
        'independently supported rather than stated by Appendix B, which is\n'
        'silent.\n'
        % (len(phantoms), listed))


def boichard_phantom_ids(pedobj):
    """
    The phantom founders Reading C mints for this pedigree, as
    ``{phantom_id: (animal_id, 'sire'|'dam')}``.

    Provided so that a caller consuming ``boichard_marginal_contributions``
    can tell a synthetic ID from a real animal and say what it stands for.
    Deterministic, and computed the same way the engine computes it.

    Returns an empty dict for a pedigree with no half-founders, which is every
    published example in the paper.
    """
    _ids, _sire, _dam, phantoms, _n_founders = _boichard_completed_arrays(pedobj)
    return phantoms


def _boichard_completed_arrays(pedobj):
    """
    Reading C's completion, as parallel arrays for the Appendix-B engine.

    Returns ``(ids, sire, dam, phantoms, n_founders)``.

    ANALYSIS-LOCAL. Nothing here touches the caller's pedigree: no NewAnimal is
    created, ``pedobj.pedigree`` is not appended to, no attribute is set, and
    ``pedobj.kw`` is not read for anything but the missing-parent sentinel. The
    completed pedigree exists only for the duration of one analysis.

    PHANTOM IDENTITY. One phantom per unknown parental SLOT, numbered
    contiguously from ``max(real animalID) + 1`` in slot order.

      * Strictly above every real ID, so no collision is possible and no
        phantom can equal the ``missing_parent`` sentinel.
      * Above every real ID also means a real animal wins a real-vs-phantom
        tie under the lowest-ID convention, so that convention keeps meaning what it meant.
      * Keyed on ``(animal_id, side)``, so two half-founders that both record
        the sentinel get DIFFERENT phantoms. Two unknown parents are two
        unknown individuals; the input contains nothing that would license
        inferring they are the same one, and this makes that structural rather
        than a convention a test has to catch.
      * The numbers carry no biological meaning. That they carry no
        MATHEMATICAL meaning either is measured, not assumed, against a second
        arbitrary encoding in tests/test_boichard_half_founder.py.

    Phantoms are emitted FIRST so parents still precede offspring, which
    Appendix B's single-pass steps 5 and 6 require.

    ``n_founders`` is the p.7 founder count of the COMPLETED pedigree, and is
    the single source of ``f`` for both public routines. It is not
    ``NewAnimal.founder == 'y'``, which describes only the caller's real
    pedigree and is 'y' solely when BOTH parents are unknown.

    On a pedigree with no half-founders this is an exact no-op: the arrays are
    the same numbers in the same order they would have been built inline.
    """
    missing = int(pedobj.kw['missing_parent'])
    animals = list(pedobj.pedigree)

    # The enumeration of which parental slots are empty is a mechanical fact
    # about the pedigree, so it is shared with the Lacy path rather than
    # written twice. Only the enumeration is shared: Boichard completion and Lacy phantoms remain
    # separate scientific items resting on separate sources, and neither
    # closes the other.
    slots = lacy_phantom_slots(pedobj)

    ids = [int(a.animalID) for a in animals]
    sire = [int(a.sireID) for a in animals]
    dam = [int(a.damID) for a in animals]

    if not slots:
        n_founders = sum(1 for i in range(len(ids))
                         if sire[i] == missing and dam[i] == missing)
        return ids, sire, dam, {}, n_founders

    base = max(ids)
    fill = {slot: base + index + 1 for index, slot in enumerate(slots)}
    phantoms = {pid: slot for slot, pid in fill.items()}

    index_of = {aid: i for i, aid in enumerate(ids)}
    for (animal_id, side), pid in fill.items():
        position = index_of[animal_id]
        if side == 'sire':
            sire[position] = pid
        else:
            dam[position] = pid

    ordered = sorted(phantoms)
    ids = ordered + ids
    sire = [missing] * len(ordered) + sire
    dam = [missing] * len(ordered) + dam

    n_founders = sum(1 for i in range(len(ids))
                     if sire[i] == missing and dam[i] == missing)
    return ids, sire, dam, phantoms, n_founders


def _boichard_require_antichain(pedobj, reference, routine):
    """
    Reference-population antichain, the detectable half.

    Whether members of the population under study may themselves be selected as
    ancestors is source-silent: the paper's figures draw the population as
    unnumbered nodes, so the question never arises there. PyPedal's convention,
    unchanged from 2.0.4, is to exclude them -- their q is zeroed before
    selection. That convention is kept, declared, and pinned by a test.

    It is only lossy in one detectable case: when the reference population is
    not an antichain, i.e. one member is an ancestor of another. There the
    zeroing discards a real contribution rather than an uninteresting one, and
    the unresolved question stops being academic. That case raises.
    """
    reference = {int(x) for x in reference}
    missing = str(pedobj.kw['missing_parent'])
    for animal in pedobj.pedigree:
        if int(animal.animalID) not in reference:
            continue
        for parent in (animal.sireID, animal.damID):
            if str(parent) != missing and int(parent) in reference:
                raise PyPedalError(
                    '%s: the reference population is not an antichain -- animal '
                    '%s is in it and so is its parent %s. PyPedal excludes '
                    'reference-population members from ancestor selection, a '
                    'convention the source does not address, and that '
                    'convention silently discards a real contribution here. '
                    'Choose a reference population in which no member is an '
                    'ancestor of another.'
                    % (routine, animal.animalID, parent))


# ---------------------------------------------------------------------------
# Appendix A/B step 1 -- "define the population under study"
#
# There are exactly two ways to obtain R, and they live side by side here so
# neither can drift from the other.
#
# SOURCE-EXPLICIT, Appendix A step 1 and Appendix B step 1, article p.22,
# identical wording in both:
#
#     "(1) define the population under study, ie, the group of N animals
#     carrying the gene pool of interest;"
#
# The verb is DEFINE. It is an instruction addressed to the analyst, not a rule
# to compute, and no generation number appears anywhere in either appendix.
# `_boichard_explicit_reference` exposes that input directly.
# `_boichard_reference_population` is the legacy route, which reaches the same
# place by matching the pedformat `g` column -- kept bit-exact for
# compatibility, and the reason the routines are unreachable on a pedigree that
# has no such column.
#
# Representing R as a Python collection of animalIDs is SOFTWARE/API DESIGN.
# Nothing here claims Boichard specifies a Python API.
# ---------------------------------------------------------------------------

def _boichard_explicit_reference(pedobj, reference, routine):
    """
    Validate and canonicalise an analyst-supplied population under study.

    Returns the same mathematical set as a sorted list of ``animalID`` ints.

    THE IDENTIFIER DOMAIN is real, current ``animalID`` integers -- exactly
    what :func:`boichard_probabilities_of_gene_origin` and
    :func:`boichard_marginal_contributions` already consume. ``originalID``,
    string IDs and names are **not** translated: on a renumbered pedigree the
    two domains overlap, so guessing which one the caller meant could silently
    analyse the wrong animals. The failure class this exists to
    remove. A caller holding original IDs maps them through ``pedobj.idmap``
    first. This is scope containment, not a claim that internal integer IDs are
    the best long-term user interface.

    EVERY REJECTION BELOW IS FAIL-LOUD, and three of them are arithmetic rather
    than fastidiousness:

    * **Duplicates.** Both engines normalise by the LENGTH of the sequence --
      ``q /= float(len(reference))`` in Appendix A step 4, ``n =
      float(len(ref_index))`` in Appendix B step 8 -- so a repeated ID divides
      by the wrong N and every contribution comes out wrong. Silently applying
      ``set()`` would change the caller's stated intent instead of reporting
      it, so the duplicate is refused and named.
    * **An empty R.** Same two divisions, with N = 0. The legacy route already
      declines an empty generation.
    * **Reading-C phantom IDs.** ``boichard_marginal_contributions`` runs on
      the COMPLETED pedigree, whose ``index_of`` does contain the phantoms, so
      a caller-supplied phantom would be indexed happily and treated as a real
      animal. Rejection is therefore structural -- membership in
      ``pedobj.pedigree``, which phantoms are never part of -- rather than a
      numeric range test that a change in phantom numbering could invalidate.

    IDs are never coerced. ``"12"`` and ``12.0`` are refused rather than
    converted, and ``bool`` is refused ahead of the integer check because it
    subclasses ``int`` and ``True`` would otherwise become animal 1.

    CANONICALISATION IS LAST, after duplicates and membership have been
    checked, so no rejection can be masked by deduplicating or reordering
    first. Sorting is safe because the order of R is arithmetically inert --
    it only sets ``q`` entries to 1 and builds an index list, both
    order-independent. Python iteration order must never acquire scientific meaning.
    """
    if isinstance(reference, (str, bytes, bytearray)):
        raise PyPedalUsageError(
            '%s: reference must be a collection of animalIDs, not %s. A string '
            'is iterable, so this would be taken apart into single characters '
            'and each one read as an ID.'
            % (routine, type(reference).__name__))
    try:
        supplied = list(reference)
    except TypeError as exc:
        raise PyPedalUsageError(
            '%s: reference must be an iterable of animalIDs; got %r.'
            % (routine, reference)) from exc

    if not supplied:
        raise PyPedalUsageError(
            '%s: the reference population is empty, so there is no population '
            'under study. Appendix A step 4 and Appendix B step 8 both divide '
            'by N, the number of animals in it.' % routine)

    real = {int(animal.animalID) for animal in pedobj.pedigree}
    missing = int(pedobj.kw['missing_parent'])

    seen = set()
    canonical = []
    for position, item in enumerate(supplied):
        if isinstance(item, bool):
            raise PyPedalUsageError(
                '%s: reference[%d] is %r. bool is a subclass of int in Python, '
                'so this would silently be read as animal %d.'
                % (routine, position, item, int(item)))
        try:
            animal_id = operator.index(item)
        except TypeError as exc:
            raise PyPedalUsageError(
                '%s: reference[%d] is %r, which is not an integer animalID. '
                'IDs are not coerced; convert it yourself if that is what you '
                'meant.' % (routine, position, item)) from exc

        if animal_id == missing:
            raise PyPedalUsageError(
                '%s: reference[%d] is %r, which is the missing-parent '
                'sentinel, not an animal.' % (routine, position, item))
        if animal_id not in real:
            raise PyPedalUsageError(
                '%s: reference[%d] is %r, which is not an animal in this '
                'pedigree. reference holds current animalID values; original '
                'or string IDs are not translated, so map them through '
                'pedobj.idmap first.' % (routine, position, item))
        if animal_id in seen:
            raise PyPedalUsageError(
                '%s: animalID %d appears more than once in reference. The '
                'population under study is a set of N animals and both '
                'appendices divide by N, so a repeated animal would change '
                'every contribution. Remove the duplicate rather than relying '
                'on PyPedal to discard it.' % (routine, animal_id))

        seen.add(animal_id)
        canonical.append(animal_id)

    return sorted(canonical)


def _boichard_select_reference(pedobj, gen, reference, routine):
    """
    Obtain R, from exactly one of the two mechanisms.

    Returns the validated explicit population, or ``None`` to mean "use the
    legacy generation route" -- which the caller then runs unchanged.

    ``gen`` and ``reference`` both name the population under study, so
    supplying both is a caller error rather than a precedence question.
    Quietly preferring one is how an analyst ends up with an analysis of a
    population they did not choose. Omitting
    ``gen`` is not "explicitly supplying" it: its default is already ``None``,
    so no sentinel object is needed to tell the two apart.
    """
    if reference is None:
        return None
    if gen is not None:
        raise PyPedalUsageError(
            '%s: gen and reference both select the population under study, so '
            'supply exactly one. gen picks it by matching the pedformat "g" '
            'column; reference states it directly.' % routine)
    return _boichard_explicit_reference(pedobj, reference, routine)


def _boichard_reference_population(pedobj, gen, routine):
    """The N animals carrying the gene pool of interest (Appendix B step 1)."""
    gens = []
    for animal in pedobj.pedigree:
        if animal.gen not in gens:
            gens.append(animal.gen)
    most_recent = _most_recent_generation(gens, routine)
    if gen is None:
        gen = most_recent
    elif str(gen) not in {str(g) for g in gens}:
        raise PyPedalError(
            '%s: generation %r was requested but the pedigree contains only '
            '%r.' % (routine, gen, sorted(str(g) for g in gens)))
    wanted = str(gen)
    reference = [int(a.animalID) for a in pedobj.pedigree if str(a.gen) == wanted]
    if not reference:
        raise PyPedalError(
            '%s: generation %r contains no animals, so there is no population '
            'under study.' % (routine, gen))
    return reference, gen, most_recent, gens


def boichard_marginal_contributions(pedobj, reference, tie_break=BOICHARD_TIE_BREAK,
                                    tol=1e-12):
    """
    Appendix B, verbatim. Yields ``(animalID, marginal contribution)`` in the
    order the algorithm selects ancestors, contributions already divided by N.

    Each round (steps 3-8):

    3. delete the pedigree information of the ancestors already found, so each
       becomes a 'pseudo founder';
    4. initialise ``q`` with 1 for the population under study and ``a`` with 1
       for the already-selected ancestors;
    5. process the pedigree YOUNGEST to OLDEST, ``q(parent) += 0.5*q(i)``;
    6. process the pedigree OLDEST to YOUNGEST, ``a(i) += 0.5*a(parent)``;
    7. ``p(i) = q(i) * (1 - a(i))``;
    8. select the highest ``p`` and divide it by N.

    The two passes run in OPPOSITE directions and both directions are explicit
    in the paper. Running either the wrong way changes the answer silently.

    The caller's pedigree is never mutated: step 3's deletions are applied to
    local parent arrays, not to the ``NewAnimal`` records.

    Half-founder completion. The pedigree the algorithm below runs on is the COMPLETED one --
    every unknown parental slot given a founder of its own. That completion is
    the whole of Reading C; steps 3 to 8 have no half-founder case, because
    after completion there are no half-founders. Appendix A step 4's halving is
    therefore NOT applied here as well: its precondition is false everywhere in
    the completed pedigree.

    A phantom may itself be selected and yielded, on a half-founder pedigree.
    That is forced rather than chosen -- a phantom is an ordinary founder, it
    is not in the population under study so R3's zeroing never reaches it, and
    withholding its contribution would leave the marginal contributions summing
    to less than one, contradicting p.8. Callers that need to tell a synthetic
    ID from a real animal have ``boichard_phantom_ids()``; the value domain of
    this generator now includes them, though its signature is unchanged.
    """
    if tie_break not in ('lowest_id', 'highest_id'):
        raise PyPedalUsageError(
            'boichard_marginal_contributions: tie_break must be "lowest_id" or '
            '"highest_id", not %r' % (tie_break,))

    if not len(pedobj.pedigree):
        return
    missing = int(pedobj.kw['missing_parent'])
    ids, base_sire, base_dam, _phantoms, _n_f = _boichard_completed_arrays(pedobj)
    lenped = len(ids)
    index_of = {aid: i for i, aid in enumerate(ids)}

    ref_index = [index_of[int(x)] for x in reference]
    n = float(len(ref_index))
    selected = []

    while len(selected) < lenped:
        # Step 3 -- pseudo founders, applied to local arrays only.
        sire, dam = list(base_sire), list(base_dam)
        for pos in selected:
            sire[pos] = missing
            dam[pos] = missing

        # Step 4
        q = np.zeros(lenped, dtype=float)
        a_vec = np.zeros(lenped, dtype=float)
        for pos in ref_index:
            q[pos] = 1.0
        for pos in selected:
            a_vec[pos] = 1.0

        # Step 5 -- q, YOUNGEST to OLDEST.
        for i in range(lenped - 1, -1, -1):
            if sire[i] != missing:
                q[index_of[sire[i]]] += 0.5 * q[i]
            if dam[i] != missing:
                q[index_of[dam[i]]] += 0.5 * q[i]

        # Step 6 -- a, OLDEST to YOUNGEST.
        chosen = set(selected)
        for i in range(lenped):
            if i in chosen:
                continue            # pinned at 1; its parents were deleted
            acc = 0.0
            if sire[i] != missing:
                acc += 0.5 * a_vec[index_of[sire[i]]]
            if dam[i] != missing:
                acc += 0.5 * a_vec[index_of[dam[i]]]
            a_vec[i] += acc

        # Step 6 postcondition. `a` is a convex combination of settled parental
        # values, so it cannot leave [0, 1]. This is asserted, NOT clamped: a
        # clamp would convert an implementation defect into a plausible number.
        if a_vec.size and (a_vec.min() < -1e-9 or a_vec.max() > 1.0 + 1e-9):
            bad = int(np.argmax(np.abs(a_vec - 0.5)))
            raise pyp_validate.PyPedalValidationError(
                'boichard_marginal_contributions: a(%d) = %r is outside [0, 1]. '
                'Appendix B step 6 is a pull, so a is a convex combination of '
                'settled parental values and cannot leave the interval. This is '
                'an implementation defect.' % (ids[bad], float(a_vec[bad])))

        # Step 7
        p = q * (1.0 - a_vec)

        # Declared convention: reference-population members are
        # not candidate ancestors. _boichard_require_antichain() has already
        # refused the one case where this discards a real contribution.
        for pos in ref_index:
            p[pos] = 0.0
        for pos in selected:
            p[pos] = 0.0

        best = float(p.max())
        if best <= tol:
            return

        # Step 8, with the lowest-ID tie-break applied explicitly rather than left to argmax.
        tied = [i for i in range(lenped) if p[i] >= best - tol]
        pos = min(tied, key=lambda i: ids[i]) if tie_break == 'lowest_id' \
            else max(tied, key=lambda i: ids[i])
        selected.append(pos)
        yield ids[pos], float(p[pos] / n)


def a_effective_ancestors_definite(pedobj, a: Optional[np.ndarray] = None, gen: Optional[int] = None,
                                   *, reference: Optional[Iterable[int]] = None, output: bool = True) -> float:
    """
    Effective number of ancestors, f_a = 1 / sum(p_k^2).

    Implements Appendix B of Boichard, Maignel & Verrier (1997), *Genet Sel
    Evol* 29:5-23, article pp.22-23, via :func:`boichard_marginal_contributions`.

    HISTORY. Until the paper was read this routine had no ``a`` vector at all:
    it zeroed selected ancestors in ``q`` and re-propagated, which is not the
    marginal-contribution step, and its update pass ran oldest-to-youngest,
    the reverse of step 5. It returned 0.0816 on valid input -- below 1, and so
    impossible for the reciprocal of a sum of squared probabilities. Both
    defects are now repaired against the source, and PyPedal 2.0.4 is no
    reference here: it carried the ``a`` vector but wrote its half-founder
    halving to a positional index of a reversed list.

    PRECONDITIONS. Generations must be defined; parents must precede offspring;
    the pedigree must contain no unresolved half-founder in the uncompleted graph
    (completion is applied internally); and the reference population must be an antichain.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : np.ndarray, optional
        Accepted for API compatibility and **ignored**. Appendix B needs no
        numerator relationship matrix -- it propagates contributions directly.
        The previous implementation built one on every call and never read it.
    gen : int, optional
        Generation to use as the population under study. Defaults to the most
        recent generation, compared numerically rather than lexicographically. Mutually
        exclusive with ``reference``.
    reference : iterable of int, keyword-only, optional
        The population under study stated directly (Appendix B step 1, article
        p.22). Same contract as in
        :func:`a_effective_founders_boichard`: real, current ``animalID``
        integers, order irrelevant, every malformed case refused with a
        :class:`~PyPedal.pyp_errors.PyPedalUsageError`.

        This routine's validation domain is unchanged, so an explicitly
        supplied R must still be an antichain. How R is supplied is independent
        of which R is structurally admissible.
    output : bool, optional, keyword-only
        If True (the default), write ``{filetag}_fa_boichard_definite_.dat``.
        If False, perform the calculation and return the same ``f_a`` without
        writing that analysis file.

    Returns
    -------
    float
        The effective number of ancestors.
    """
    routine = 'a_effective_ancestors_definite'
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered a_effective_ancestors_definite()')

    explicit = _boichard_select_reference(pedobj, gen, reference, routine)
    if explicit is None:
        pyp_validate.require_generations(pedobj, routine)
    _boichard_require_topological(pedobj, routine)

    # `reference` is what the caller passed and may be None; `population` is
    # the R actually analysed.
    if explicit is None:
        population, gen, most_recent, gens = _boichard_reference_population(
            pedobj, gen, routine)
    else:
        population = explicit
    _boichard_require_antichain(pedobj, population, routine)

    _ids, _s, _d, phantoms, n_founders = _boichard_completed_arrays(pedobj)
    order = list(boichard_marginal_contributions(pedobj, population))
    contribs = {animal_id: value for animal_id, value in order}
    ancestors = [animal_id for animal_id, _ in order]

    sum_p_sq = sum(value ** 2 for value in contribs.values())
    f_a = 1.0 / sum_p_sq if sum_p_sq > 0 else 0.0

    # p.8: "the marginal contributions over all ancestors sum to one". This is
    # not enforced by construction -- nothing normalises the contributions --
    # so it is a real check on the propagation and on step 3's deletions.
    pyp_validate.check_contribution_vector(
        pedobj, list(contribs.values()), '%s: marginal contributions' % routine)

    n_f = len([x for x in pedobj.pedigree if int(x.animalID) not in set(population)])
    n_d = len(population)

    if pedobj.kw.get('messages') == 'verbose':
        print('=' * 60)
        print(f'animals:\t{len(pedobj.pedigree)}')
        print(f'ancestors:\t{len(ancestors)}')
        print(f'descendants:\t{n_d}')
        print(f'f_a:\t\t{f_a:.3f}')
        print('=' * 60)

    if output:
        _write_boichard_definite_output(
            pedobj, n_f, n_d, population, explicit,
            gens if explicit is None else None, gen if explicit is None else None,
            f_a, ancestors, phantoms, order)

    if pedobj.kw.get('debug_messages'):
        logger.info('Exited a_effective_ancestors_definite()')

    # p.8: "The number of ancestors with a positive contribution is less than or
    # equal to the total number of founders", so the founder count is the k in
    # f_a <= k. Founders follow the paper's own definition (p.7) on the pedigree
    # the engine actually ran on, which after half-founder completion is the COMPLETED one.
    pyp_validate.check_effective_number(
        pedobj, float(f_a), '%s: f_a' % routine,
        n_contributors=max(n_founders, 1))

    return float(f_a)


def a_effective_ancestors_indefinite(pedobj, a: Optional[np.ndarray] = None, gen: Optional[int] = None, n: int = 25,
                                     *, reference: Optional[Iterable[int]] = None, output: bool = True) -> Tuple[float, float]:
    """
    Lower and upper bounds on the effective number of ancestors, ``(f_l, f_u)``.

    Implements the scheme on article pp.9-10 of Boichard, Maignel & Verrier
    (1997). With ``c`` the contribution explained by the first ``n`` ancestors
    and ``f`` the total number of founders:

        f_u = 1 / [ sum(p_i^2) + (1-c)^2 / (f - n) ]
        f_l = 1 / [ sum(p_i^2) + m * p_n^2 ],    m = (1 - c) / p_n

    The upper bound assumes the unexplained mass is spread equally over the
    remaining founders; the lower bound assumes it is concentrated on ``m``
    founders each contributing ``p_n``. ``m`` is generally not an integer and
    the paper does not round it.

    HISTORY. This routine was made to refuse, because what
    stood here computed ``f_l`` and ``f_u`` from a single expression evaluated
    twice -- no input could make them differ -- and returned values below 1,
    which is impossible. PyPedal 2.0.4 is no better: it returns ``NaN`` for the
    lower bound. The refusal is now retired, the scheme having been read out of
    the paper rather than reconstructed.

    SHARES THE APPENDIX-B ENGINE. The bounds are a truncation of exactly the sequence
    Appendix B produces, not a second method, so this consumes
    :func:`boichard_marginal_contributions` and inherits its conventions and
    refusals unchanged. In particular **this routine refuses precisely when
    ``a_effective_ancestors_definite`` refuses** -- a bound computed on
    unresolved semantics is the invisibly-wrong answer the refusal existed to
    prevent.

    PRECONDITION. Both formulas are built on ``1 - c``, which is the mass left
    unexplained by the first ``n`` ancestors **only because article p.8 makes
    the contributions over all ancestors sum to one**. If the full sequence
    sums to ``S != 1`` the unexplained mass is ``S - c``, so neither bound is
    defined, and ``f_l <= f_a <= f_u`` has no middle term because ``f_a`` is
    just what :func:`a_effective_ancestors_definite` refuses to produce there.
    The full sequence is therefore validated with the same
    ``check_contribution_vector`` criterion the exact routine uses, before any
    of the arithmetic below. Nothing is renormalised or clipped: a reference
    population for which this implementation cannot produce a valid
    contribution sequence is refused, not quietly repaired.

    ENDPOINT AND DOMAIN HANDLING. Both formulas divide, so the degenerate cases
    are handled explicitly. These are algebraic safeguards derived from the
    published formulas, not new semantics:

    * when the unexplained mass is zero within tolerance the truncation has
      reached the exact answer, so ``f_l = f_u = f_a`` is returned and neither
      singular residual term is evaluated;
    * ``(f - n)`` is never divided by when ``n == f``;
    * ``p_n == 0`` with mass still unexplained is an inconsistent internal
      state and raises, rather than inventing a bound;
    * the mass test uses an explicit tolerance, never ``== 0.0``, because ``c``
      is a sum of many floating-point contributions.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : np.ndarray, optional
        Accepted for API compatibility and ignored; see
        :func:`a_effective_ancestors_definite`.
    gen : int, optional
        Generation to use as the population under study. Mutually exclusive
        with ``reference``.
    n : int, optional
        Number of ancestors to take before bounding. Truncated to the number of
        ancestors that actually have a positive contribution, since taking more
        is not defined. Default 25, as before.
    reference : iterable of int, keyword-only, optional
        The population under study stated directly. Same contract as in
        :func:`a_effective_founders_boichard`, and -- because this routine
        shares the Appendix-B engine and must refuse precisely where
        :func:`a_effective_ancestors_definite` refuses -- an explicitly
        supplied R is still subject to the antichain requirement.
    output : bool, optional, keyword-only
        If True (the default), write ``{filetag}_fa_boichard_indefinite_.dat``.
        If False, perform the calculation and return the same ``(f_l, f_u)``
        without writing that analysis file.

    Returns
    -------
    tuple of float
        ``(f_l, f_u)``, satisfying ``f_l <= f_a <= f_u``.
    """
    routine = 'a_effective_ancestors_indefinite'
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered a_effective_ancestors_indefinite()')

    if n < 1:
        raise PyPedalUsageError(
            '%s: n must be at least 1, not %r. The bounds are defined on the '
            'first n most important ancestors.' % (routine, n))

    explicit = _boichard_select_reference(pedobj, gen, reference, routine)
    if explicit is None:
        pyp_validate.require_generations(pedobj, routine)
    _boichard_require_topological(pedobj, routine)

    # `reference` is what the caller passed and may be None; `population` is
    # the R actually analysed.
    if explicit is None:
        population, gen, most_recent, gens = _boichard_reference_population(
            pedobj, gen, routine)
    else:
        population = explicit
    _boichard_require_antichain(pedobj, population, routine)

    _ids, _s, _d, phantoms, n_founders = _boichard_completed_arrays(pedobj)
    order = list(boichard_marginal_contributions(pedobj, population))
    if not order:
        raise PyPedalError(
            '%s: no ancestor has a positive contribution, so there is nothing '
            'to bound.' % routine)

    # p.8: "the marginal contributions over all ancestors sum to one". The
    # residual below is written `1.0 - c`, so that 1.0 IS this invariant,
    # hard-coded into the arithmetic. Unless it holds, (1-c) is not the
    # unexplained mass and neither published bound formula is defined -- and
    # f_l <= f_a <= f_u has no middle term, because f_a is exactly what
    # a_effective_ancestors_definite refuses on such a sequence.
    #
    # Checked on the FULL sequence, before truncation: the first n terms
    # legitimately sum to c < 1, which is the entire point of a bound. Same
    # validator, same criterion, same sequence as the exact routine, so the
    # two cannot disagree about which pedigrees they answer for.
    pyp_validate.check_contribution_vector(
        pedobj, [value for _, value in order],
        '%s: marginal contributions' % routine)

    # Taking more ancestors than exist is not defined; the paper's n indexes
    # into the selection sequence.
    taken = order[:min(n, len(order))]
    contributions = [value for _, value in taken]
    ssq = sum(value ** 2 for value in contributions)
    c = sum(contributions)
    p_n = contributions[-1]
    residual = 1.0 - c
    n_taken = len(taken)

    # `f` is the total number of founders, and comes from the same completion
    # the engine ran on rather than from a second, independent count. Under
    # After completion that distinction has teeth: NewAnimal.founder is 'y' only when
    # BOTH parents are unknown, so on a half-founder pedigree it would feed a
    # different f into f_u = 1/(ssq + (1-c)^2/(f-n)) than the sequence being
    # bounded was generated from. That is arithmetic, not just a guard.

    tol = 1e-12
    if residual <= tol:
        # The truncation has reached the exact answer. Neither residual term is
        # evaluated: (f - n) may be zero here and (1-c)/p_n is vacuous.
        exact = 1.0 / ssq if ssq > 0 else 0.0
        f_l = f_u = exact
    else:
        if n_taken >= n_founders:
            raise pyp_validate.PyPedalValidationError(
                '%s: %d ancestors have been taken but the pedigree has only %d '
                'founders, and %r of the contribution is still unexplained. The '
                'upper bound would divide by zero. The selection sequence is '
                'inconsistent with the founder count.'
                % (routine, n_taken, n_founders, residual))
        if p_n <= tol:
            raise pyp_validate.PyPedalValidationError(
                '%s: the %dth marginal contribution is %r while %r of the '
                'contribution is still unexplained, so m = (1-c)/p_n is '
                'undefined. The selection sequence is inconsistent.'
                % (routine, n_taken, p_n, residual))
        f_u = 1.0 / (ssq + residual * residual / float(n_founders - n_taken))
        m = residual / p_n
        f_l = 1.0 / (ssq + m * p_n * p_n)

    if pedobj.kw.get('messages') == 'verbose':
        print('=' * 60)
        print(f'animals:\t{len(pedobj.pedigree)}')
        print(f'ancestors taken:\t{n_taken}')
        print(f'f_l:\t\t{f_l:.3f}')
        print(f'f_u:\t\t{f_u:.3f}')
        print('=' * 60)

    if output:
        _write_boichard_indefinite_output(
            pedobj, n_taken, order, n_founders, c, taken, phantoms, f_l, f_u)

    # Both bounds are effective numbers of ancestors and obey the same
    # invariant; and the interval must not be inverted.
    pyp_validate.check_effective_number(
        pedobj, float(f_l), '%s: f_l' % routine,
        n_contributors=max(n_founders, 1))
    pyp_validate.check_effective_number(
        pedobj, float(f_u), '%s: f_u' % routine,
        n_contributors=max(n_founders, 1))
    if f_l > f_u + 1e-9:
        raise pyp_validate.PyPedalValidationError(
            '%s: the lower bound %r exceeds the upper bound %r.'
            % (routine, f_l, f_u))

    if pedobj.kw.get('debug_messages'):
        logger.info('Exited a_effective_ancestors_indefinite()')

    return float(f_l), float(f_u)


def a_coefficients(pedobj, a: Optional[np.ndarray] = None, method: str = 'nrm', output: bool = True) -> Dict[str, float]:
    """
    Write population average coefficients of inbreeding and relationship to a file, 
    as well as individual animal IDs and coefficients of inbreeding. For large pedigrees 
    that cannot be allocated due to memory restrictions, outputs -999.9 for all outputs.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : np.ndarray, optional
        A numerator relationship matrix.
    method : str, optional
        Determines which procedure should be called to build a relationship matrix (nrm|frm). Default is 'nrm'.
    output : bool, optional
        If True (the default), write the historical coefficient ``.dat`` files.
        If False, perform the calculation and return the same dictionary
        without writing those analysis files.

    Returns
    -------
    Dict[str, float]
        A dictionary of non-zero individual inbreeding coefficients.
    """
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered a_coefficients()')

    # Validate method
    if method not in ['nrm', 'frm']:
        method = 'nrm'

    # Retrieve or create the numerator relationship matrix
    if pedobj.kw.get('form_nrm'):
        a = pedobj.nrm.nrm
    elif a is None:
        try:
            if method == 'nrm':
                a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
            else:
                a = pyp_nrm.fast_a_matrix_r(pedobj.pedigree, pedobj.kw)
        except Exception:
            return {}

    lenped = len(pedobj.pedigree)
    f_sum = f_n = fnz_sum = fnz_n = r_sum = r_n = rnz_sum = rnz_n = 0.0
    individual_coi = {}

    # Calculate coefficients of inbreeding
    for row in range(lenped):
        coi = pyp_nrm._coi_from_matrix(a, row)
        f_sum += coi
        f_n += 1.0
        if coi > 0.0:
            fnz_sum += coi
            fnz_n += 1.0
            individual_coi[pedobj.pedigree[row].animalID] = coi

    f_avg = f_sum / f_n
    fnz_avg = fnz_sum / fnz_n if fnz_n > 0 else 0.0

    # Calculate coefficients of relationship
    for row in range(lenped):
        for col in range(row):
            rel_value = pyp_nrm._matrix_value(a, row, col)
            r_sum += rel_value
            r_n += 1.0
            if rel_value > 0.0:
                rnz_sum += rel_value
                rnz_n += 1.0

    r_avg = r_sum / r_n
    rnz_avg = rnz_sum / rnz_n if rnz_n > 0 else 0.0

    if output:
        _write_a_coefficients_output(
            pedobj, a, lenped, f_n, f_sum, f_avg, fnz_n, fnz_sum, fnz_avg,
            r_n, r_sum, r_avg, rnz_n, rnz_sum, rnz_avg)

    if pedobj.kw.get('debug_messages'):
        logger.info('Exited a_coefficients()')

    return individual_coi


def fast_a_coefficients(
    pedobj, 
    a: Optional[np.ndarray] = None, 
    method: str = 'nrm', 
    debug: bool = False, 
    storage: str = 'dense',
    output: bool = True,
) -> Dict[str, float]:
    """
    Writes population average coefficients of inbreeding and relationship to a file, 
    as well as individual animal IDs and coefficients of inbreeding. Returns a dictionary 
    of non-zero individual CoI.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    a : np.ndarray, optional
        A numerator relationship matrix. Default is None.
    method : str, optional
        Determines which procedure should be called to build a relationship matrix (nrm|frm). Default is 'nrm'.
    debug : bool, optional
        Print debugging messages if True, don't print otherwise. Default is False.
    storage : str, optional
        Use dense or sparse matrix storage. Default is 'dense'.
    output : bool, optional
        Combined with ``pedobj.kw['file_io']``. Files are written only when
        both ``output`` is True (the default) and ``file_io`` is True
        (``file_io`` itself defaults to True on a loaded pedigree).
        ``output=False`` suppresses the files even when ``file_io`` is True.
        Setting ``file_io`` to False still suppresses them even when
        ``output`` is True.

    Returns
    -------
    Dict[str, float]
        A dictionary of non-zero individual coefficients of inbreeding.
    """
    if debug or pedobj.kw.get('debug_messages', False):
        logger.info('Entered fast_a_coefficients()')

    # Validate method and storage options
    if method not in ['nrm', 'frm']:
        method = 'nrm'
    if storage not in ['dense', 'sparse']:
        storage = 'dense'

    # Retrieve or create the numerator relationship matrix
    if pedobj.kw.get('form_nrm', False):
        a = pedobj.nrm.nrm
    elif a is None:
        try:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw, method=storage)
        except Exception:
            logger.error("Matrix creation failed.")
            return {}

    lenped = len(pedobj.pedigree)
    f_sum = f_n = fnz_sum = fnz_n = r_sum = r_n = rnz_sum = rnz_n = 0.0
    individual_coi = {}

    # Calculate coefficients of inbreeding and relationship
    for row in range(lenped):
        coi = pyp_nrm._coi_from_matrix(a, row)
        f_sum += coi
        f_n += 1.0
        if coi > 0.0:
            fnz_sum += coi
            fnz_n += 1.0
            individual_coi[pedobj.pedigree[row].animalID] = coi

        for col in range(row):
            rel_value = pyp_nrm._matrix_value(a, row, col)
            r_sum += rel_value
            r_n += 1.0
            if rel_value > 0.0:
                rnz_sum += rel_value
                rnz_n += 1.0

    f_avg = f_sum / f_n if f_n > 0 else 0.0
    fnz_avg = fnz_sum / fnz_n if fnz_n > 0 else 0.0
    r_avg = r_sum / r_n if r_n > 0 else 0.0
    rnz_avg = rnz_sum / rnz_n if rnz_n > 0 else 0.0

    if output and pedobj.kw.get('file_io', False):
        _write_fast_a_coefficients_output(
            pedobj, a, lenped, f_n, f_sum, f_avg, fnz_n, fnz_sum, fnz_avg,
            r_n, r_sum, r_avg, rnz_n, rnz_sum, rnz_avg)

    if debug or pedobj.kw.get('debug_messages', False):
        logger.info('Exited fast_a_coefficients()')

    return individual_coi


def theoretical_ne_from_metadata(pedobj, output: bool = True) -> bool:
    """
    Computes the theoretical effective population size (N_e) based on the number 
    of sires and dams in a pedigree metadata object. Writes results to an output file
    when ``output`` is True.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    output : bool, optional
        If True (the default), write ``{filetag}_ne_from_metadata_.dat``. If
        False, perform the calculation without writing that analysis file.
        The return value remains True on success and False on failure.

    Returns
    -------
    bool
        True on success, False on failure.
    """
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered theoretical_ne_from_metadata()')

    try:
        # Retrieve the number of unique sires and dams
        ns = float(pedobj.metadata.num_unique_sires)
        nd = float(pedobj.metadata.num_unique_dams)

        # Calculate the theoretical effective population size
        ne = 1.0 / ((1.0 / (4.0 * ns)) + (1.0 / (4.0 * nd)))

        if output:
            _write_theoretical_ne_output(pedobj, ns, nd, ne)

        if pedobj.kw.get('debug_messages'):
            logger.info('Exited theoretical_ne_from_metadata()')

        return True

    except Exception as e:
        logger.error(f"Error in theoretical_ne_from_metadata: {e}")
        if pedobj.kw.get('debug_messages'):
            logger.info('Exited theoretical_ne_from_metadata() with failure')

        return False


def pedigree_completeness(pedobj, gens: int = 4) -> Dict[str, float]:
    """
    Computes the proportion of known ancestors in the pedigree of each animal in the population
    for a user-determined number of generations. Also computes mean completeness for all animals 
    and non-founders as summary statistics.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    gens : int, optional
        The number of generations to trace for completeness. Default is 4.

    Returns
    -------
    Dict[str, float]
        Dictionary of summary statistics.
    """
    if pedobj.kw.get('debug_messages'):
        logger.info('Entered pedigree_completeness()')

    l = len(pedobj.pedigree)
    c_summary = {}
    mp = str(pedobj.kw['missing_parent'])

    # Initialize summary statistics
    c_max, nf_c_max = 0.0, 0.0
    c_min, nf_c_min = 1.0, 1.0
    c_sum, nf_c_sum = 0.0, 0.0
    c_cnt, nf_c_cnt = 0, 0

    for i in range(l):
        animal = pedobj.pedigree[i]
        animalid = int(animal.animalID)
        sireid = int(animal.sireID)
        damid = int(animal.damID)

        # Debug information
        if pedobj.kw.get('debug_messages'):
            logger.debug(f"Animal: {animalid}")
            if str(sireid) != mp:
                logger.debug(f"\tSire: {sireid}")
            if str(damid) != mp:
                logger.debug(f"\tDam: {damid}")

        # Founders and missing parents
        if animal.founder == 'y' or (str(sireid) == mp and str(damid) == mp):
            _compl = 0.0
        else:
            # Calculate completeness using recursive function
            _sire_ped = pyp_nrm.recurse_pedigree_n(pedobj, sireid, [], gens - 1)
            _dam_ped = pyp_nrm.recurse_pedigree_n(pedobj, damid, [], gens - 1)
            n_max_ancestors = 2 * ((2 ** gens) - 1)
            _compl = (len(_sire_ped) + len(_dam_ped)) / float(n_max_ancestors)

        animal.pedcomp = _compl

        # All-animal summary statistics
        c_sum += _compl
        c_cnt += 1
        c_max = max(c_max, _compl)
        c_min = min(c_min, _compl)

        # Non-founder summary statistics
        if animal.founder != 'y':
            nf_c_sum += _compl
            nf_c_cnt += 1
            nf_c_max = max(nf_c_max, _compl)
            nf_c_min = min(nf_c_min, _compl)

    # Calculate summary metrics
    c_summary['sum'] = c_sum
    c_summary['n'] = c_cnt
    c_summary['min'] = c_min
    c_summary['max'] = c_max
    c_summary['range'] = c_max - c_min
    c_summary['average'] = c_sum / c_cnt if c_cnt > 0 else 0.0

    c_summary['nonfounder_sum'] = nf_c_sum
    c_summary['nonfounder_n'] = nf_c_cnt
    c_summary['nonfounder_min'] = nf_c_min
    c_summary['nonfounder_max'] = nf_c_max
    c_summary['nonfounder_range'] = nf_c_max - nf_c_min
    c_summary['nonfounder_average'] = nf_c_sum / nf_c_cnt if nf_c_cnt > 0 else 0.0

    # Write verbose output if enabled
    if pedobj.kw.get('messages') == 'verbose':
        verbose_summary = [
            '-' * 80,
            'Pedigree completeness summary statistics (all animals):',
            f"\tN: {c_cnt}",
            f"\tSum: {c_sum:.4f}",
            f"\tMean: {c_summary['average']:.4f}",
            f"\tMin: {c_min:.4f}",
            f"\tMax: {c_max:.4f}",
            f"\tRange: {c_summary['range']:.4f}",
            '-' * 80,
            'Pedigree completeness summary statistics (non-founders):',
            f"\tN: {nf_c_cnt}",
            f"\tSum: {nf_c_sum:.4f}",
            f"\tMean: {c_summary['nonfounder_average']:.4f}",
            f"\tMin: {nf_c_min:.4f}",
            f"\tMax: {nf_c_max:.4f}",
            f"\tRange: {c_summary['nonfounder_range']:.4f}",
            '-' * 80
        ]
        print("\n".join(verbose_summary))

    if pedobj.kw.get('debug_messages'):
        logger.info('Exited pedigree_completeness()')

    return c_summary


def common_ancestors(anim_a: int, anim_b: int, pedobj) -> List[int]:
    """
    Returns a list of the ancestors that two animals share in common.

    Parameters
    ----------
    anim_a : int
        The renumbered ID of the first animal.
    anim_b : int
        The renumbered ID of the second animal.
    pedobj : object
        A PyPedal pedigree object.

    Returns
    -------
    List[int]
        A list of animals related to both anim_a and anim_b.
    """
    if pedobj.kw.get('debug_messages'):
        logger.info(f'Entered common_ancestors() for animals {anim_a} and {anim_b}')

    try:
        # Get the lists of related animals for both input animals
        ped_a = related_animals(anim_a, pedobj)
        ped_b = related_animals(anim_b, pedobj)

        # Find the intersection of the two lists
        shared = list(set(ped_a) & set(ped_b))
        shared.sort()

        if pedobj.kw.get('debug_messages'):
            logger.debug(f"Ancestors of {anim_a}: {ped_a}")
            logger.debug(f"Ancestors of {anim_b}: {ped_b}")
            logger.debug(f"Shared ancestors: {shared}")

        return shared

    except Exception as e:
        logger.error(f"Error in common_ancestors: {e}")
        return []


def related_animals(anim: int, pedobj) -> List[int]:
    """
    Returns a list of the ancestors of an animal.

    Parameters
    ----------
    anim : int
        The renumbered ID of an animal.
    pedobj : object
        A PyPedal pedigree object.

    Returns
    -------
    List[int]
        A list of ancestors of the given animal.
    """
    if pedobj.kw.get('debug_messages'):
        logger.info(f"Entered related_animals() for animal {anim}")

    _ped = []
    try:
        _ped = pyp_nrm.recurse_pedigree_idonly(pedobj, anim, _ped)
    except Exception as e:
        logger.error(f"Error in related_animals: {e}")

    if pedobj.kw.get('debug_messages'):
        logger.info(f"Exited related_animals() with ancestors: {_ped}")

    return _ped


def relationship(anim_a, anim_b, pedobj, renumber=False):
    """
    Returns the coefficient of relationship for two animals, anim_a and anim_b.

    Parameters:
    ----------
    anim_a : int
        The renumbered ID of an animal, a.
    anim_b : int
        The renumbered ID of an animal, b.
    pedobj : object
        A PyPedal pedigree object.
    renumber : bool, optional
        Whether to renumber the pedigree if it has not been renumbered (default is False).

    Returns:
    -------
    float
        The coefficient of relationship between anim_a and anim_b.

    Raises:
    ------
    PyPedalUsageError
        If either ID is not a current/renumbered animalID in this pedigree.
        Unresolved IDs are not returned as 0.0; 0.0 remains the coefficient
        for a genuinely unrelated existing pair.
    """
    if not pedobj.kw['pedigree_is_renumbered'] and not renumber:
        if pedobj.kw['messages'] != 'quiet':
            print('[WARNING]: The pedigree you passed to pyp_metrics/relationship() is not renumbered; '
                  'this may result in incorrect calculations!')
        logger.warning('The pedigree you passed to pyp_metrics/relationship() is not renumbered; '
                        'this may result in incorrect calculations!')

    elif not pedobj.kw['pedigree_is_renumbered'] and renumber:
        if pedobj.kw['messages'] != 'quiet':
            print('[INFO]: Renumbering the pedigree in pyp_metrics/relationship().')
        logger.info('Renumbering the pedigree in pyp_metrics/relationship().')
        pedobj.kw['renumber'] = True
        pedobj.renumber()

    # Current/renumbered animalID domain. originalID is not translated:
    # 0.0 is a valid unrelated coefficient, so an unresolved ID must not
    # collapse to the same number (post-rc1 API safety candidate A).
    current_ids = {int(animal.animalID) for animal in pedobj.pedigree}

    def _require_current_id(label, anim):
        try:
            aid = int(anim)
        except (TypeError, ValueError) as exc:
            raise PyPedalUsageError(
                'relationship: %s animal ID %r is not an animal in this '
                'pedigree. relationship() requires existing current/'
                'renumbered animalID values.' % (label, anim)
            ) from exc
        if aid not in current_ids:
            raise PyPedalUsageError(
                'relationship: %s animal ID %r is not an animal in this '
                'pedigree. relationship() requires existing current/'
                'renumbered animalID values.' % (label, anim)
            )

    _require_current_id('first', anim_a)
    _require_current_id('second', anim_b)

    _r = 0.0  # Default relationship
    if anim_a == anim_b:
        return 1.0

    try:
        if pedobj.kw.get('form_nrm') and getattr(pedobj, 'nrm', None) is not None:
            if pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records:
                return pyp_nrm._matrix_value(pedobj.nrm.nrm, int(anim_a) - 1, int(anim_b) - 1)
    except Exception:
        pass

    try:
        _ped_a = pyp_nrm.recurse_pedigree(pedobj, anim_a, [])
        _ped_b = pyp_nrm.recurse_pedigree(pedobj, anim_b, [])
        _ped = []
        _seen = {}

        for _a in _ped_a:
            if _a.animalID not in _seen:
                _ped.append(_a)
                _seen[_a.animalID] = _a.animalID

        for _b in _ped_b:
            if _b.animalID not in _seen:
                _ped.append(_b)
                _seen[_b.animalID] = _b.animalID

        _tag = f"{pedobj.kw['filetag']}"
        _reord = [copy.deepcopy(animal) for animal in _ped]

        if pedobj.kw['slow_reorder']:
            _reord = pyp_utils.reorder(_reord, _tag, debug=pedobj.kw['debug_messages'],
                                       missingparent=pedobj.kw['missing_parent'])
        else:
            _reord = pyp_utils.fast_reorder(_reord, _tag,
                                            missingparent=pedobj.kw['missing_parent'])

        _s, _map = pyp_utils.renumber(_reord, _tag, returnmap=True,
                                      debug=pedobj.kw['debug_messages'],
                                      missingparent=pedobj.kw['missing_parent'],
                                      animaltype=pedobj.kw['animal_type'])

        _opts = copy.deepcopy(pedobj.kw)
        _opts['filetag'] = _tag

        if pedobj.kw['nrm_method'] == 'nrm':
            _a = pyp_nrm.fast_a_matrix(_s, _opts)
        else:
            _a = pyp_nrm.fast_a_matrix_r(_s, _opts)

        if _a is False:
            return 0.0
        _r = pyp_nrm._matrix_value(_a, _map[anim_a] - 1, _map[anim_b] - 1)
    except Exception as e:
        logger.warning(
            'Could not compute the relationship between animals %s and %s; defaulting to 0.0. Error: %s',
            anim_a, anim_b, str(e)
        )
        _r = 0.0
    return _r


def _mating_require_gens(gens, routine):
    """Accept only full-pedigree gens values (0 and -1)."""
    if gens not in (0, -1):
        raise PyPedalUsageError(
            '%s: gens=%r is not supported. PyPedal 4.0 accepts gens=0 or '
            'gens=-1 (both mean the full available pedigree, read-only). '
            'Truncated-generation approximations are not implemented.'
            % (routine, gens)
        )


def _mating_require_renumbered(pedobj, routine):
    """Mating COI uses the current/renumbered animalID domain."""
    if not pedobj.kw.get('pedigree_is_renumbered'):
        raise PyPedalUsageError(
            '%s requires a supported renumbered pedigree '
            '(pedigree_is_renumbered is false). Renumber first; this '
            'call did not mutate.' % routine
        )


def _mating_current_ids(pedobj):
    return {int(animal.animalID) for animal in pedobj.pedigree}


def _mating_require_current_id(value, current_ids, routine, label):
    """Require a current animalID. originalID and call names are not translated."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise PyPedalUsageError(
            '%s: %s ID %r is not a current/renumbered animalID in this '
            'pedigree. Mating COI does not translate originalID or call names.'
            % (routine, label, value)
        )
    aid = int(value)
    if aid not in current_ids:
        raise PyPedalUsageError(
            '%s: %s ID %r is not a current/renumbered animalID in this '
            'pedigree. Mating COI does not translate originalID or call names.'
            % (routine, label, value)
        )
    return aid


def _mating_attached_nrm(pedobj):
    """Return a usable attached NRM, or None.

    Non-None ``pedobj.nrm`` is the cache-validity signal. Shape is checked
    only after that: a matching shape does not revive an invalidated cache.
    """
    nrm_obj = getattr(pedobj, 'nrm', None)
    if nrm_obj is None:
        return None
    matrix = getattr(nrm_obj, 'nrm', None)
    if matrix is None:
        return None
    n = len(pedobj.pedigree)
    try:
        shape = matrix.shape
    except AttributeError:
        return None
    if len(shape) < 2 or int(shape[0]) != n or int(shape[1]) != n:
        return None
    metadata = getattr(pedobj, 'metadata', None)
    if metadata is not None:
        num_records = getattr(metadata, 'num_records', n)
        if int(num_records) != n:
            return None
    return matrix


def _local_additive_relationship(pedobj, anim_a, anim_b):
    """A_ij from an ancestor sub-pedigree, including a true A_ii diagonal.

    Uses the same copy/reorder/renumber/fast_a_matrix construction as
    :func:`relationship` for distinct animals, but never the self shortcut
    that returns 1.0. Operates only on copies; does not attach an NRM to
    ``pedobj`` and does not mutate the caller's pedigree.
    """
    _ped_a = pyp_nrm.recurse_pedigree(pedobj, anim_a, [])
    if anim_a == anim_b:
        _ped_b = []
    else:
        _ped_b = pyp_nrm.recurse_pedigree(pedobj, anim_b, [])
    _ped = []
    _seen = {}
    for _a in _ped_a:
        if _a.animalID not in _seen:
            _ped.append(_a)
            _seen[_a.animalID] = _a.animalID
    for _b in _ped_b:
        if _b.animalID not in _seen:
            _ped.append(_b)
            _seen[_b.animalID] = _b.animalID

    _tag = '%s' % (pedobj.kw['filetag'],)
    _reord = [copy.deepcopy(animal) for animal in _ped]
    _missing = pedobj.kw.get('missing_parent', 0)
    if pedobj.kw.get('slow_reorder'):
        _reord = pyp_utils.reorder(
            _reord, _tag, debug=pedobj.kw.get('debug_messages', False),
            missingparent=_missing,
        )
    else:
        _reord = pyp_utils.fast_reorder(_reord, _tag, missingparent=_missing)

    _s, _map = pyp_utils.renumber(
        _reord, _tag, returnmap=True,
        debug=pedobj.kw.get('debug_messages', False),
        missingparent=_missing,
        animaltype=pedobj.kw.get('animal_type', 'default'),
    )
    _opts = copy.deepcopy(pedobj.kw)
    _opts['filetag'] = _tag
    if pedobj.kw.get('nrm_method', 'nrm') == 'nrm':
        _a = pyp_nrm.fast_a_matrix(_s, _opts)
    else:
        _a = pyp_nrm.fast_a_matrix_r(_s, _opts)
    if _a is False:
        raise PyPedalUsageError(
            'mating_coi: could not form the additive relationship for '
            'animals %s and %s.' % (anim_a, anim_b)
        )
    return float(pyp_nrm._matrix_value(_a, _map[anim_a] - 1, _map[anim_b] - 1))


def mating_coi(anim_a, anim_b, pedobj, gens=0):
    """
    Prospective inbreeding coefficient of the offspring of two animals.

    This is a read-only calculation. It does not insert a hypothetical
    child, does not call ``addanimal`` / ``delanimal``, and does not
    form or attach a pedigree-wide NRM.

    IDs are current / renumbered ``animalID`` values, the same domain as
    :func:`relationship`. ``originalID`` is not translated. Call/display
    names are not identities.

    For distinct parents i and j the numerator relationship A_ij is halved:

        F_offspring(i x j) = A_ij / 2

    For self-mating, :func:`relationship` returns 1.0 rather than
    A_ii = 1 + F_i. This function uses the true diagonal instead:

        F_offspring(i x i) = A_ii / 2 = (1 + F_i) / 2

    ``gens`` may be 0 or -1 (both: full available pedigree). Any other
    value is a usage error; truncated-generation approximations are not
    implemented.

    Parameters
    ----------
    anim_a : int
        Current/renumbered animalID of the first prospective parent.
    anim_b : int
        Current/renumbered animalID of the second prospective parent.
    pedobj : object
        A PyPedal pedigree object. Must be in the supported renumbered state.
    gens : int, optional
        0 or -1 for a full-pedigree calculation. Other values are refused.

    Returns
    -------
    float
        Inbreeding coefficient of the prospective offspring.

    Raises
    ------
    PyPedalUsageError
        If the pedigree is not renumbered, an ID is not a current animalID,
        or ``gens`` is outside {0, -1}.
    """
    _mating_require_gens(gens, 'mating_coi')
    _mating_require_renumbered(pedobj, 'mating_coi')
    current_ids = _mating_current_ids(pedobj)
    a = _mating_require_current_id(anim_a, current_ids, 'mating_coi', 'first')
    b = _mating_require_current_id(anim_b, current_ids, 'mating_coi', 'second')

    matrix = _mating_attached_nrm(pedobj)
    if a == b:
        if matrix is not None:
            a_ii = pyp_nrm._matrix_value(matrix, a - 1, a - 1)
        else:
            a_ii = _local_additive_relationship(pedobj, a, a)
        return float(a_ii) / 2.0

    if matrix is not None:
        a_ij = pyp_nrm._matrix_value(matrix, a - 1, b - 1)
    else:
        a_ij = relationship(a, b, pedobj, renumber=False)
    return float(a_ij) / 2.0


def _mating_group_stats(values):
    """count/sum/min/max/range/mean, with None extrema when the set is empty."""
    if not values:
        return {
            'count': 0,
            'sum': 0.0,
            'min': None,
            'max': None,
            'range': None,
            'mean': None,
        }
    total = float(sum(values))
    lo = float(min(values))
    hi = float(max(values))
    return {
        'count': len(values),
        'sum': total,
        'min': lo,
        'max': hi,
        'range': hi - lo,
        'mean': total / float(len(values)),
    }


def _is_signed_int_token(token):
    if not isinstance(token, str) or not token:
        return False
    if token[0] in '+-':
        return token[1:].isdigit()
    return token.isdigit()


def _parse_mating_pair_item(item, names_flag, routine):
    """Return a two-element tuple from a modern pair or a legacy 'a_b' string."""
    if isinstance(item, bytes):
        raise PyPedalUsageError(
            '%s: each mating must be a two-ID pair or a numeric '
            "'sire_dam' string, not %r." % (routine, item)
        )
    if isinstance(item, str):
        parts = item.split('_')
        if len(parts) != 2 or parts[0] == '' or parts[1] == '':
            raise PyPedalUsageError(
                '%s: legacy mating string %r is not an unambiguous '
                "two-ID 'a_b' pair." % (routine, item)
            )
        if not names_flag:
            if not _is_signed_int_token(parts[0]) or not _is_signed_int_token(parts[1]):
                raise PyPedalUsageError(
                    '%s: legacy mating string %r is not an unambiguous '
                    "numeric 'a_b' pair." % (routine, item)
                )
            return int(parts[0]), int(parts[1])
        return parts[0], parts[1]
    try:
        seq = tuple(item)
    except TypeError:
        raise PyPedalUsageError(
            '%s: each mating must be a two-ID pair or a numeric '
            "'sire_dam' string, not %r." % (routine, item)
        )
    if len(seq) != 2:
        raise PyPedalUsageError(
            '%s: each mating must contain exactly two IDs, not %r.'
            % (routine, item)
        )
    return seq[0], seq[1]


def _resolve_string_identity(token, pedobj, current_ids, routine, label):
    """Resolve a unique string identity through namemap, never call names.

    On ASD pedigrees ``namemap`` stores the hashed originalID. The names=1
    adapter then maps that unique string identity to the current animalID
    through ``idmap``. Integer ``n`` call names are not identities.
    """
    mapped = None
    if token in pedobj.namemap:
        mapped = pedobj.namemap[token]
    elif (
        isinstance(token, numbers.Integral)
        and not isinstance(token, bool)
        and str(int(token)) in pedobj.namemap
    ):
        mapped = pedobj.namemap[str(int(token))]
    if mapped is None:
        raise PyPedalUsageError(
            '%s: %s identity %r is not a unique string identity in this '
            'pedigree string-ID map. names=1 does not search call names.'
            % (routine, label, token)
        )
    if mapped in current_ids:
        return int(mapped)
    try:
        mapped_key = mapped
        if mapped_key in pedobj.idmap:
            aid = int(pedobj.idmap[mapped_key])
            if aid in current_ids:
                return aid
    except (TypeError, ValueError):
        pass
    raise PyPedalUsageError(
        '%s: %s identity %r did not resolve to a current animalID.'
        % (routine, label, token)
    )


def mating_coi_group(matings, pedobj, names=0, gens=0):
    """
    Prospective offspring inbreeding for an explicit list of matings.

    Evaluates exactly the pairs supplied by the caller. It does not form a
    Cartesian product, does not select mates, and does not mutate the
    pedigree. Duplicate exact pairs are evaluated once, in first-occurrence
    order. Reversed pairs (A, B) and (B, A) remain distinct keys; their
    coefficients are equal.

    Modern input is an iterable of two-ID pairs::

        [(1, 2), (1, 6)]

    Legacy numeric strings such as ``["1_2", "1_6"]`` are accepted when
    parsing is unambiguous. Result keys are always current-animalID tuples,
    never underscore-delimited strings.

    ``names=1`` is deprecated compatibility behaviour. It resolves unique
    string identities through the pedigree string-ID map (``namemap``) and
    does **not** search the non-unique call/display name field ``n``.

    Parameters
    ----------
    matings : iterable
        Explicit proposed matings as two-ID sequences or legacy ``'a_b'``
        strings.
    pedobj : object
        A PyPedal pedigree object. Must be in the supported renumbered state.
    names : int, optional
        0 (default): current animalIDs. 1: deprecated unique string-ID lookup.
    gens : int, optional
        0 or -1 for a full-pedigree calculation. Other values are refused.

    Returns
    -------
    dict
        ``{'matings': {(a, b): F, ...}, 'metadata': {'all': ..., 'nonzero': ...}}``
        Metadata fields are count, sum, min, max, range, and mean. An empty
        set (including no nonzero matings) uses count=0, sum=0.0, and None
        for min, max, range, and mean.

    Raises
    ------
    PyPedalUsageError
        If input is malformed, an ID cannot be resolved, the pedigree is
        not renumbered, or ``gens`` is outside {0, -1}.
    """
    _mating_require_gens(gens, 'mating_coi_group')
    _mating_require_renumbered(pedobj, 'mating_coi_group')
    if names not in (0, 1, False, True):
        raise PyPedalUsageError(
            'mating_coi_group: names=%r is not supported; use 0 (current '
            'animalID) or 1 (deprecated unique string identity).' % (names,)
        )
    names_flag = 1 if names in (1, True) else 0
    if names_flag:
        pedformat = pedobj.kw.get('pedformat', '') or ''
        if 'A' not in pedformat:
            raise PyPedalUsageError(
                'mating_coi_group: names=1 resolves unique string identities '
                '(ASD) through the pedigree string-ID map. It does not search '
                'the non-unique call/display name field n.'
            )
        warnings.warn(
            'mating_coi_group(names=1) is deprecated. It resolves unique '
            'string identities through the pedigree string-ID map, not call '
            'names. Prefer current animalID pairs.',
            DeprecationWarning,
            stacklevel=2,
        )

    if matings is None:
        raise PyPedalUsageError(
            'mating_coi_group: matings must be an iterable of explicit pairs.'
        )
    try:
        items = list(matings)
    except TypeError:
        raise PyPedalUsageError(
            'mating_coi_group: matings must be an iterable of explicit pairs, '
            'not %r.' % (type(matings).__name__,)
        )

    current_ids = _mating_current_ids(pedobj)
    resolved = []
    seen = set()
    for item in items:
        raw_a, raw_b = _parse_mating_pair_item(item, names_flag, 'mating_coi_group')
        if names_flag:
            a = _resolve_string_identity(
                raw_a, pedobj, current_ids, 'mating_coi_group', 'first'
            )
            b = _resolve_string_identity(
                raw_b, pedobj, current_ids, 'mating_coi_group', 'second'
            )
        else:
            a = _mating_require_current_id(
                raw_a, current_ids, 'mating_coi_group', 'first'
            )
            b = _mating_require_current_id(
                raw_b, current_ids, 'mating_coi_group', 'second'
            )
        key = (a, b)
        if key in seen:
            continue
        seen.add(key)
        resolved.append(key)

    results = {}
    for a, b in resolved:
        results[(a, b)] = mating_coi(a, b, pedobj, gens=gens)

    values = list(results.values())
    nonzero = [value for value in values if value != 0.0]
    return {
        'matings': results,
        'metadata': {
            'all': _mating_group_stats(values),
            'nonzero': _mating_group_stats(nonzero),
        },
    }


class _GeneDropPlan:
    """
    Everything ``effective_founder_genomes`` needs to drop genes, resolved once
    before any replicate runs.

    Built as local state on purpose. The routine used to
    write each replicate's sampled alleles onto ``NewAnimal.alleles`` on the
    caller's own records, which made a second call start from the first call's
    leftovers and left simulation scratch visible to every downstream reader of
    that field.

    Attributes
    ----------
    parents : list of (int or None, int or None)
        Parent POSITIONS, not IDs. The routine used to
        index ``pedigree[int(animal.sireID) - 1]``, which silently became
        ``pedigree[-1]`` for the missing-parent sentinel 0 and adopted the last
        animal in the list -- sometimes the animal itself -- as a parent.
    founder_genes : list of (str, str) or None
        The two founder genes of an animal with both parents unknown.
    slot_genes : list of (str or None, str or None)
        The single unique gene standing in for an unknown parental slot.
    n_founder_genes : int
        The ``2f`` of equation 2, counting **two** conceptual genes per phantom
        founder even though only one is ever materialised -- see
        :func:`effective_founder_genomes`.
    """

    __slots__ = ('parents', 'founder_genes', 'slot_genes', 'n_founder_genes',
                 'n_true_founders', 'n_slots', 'reference')

    def __init__(self, parents, founder_genes, slot_genes, n_true_founders,
                 n_slots, reference):
        self.parents = parents
        self.founder_genes = founder_genes
        self.slot_genes = slot_genes
        self.n_true_founders = n_true_founders
        self.n_slots = n_slots
        self.n_founder_genes = 2 * (n_true_founders + n_slots)
        self.reference = reference


def _build_gene_drop_plan(pedobj, most_recent, routine):
    """
    Resolve parents to positions, name the founder genes, and refuse a pedigree
    the simulation cannot legally traverse.

    Validation happens ONCE, before the first replicate and before any report
    file is opened, so a refused calculation never leaves a plausible-looking
    partial report behind.

    ``NewPedigree.load()`` guarantees a canonically ordered, renumbered pedigree,
    but a caller can hand this routine anything -- ``renumber``
    and ``reorder`` are options, and ``pedigree_is_renumbered=True`` disables
    the checks that would otherwise catch a hand-built list. The structural
    conditions are therefore verified here rather than assumed.
    """
    pedigree = pedobj.pedigree
    if not pedigree:
        raise PyPedalValidationError(
            '%s: the pedigree is empty, so there are no founder genes to drop '
            'and N_g is undefined.' % routine)

    missing = str(pedobj.kw.get('missing_parent', 0))
    position = {}
    for index, animal in enumerate(pedigree):
        key = str(animal.animalID)
        if key in position:
            raise PyPedalPedigreeStructureError(
                '%s: animal ID %s appears more than once, so a parent '
                'reference to it is ambiguous and gene dropping cannot '
                'proceed.' % (routine, key), animals=[animal.animalID])
        position[key] = index

    parents, founder_genes, slot_genes = [], [], []
    labels = set()
    n_true_founders = n_slots = 0

    for index, animal in enumerate(pedigree):
        padded = str(getattr(animal, 'paddedID', None) or animal.animalID)
        resolved, slots = [], []
        for side, parent_id in (('sire', animal.sireID), ('dam', animal.damID)):
            if str(parent_id) == missing:
                resolved.append(None)
                slots.append(None)
                continue
            if str(parent_id) == str(animal.animalID):
                raise PyPedalPedigreeStructureError(
                    '%s: animal %s is recorded as its own %s.'
                    % (routine, animal.animalID, side), animals=[animal.animalID])
            parent_index = position.get(str(parent_id))
            if parent_index is None:
                raise PyPedalPedigreeStructureError(
                    '%s: animal %s names %s %s, which has no record in this '
                    'pedigree. A known parent is not the same thing as a '
                    'missing one and will not be treated as one.'
                    % (routine, animal.animalID, side, parent_id),
                    animals=[animal.animalID])
            if parent_index >= index:
                raise PyPedalPedigreeStructureError(
                    '%s: animal %s precedes its %s %s, so the parent has no '
                    'genotype yet when the offspring samples from it. Gene '
                    'dropping needs parents before offspring; reorder the '
                    'pedigree.' % (routine, animal.animalID, side, parent_id),
                    animals=[animal.animalID])
            resolved.append(parent_index)
            slots.append(None)

        if resolved[0] is None and resolved[1] is None:
            # An ordinary founder: two unique founder genes, no transmission.
            genes = ('%s__1' % padded, '%s__2' % padded)
            founder_genes.append(genes)
            labels.update(genes)
            n_true_founders += 1
            slot_genes.append((None, None))
        else:
            founder_genes.append(None)
            # Each unknown parental SLOT is its own
            # phantom founder -- Lacy p.113, Boichard Appendix A step 4, and
            # Baumung et al. (2015) p.102, which states the gene-drop rule
            # outright. FG-10 proved that materialising ONE unique gene in the
            # slot is exactly distributionally equivalent to Baumung's dummy
            # founder carrying two: the phantom has exactly one offspring, so
            # its untransmitted gene has frequency zero in R and contributes
            # nothing to SUM f_k^2. The conceptual founder count still counts
            # both -- see _GeneDropPlan.n_founder_genes.
            slot = []
            for side, parent_index in zip(('s', 'd'), resolved):
                if parent_index is None:
                    label = '%s__%s' % (padded, side)
                    slot.append(label)
                    labels.add(label)
                    n_slots += 1
                else:
                    slot.append(None)
            slot_genes.append(tuple(slot))
        parents.append((resolved[0], resolved[1]))

    if len(labels) != 2 * n_true_founders + n_slots:
        raise PyPedalInternalError(
            '%s: founder gene labels are not unique. Two animals produced the '
            'same padded ID, so distinct founder genes would be counted as '
            'one.' % routine)

    reference = [i for i, animal in enumerate(pedigree) if animal.gen == most_recent]
    if not reference:
        raise PyPedalValidationError(
            '%s: no animal belongs to generation %r, so the population under '
            'study is empty and N_g is undefined.' % (routine, most_recent))

    return _GeneDropPlan(parents, founder_genes, slot_genes, n_true_founders,
                         n_slots, reference)


def effective_founder_genomes(pedobj, rounds=10, chrometype='autosome',
                              heterogametic='m', quiet=False, *, seed=None,
                              output=True):
    """
    Estimate N_g, the effective number of founder genomes still present in the
    population under study, by replicated Mendelian gene dropping.

    Implements **equation 2 of Boichard, Maignel & Verrier (1997)**, *Genet Sel
    Evol* 29:5-23, article p.11::

        N_g = 1 / (2 * SUM_{k=1..2f} f_k^2)

    where ``f_k`` is the realised frequency of founder gene *k*, obtained by
    gene counting in the population under study after simulating segregation
    through the pedigree (MacCluer *et al.*, 1986). Each replicate produces one
    N_g and the routine returns their **arithmetic mean**, so ``rounds``
    genuinely buys precision.

    **This is not Lacy's founder genome equivalent.** Lacy (1989) p.115 defines

        f_g = 1 / SUM (p_i^2 / r_i)

    devaluing each founder by ``r_i``, the expected proportion of its alleles
    *retained* -- a survival probability, not a frequency. The two are different
    functionals and do not agree: on ``new_lacy.ped``, which is Lacy's own
    Appendix A worked example, the published f_g is 2.18 while the exact N_g is
    1.842. Boichard's parenthetical identification of the two is a naming remark
    rather than a derivation. PyPedal provides no f_g routine.

    The population under study is selected, as it always has been, by the
    ``g`` pedigree-format column: the animals whose ``gen`` equals the
    numerically largest generation label. Founders belonging to it contribute
    their gene copies like any other member. An animal with one unknown parent
    is a descendant of its known parent and of a distinct phantom founder
    contributing one unique gene to that slot.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object. It is **not** modified.
    rounds : int
        Number of Monte Carlo replicates. Must be >= 1.
    chrometype : str
        ``'autosome'`` is the supported domain. Any other value, including
        ``'sex'``, raises :class:`~PyPedal.pyp_errors.PyPedalUsageError`.
    heterogametic : str
        Accepted for API compatibility and irrelevant to an autosomal
        calculation, which is the only kind implemented.
    quiet : bool
        If False, print the summary when ``kw['messages'] == 'verbose'``.
    seed : int, optional
        Seed for this call's own local generator. With a seed the result is
        reproducible; without one it is not. Either way the process-global
        ``random`` and ``numpy.random`` states are left untouched.
    output : bool
        Write ``<filetag>_gene_drop.out``. Default True, as it has always been.

    Returns
    -------
    float
        Mean N_g over ``rounds`` replicates. Bounded below by 0.5, which is
        attained exactly when one founder gene is fixed in the population under
        study -- a legitimate result, not an error.

    Raises
    ------
    PyPedalUsageError
        ``rounds`` < 1, or ``chrometype`` is not ``'autosome'``.
    PyPedalPedigreeStructureError
        Duplicate animal ID, self-parent, a known parent with no record, or an
        offspring preceding its parent.
    PyPedalValidationError
        Empty pedigree, empty population under study, or a computed value
        outside its mathematically possible range.
    """
    logger.info('Entered effective_founder_genomes()')

    if rounds < 1:
        raise PyPedalUsageError(
            "effective_founder_genomes: rounds must be a positive integer, got "
            "%r. Each round is one Monte Carlo replicate of the gene drop, so "
            "there is no meaningful calculation with fewer than one. "
            "(This used to be coerced to 1 silently.)" % (rounds,))

    if chrometype != 'autosome':
        # 'sex' was accepted and then never read, so the
        # routine returned the AUTOSOMAL N_g to a caller who had asked for a
        # sex chromosome. Autosome is the supported 4.0 domain; any other
        # value is outside that domain (not a stub feature).
        raise PyPedalUsageError(
            "effective_founder_genomes supports chrometype='autosome' only "
            f"(got {chrometype!r}). Sex-chromosome gene dropping is outside "
            "the PyPedal 4.0 domain."
        )

    if heterogametic not in ['m', 'f']:
        logger.warning(f"Unrecognized heterogametic value '{heterogametic}' in effective_founder_genomes(). Defaulting to 'm'.")
        heterogametic = 'm'

    gens = [animal.gen for animal in pedobj.pedigree]
    # This was `sorted(set(gens))[-1]`, which orders
    # `NewAnimal.gen` -- a string -- lexicographically, so a pedigree declaring
    # ten or more generations analysed generation '9' rather than '12'. That is
    # already repaired for the three Boichard routines by this
    # helper and missed here.
    most_recent = _most_recent_generation(gens, 'effective_founder_genomes')
    plan = _build_gene_drop_plan(pedobj, most_recent, 'effective_founder_genomes')

    summary_freqs = {
        'allele_count': 0,
        'distinct_allele_count': 0,
        'distinct_alleles': {},
        'n_g': 0.0
    }

    # This was an argument-less `random.seed()` once per
    # round, which reseeded the PROCESS-GLOBAL generator from OS entropy: the
    # result was irreproducible and any caller-side seeding was destroyed.
    # PyPedal 2.0.4 had the same defect against a different generator -- its
    # `from numpy import random` shadowed the standard library import, so the
    # legacy routine consumed and reset NumPy's global state instead (D6a). A
    # local generator supersedes both: neither global stream is touched.
    rng = random.Random(seed)

    outputfile = f"{pedobj.kw['filetag']}_gene_drop.out"
    myline = "=" * 80
    myline2 = "*" * 80

    for r in range(rounds):
        allele_freqs = {}

        # One replicate's genotypes, LOCAL to this call.
        drawn = [None] * len(pedobj.pedigree)
        for index in range(len(pedobj.pedigree)):
            genes = plan.founder_genes[index]
            if genes is None:
                sire_index, dam_index = plan.parents[index]
                sire_slot, dam_slot = plan.slot_genes[index]
                sire_gene = (sire_slot if sire_index is None
                             else rng.choice(drawn[sire_index]))
                dam_gene = (dam_slot if dam_index is None
                            else rng.choice(drawn[dam_index]))
                genes = (sire_gene, dam_gene)
            drawn[index] = genes

        # The `continue` that skipped a founder's
        # TRANSMISSION also skipped its TALLY, so a founder belonging to the
        # population under study contributed none of its gene copies to the
        # eq. 2 frequencies. Those are two different questions: a founder needs
        # no transmission because it already carries two unique founder genes,
        # but if it is a member of R both copies must be counted. Boichard
        # App. A/B step 1 lets the analyst put anything in R, and Lacy's okapi
        # and Goeldi's-monkey populations contain wild-caught founders outright.
        for index in plan.reference:
            for allele in drawn[index]:
                allele_freqs[allele] = allele_freqs.get(allele, 0) + 1

        nalleles = len(allele_freqs)
        allelecount = sum(allele_freqs.values())
        if not allelecount:
            raise PyPedalValidationError(
                'effective_founder_genomes: the population under study '
                '(generation %r) contains no countable gene copies, so '
                'SUM f_k^2 and therefore N_g are undefined.' % (most_recent,))
        _ngs = sum((freq / allelecount) ** 2 for freq in allele_freqs.values())
        _ng = 1.0 / (2.0 * _ngs)

        summary_freqs['allele_count'] += allelecount
        summary_freqs['distinct_allele_count'] += nalleles
        for allele, freq in allele_freqs.items():
            summary_freqs['distinct_alleles'][allele] = summary_freqs['distinct_alleles'].get(allele, 0) + freq
        # This was `(_ng + n_g) / 2.0`, an exponentially
        # weighted moving average whose weights are 2**-(R-r). They sum to 1, so
        # the estimator was unbiased, but SUM w^2 -> 1/3: the effective sample
        # size was 3 no matter how large `rounds` was, and the documented
        # parameter bought no precision past about three replicates. Boichard
        # p.11 replicates "to obtain an accurate estimate", and his published
        # Table II Family 2 value of 1.1 is the mean of the PER-REPLICATE N_g
        # (1.10625) rather than 1/(2 * mean SUM f_k^2), which is 1.0 exactly.
        # So the per-replicate inversion above is right and only the averaging
        # was wrong. Online arithmetic mean, no list of every replicate:
        summary_freqs['n_g'] += (_ng - summary_freqs['n_g']) / (r + 1)

        if not output:
            continue
        with open(outputfile, 'a' if r > 0 else 'w') as dout:
            if r == 0:
                dout.write(f"# FILE: {outputfile}\n")
                dout.write(f"# Results from {rounds}-round PyPedal gene-drop simulation.\n")
            dout.write(f"{myline}\n")
            dout.write(f"Allele frequency data from gene drop simulation, round {r + 1}\n")
            dout.write(f"\tNumber of distinct alleles: {nalleles}\n")
            dout.write(f"\tNumber of alleles in latest generation: {allelecount}\n")
            for allele, freq in allele_freqs.items():
                freq_ratio = freq / allelecount
                dout.write(f"\t\tAllele {allele}: {freq_ratio:.4f} ({freq_ratio ** 2:.4f})\n")
            dout.write(f"\tEffective number of founder genomes: {_ng:.4f}\n")

    summary_stats = {
        'allele_count': summary_freqs['allele_count'] / rounds,
        'distinct_allele_count': summary_freqs['distinct_allele_count'] / rounds,
        'distinct_alleles': {
            allele: freq / sum(summary_freqs['distinct_alleles'].values())
            for allele, freq in summary_freqs['distinct_alleles'].items()
        },
        'n_g': summary_freqs['n_g']
    }

    if pedobj.kw['messages'] == 'verbose' and not quiet:
        print(myline2)
        print(f"Summary statistics from {rounds}-round gene-drop simulation")
        print(f"\tNumber of distinct founder alleles: {plan.n_founder_genes}")
        print(f"\tMean allele count in latest generation: {summary_stats['allele_count']:.4f}")
        print(f"\tMean number of distinct alleles in latest generation: {summary_stats['distinct_allele_count']:.4f}")
        print("\tFrequency of distinct alleles sampled:")
        for allele, freq in summary_stats['distinct_alleles'].items():
            print(f"\t\tAllele {allele}: {freq:.4f} ({freq ** 2:.4f})")
        print(f"\tMean effective number of founder genomes: {summary_stats['n_g']:.4f}")
        print(myline2)

    if output:
        _write_gene_drop_summary(outputfile, myline2, rounds, plan, summary_stats)

    # SUM f_k = 1 forces SUM f_k^2 <= 1, so N_g = 1/(2 SUM f_k^2) >= 1/2, and
    # the mean of values each >= 1/2 is >= 1/2. Equality is FIXATION of a single
    # founder gene, which is a legitimate gene-drop outcome -- this is
    # deliberately not the `>= 1` a reader might expect. A violation could only
    # come from a defect in this routine, hence PyPedalValidationError.
    if summary_stats['n_g'] < 0.5 - 1e-9:
        raise PyPedalValidationError(
            'effective_founder_genomes: computed N_g = %r, below the '
            'mathematically attainable minimum of 0.5. This is a PyPedal '
            'defect, not a property of the data.' % (summary_stats['n_g'],))

    logger.info('Exited effective_founder_genomes()')
    return summary_stats['n_g']


def _write_gene_drop_summary(outputfile, myline2, rounds, plan, summary_stats):
    """The trailing summary block of ``<filetag>_gene_drop.out``."""
    with open(outputfile, 'a') as dout:
        dout.write(f"{myline2}\n")
        dout.write(f"Summary statistics from {rounds}-round gene-drop simulation\n")
        dout.write(f"\tNumber of distinct founder alleles: {plan.n_founder_genes}\n")
        dout.write(f"\tMean allele count in latest generation: {summary_stats['allele_count']:.4f}\n")
        dout.write(f"\tMean number of distinct alleles in latest generation: {summary_stats['distinct_allele_count']:.4f}\n")
        dout.write("\tFrequency of distinct alleles sampled:\n")
        for allele, freq in summary_stats['distinct_alleles'].items():
            dout.write(f"\t\tAllele {allele}: {freq:.4f} ({freq ** 2:.4f})\n")
        dout.write(f"\tMean effective number of founder genomes: {summary_stats['n_g']:.4f}\n")


def generation_intervals(pedobj, units='y'):
    """
    Computes the average age of parents at the time of birth of their first (oldest) offspring.
    Average ages are computed for each of four paths: sire-son, sire-daughter, dam-son, and dam-daughter.
    An overall mean is also computed.

    Parameters:
    ----------
    pedobj : object
        A PyPedal pedigree object.
    units : str
        Units in which the generation lengths should be returned ('y' for years).

    Returns:
    -------
    dict
        A dictionary containing the five average ages (for each path and overall).
    """
    logger.info('Entered generation_intervals()')

    sire_son = {}
    sire_dau = {}
    dam_son = {}
    dam_dau = {}
    n_unks = 0

    if not pedobj.kw['set_offspring']:
        pyp_utils.assign_offspring(pedobj)
        pedobj.kw['set_offspring'] = 1

    for m in pedobj.pedigree:
        if pedobj.kw['debug_messages']:
            if m.sons or m.daus:
                print(f"Animal {m.animalID} has sex {m.sex}")

        if m.sex.lower() == 'u':
            n_unks += 1

        oldest_son = None
        oldest_dau = None

        # Find the oldest son with a known recorded year
        for s in m.sons:
            s = int(s)
            child_year = pyp_chronology.recorded_year(pedobj.pedigree[s - 1])
            if child_year is None:
                continue
            if oldest_son is None or child_year < pyp_chronology.recorded_year(
                pedobj.pedigree[oldest_son - 1]
            ):
                oldest_son = s

        # Find the oldest daughter with a known recorded year
        for d in m.daus:
            d = int(d)
            child_year = pyp_chronology.recorded_year(pedobj.pedigree[d - 1])
            if child_year is None:
                continue
            if oldest_dau is None or child_year < pyp_chronology.recorded_year(
                pedobj.pedigree[oldest_dau - 1]
            ):
                oldest_dau = d

        # Assign paths to dictionaries
        if m.sex.lower() == 'm':  # Male
            if oldest_son:
                sire_son[m.animalID] = oldest_son
            if oldest_dau:
                sire_dau[m.animalID] = oldest_dau
        elif m.sex.lower() == 'f':  # Female
            if oldest_son:
                dam_son[m.animalID] = oldest_son
            if oldest_dau:
                dam_dau[m.animalID] = oldest_dau

    if pedobj.kw['messages'] == 'verbose' and n_unks > 0:
        print(f"[MESSAGE]: {n_unks} of {len(pedobj.pedigree)} animals in the pedigree were of unknown sex and excluded.")
        print("Paths:")
        print(f"\tSire-Son: {sire_son}")
        print(f"\tSire-Dau: {sire_dau}")
        print(f"\tDam-Son: {dam_son}")
        print(f"\tDam-Dau: {dam_dau}")

    # Compute generation lengths from pairs with known recorded years only.
    def compute_mean_age(path_dict):
        intervals = []
        for k, v in path_dict.items():
            parent_year = pyp_chronology.recorded_year(pedobj.pedigree[k - 1])
            child_year = pyp_chronology.recorded_year(pedobj.pedigree[v - 1])
            if parent_year is None or child_year is None:
                continue
            intervals.append(child_year - parent_year)
        if not intervals:
            return None
        return sum(intervals) / len(intervals)

    ss_mean = compute_mean_age(sire_son)
    sd_mean = compute_mean_age(sire_dau)
    ds_mean = compute_mean_age(dam_son)
    dd_mean = compute_mean_age(dam_dau)
    defined = [value for value in (ss_mean, sd_mean, ds_mean, dd_mean) if value is not None]
    overall_mean = (sum(defined) / len(defined)) if defined else None

    genlens = {
        'ss': ss_mean,
        'sd': sd_mean,
        'ds': ds_mean,
        'dd': dd_mean,
        'mean': overall_mean
    }

    if pedobj.kw['messages'] == 'verbose':
        print("Means:")
        print(f"\tSire-Son: {ss_mean}")
        print(f"\tSire-Dau: {sd_mean}")
        print(f"\tDam-Son: {ds_mean}")
        print(f"\tDam-Dau: {dd_mean}")
        print(f"\tOverall: {overall_mean}")

    logger.info('Exited generation_intervals()')
    return genlens


def generation_intervals_all(pedobj, units='y'):
    """
    Computes the average age of parents at the time of birth of their offspring. 
    The computation is made using birth years for all known offspring of sires and dams, 
    which implies discrete generations. Average ages are computed for each of four paths: 
    sire-son, sire-daughter, dam-son, and dam-daughter. An overall mean is computed as well.

    Parameters:
    ----------
    pedobj : object
        A PyPedal pedigree object.
    units : str
        Units in which the generation lengths should be returned ('y' for years).

    Returns:
    -------
    dict
        A dictionary containing the average ages for each path and overall mean.
    """
    logger.info('Entered generation_intervals_all')

    sire_son = {}
    sire_dau = {}
    dam_son = {}
    dam_dau = {}
    n_unks = 0

    if not pedobj.kw['set_offspring']:
        pyp_utils.assign_offspring(pedobj)
        pedobj.kw['set_offspring'] = 1

    for m in pedobj.pedigree:
        if pedobj.kw['debug_messages']:
            print(f"Animal {m.animalID} has sex {m.sex}")

        if m.sex.lower() == 'u':
            n_unks += 1

        if pedobj.kw['debug_messages']:
            print(f"\tAnimal: {m.animalID}")
            print(f"\t\tsons: {m.sons}")
            print(f"\t\tdaus: {m.daus}")

        # Add sire-offspring pairs with known recorded years
        if m.sex.lower() == 'm':
            for s in m.sons:
                s = int(s)
                if (
                    pyp_chronology.recorded_year(m) is None
                    or pyp_chronology.recorded_year(pedobj.pedigree[s - 1]) is None
                ):
                    continue
                if pedobj.kw['debug_messages']:
                    print(f"\tAdding sire-son pair {m.animalID}-{s} to sire_son")
                sire_son[m.animalID] = s

            for d in m.daus:
                d = int(d)
                if (
                    pyp_chronology.recorded_year(m) is None
                    or pyp_chronology.recorded_year(pedobj.pedigree[d - 1]) is None
                ):
                    continue
                if pedobj.kw['debug_messages']:
                    print(f"\tAdding sire-daughter pair {m.animalID}-{d} to sire_dau")
                sire_dau[m.animalID] = d

        # Add dam-offspring pairs with known recorded years
        if m.sex.lower() == 'f':
            for s in m.sons:
                s = int(s)
                if (
                    pyp_chronology.recorded_year(m) is None
                    or pyp_chronology.recorded_year(pedobj.pedigree[s - 1]) is None
                ):
                    continue
                if pedobj.kw['debug_messages']:
                    print(f"\tAdding dam-son pair {m.animalID}-{s} to dam_son")
                dam_son[m.animalID] = s

            for d in m.daus:
                d = int(d)
                if (
                    pyp_chronology.recorded_year(m) is None
                    or pyp_chronology.recorded_year(pedobj.pedigree[d - 1]) is None
                ):
                    continue
                if pedobj.kw['debug_messages']:
                    print(f"\tAdding dam-daughter pair {m.animalID}-{d} to dam_dau")
                dam_dau[m.animalID] = d

    if pedobj.kw['messages'] == 'verbose':
        if n_unks > 0:
            print(
                f"[MESSAGE]: {n_unks} of {len(pedobj.pedigree)} animals in the pedigree were of unknown sex and excluded from calculations."
            )
        print("Paths:")
        print(f"\tSire-Son: {sire_son}")
        print(f"\tSire-Dau: {sire_dau}")
        print(f"\tDam-Son: {dam_son}")
        print(f"\tDam-Dau: {dam_dau}")

    # Compute generation lengths from pairs with known recorded years only.
    def compute_mean_age(path_dict):
        intervals = []
        for k, v in path_dict.items():
            parent_year = pyp_chronology.recorded_year(pedobj.pedigree[k - 1])
            child_year = pyp_chronology.recorded_year(pedobj.pedigree[v - 1])
            if parent_year is None or child_year is None:
                continue
            intervals.append(child_year - parent_year)
        if not intervals:
            return None
        return sum(intervals) / len(intervals)

    ss_mean = compute_mean_age(sire_son)
    sd_mean = compute_mean_age(sire_dau)
    ds_mean = compute_mean_age(dam_son)
    dd_mean = compute_mean_age(dam_dau)
    defined = [value for value in (ss_mean, sd_mean, ds_mean, dd_mean) if value is not None]
    overall_mean = (sum(defined) / len(defined)) if defined else None

    genlens = {
        'ss': ss_mean,
        'sd': sd_mean,
        'ds': ds_mean,
        'dd': dd_mean,
        'mean': overall_mean,
    }

    if pedobj.kw['messages'] == 'verbose':
        print("Means:")
        print(f"\tSire-Son: {ss_mean}")
        print(f"\tSire-Dau: {sd_mean}")
        print(f"\tDam-Son: {ds_mean}")
        print(f"\tDam-Dau: {dd_mean}")
        print(f"\tOverall: {overall_mean}")

    logger.info('Exited generation_intervals_all')
    return genlens


def founder_descendants(pedobj):
    """
    founder_descendants() returns a dictionary containing a list of descendants of
    each founder in the pedigree.

    Parameters:
    ----------
    pedobj : object
        An instance of a PyPedal NewPedigree object.

    Returns:
    -------
    dict
        A dictionary containing a list of descendants for each founder in the pedigree.
    """
    logger.info('Entered founder_descendants()')

    founder_peds = {}
    for f in pedobj.metadata.unique_founder_list:
        _desc = descendants(_current_animal_id(pedobj, f), pedobj, {})
        founder_peds[f] = _desc

    logger.info('Exited founder_descendants()')
    return founder_peds


def descendants(anid, pedobj, _desc):
    """
    descendants() uses pedigree metadata to walk a pedigree and return a list of all
    of the descendants of a given animal.

    Parameters:
    ----------
    anid : int
        An animal ID.
    pedobj : object
        A PyPedal pedigree object containing Animal objects.
    _desc : dict
        A dictionary of descendants of the given animal ID.

    Returns:
    -------
    dict
        A dictionary containing all descendants of the given animal ID.
    """
    logger.info('Entered descendants()')

    stack = [int(anid)]
    while stack:
        current = stack.pop()
        animal = pedobj.pedigree[current - 1]
        for child_id in list(animal.unks.keys()) + list(animal.sons.keys()) + list(animal.daus.keys()):
            if child_id not in _desc:
                _desc[child_id] = child_id
                stack.append(int(child_id))

    logger.info('Exited descendants()')
    return _desc


# Probability that a given one of a parent's two alleles is the one transmitted.
# This is Mendelian segregation and is 0.5 by definition -- it is not a tunable
# parameter and not an allele frequency.
MENDELIAN_TRANSMISSION_P = 0.5


def dropped_ancestral_inbreeding(pedobj, rounds=100, loci=100, frequency=None, seed=5048665):
    """
    dropped_ancestral_inbreeding() uses a gene dropping approach to calculate
    ancestral inbreeding, the probability of an individual inheriting an allele
    that has undergone inbreeding in the past at least once.

    Parameters:
    ----------
    pedobj : object
        A PyPedal pedigree object.
    rounds : int, optional
        Number of times to simulate segregation through the entire pedigree (default is 100).
    loci : int, optional
        Number of biallelic, unlinked loci to simulate (default is 100).
    frequency : float, optional
        DEPRECATED and ignored. It is accepted so that existing calls keep
        working, and passing it raises a DeprecationWarning.

        It was documented as a minor allele frequency and then used as the
        probability of transmitting a parent's first allele. That is not what a
        minor allele frequency is, and in this routine an allele frequency has
        no mathematical role at all: founder alleles are unique origin labels
        (paddedID_locus_copy), so this is IBD gene dropping, which tracks
        descent rather than allelic state. The transmission probability is
        Mendelian segregation, 0.5 by definition -- see
        MENDELIAN_TRANSMISSION_P.
    seed : int, optional
        Seed for the random number generator (default is 5048665).

    Returns:
    -------
    dict
        A dictionary of ancestral inbreeding coefficients keyed to animal IDs.
    """
    logger.info('Entered dropped_ancestral_inbreeding()')

    # NO DOMAIN PRECONDITION. A half-founder used to be refused here under
    # because Suwanlee et al. (2007) p.490 does not define the gene
    # source of a missing parental side. Baumung et al. (2015) p.102 does, and
    # the drop loop below implements it, so all three parentage cases -- zero,
    # one and two known parents -- are now in the supported domain.
    #
    # The mechanism repair does NOT rest on that removed refusal. The
    # sentinel is never used as an index: parents are resolved by ID, the
    # missing-parent path is an explicit construction rather than a fall-through,
    # and allele identity is int64 so no placeholder can compare equal to itself
    # as false autozygosity.

    # Validate parameters
    if rounds < 1:
        logger.error("Rounds must be greater than 0. Defaulting to 100.")
        rounds = 100
    if loci < 1:
        logger.error("Loci must be greater than 0. Defaulting to 100.")
        loci = 100
    # `frequency` is accepted for call compatibility and ignored. Warn only when
    # a caller actually passes one, so that ordinary use produces no noise while
    # the deprecation stays visible to the callers it concerns.
    if frequency is not None:
        warnings.warn(
            "dropped_ancestral_inbreeding(frequency=...) is deprecated and "
            "ignored. It was used as the allele-transmission probability, but "
            "this routine performs IBD gene dropping, in which transmission is "
            "Mendelian segregation (0.5) and an allele frequency has no role. "
            "Results no longer depend on this argument.",
            DeprecationWarning, stacklevel=2)

    try:
        seed = int(seed)
    except ValueError:
        logger.error("Seed must be an integer. Defaulting to 5048665.")
        seed = 5048665

    # Simulation-local RNG. Previously np.random.seed(seed) + np.random.rand(),
    # which reseeded numpy's process-global state for every other caller. The
    # public `seed` argument is unchanged, and no `rng` parameter is added; only
    # where the state lives has changed. This alters the draw sequence, which is
    # an intentional reproducibility change and not a scientific one -- see
    # a dedicated gene-drop verification.
    rng = np.random.default_rng(seed)

    missing = str(pedobj.kw['missing_parent'])
    id2aic = {p.animalID: 0.0 for p in pedobj.pedigree}

    for _replicate in range(rounds):
        # Simulation state, all local to this call. Keyed by
        # animalID and populated only as an animal is processed, so a parent
        # that has not been reached yet is ABSENT rather than holding a
        # placeholder -- which is what made the old `["", ""]` initialisation
        # able to read as autozygous.
        #
        #   labels[id] = (int64[loci], int64[loci])   founder-allele identity
        #   flags[id]  = (bool[loci],  bool[loci])    persistent IBD flags
        #
        # TWO flag arrays per animal, not one indicator per locus. D1 scores an
        # individual carrying exactly one flagged allele at a locus as one half
        # there, which a per-locus boolean cannot express -- collapsing them
        # would silently implement D2.
        labels = {}
        flags = {}
        next_label = 0

        for animal in pedobj.pedigree:
            animal_id = animal.animalID

            if animal.founder == 'y':
                # "Two unique alleles are assigned to each founder"
                # (Suwanlee et al. 2007, p.490). The label is one int64 per
                # (founder, allele copy), BROADCAST across the locus axis: the
                # same integer appears at every locus of one copy. That is
                # correct because identity is only ever compared elementwise,
                # within a locus -- the locus coordinate is part of the identity
                # domain by construction, and no cross-locus comparison is
                # performed anywhere. This is deliberately NOT a claim that
                # labels are unique across loci; they are not.
                labels[animal_id] = (np.full(loci, next_label, dtype=np.int64),
                                     np.full(loci, next_label + 1, dtype=np.int64))
                flags[animal_id] = (np.zeros(loci, dtype=bool),
                                    np.zeros(loci, dtype=bool))
                next_label += 2
                # Nothing can arrive flagged at a founder.
                continue

            # Parents are resolved BY ID, never by position: pedigree[id - 1] is
            # what turned the sentinel 0 into pedigree[-1], the last animal in
            # the pedigree, silently and with no IndexError.
            inherited = []
            for parent_id in (animal.sireID, animal.damID):
                if str(parent_id) == missing:
                    # THE HALF-FOUNDER RULE. Baumung et al. (2015)
                    # p.102, "Computational strategy", source-explicit:
                    #
                    #   "For each half-founder, that is an animal with just one
                    #   parent known, a dummy founder is created and the unknown
                    #   second parent is assigned an artificially created new
                    #   identification number and also provided with two unique
                    #   alleles."
                    #
                    # So the unknown parental side is an ORDINARY FOUNDER: two
                    # unique alleles, nothing flagged, one gamete transmitted by
                    # the same Mendelian draw as any other parent. Suwanlee et
                    # al. (2007) p.490 is silent on this case; the authority here
                    # is GRAIN.
                    #
                    # Minted here and deliberately NOT stored in labels/flags. A
                    # dummy serves exactly one offspring, so it is never looked
                    # up again -- and keeping it out of the dict means there is
                    # no synthetic key that could collide with a real animalID,
                    # and no way for two unknown parental slots to become one.
                    # Two half-founders that both record the sentinel are NOT
                    # the same individual: the pedigree carries no information
                    # that they are, and merging them would invent a
                    # relationship. This is a structural guarantee, not a
                    # convention the tests have to catch.
                    #
                    # next_label is the allocator the founder branch above uses,
                    # so a dummy is unrelated to every recorded founder and to
                    # every other dummy by construction.
                    parent_labels = (np.full(loci, next_label, dtype=np.int64),
                                     np.full(loci, next_label + 1, dtype=np.int64))
                    parent_flags = (np.zeros(loci, dtype=bool),
                                    np.zeros(loci, dtype=bool))
                    next_label += 2
                elif parent_id not in labels:
                    # A KNOWN parent that has not been reached yet. This is not a
                    # data problem -- the sentinel was handled above -- it is a
                    # broken invariant: the pedigree is not ordered
                    # parents-before-offspring, or this routine has a defect.
                    # Raised explicitly rather than surfacing a bare KeyError,
                    # which would leak an implementation detail as the contract.
                    raise PyPedalInternalError(
                        'dropped_ancestral_inbreeding: animal %s has parent %s, '
                        'which has not been dropped through yet. Suwanlee\'s '
                        'drop is a single pass and requires parents to precede '
                        'offspring in pedobj.pedigree; load the pedigree with '
                        'reorder enabled.' % (animal_id, parent_id))
                else:
                    parent_labels, parent_flags = labels[parent_id], flags[parent_id]
                # One Mendelian draw picks the allele AND its flag together, so
                # a flag can never separate from the allele it belongs to.
                take_first = rng.random(loci) < MENDELIAN_TRANSMISSION_P
                inherited.append(
                    (np.where(take_first, parent_labels[0], parent_labels[1]),
                     np.where(take_first, parent_flags[0], parent_flags[1])))

            (label_1, flag_1), (label_2, flag_2) = inherited

            # This animal's coefficient counts the alleles that arrived ALREADY
            # flagged, and is taken BEFORE any new flag is raised below.
            # Autozygosity arising here contributes to its OFFSPRING, not to
            # itself -- "alleles which are identical by descent for the first
            # time are flagged and contribute to ancestral inbreeding
            # coefficient of offspring" (p.490). Same convention as Ballou's
            # recursion, where f_a is built from the parents.
            #
            # Denominator: flagged alleles / (2 * loci),
            # taken from the published experiment rather than a source formula.
            id2aic[animal_id] += float(
                (flag_1.sum() + flag_2.sum()) / (2.0 * loci))

            # "Flagging alleles once they are in IBD state for the first time."
            # The flag is a persistent property of the ALLELE and travels with
            # every copy of it, so it is stored and then inherited by
            # descendants -- this is what makes inbreeding more than one
            # generation back visible.
            autozygous = label_1 == label_2
            labels[animal_id] = (label_1, label_2)
            flags[animal_id] = (flag_1 | autozygous, flag_2 | autozygous)

    # Mean over replicates. Each replicate contributes one already-per-locus
    # estimate, so with equal-size rounds this is the mean over all
    # rounds * loci independent loci (the missing
    # outer division).
    for animal_id in id2aic:
        id2aic[animal_id] /= float(rounds)

    logger.info('Exited dropped_ancestral_inbreeding()')

    # Postcondition. F_a is the probability of having
    # inherited an allele that was inbred at some point, so it lies in [0, 1].
    # O(n) over values already accumulated.
    #
    # This no longer fires. An earlier revision of this comment said it did,
    # because the routine returned rounds x the coefficient;
    # that was repaired by the division immediately above and the comment was
    # left behind describing the defect as live.
    pyp_validate.check_ancestral_inbreeding(
        pedobj, id2aic, 'dropped_ancestral_inbreeding')

    return id2aic


def ballou_ancestral_inbreeding(pedobj):
    """
    ballou_ancestral_inbreeding() calculates ancestral inbreeding,
    the probability of an individual inheriting an allele that has
    undergone inbreeding in the past at least once, using the method
    of Ballou (1997).

    Parameters:
    ----------
    pedobj : object
        A PyPedal pedigree object.

    Returns:
    -------
    dict
        A dictionary of ancestral inbreeding coefficients keyed to animal IDs.
    """
    logger.info("Entered ballou_ancestral_inbreeding()")

    # Initialize the dictionary mapping animal IDs to ancestral inbreeding coefficients
    id2aic = {p.animalID: 0.0 for p in pedobj.pedigree}

    # Calculate coefficients of inbreeding if they are not already in the pedigree
    if not pedobj.kw.get("f_computed", False):
        pyp_nrm.inbreeding(pedobj, output=False)

    # Calculate ancestral inbreeding
    for p in pedobj.pedigree:
        if str(p.sireID) == str(pedobj.kw["missing_parent"]):
            f_s = 0.0
            f_as = 0.0
        else:
            f_s = pedobj.pedigree[p.sireID - 1].fa
            f_as = id2aic[p.sireID]

        if str(p.damID) == str(pedobj.kw["missing_parent"]):
            f_d = 0.0
            f_ad = 0.0
        else:
            f_d = pedobj.pedigree[p.damID - 1].fa
            f_ad = id2aic[p.damID]

        id2aic[p.animalID] = (
            f_as + (1.0 - f_as) * f_s + f_ad + (1.0 - f_ad) * f_d
        ) / 2.0

    # Ballou's f_a is a probability, so the [0, 1] bound is legitimate here.
    # It is a postcondition on values already computed -- O(n), no recomputation.
    # This routine now agrees with the independent Ballou oracle on all eight corpus
    # pedigrees plus three purpose-built hand-derivable fixtures.
    pyp_validate.check_ancestral_inbreeding(
        pedobj, id2aic, "ballou_ancestral_inbreeding")

    logger.info("Exited ballou_ancestral_inbreeding()")
    return id2aic
