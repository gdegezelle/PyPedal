#!/usr/bin/env python3

"""
pyp_nrm.py - A module for computing numerator relationship matrices (NRM).

Version: see PyPedal.__version__
Author: John B. Cole, PhD (john.cole@ars.usda.gov)
License: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
"""

import sys
import os
import shutil
import heapq
import numpy as np
import logging
from scipy.sparse import lil_matrix, coo_matrix, csr_matrix
import math
import copy
import datetime
import subprocess
import time
from math import sqrt
from . import pyp_errors, pyp_network, pyp_utils
from . import pyp_validate
from typing import List


def _matrix_value(matrix, row, col):
    """Return A[row, col] as a float for dense or sparse relationship matrices."""
    try:
        value = matrix[row, col]
    except Exception:
        value = matrix[row][col]
    if hasattr(value, "toarray"):
        value = value.toarray().item()
    return float(value)


def _coi_from_matrix(matrix, row):
    """Return A[ii] - 1 for dense or sparse relationship matrices."""
    return _matrix_value(matrix, row, row) - 1.0


def _unpack_inbreeding(result, rels):
    """Normalize inbreeding helper return values to (fx, reldict-or-None)."""
    if isinstance(result, tuple):
        fx = result[0]
        reldict = result[1] if len(result) > 1 else None
    else:
        fx = result
        reldict = None
    if fx is False or fx is None:
        fx = {}
    if rels and reldict is None:
        reldict = {}
    return fx, reldict


def a_matrix(pedobj, save=False):
    """
    Form a numerator relationship matrix (NRM) from a pedigree. DEPRECATED.
    :param pedobj: A PyPedal pedigree object.
    :param save: Flag to indicate whether or not the relationship matrix is written to a file.
    :return: The NRM as a NumPy matrix.
    """
    logging.info('Entered a_matrix()')

    try:
        l = pedobj.metadata.num_records
        a = np.zeros((l, l), dtype=float)  # Initialize a matrix of zeros of appropriate dimension

        for row in range(l):
            for col in range(row, l):
                # Cast IDs to integers (in case they are read as strings)
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)

                if (
                    str(pedobj.pedigree[col].sireID) == str(pedobj.kw['missing_parent'])
                    and str(pedobj.pedigree[col].damID) == str(pedobj.kw['missing_parent'])
                ):
                    # Both parents unknown and assumed unrelated
                    a[row, col] = 1.0 if row == col else 0.0
                elif str(pedobj.pedigree[col].sireID) == str(pedobj.kw['missing_parent']):
                    # Sire unknown, dam known
                    a[row, col] = 1.0 if row == col else 0.5 * a[row, pedobj.pedigree[col].damID - 1]
                elif str(pedobj.pedigree[col].damID) == str(pedobj.kw['missing_parent']):
                    # Sire known, dam unknown
                    a[row, col] = 1.0 if row == col else 0.5 * a[row, pedobj.pedigree[col].sireID - 1]
                else:
                    # Both parents known
                    if row == col:
                        a[row, col] = 1.0 + 0.5 * a[
                            pedobj.pedigree[col].sireID - 1, pedobj.pedigree[col].damID - 1
                        ]
                    else:
                        intermediate = (
                            a[row, pedobj.pedigree[col].sireID - 1]
                            + a[row, pedobj.pedigree[col].damID - 1]
                        )
                        a[row, col] = 0.5 * intermediate

                # Symmetrize the matrix
                a[col, row] = a[row, col]

        if save:
            a_outputfile = f"{pedobj.kw['filetag']}_a_matrix_.dat"
            with open(a_outputfile, 'w') as aout:
                for row in range(l):
                    line = ",".join(f"{a[row, col]:.5f}" for col in range(l)) + "\n"
                    aout.write(line)

        logging.info('Exited a_matrix()')
        return a

    except Exception as e:
        logging.error(f"Error in a_matrix(): {e}")
        return np.zeros((1, 1), dtype=float)


def fast_a_matrix(pedigree, pedopts, save=False, method='sparse', debug=False, fill=1):
    """
    Form a numerator relationship matrix from a pedigree.

    Args:
        pedigree (PyPedal pedigree object)
        pedopts (Dict[str, Any]): PyPedal options.
        save (bool, optional): Flag to indicate whether or not the relationship matrix is written to a file.
        method (str, optional): Use dense or sparse matrix storage. Defaults to 'sparse' for better memory efficiency.
        debug (bool, optional): Flag turning debugging messages on (True) and off (False)
        fill (bool, optional): Fill both the upper and lower off-diagonals when True, only the upper otherwise.

    Returns:
        The NRM as a NumPy matrix or sparse matrix, depending on method.
    """
    logging.info('Entered fast_a_matrix()')

    foundercoi = int(pedopts.get("foundercoi", 0))
    if foundercoi not in [0, 1]:
        foundercoi = 0

    ped_length = len(pedigree)
    
    # Priority 3: Pre-extract ID arrays for better performance
    # Convert missing_parent to integer once
    missing_parent_int = int(pedopts["missing_parent"])
    
    # Extract all IDs as numpy arrays (faster access)
    animal_ids = np.array([int(a.animalID) for a in pedigree], dtype=np.int32)
    sire_ids = np.array([int(a.sireID) for a in pedigree], dtype=np.int32)
    dam_ids = np.array([int(a.damID) for a in pedigree], dtype=np.int32)
    
    # Auto-select method for large pedigrees if not specified
    if method not in ["dense", "sparse"]:
        method = "sparse" if ped_length > 10000 else "dense"

    # Priority 1: Use sparse matrices by default, with better format
    use_sparse = (method == "sparse")
    
    try:
        if use_sparse:
            # Use LIL for efficient incremental construction, then convert to CSR
            a = lil_matrix((ped_length, ped_length), dtype=np.float64)
        else:
            a = np.zeros([ped_length, ped_length], dtype="float64")
    except MemoryError:
        # Fallback to memmap if even sparse fails
        if use_sparse:
            logging.warning("Sparse matrix allocation failed, falling back to memmap")
            use_sparse = False
        a = np.memmap(
            "fast_a_matrix_mmap.bin", dtype="float32", mode="w+", shape=(ped_length, ped_length)
        )
    except Exception as e:
        logging.error(f"Unable to allocate a matrix of rank {ped_length} in fast_a_matrix(): {e}")
        return False

    # Priority 2: Eliminate string conversions - use integer comparisons directly
    # Initialize diagonal with 1.0
    for i in range(ped_length):
        a[i, i] = 1.0
        if foundercoi == 1 and sire_ids[i] == missing_parent_int and dam_ids[i] == missing_parent_int:
            a[i, i] = 1.0 + pedigree[i].fa

    # Build relationship matrix using integer comparisons (no string conversions)
    for row in range(ped_length):
        for col in range(row, ped_length):
            sire_col = sire_ids[col]
            dam_col = dam_ids[col]
            
            # Use integer comparisons instead of string conversions
            if sire_col != missing_parent_int and dam_col != missing_parent_int:
                if row == col:
                    a[row, col] += 0.5 * float(a[sire_col - 1, dam_col - 1])
                else:
                    a[row, col] = 0.5 * float(a[row, sire_col - 1] + a[row, dam_col - 1])
                    if fill:
                        a[col, row] = a[row, col]
            elif sire_col == missing_parent_int and dam_col != missing_parent_int:
                if row != col:
                    a[row, col] = 0.5 * float(a[row, dam_col - 1])
                    if fill:
                        a[col, row] = a[row, col]
            elif sire_col != missing_parent_int and dam_col == missing_parent_int:
                if row != col:
                    a[row, col] = 0.5 * float(a[row, sire_col - 1])
                    if fill:
                        a[col, row] = a[row, col]

    # Convert sparse matrix to CSR format for better performance
    if use_sparse:
        a = a.tocsr()

    if save:
        a_outputfile = f"{pedopts['filetag']}_new_a_matrix_.dat"
        with open(a_outputfile, "w") as aout:
            aout.write("Produced by pyp_nrm/fast_a_matrix()\n")
            # Handle both dense and sparse matrices
            if use_sparse:
                a_dense = a.toarray() if hasattr(a, 'toarray') else a
            else:
                a_dense = a
            for row in range(ped_length):
                line = ",".join(f"{float(a_dense[row, col]):.5f}" for col in range(ped_length)) + "\n"
                aout.write(line)

    logging.info('Exited fast_a_matrix()')
    
    # Return sparse matrix directly if sparse, or convert to numpy array if dense
    if use_sparse and not save:
        return a
    else:
        return np.array(a, dtype=float) if not use_sparse else np.array(a.toarray(), dtype=float)



def fast_a_matrix_r(pedigree, pedopts, save=False, method="sparse"):
    """
    Form a relationship matrix from a pedigree. fast_a_matrix_r() differs from
    fast_a_matrix() in that the coefficients of relationship are corrected for the
    inbreeding of the parents.

    :param pedigree: A PyPedal pedigree.
    :param pedopts: PyPedal options.
    :param save: Flag to indicate whether or not the relationship matrix is written to a file.
    :param method: Use dense or sparse matrix storage. Defaults to 'sparse' for better memory efficiency.
    """
    logging.info("Entered fast_a_matrix_r()")

    ped_length = len(pedigree)
    
    # Priority 3: Pre-extract ID arrays for better performance
    # Convert missing_parent to integer once
    missing_parent_int = int(pedopts["missing_parent"])
    
    # Extract all IDs as numpy arrays (faster access)
    animal_ids = np.array([int(a.animalID) for a in pedigree], dtype=np.int32)
    sire_ids = np.array([int(a.sireID) for a in pedigree], dtype=np.int32)
    dam_ids = np.array([int(a.damID) for a in pedigree], dtype=np.int32)
    
    # Auto-select method for large pedigrees if not specified
    if method not in ["dense", "sparse"]:
        method = "sparse" if ped_length > 10000 else "dense"

    # Priority 1: Use sparse matrices by default
    use_sparse = (method == "sparse")
    
    try:
        if use_sparse:
            a = lil_matrix((ped_length, ped_length), dtype=np.float64)
        else:
            a = np.zeros((ped_length, ped_length), dtype="float64")
    except MemoryError:
        if use_sparse:
            logging.warning("Sparse matrix allocation failed, falling back to memmap")
            use_sparse = False
        a = np.memmap(
            "fast_a_matrix_r_mmap.bin", dtype="float32", mode="w+", shape=(ped_length, ped_length)
        )
    except Exception as e:
        logging.error(f"Unable to allocate a matrix of rank {ped_length} in fast_a_matrix_r(): {e}")
        return False

    # Priority 2: Use integer comparisons directly (no string conversions)
    for row in range(ped_length):
        for col in range(row, ped_length):
            sire_col = sire_ids[col]
            dam_col = dam_ids[col]
            
            if sire_col == missing_parent_int and dam_col == missing_parent_int:
                if row == col:
                    a[row, col] = 1.0
            elif sire_col == missing_parent_int:
                if row == col:
                    a[row, col] = 1.0
                else:
                    a[row, col] = 0.5 * float(a[row, dam_col - 1])
                    a[col, row] = a[row, col]
            elif dam_col == missing_parent_int:
                if row == col:
                    a[row, col] = 1.0
                else:
                    a[row, col] = 0.5 * float(a[row, sire_col - 1])
                    a[col, row] = a[row, col]
            elif sire_col != missing_parent_int and dam_col != missing_parent_int:
                if row == col:
                    a[row, col] = 1.0 + 0.5 * float(a[sire_col - 1, dam_col - 1])
                else:
                    intermediate = a[row, sire_col - 1] + a[row, dam_col - 1]
                    a[row, col] = 0.5 * float(intermediate)
                    a[col, row] = a[row, col]

    # Convert sparse matrix to CSR format for better performance
    if use_sparse:
        a = a.tocsr()

    if save:
        a_outputfile = f"{pedopts['filetag']}_a_matrix_r_.dat"
        with open(a_outputfile, "w") as aout:
            aout.write("Produced by pyp_nrm/fast_a_matrix_r()\n")
            # Handle both dense and sparse matrices
            if use_sparse:
                a_dense = a.toarray() if hasattr(a, 'toarray') else a
            else:
                a_dense = a
            for row in range(ped_length):
                line = ",".join(f"{float(a_dense[row, col]):.5f}" for col in range(ped_length)) + "\n"
                aout.write(line)

    logging.info("Exited fast_a_matrix_r()")
    
    # Return sparse matrix directly if sparse, or convert to numpy array if dense
    if use_sparse and not save:
        return a
    else:
        return np.array(a, dtype=float) if not use_sparse else np.array(a.toarray(), dtype=float)


def inbreeding(pedobj, method='tabular', gens=0, rels=0, output=True, force=False, amethod=3):
    """
    Dispatch pedigrees to the appropriate function for computing CoI.

    The default ``method='tabular'`` uses the tabular/full-NRM-style path and
    is intended for small pedigrees. There is no automatic switch at 10,000
    animals or at any other size. For large pedigrees, including tens of
    thousands of animals, use ``method='meu_luo'``. Forming a full NRM or
    calling other dense relationship-matrix routines is not the intended
    large-pedigree path.

    :param pedobj: A PyPedal pedigree object.
    :param method: Method to compute CoI ('tabular'|'vanraden'|'meu_luo'|'mod_meu_luo'|'aguilar').
    :param gens: Number of generations from the pedigree to use for CoI. Default: 0 (complete pedigree).
    :param rels: Flag to compute summary statistics for relationships (0: no, 1: yes).
    :param output: Flag to write output files (True: yes, False: no).
    :param force: Override use of NRM for CoI (0: use NRM, 1: ignore NRM).
    :param amethod: Method parameter for Aguilar's INBUPGF90 program.
    :return: Dictionary of CoI keyed to renumbered animal IDs.
    """
    logging.info('Entered inbreeding()')
    
    fx = {}
    reldict = None
    
    valid_methods = ['vanraden', 'tabular', 'meu_luo', 'mod_meu_luo', 'aguilar']
    if method not in valid_methods:
        logging.warning(f"Unrecognized method '{method}' provided to inbreeding(); defaulting to 'tabular'.")
        method = 'tabular'

    if gens < 0:
        logging.warning(f"Invalid gens value '{gens}' provided to inbreeding(); defaulting to 0.")
        gens = 0

    if rels:
        rel_dict = {
            'r_count': (pedobj.metadata.num_records * (pedobj.metadata.num_records + 1)) // 2,
            'r_nonzero_count': 0,
            'r_min': 0.0,
            'r_max': 0.0,
            'r_rng': 0.0,
            'r_avg': 0.0,
            'r_nonzero_avg': 0.0,
            'r_sum': 0.0,
            'r_nonzero_sum': 0.0,
        }

    # Use precomputed NRM if available
    if pedobj.kw['form_nrm'] and pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records and not force:
        logging.info("Using precomputed NRM for CoI.")
        for _i in range(pedobj.metadata.num_records):
            fx[pedobj.pedigree[_i].animalID] = _coi_from_matrix(pedobj.nrm.nrm, _i)

        if rels:
            reldict = compute_relationship_stats(pedobj.nrm.nrm)
    else:
        if method == 'vanraden':
            fx, reldict = _unpack_inbreeding(
                inbreeding_vanraden(pedobj, gens=gens, rels=rels), rels
            )
        elif method == 'meu_luo':
            warn_no_relationships(method, rels)
            fx, reldict = _unpack_inbreeding(inbreeding_meuwissen_luo(pedobj, gens=gens), rels)
        elif method == 'mod_meu_luo':
            warn_no_relationships(method, rels)
            fx, reldict = _unpack_inbreeding(
                inbreeding_modified_meuwissen_luo(pedobj, gens=gens), rels
            )
        elif method == 'aguilar':
            logging.info(f"Using INBUPGF90 program with method {amethod}.")
            fx, reldict = _unpack_inbreeding(inbreeding_aguilar(pedobj, amethod), rels)
        else:
            fx, reldict = _unpack_inbreeding(
                inbreeding_tabular(pedobj, gens=gens, rels=rels), rels
            )

    # Postcondition. F is the probability that an animal's
    # two alleles at a locus are identical by descent, so it must lie in [0, 1].
    # O(n) over coefficients already computed -- nothing is recomputed here and
    # the NRM is not touched, so this does not change any complexity class.
    pyp_validate.check_inbreeding_coefficients(pedobj, fx, "inbreeding(method=%s)" % method)

    if output:
        write_inbreeding_output(pedobj, fx)

    update_pedigree_metadata(fx, pedobj)

    out_dict = {'metadata': compute_inbreeding_stats(fx), 'fx': fx}
    if rels:
        out_dict['rel_dict'] = reldict if reldict is not None else rel_dict
    return out_dict


def warn_no_relationships(method, rels):
    """Log a warning if relationships were requested for a method that does not provide them."""
    if rels:
        logging.warning(
            f"Method '{method}' does not compute relationships. Only coefficients of inbreeding will be returned."
        )


def compute_relationship_stats(nrm_matrix):
    """Compute relationship statistics from the NRM matrix."""
    n = nrm_matrix.shape[0]
    reldict = {
        'r_count': (n * (n + 1)) // 2,
        'r_nonzero_count': 0,
        'r_min': 1.0,
        'r_max': 0.0,
        'r_sum': 0.0,
        'r_nonzero_sum': 0.0,
    }

    for i in range(n):
        for j in range(i, n):
            val = _matrix_value(nrm_matrix, i, j)
            reldict['r_sum'] += val
            if val > 0.0:
                reldict['r_nonzero_count'] += 1
                reldict['r_nonzero_sum'] += val
                reldict['r_max'] = max(reldict['r_max'], val)
                reldict['r_min'] = min(reldict['r_min'], val)

    reldict['r_avg'] = reldict['r_sum'] / reldict['r_count'] if reldict['r_count'] else 0.0
    reldict['r_nonzero_avg'] = (
        reldict['r_nonzero_sum'] / reldict['r_nonzero_count'] if reldict['r_nonzero_count'] else 0.0
    )
    return reldict


def compute_inbreeding_stats(fx):
    """Compute inbreeding statistics for all animals and for those with CoI > 0."""

    def _stats(values):
        if not values:
            return {
                "f_count": 0,
                "f_sum": 0.0,
                "f_min": 0.0,
                "f_max": 0.0,
                "f_avg": 0.0,
                "f_rng": 0.0,
            }
        f_count = len(values)
        f_sum = sum(values)
        f_min = min(values)
        f_max = max(values)
        return {
            "f_count": f_count,
            "f_sum": f_sum,
            "f_min": f_min,
            "f_max": f_max,
            "f_avg": f_sum / f_count,
            "f_rng": f_max - f_min,
        }

    values = list(fx.values()) if fx else []
    nonzero = [value for value in values if value > 0.0]
    return {"all": _stats(values), "nonzero": _stats(nonzero)}


def update_pedigree_metadata(fx, pedobj):
    """Update the metadata of the pedigree object with inbreeding coefficients."""
    if not isinstance(fx, dict) or not fx:
        return
    for k, v in fx.items():
        try:
            pedobj.pedigree[int(k) - 1].fa = v
        except (ValueError, TypeError, IndexError, AttributeError):
            continue
    pedobj.kw['f_computed'] = True


def write_inbreeding_output(pedobj, fx):
    """Write inbreeding coefficients to an output file."""
    if not isinstance(fx, dict) or not fx:
        return
    output_file = f"{pedobj.kw['filetag']}_inbreeding.dat"
    try:
        with open(output_file, 'w') as f:
            f.write("# Inbreeding coefficients\n")
            if 'ASD' in pedobj.kw['pedformat']:
                f.write("# Name\tRenum ID\tf_x\n")
            else:
                f.write("# Orig ID\tRenum ID\tf_x\n")

            for k, v in fx.items():
                try:
                    animal = pedobj.pedigree[int(k) - 1]
                except (ValueError, TypeError, IndexError):
                    continue
                if 'ASD' in pedobj.kw['pedformat']:
                    f.write(f"{animal.name}\t{k}\t{v}\n")
                else:
                    f.write(f"{animal.originalID}\t{k}\t{v}\n")
    except OSError as exc:
        logging.error("Unable to write inbreeding output %s: %s", output_file, exc)


def inbreeding_vanraden(pedobj, cleanmaps=True, gens=0, rels=False):
    """
    inbreeding_vanraden() uses VanRaden's (1992) method for computing coefficients of
    inbreeding in a large pedigree.
    
    :param pedobj: A PyPedal pedigree object.
    :param cleanmaps: Whether to delete subpedigree ID maps after use.
    :param gens: Number of generations from the pedigree to use for calculating CoI (default is 0, meaning the full pedigree).
    :param rels: Whether to compute summary statistics for coefficients of relationship.
    :return: A dictionary of CoI keyed to renumbered animal IDs (and optionally the relationship dictionary).
    """
    logging.info("Entered inbreeding_vanraden()")
    from . import pyp_network

    ng = pyp_network.ped_to_graph(pedobj)

    _ped = []
    top_ped = []

    if gens > 0:
        top_peddict = pyp_network.find_ancestors_g(ng, len(pedobj.idmap), {}, gens)
        top_peddict[len(pedobj.idmap)] = 1
        top_ped = list(top_peddict.keys())
        top_r = []
        _anids = []
        for _j in top_ped:
            if top_peddict[_j] <= gens:
                top_r.append(copy.deepcopy(pedobj.pedigree[int(_j) - 1]))
                if top_peddict[_j] == 1:
                    top_r[-1].sireID = 0
                    top_r[-1].damID = 0
                _anids.append(top_r[-1].animalID)
    else:
        _anids = list(pedobj.backmap.keys())

    fx = {}
    _parents = {}
    _anids.sort(reverse=True)
    _cum_pct_proc = 0.0
    _related = {}

    if rels:
        reldict = {
            'r_count': 0,
            'r_nonzero_count': 0,
            'r_nonzero_sum': 0.0,
            'r_max': 0.0,
            'r_min': 1.0,
            'r_sum': 0.0,
        }

    for i in _anids:
        if gens == 0:
            _parent_key = f"{pedobj.pedigree[int(i) - 1].sireID}_{pedobj.pedigree[int(i) - 1].damID}"
        else:
            _animal = next((animal for animal in top_r if int(animal.animalID) == int(i)), None)
            if _animal is None:
                _animal = pedobj.pedigree[int(i) - 1]
            _parent_key = f"{_animal.sireID}_{_animal.damID}"
        
        if i not in fx:
            try:
                fx[i] = fx[_parents[_parent_key]]
            except KeyError:
                _ped = top_peddict if gens > 0 else pyp_network.find_ancestors(ng, i, [])
                if gens == 0:
                    _ped.append(i)
                _r = [copy.deepcopy(pedobj.pedigree[int(j) - 1]) for j in _ped]
                
                _missing = pedobj.kw.get('missing_parent', 0)
                if pedobj.kw.get('slow_reorder', False):
                    _r = pyp_utils.reorder(_r, f"{pedobj.kw['filetag']}_{i}",
                                           missingparent=_missing)
                else:
                    _r = pyp_utils.fast_reorder(_r, f"{pedobj.kw['filetag']}_{i}",
                                                missingparent=_missing)

                _s, _map = pyp_utils.renumber(_r, f"{pedobj.kw['filetag']}_{i}", returnmap=True, 
                                              debug=pedobj.kw.get('debug_messages', False), 
                                              missingparent=_missing,
                                              animaltype=pedobj.kw.get('animal_type', 'new'))
                
                _backmap = {v: k for k, v in _map.items()}
                _opts = copy.deepcopy(pedobj.kw)
                _opts['filetag'] = f"{pedobj.kw['filetag']}_{i}"

                if pedobj.kw.get('nrm_method') == 'nrm':
                    _a = fast_a_matrix(_s, _opts, method=pedobj.kw.get('matrix_type', 'sparse'))
                else:
                    _a = fast_a_matrix_r(_s, _opts, method=pedobj.kw.get('matrix_type', 'sparse'))

                if _a is False:
                    continue

                for j in range(len(_s)):
                    _orig_id = _backmap[_s[j].animalID]
                    fx.setdefault(_orig_id, _coi_from_matrix(_a, j))
                    
                    if rels:
                        for k in range(j, len(_s)):
                            if j != k:
                                rel_value = _matrix_value(_a, j, k)
                                _rxykey = f"{_backmap[_s[j].animalID]}_{_backmap[_s[k].animalID]}"
                                if rel_value > 0.0:
                                    _related.setdefault(_rxykey, rel_value)
                                    reldict['r_nonzero_count'] += 1
                                    reldict['r_nonzero_sum'] += rel_value
                                    reldict['r_max'] = max(reldict['r_max'], rel_value)
                                    reldict['r_min'] = min(reldict['r_min'], rel_value)
                                reldict['r_count'] += 1
                                reldict['r_sum'] += rel_value

                _parents[_parent_key] = i

                if cleanmaps:
                    pyp_utils.delete_id_map(f"{pedobj.kw['filetag']}_{i}")

    logging.info("Exited inbreeding_vanraden()")
    return (fx, reldict) if rels else fx


def inbreeding_aguilar(pedobj, amethod=3):
    """
    inbreeding_aguilar() uses Ignacio Aguilar's INBUPGF90 program to compute coefficients of
    inbreeding in large pedigrees.
    :param pedobj: A PyPedal pedigree object.
    :param amethod: The method for computing inbreeding (1, 2, or 3).
    :return: A dictionary of CoI keyed to renumbered animal IDs.
    """
    # Define log and file paths
    logfile = f"{pedobj.kw['pedname']}_aguilar.log"
    pedfile = f"aguilar_pedigree_{pedobj.kw['pedname']}.txt"
    coifile = f"{pedfile}.solinb"

    # Validate `amethod`
    if amethod not in [1, 2, 3]:
        amethod = 3

    # Verbose logging
    if pedobj.kw.get('messages') == 'verbose':
        print(f"[inbreeding_aguilar]: Started INBUPGF90 to calculate COI at {datetime.datetime.now():%Y-%m-%d %H:%M}")
        logging.info("[inbreeding_aguilar]: Started INBUPGF90 to calculate COI.")

    # Check for INBUPGF90 binary in the system's PATH
    if not shutil.which("inbupgf90"):
        raise FileNotFoundError("The INBUPGF90 binary is not found in your system's PATH.")

    with open(pedfile, "w", encoding="utf-8") as handle:
        for animal in pedobj.pedigree:
            handle.write(f"{animal.animalID} {animal.sireID} {animal.damID}\n")

    # Command to execute INBUPGF90
    callinbupgf90 = [
        "inbupgf90", 
        "--pedfile", pedfile, 
        "--method", str(amethod)
    ]

    # Execute INBUPGF90
    process = subprocess.Popen(
        callinbupgf90, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    # Poll process and wait for completion
    time_waited = 0
    while process.poll() is None:
        time.sleep(10)
        time_waited += 10
        if time_waited % 60 == 0 and pedobj.kw.get('messages') == 'verbose':
            print(f"\t[inbreeding_aguilar]: Waiting for INBUPGF90 to finish -- {time_waited // 60} minutes so far...")
            logging.info("[inbreeding_aguilar]: Waiting for INBUPGF90 to finish -- %d minutes so far...", time_waited // 60)

    # Capture output and errors
    results, errors = process.communicate()

    if errors:
        logging.error("[inbreeding_aguilar]: INBUPGF90 finished with errors: %s", errors.decode("utf-8"))
        if pedobj.kw.get('messages') == 'verbose':
            print(f"\t[inbreeding_aguilar]: INBUPGF90 finished with errors: {errors.decode('utf-8')}")
        raise RuntimeError(f"INBUPGF90 finished with errors: {errors.decode('utf-8')}")
    else:
        logging.info("[inbreeding_aguilar]: INBUPGF90 finished successfully.")
        if pedobj.kw.get('messages') == 'verbose':
            print(f"\t[inbreeding_aguilar]: INBUPGF90 finished successfully at {datetime.datetime.now():%Y-%m-%d %H:%M}")

    # Load coefficients of inbreeding into a dictionary
    if pedobj.kw.get('messages') == 'verbose':
        print(f"[inbreeding_aguilar]: Loading coefficients of inbreeding from {coifile} at {datetime.datetime.now():%Y-%m-%d %H:%M}")
        logging.info("[inbreeding_aguilar]: Loading coefficients of inbreeding from %s.", coifile)

    inbr = {}
    if not os.path.exists(coifile):
        raise FileNotFoundError(f"Expected output file {coifile} not found.")

    with open(coifile, "r") as ifh:
        for line in ifh:
            pieces = line.split()
            inbr[pieces[0]] = float(pieces[1])

    # Cleanup temporary files if necessary
    os.remove(pedfile)
    os.remove(coifile)
    os.remove(logfile)

    # Return the coefficients of inbreeding
    return inbr


def recurse_pedigree(pedobj, anid: int, _ped: List) -> List:
    """
    Recursively builds the subpedigrees used by inbreeding calculations.

    For the animal with animalID `anid`, recursively traverses the pedigree 
    and adds references to relatives to the temporary pedigree `_ped`.

    Parameters
    ----------
    pedobj : object
        A PyPedal pedigree object.
    anid : int
        The ID of the animal whose relatives are being located.
    _ped : List
        A temporary list storing references to relatives of `anid`.

    Returns
    -------
    List
        A list of references to the relatives of `anid` contained in the pedigree.
    """
    try:
        anid = int(anid)
        if anid != 0:
            if pedobj.pedigree[anid - 1] not in _ped:
                _ped.append(pedobj.pedigree[anid - 1])

        sire_id = pedobj.pedigree[anid - 1].sireID
        dam_id = pedobj.pedigree[anid - 1].damID

        if sire_id != pedobj.kw['missing_parent']:
            recurse_pedigree(pedobj, sire_id, _ped)

        if dam_id != pedobj.kw['missing_parent']:
            recurse_pedigree(pedobj, dam_id, _ped)

    except RecursionError:
        # Do NOT swallow this. A blanket `except Exception`
        # caught RecursionError and returned the ancestors collected so far, so
        # a pedigree deeper than the interpreter's recursion limit produced a
        # SILENTLY TRUNCATED ancestor list and every statistic derived from it
        # was computed on part of the pedigree. Failing loudly is the only
        # honest outcome: the caller must raise the recursion limit or use an
        # iterative path.
        logging.error(
            "recurse_pedigree() exceeded the recursion limit (%d) while walking "
            "the ancestors of animal %s. The ancestor list would have been "
            "silently truncated, so the traversal is being failed instead. "
            "Raise sys.setrecursionlimit() or use a shallower pedigree.",
            sys.getrecursionlimit(), anid)
        raise
    except (IndexError, KeyError, ValueError, TypeError, AttributeError) as e:
        # Genuine data problems: an ID with no record, or an unparseable one.
        logging.error(f"Error in recurse_pedigree for animal {anid}: {e}")

    return _ped


def recurse_pedigree_n(pedobj, anid, _ped, depth=3):
    """
    recurse_pedigree_n() recurses to build a pedigree of depth `n`. A depth
    less than 1 returns the animal whose relatives were to be identified.

    :param pedobj: A PyPedal pedigree object.
    :param anid: The ID of the animal whose relatives are being located.
    :param _ped: A temporary PyPedal pedigree that stores references to relatives of `anid`.
    :param depth: The depth of the pedigree to return.
    :return: A list of references to the relatives of `anid` contained in `pedobj`.
    """
    try:
        anid = int(anid)
        if anid != pedobj.kw['missing_parent']:
            if pedobj.pedigree[anid - 1] not in _ped:
                _ped.append(pedobj.pedigree[anid - 1])

        if depth > 0:
            sire_id = pedobj.pedigree[anid - 1].sireID
            dam_id = pedobj.pedigree[anid - 1].damID

            if sire_id != pedobj.kw['missing_parent']:
                recurse_pedigree_n(pedobj, sire_id, _ped, depth - 1)
            if dam_id != pedobj.kw['missing_parent']:
                recurse_pedigree_n(pedobj, dam_id, _ped, depth - 1)
    except Exception as e:
        # Log the exception or handle it if necessary
        logging.error(f"An error occurred in recurse_pedigree_n: {e}")
        pass

    return _ped


def recurse_pedigree_onesided(pedobj, anid, _ped, side):
    """
    recurse_pedigree_onesided() recurses to build a subpedigree from either the sire
    or dam side of a pedigree.

    :param pedobj: A PyPedal pedigree object.
    :param anid: The ID of the animal whose relatives are being located.
    :param _ped: A temporary PyPedal pedigree that stores references to relatives of `anid`.
    :param side: The side to build: 's' for sire and 'd' for dam.
    :return: A list of references to the relatives of `anid` contained in `pedobj`.
    """
    try:
        anid = int(anid)
        if anid != 0:
            if pedobj.pedigree[anid - 1] not in _ped:
                _ped.append(pedobj.pedigree[anid - 1])

        if side == 's':
            sire_id = pedobj.pedigree[anid - 1].sireID
            if sire_id != pedobj.kw['missing_parent']:
                recurse_pedigree(pedobj, sire_id, _ped)
        elif side == 'd':
            dam_id = pedobj.pedigree[anid - 1].damID
            if dam_id != pedobj.kw['missing_parent']:
                recurse_pedigree(pedobj, dam_id, _ped)
    except Exception as e:
        # Log the exception or handle it if necessary
        logging.error(f"An error occurred in recurse_pedigree_onesided: {e}")
        pass

    return _ped


def recurse_pedigree_idonly(pedobj, anid, _ped):
    """
    recurse_pedigree_idonly() performs the recursion needed to build subpedigrees.

    :param pedobj: A PyPedal pedigree object.
    :param anid: The ID of the animal whose relatives are being located.
    :param _ped: A list that stores the animalIDs of relatives of `anid`.
    :return: A list of animalIDs of the relatives of `anid` contained in `pedobj`.
    """
    try:
        anid = int(anid)
        if anid != 0:
            # Check if the animal ID is not already in the list
            if pedobj.pedigree[anid - 1].animalID not in _ped:
                _ped.append(pedobj.pedigree[anid - 1].animalID)

        # Get the sire and dam IDs
        sire_id = pedobj.pedigree[anid - 1].sireID
        dam_id = pedobj.pedigree[anid - 1].damID

        # Recurse through the sire's side
        if sire_id != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly(pedobj, sire_id, _ped)

        # Recurse through the dam's side
        if dam_id != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly(pedobj, dam_id, _ped)

    except Exception as e:
        # Log the exception or handle it if necessary
        logging.error(f"An error occurred in recurse_pedigree_idonly: {e}")
        pass

    return _ped


def recurse_pedigree_idonly_side(pedobj, anid, _ped, side='s'):
    """
    recurse_pedigree_idonly_side() performs the recursion needed to build
    a subpedigree containing only animal IDs for either all sires or all
    dams. That is, a pedigree would go sire-paternal grandsire-paternal
    great-grandsire, etc.

    :param pedobj: A PyPedal pedigree object.
    :param anid: The ID of the animal whose relatives are being located.
    :param _ped: A list that stores the animalIDs of relatives of `anid`.
    :param side: The side of the pedigree to follow ('s' for sires, 'd' for dams).
    :return: A list of animalIDs of the relatives of `anid` contained in `pedobj`.
    """
    if side not in ['s', 'd']:
        side = 's'  # Default to 's' if invalid input

    try:
        anid = int(anid)
        if anid != 0:
            # Add the current animal ID to the pedigree if it's not already included
            if pedobj.pedigree[anid - 1].animalID not in _ped:
                _ped.append(pedobj.pedigree[anid - 1].animalID)

        # Determine the sire or dam based on the selected side
        sire_id = pedobj.pedigree[anid - 1].sireID
        dam_id = pedobj.pedigree[anid - 1].damID

        # Recurse through the chosen side
        if side == 's' and sire_id != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly_side(pedobj, sire_id, _ped, side='s')
        elif side == 'd' and dam_id != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly_side(pedobj, dam_id, _ped, side='d')

    except Exception as e:
        # Log the exception or handle it if necessary
        logging.error(f"Error in recurse_pedigree_idonly_side: {e}")
        pass

    return _ped


def inbreeding_tabular(pedobj, gens=0, rels=0):
    """
    inbreeding_tabular() computes CoI using the tabular method by calling
    fast_a_matrix() to form the NRM directly. In order for this routine
    to return successfully, you must be able to allocate a matrix of floats
    of dimension len(myped)**2.

    :param pedobj: A PyPedal pedigree object.
    :param gens: The number of generations from the pedigree to be used for calculating CoI. Default is 0, which uses the complete pedigree.
    :param rels: Flag indicating whether to compute summary statistics for coefficients of relationship.
    :return: A dictionary of CoI keyed to renumbered animal IDs, and optionally a dictionary of relationship statistics if `rels` is True.
    """
    logging.info('Entered inbreeding_tabular()')

    if rels:
        reldict = {
            'r_count': 0,
            'r_nonzero_count': 0,
            'r_nonzero_sum': 0.0,
            'r_max': 0.0,
            'r_min': 1.0,
            'r_sum': 0.0,
        }

    try:
        # If generations are specified, handle subpedigree creation
        if int(gens) > 0:
            ng = pyp_network.ped_to_graph(pedobj)
            _ped = pyp_network.find_ancestors_g(ng, pedobj.metadata.num_records, [], gens)
            _a, _s, _r = [], [], []
            _map = {}

            for j in _ped:
                _r.append(copy.deepcopy(pedobj.pedigree[int(j) - 1]))

            _missing = pedobj.kw.get('missing_parent', 0)
            if pedobj.kw['slow_reorder']:
                _r = pyp_utils.reorder(_r, pedobj.kw['filetag'],
                                       missingparent=_missing)
            else:
                _r = pyp_utils.fast_reorder(_r, pedobj.kw['filetag'],
                                            missingparent=_missing)

            _s, _map = pyp_utils.renumber(
                _r,
                pedobj.kw['filetag'],
                returnmap=True,
                debug=pedobj.kw.get('debug_messages', False),
                missingparent=_missing,
                animaltype=pedobj.kw.get('animal_type', 'default')
            )

            _backmap = {v: k for k, v in _map.items()}
            _opts = copy.deepcopy(pedobj.kw)
            _opts['filetag'] = pedobj.kw['filetag']

            if pedobj.kw['nrm_method'] == 'nrm':
                _a = fast_a_matrix(_s, _opts, method=pedobj.kw.get('matrix_type', 'dense'))
            else:
                _a = fast_a_matrix_r(_s, _opts, method=pedobj.kw.get('matrix_type', 'dense'))

            if _a is False:
                raise RuntimeError("Unable to allocate the relationship matrix")

            fx = {}
            for j in range(len(_s)):
                fx[_backmap[_s[j].animalID]] = _coi_from_matrix(_a, j)

        else:
            if pedobj.kw['nrm_method'] == 'nrm':
                _a = fast_a_matrix(pedobj.pedigree, pedobj.kw, method=pedobj.kw.get('matrix_type', 'sparse'))
            else:
                _a = fast_a_matrix_r(pedobj.pedigree, pedobj.kw, method=pedobj.kw.get('matrix_type', 'sparse'))

            if _a is False:
                raise RuntimeError("Unable to allocate the relationship matrix")

            fx = {}
            for i in range(pedobj.metadata.num_records):
                fx[pedobj.pedigree[i].animalID] = _coi_from_matrix(_a, i)

                if rels:
                    for j in range(i, pedobj.metadata.num_records):
                        if i != j:
                            rel_value = _matrix_value(_a, i, j)
                            if rel_value > 0.0:
                                reldict['r_nonzero_count'] += 1
                                reldict['r_nonzero_sum'] += rel_value
                                reldict['r_max'] = max(reldict['r_max'], rel_value)
                                reldict['r_min'] = min(reldict['r_min'], rel_value)

                            reldict['r_count'] += 1
                            reldict['r_sum'] += rel_value

        logging.info('Exited inbreeding_tabular()')

        if rels:
            return fx, reldict
        else:
            return fx

    except pyp_errors.PyPedalError:
        # A structural refusal must reach the caller. Returning {} here is what
        # made `inbreeding(gens=N)` look like it had simply found nothing to
        # report: find_ancestors_g truncates at a generation boundary, the
        # animals on that boundary keep parent references to animals that were
        # cut, reorder could not resolve them, and this handler turned the
        # failure into an empty result. Every generation limit that produced a
        # dangling parent reference would otherwise look like "no inbreeding".
        logging.error("inbreeding_tabular() refused", exc_info=True)
        raise
    except Exception as e:
        logging.error(f"Error in inbreeding_tabular(): {e}")
        if rels:
            return {}, {}
        else:
            return {}


def _require_meuwissen_luo_numbering(pedigree, missing_parent):
    """
    Meuwissen and Luo (1992) index parents as ``animalID - 1`` and process
    animals in that order so a parent's *F* is already final. This call does
    not renumber. An unnumbered or parent-after-child list would produce
    silently wrong coefficients, so it is refused.
    """
    n = len(pedigree)
    missing = str(missing_parent)
    for i, animal in enumerate(pedigree):
        try:
            animal_id = int(animal.animalID)
        except (TypeError, ValueError) as exc:
            raise pyp_errors.PyPedalUsageError(
                "inbreeding_meuwissen_luo requires a supported renumbered "
                "pedigree (animalID must be integers 1..n with parents before "
                "offspring). Renumber first; this call did not compute "
                "coefficients."
            ) from exc
        if animal_id != i + 1:
            raise pyp_errors.PyPedalUsageError(
                "inbreeding_meuwissen_luo requires a supported renumbered "
                "pedigree (animalID must be integers 1..n with parents before "
                "offspring). Renumber first; this call did not compute "
                "coefficients."
            )
        for parent_id in (animal.sireID, animal.damID):
            if parent_id == missing_parent or str(parent_id) == missing:
                continue
            try:
                parent = int(parent_id)
            except (TypeError, ValueError) as exc:
                raise pyp_errors.PyPedalUsageError(
                    "inbreeding_meuwissen_luo requires a supported renumbered "
                    "pedigree (animalID must be integers 1..n with parents "
                    "before offspring). Renumber first; this call did not "
                    "compute coefficients."
                ) from exc
            if parent < 1 or parent > n or parent >= animal_id:
                raise pyp_errors.PyPedalUsageError(
                    "inbreeding_meuwissen_luo requires a supported renumbered "
                    "pedigree (animalID must be integers 1..n with parents "
                    "before offspring). Renumber first; this call did not "
                    "compute coefficients."
                )


def inbreeding_meuwissen_luo(pedobj, gens=0, **kw):
    """
    inbreeding_meuwissen_luo() computes CoI using the method of Meuwissen and
    Luo (1992). It calculates only inbreeding coefficients, not relationships.
    This code is a direct implementation of the algorithm presented on
    pp. 311-312 of Meuwissen, T.H.E., and Z. Luo. 1992. Computing inbreeding
    coefficients in large populations. Genet. Sel. Evol. 24:305-313.

    The ancestor set is a max-heap plus a membership set so each visit is
    logarithmic in the frontier size, not linear. ``L`` is reset only on
    indices touched for the current animal.

    The pedigree must already be numbered 1..n with parents before offspring.
    This routine does not renumber. An unnumbered list raises
    :class:`~PyPedal.pyp_errors.PyPedalUsageError`.
    """
    try:
        logging.info("Entered inbreeding_meuwissen_luo()")
    except Exception:
        pass

    n = len(pedobj.pedigree)
    missing_parent = pedobj.kw["missing_parent"]
    pedigree = pedobj.pedigree
    _require_meuwissen_luo_numbering(pedigree, missing_parent)

    sire_idx = [-1] * n
    dam_idx = [-1] * n
    animal_id = [0] * n
    for i, animal in enumerate(pedigree):
        animal_id[i] = animal.animalID
        sire_id = animal.sireID
        dam_id = animal.damID
        if sire_id != missing_parent and str(sire_id) != str(missing_parent):
            sire_idx[i] = int(sire_id) - 1
        if dam_id != missing_parent and str(dam_id) != str(missing_parent):
            dam_idx[i] = int(dam_id) - 1

    try:
        logging.info("Allocating vectors in inbreeding_meuwissen_luo().")
        lvec = [0.0] * n
        dvec = [0.0] * n
        fvec = [0.0] * n
    except Exception as e:
        logging.error(
            f"Unable to allocate vectors in inbreeding_meuwissen_luo(): {str(e)}"
        )
        return False

    if pedobj.kw.get("debug_messages", False):
        print("[DEBUG]: Starting loop over pedigree with", n, "animals")

    heappush = heapq.heappush
    heappop = heapq.heappop

    for i in range(n):
        if pedobj.kw.get("debug_messages", False):
            print(f"\t[DEBUG]: Initializing for animal {animal_id[i]} (idx: {i})")

        si = sire_idx[i]
        di = dam_idx[i]
        if si >= 0 and di >= 0:
            dvec[i] = 0.5 - 0.25 * (fvec[si] + fvec[di])
        elif si < 0 and di < 0:
            dvec[i] = 1.0
        elif si >= 0:
            dvec[i] = 0.5 - 0.25 * (fvec[si] - 1.0)
        else:
            dvec[i] = 0.5 - 0.25 * (fvec[di] - 1.0)

        lvec[i] = 1.0
        # Each ancestor is processed once. If already queued, only L is updated
        # (Meuwissen and Luo 1992).
        queued = {i}
        heap = [-i]
        touched = [i]
        aii = 0.0

        while heap:
            jidx = -heappop(heap)
            if jidx not in queued:
                continue
            queued.remove(jidx)

            lj = lvec[jidx]
            parent_sire = sire_idx[jidx]
            if parent_sire >= 0:
                if lvec[parent_sire] == 0.0:
                    touched.append(parent_sire)
                lvec[parent_sire] += 0.5 * lj
                if parent_sire not in queued:
                    queued.add(parent_sire)
                    heappush(heap, -parent_sire)
            parent_dam = dam_idx[jidx]
            if parent_dam >= 0:
                if lvec[parent_dam] == 0.0:
                    touched.append(parent_dam)
                lvec[parent_dam] += 0.5 * lj
                if parent_dam not in queued:
                    queued.add(parent_dam)
                    heappush(heap, -parent_dam)

            aii += lj * lj * dvec[jidx]

        fvec[i] = aii - 1.0
        for idx in touched:
            lvec[idx] = 0.0

    fx = {animal_id[i]: fvec[i] for i in range(n)}

    logging.info("Exited inbreeding_meuwissen_luo()")
    return fx


def inbreeding_modified_meuwissen_luo(pedobj, gens=0, **kw):
    """
    inbreeding_modified_meuwissen_luo() computes CoI using the method of Meuwissen
    and Luo (1992) as modified by Quaas (1995). It calculates only inbreeding coefficients,
    not relationships. This code is a direct implementation of the algorithm presented
    in Appendix B.2 of Mrode (2005). Mrode cites Quaas's method as: Quaas, R. L. 1995. Fx
    algorithms. An unpublished note.
    """
    try:
        logging.info("Entered inbreeding_modified_meuwissen_luo()")
    except Exception:
        pass

    # Setup dictionary to accumulate coefficients of inbreeding
    fx = {p.animalID: 0.0 for p in pedobj.pedigree}

    # Allocate memory for vectors
    try:
        logging.info("Allocating vectors in inbreeding_modified_meuwissen_luo().")
        lvecs = np.zeros(len(pedobj.pedigree), dtype=np.float64)
        lvecd = np.zeros(len(pedobj.pedigree), dtype=np.float64)
        avec = np.zeros(len(pedobj.pedigree), dtype=np.float64)
        dvec = np.zeros(len(pedobj.pedigree), dtype=np.float64)
    except MemoryError:
        logging.info(
            "Unable to allocate vectors in RAM, trying to allocate memory-mapped files."
        )
        lvecs = np.memmap("lvecs_memmap.bin", dtype="float64", mode="w+", shape=(len(pedobj.pedigree),))
        lvecd = np.memmap("lvecd_memmap.bin", dtype="float64", mode="w+", shape=(len(pedobj.pedigree),))
        avec = np.memmap("avec_memmap.bin", dtype="float64", mode="w+", shape=(len(pedobj.pedigree),))
        dvec = np.memmap("dvec_memmap.bin", dtype="float64", mode="w+", shape=(len(pedobj.pedigree),))
    except Exception as e:
        logging.error(
            f"Unable to allocate vectors in inbreeding_modified_meuwissen_luo(): {str(e)}"
        )
        return False

    if pedobj.kw.get("debug_messages", False):
        print(f"[DEBUG]: Starting loop over pedigree with {len(pedobj.pedigree)} animals")

    for i in range(len(pedobj.pedigree)):
        if pedobj.kw.get("debug_messages", False):
            print(
                f"\t[DEBUG]: Initializing local data structures for animal {pedobj.pedigree[i].animalID} (idx: {i})"
            )

        ancs = []
        ancd = []
        lvecs.fill(0.0)
        lvecd.fill(0.0)

        # Compute dvec based on parent IDs
        sire_id = pedobj.pedigree[i].sireID
        dam_id = pedobj.pedigree[i].damID
        missing_parent = pedobj.kw["missing_parent"]

        if sire_id != missing_parent and dam_id != missing_parent:
            dvec[i] = 0.5 - 0.25 * (fx[sire_id] + fx[dam_id])
        elif sire_id == missing_parent and dam_id == missing_parent:
            dvec[i] = 1.0
        elif sire_id != missing_parent:
            dvec[i] = 0.5 - 0.25 * (fx[sire_id] - 1.0)
        else:
            dvec[i] = 0.5 - 0.25 * (fx[dam_id] - 1.0)

        if sire_id != missing_parent and sire_id - 1 not in ancs:
            ancs.append(sire_id - 1)
            lvecs[sire_id - 1] = 1.0

        if dam_id != missing_parent and dam_id - 1 not in ancd:
            ancd.append(dam_id - 1)
            lvecd[dam_id - 1] = 1.0

        # Process ancestors
        while ancs and ancd:
            j = max(ancs)
            k = max(ancd)

            # An ancestor reachable by more than one path must appear in the
            # ancestor list ONCE. The migration dropped the "not in" membership
            # guard, so such an ancestor was pushed repeatedly and its
            # contribution counted once per path -- overstating F by 43% on
            # hartlandclark.ped and producing frankly impossible values (4.72,
            # 1.0) elsewhere. The L accumulation below stays unconditional, as
            # in PyPedal 2.0.4: it is the list membership that must be unique,
            # not the contribution.
            if j > k:
                if pedobj.pedigree[j].sireID != missing_parent:
                    if pedobj.pedigree[j].sireID - 1 not in ancs:
                        ancs.append(pedobj.pedigree[j].sireID - 1)
                    lvecs[pedobj.pedigree[j].sireID - 1] += 0.5 * lvecs[j]

                if pedobj.pedigree[j].damID != missing_parent:
                    if pedobj.pedigree[j].damID - 1 not in ancs:
                        ancs.append(pedobj.pedigree[j].damID - 1)
                    lvecs[pedobj.pedigree[j].damID - 1] += 0.5 * lvecs[j]

                ancs.remove(j)

            elif k > j:
                if pedobj.pedigree[k].sireID != missing_parent:
                    if pedobj.pedigree[k].sireID - 1 not in ancd:
                        ancd.append(pedobj.pedigree[k].sireID - 1)
                    lvecd[pedobj.pedigree[k].sireID - 1] += 0.5 * lvecd[k]

                if pedobj.pedigree[k].damID != missing_parent:
                    if pedobj.pedigree[k].damID - 1 not in ancd:
                        ancd.append(pedobj.pedigree[k].damID - 1)
                    lvecd[pedobj.pedigree[k].damID - 1] += 0.5 * lvecd[k]

                ancd.remove(k)

            else:
                if pedobj.pedigree[j].sireID != missing_parent:
                    if pedobj.pedigree[j].sireID - 1 not in ancs:
                        ancs.append(pedobj.pedigree[j].sireID - 1)
                    lvecs[pedobj.pedigree[j].sireID - 1] += 0.5 * lvecs[j]

                if pedobj.pedigree[j].damID != missing_parent:
                    if pedobj.pedigree[j].damID - 1 not in ancs:
                        ancs.append(pedobj.pedigree[j].damID - 1)
                    lvecs[pedobj.pedigree[j].damID - 1] += 0.5 * lvecs[j]

                if pedobj.pedigree[k].sireID != missing_parent:
                    if pedobj.pedigree[k].sireID - 1 not in ancd:
                        ancd.append(pedobj.pedigree[k].sireID - 1)
                    lvecd[pedobj.pedigree[k].sireID - 1] += 0.5 * lvecd[k]

                if pedobj.pedigree[k].damID != missing_parent:
                    if pedobj.pedigree[k].damID - 1 not in ancd:
                        ancd.append(pedobj.pedigree[k].damID - 1)
                    lvecd[pedobj.pedigree[k].damID - 1] += 0.5 * lvecd[k]

                fx[pedobj.pedigree[i].animalID] += lvecs[j] * lvecd[k] * 0.5 * dvec[j]

                ancs.remove(j)
                ancd.remove(k)

    # Clean-up allocated resources
    del lvecs, lvecd, avec, dvec
    logging.info("Exited inbreeding_modified_meuwissen_luo()")
    return fx


def a_decompose(pedobj):
    """
    Form the decomposed form of A, TDT', directly from a pedigree (after
    Henderson, 1976; Thompson, 1977; Mrode, 1996).  Return D, a diagonal
    matrix, and T, a lower triangular matrix such that A = TDT'.
    """
    logging.info("Entered a_decompose()")

    l = pedobj.metadata.num_records

    if not (pedobj.kw['form_nrm'] and pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records):
        if pedobj.kw['nrm_method'] == 'nrm':
            a = fast_a_matrix(pedobj.pedigree, pedobj.kw, method=pedobj.kw['matrix_type'])
        else:
            a = fast_a_matrix_r(pedobj.pedigree, pedobj.kw, method=pedobj.kw['matrix_type'])
    else:
        a = pedobj.nrm

    try:
        T = np.identity(l, dtype=float)
        D = np.identity(l, dtype=float)
        for row in range(l):
            for col in range(row + 1):
                # Cast these because items are read from the pedigree file as characters, not integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)

                if (
                    pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']
                    and pedobj.pedigree[row].damID == pedobj.kw['missing_parent']
                ):
                    if row == col:
                        # Both parents unknown and assumed unrelated
                        T[row, col] = 1.0
                        D[row, col] = 1.0
                    else:
                        T[row, col] = 0.0
                elif pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']:
                    # Sire unknown, dam known
                    if row == col:
                        T[row, col] = 1.0
                        fd = a[pedobj.pedigree[row].damID - 1, pedobj.pedigree[row].damID - 1] - 1.0
                        D[row, col] = 0.75 - (0.5 * fd)
                    else:
                        T[row, col] = 0.5 * T[pedobj.pedigree[row].damID - 1, col]
                elif pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    # Sire known, dam unknown
                    if row == col:
                        T[row, col] = 1.0
                        fs = a[pedobj.pedigree[row].sireID - 1, pedobj.pedigree[row].sireID - 1] - 1.0
                        D[row, col] = 0.75 - (0.5 * fs)
                    else:
                        T[row, col] = 0.5 * T[pedobj.pedigree[row].sireID - 1, col]
                elif (
                    pedobj.pedigree[row].sireID != pedobj.kw['missing_parent']
                    and pedobj.pedigree[row].damID != pedobj.kw['missing_parent']
                ):
                    # Both parents known
                    if row == col:
                        T[row, col] = 1.0
                        fs = a[pedobj.pedigree[row].sireID - 1, pedobj.pedigree[row].sireID - 1] - 1.0
                        fd = a[pedobj.pedigree[row].damID - 1, pedobj.pedigree[row].damID - 1] - 1.0
                        D[row, col] = 0.5 - (0.25 * (fs + fd))
                    else:
                        T[row, col] = 0.5 * (
                            T[pedobj.pedigree[row].sireID - 1, col]
                            + T[pedobj.pedigree[row].damID - 1, col]
                        )
                else:
                    logging.error(
                        f"[ERROR]: There is a problem with the sire (ID {pedobj.pedigree[col].sireID}) and/or "
                        f"dam (ID {pedobj.pedigree[col].damID}) of animal {pedobj.pedigree[col].animalID}"
                    )
                    break
    except Exception as e:
        logging.error(f"An error occurred in a_decompose: {e}")
        D = np.identity(1, dtype=float)
        T = np.identity(1, dtype=float)

    # Save D matrix to file
    outputfile = f"{pedobj.kw['filetag']}_a_decompose_d_.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ",".join(f"{D[row, col]:7.5f}" for col in range(l))
            aout.write(f"{line}\n")

    # Save T matrix to file
    outputfile = f"{pedobj.kw['filetag']}_a_decompose_t_.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ",".join(f"{T[row, col]:7.5f}" for col in range(l))
            aout.write(f"{line}\n")

    logging.info("Exited a_decompose()")
    return D, T


def form_d_nof(pedobj):
    """
    Form the diagonal matrix, D, used in decomposing A and forming the direct
    inverse of A. This function does not write output to a file - if you need D in
    a file, use the a_decompose() function. form_d() is a convenience function
    used by other functions. Note that inbreeding is not considered in the
    formation of D.
    """
    try:
        logging.info("Entered form_d_nof()")
    except Exception:
        pass

    try:
        l = pedobj.metadata.num_records
        D = np.identity(l, dtype=float)
        for row in range(l):
            for col in range(row + 1):
                # Cast these because items are read from the pedigree file as strings, not integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)

                if (
                    pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']
                    and pedobj.pedigree[row].damID == pedobj.kw['missing_parent']
                ):
                    if row == col:
                        # Both parents unknown and assumed unrelated
                        D[row, col] = 1.0
                elif pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']:
                    # Sire unknown, dam known
                    if row == col:
                        D[row, col] = 0.75
                elif pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    # Sire known, dam unknown
                    if row == col:
                        D[row, col] = 0.75
                elif (
                    pedobj.pedigree[row].sireID != pedobj.kw['missing_parent']
                    and pedobj.pedigree[row].damID != pedobj.kw['missing_parent']
                ):
                    # Both parents known
                    if row == col:
                        D[row, col] = 0.5
                else:
                    logging.error(
                        f"[ERROR]: There is a problem with the sire (ID {pedobj.pedigree[col].sireID}) and/or "
                        f"dam (ID {pedobj.pedigree[col].damID}) of animal {pedobj.pedigree[col].animalID}"
                    )
                    break
    except Exception as e:
        logging.error(f"An error occurred in form_d_nof: {e}")
        D = np.identity(1, dtype=float)

    logging.info("Exited form_d_nof()")
    return D


def a_inverse_dnf(pedobj, filetag='_a_inverse_dnf_'):
    """
    Form the inverse of A directly using the method of Henderson (1976) which
    does not account for inbreeding.
    """
    try:
        logging.info('Entered a_inverse_dnf()')
    except Exception:
        pass

    l = pedobj.metadata.num_records
    try:
        # Grab the diagonal matrix, d, and form its inverse
        d_inv = form_d_nof(pedobj)
        for i in range(l):
            d_inv[i, i] = 1.0 / d_inv[i, i]
        a_inv = np.zeros((l, l), dtype=float)
        for i in range(l):
            # Cast these because items are read from the pedigree file as strings, not integers
            pedobj.pedigree[i].animalID = int(pedobj.pedigree[i].animalID)
            pedobj.pedigree[i].sireID = int(pedobj.pedigree[i].sireID)
            pedobj.pedigree[i].damID = int(pedobj.pedigree[i].damID)
            s = pedobj.pedigree[i].sireID - 1
            d = pedobj.pedigree[i].damID - 1

            if pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # Both parents unknown and assumed unrelated
                a_inv[i, i] += d_inv[i, i]
            elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent']:
                # Sire unknown, dam known
                a_inv[i, i] += d_inv[i, i]
                a_inv[d, i] += (-0.5) * d_inv[i, i]
                a_inv[i, d] += (-0.5) * d_inv[i, i]
                a_inv[d, d] += 0.25 * d_inv[i, i]
            elif pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # Sire known, dam unknown
                a_inv[i, i] += d_inv[i, i]
                a_inv[s, i] += (-0.5) * d_inv[i, i]
                a_inv[i, s] += (-0.5) * d_inv[i, i]
                a_inv[s, s] += 0.25 * d_inv[i, i]
            elif pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
                # Both parents known
                a_inv[i, i] += d_inv[i, i]
                a_inv[s, i] += (-0.5) * d_inv[i, i]
                a_inv[i, s] += (-0.5) * d_inv[i, i]
                a_inv[d, i] += (-0.5) * d_inv[i, i]
                a_inv[i, d] += (-0.5) * d_inv[i, i]
                a_inv[s, s] += 0.25 * d_inv[i, i]
                a_inv[s, d] += 0.25 * d_inv[i, i]
                a_inv[d, s] += 0.25 * d_inv[i, i]
                a_inv[d, d] += 0.25 * d_inv[i, i]
            else:
                logging.error(
                    f"[ERROR]: There is a problem with the sire (ID {pedobj.pedigree[i].sireID}) "
                    f"and/or dam (ID {pedobj.pedigree[i].damID}) of animal {pedobj.pedigree[i].animalID}"
                )
                break
    except Exception as e:
        logging.error(f"An error occurred in a_inverse_dnf: {e}")
        a_inv = np.zeros((1, 1), dtype=float)

    # Write the inverse matrix to a file
    outputfile = f"{pedobj.kw['filetag']}_a_inverse_dnf_a_inv.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ','.join(f"{a_inv[row, col]:7.5f}" for col in range(l)) + '\n'
            aout.write(line)

    # Write the inverse of D to a file
    outputfile = f"{pedobj.kw['filetag']}_a_inverse_dnf_d_inv.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ','.join(f"{d_inv[row, col]:7.5f}" for col in range(l)) + '\n'
            aout.write(line)

    logging.info('Exited a_inverse_dnf()')
    return a_inv


def a_inverse_df(pedobj):
    """
    Directly form the inverse of A from the pedigree file - accounts for
    inbreeding - using the method of Quaas (1976).
    """
    try:
        logging.info('Entered a_inverse_df()')
    except Exception:
        pass

    l = pedobj.metadata.num_records

    try:
        # Initialize matrices
        d_inv = np.zeros((l, l), dtype=float)
        a_inv = np.zeros((l, l), dtype=float)
        LL = np.zeros((l, l), dtype=float)

        # Form L and D-inverse
        for row in range(l):
            for col in range(row + 1):
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)

                s = pedobj.pedigree[row].sireID - 1
                d = pedobj.pedigree[row].damID - 1
                s_sq = 0.0
                d_sq = 0.0

                if row == col:
                    for m in range(s + 1):
                        s_sq += LL[s, m] ** 2
                    s_sq *= 0.25
                    for m in range(d + 1):
                        d_sq += LL[d, m] ** 2
                    d_sq *= 0.25
                    LL[row, col] = sqrt(1.0 - s_sq - d_sq)
                    d_inv[row, col] = 1.0 / (LL[row, col] ** 2)
                else:
                    LL[row, col] = 0.5 * (LL[s, col] + LL[d, col])

        # Use D-inverse to compute A-inverse
        for i in range(l):
            s = pedobj.pedigree[i].sireID - 1
            d = pedobj.pedigree[i].damID - 1

            if pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # Both parents unknown and assumed unrelated
                a_inv[i, i] += d_inv[i, i]
            elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent']:
                # Sire unknown, dam known
                a_inv[i, i] += d_inv[i, i]
                a_inv[d, i] += -0.5 * d_inv[i, i]
                a_inv[i, d] += -0.5 * d_inv[i, i]
                a_inv[d, d] += 0.25 * d_inv[i, i]
            elif pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # Sire known, dam unknown
                a_inv[i, i] += d_inv[i, i]
                a_inv[s, i] += -0.5 * d_inv[i, i]
                a_inv[i, s] += -0.5 * d_inv[i, i]
                a_inv[s, s] += 0.25 * d_inv[i, i]
            elif pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
                # Both parents known
                a_inv[i, i] += d_inv[i, i]
                a_inv[s, i] += -0.5 * d_inv[i, i]
                a_inv[i, s] += -0.5 * d_inv[i, i]
                a_inv[d, i] += -0.5 * d_inv[i, i]
                a_inv[i, d] += -0.5 * d_inv[i, i]
                a_inv[s, s] += 0.25 * d_inv[i, i]
                a_inv[s, d] += 0.25 * d_inv[i, i]
                a_inv[d, s] += 0.25 * d_inv[i, i]
                a_inv[d, d] += 0.25 * d_inv[i, i]
    except Exception as e:
        logging.error(f"An error occurred in a_inverse_df: {e}")
        a_inv = np.zeros((1, 1), dtype=float)

    # Write the inverse matrix to a file
    outputfile = f"{pedobj.kw['filetag']}_a_inverse_df_a_inv.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ','.join(f"{a_inv[row, col]:7.5f}" for col in range(l)) + '\n'
            aout.write(line)

    # Write LL matrix to a file
    outputfile = f"{pedobj.kw['filetag']}_a_inverse_df_l.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ','.join(f"{LL[row, col]:7.5f}" for col in range(l)) + '\n'
            aout.write(line)

    # Write D-inverse matrix to a file
    outputfile = f"{pedobj.kw['filetag']}_a_inverse_df_d_inv.dat"
    with open(outputfile, 'w') as aout:
        for row in range(l):
            line = ','.join(f"{d_inv[row, col]:7.5f}" for col in range(l)) + '\n'
            aout.write(line)

    logging.info('Exited a_inverse_df()')
    return a_inv


def partial_inbreeding(pedobj, animals=None, gens=0, rels=1, cleanmaps=1):
    """
    partial_inbreeding() computes coefficients of partial inbreeding,
    which is the probability that an individual, i, is IDB at a locus
    and that the alleles were derived from ancestor j.
    Input: A PyPedal pedigree object.
    Output: A dictionary of partial CoI keyed to renumbered animal IDs.
    """
    if animals is None:
        animals = []

    try:
        logging.info('Entered partial_inbreeding()')
    except Exception:
        pass

    ng = pyp_network.ped_to_graph(pedobj)

    _ped = []  # Temporary pedigree
    top_ped = []

    if int(gens) > 0:
        top_peddict = pyp_network.find_ancestors_g(ng, len(pedobj.idmap), {}, gens)
        top_peddict[len(pedobj.idmap)] = 1
        top_ped = list(top_peddict.keys())
        top_r = []
        _anids = []
        for _j in top_ped:
            if top_peddict[_j] <= gens:
                top_r.append(copy.deepcopy(pedobj.pedigree[int(_j) - 1]))
                if top_peddict[_j] == 1:
                    top_r[-1].sireID = 0
                    top_r[-1].damID = 0
                _anids.append(top_r[-1].animalID)
    else:
        _anids = list(pedobj.backmap.keys())

    fx = {}  # Coefficients of inbreeding
    _parents = {}  # Stores sire-dam pairs with the youngest offspring
    _anids.sort(reverse=True)  # Youngest animals first

    # Relationship statistics setup
    if rels:
        reldict = {
            'r_count': 0,
            'r_nonzero_count': 0,
            'r_sum': 0.0,
            'r_max': 0.0,
            'r_min': 1.0
        }

    partial_inbreeding_dict = {}

    for i in _anids:
        if int(gens) == 0:
            _parent_key = f"{pedobj.pedigree[int(i) - 1].sireID}_{pedobj.pedigree[int(i) - 1].damID}"
        else:
            _animal = next((animal for animal in top_r if int(animal.animalID) == int(i)), None)
            if _animal is None:
                _animal = pedobj.pedigree[int(i) - 1]
            _parent_key = f"{_animal.sireID}_{_animal.damID}"

        if i not in fx:
            try:
                fx[i] = fx[_parents[_parent_key]]
            except KeyError:
                _tag = f"{pedobj.kw['filetag']}_{i}"
                if int(gens) > 0:
                    _ped = top_peddict
                else:
                    _ped = pyp_network.find_ancestors(ng, i, [])
                    _ped.append(i)

                _r = top_r if int(gens) > 0 else []
                if not _r:
                    for j in _ped:
                        _r.append(copy.deepcopy(pedobj.pedigree[int(j) - 1]))

                _r = pyp_utils.reorder(_r, _tag,
                                       missingparent=pedobj.kw['missing_parent'])
                _s, _map = pyp_utils.renumber(
                    _r, _tag, returnmap=True, debug=pedobj.kw['debug_messages'],
                    missingparent=pedobj.kw['missing_parent'],
                    animaltype=pedobj.kw['animal_type']
                )
                _flist = pyp_utils.founders_from_list(_r, pedobj.kw['missing_parent'])

                _backmap = {v: k for k, v in _map.items()}
                _opts = copy.deepcopy(pedobj.kw)

                for _f in _flist:
                    _opts['filetag'] = f"{_tag}_{_f}"
                    _f_dict = fast_partial_a_matrix(_s, _f, _flist, _opts, method=pedobj.kw['matrix_type'])
                    for k, v in _f_dict.items():
                        if k not in partial_inbreeding_dict:
                            partial_inbreeding_dict[k] = {}
                        partial_inbreeding_dict[k].update(v)

                if cleanmaps:
                    pyp_utils.delete_id_map(_tag)

    try:
        logging.info('Exited partial_inbreeding()')
    except Exception:
        pass

    return partial_inbreeding_dict


def fast_partial_a_matrix(pedigree, founder, founderlist, pedopts, method='dense', debug=0):
    """
    fast_partial_a_matrix() calculates a partial kinship matrix for a given
    founder in a pedigree and returns a dictionary of partial inbreeding
    coefficients between that founder and descendants (non-founders) in
    the pedigree.
    """
    _animals = {}
    _sires = {}
    _dams = {}
    l = len(pedigree)

    if method not in ['dense', 'sparse']:
        method = 'dense'

    # Use PySparse for large matrices, otherwise use NumPy
    if method == 'sparse':
        try:
            a = lil_matrix((l, l))
        except Exception:
            if pedopts['messages'] != 'quiet':
                print('[WARNING]: Could not create a sparse identity matrix in fast_partial_a_matrix()!')
            logging.warning('Could not create a sparse identity matrix in fast_partial_a_matrix()!.')
            return fast_partial_a_matrix(pedigree, founder, founderlist, pedopts, method='dense', debug=debug)
    else:
        try:
            a = np.zeros((l, l), dtype=np.float64)
        except MemoryError:
            a = np.memmap('fast_partial_a_matrix_mmap.bin', dtype='float32', mode='w+', shape=(l, l))
        except Exception:
            print(f'[ERROR]: Unable to allocate a matrix of rank {l} in fast_partial_a_matrix()!')
            logging.error(f'Unable to allocate a matrix of rank {l} in fast_partial_a_matrix()!')
            return False

    # Initialize animal, sire, and dam lists
    if pedopts.get('debug_messages') and pedopts['messages'] != 'quiet':
        print(f'\t\t[pyp_nrm/fast_partial_a_matrix()] Started forming animal, sire, and dam lists at {pyp_utils.pyp_nice_time()}')
    
    for i in range(l):
        _animals[i] = int(pedigree[i].animalID)
        _sires[i] = int(pedigree[i].sireID)
        _dams[i] = int(pedigree[i].damID)

    if pedopts.get('debug_messages') and pedopts['messages'] != 'quiet':
        print(f'\t\t[pyp_nrm/fast_partial_a_matrix()] Finished forming animal, sire, and dam lists at {pyp_utils.pyp_nice_time()}')
        print(f'\t\t[pyp_nrm/fast_partial_a_matrix()] Started computing A at {pyp_utils.pyp_nice_time()}')

    partial_f = {}

    if debug:
        print('n_founders:', len(founderlist))
        print('founder:', founder)

    # Step 1: Initialize the matrix
    fidx = founderlist.index(founder)
    a[fidx, fidx] = 0.5

    # Step 2: Fill the matrix
    for row in range(len(founderlist), l):
        for col in range(row, l):
            sire = _sires[col]
            dam = _dams[col]
            if sire == pedopts['missing_parent'] and dam == pedopts['missing_parent']:
                continue
            elif sire == pedopts['missing_parent']:
                if row != col:
                    a[row, col] = 0.5 * a[row, dam - 1]
                    a[col, row] = a[row, col]
            elif dam == pedopts['missing_parent']:
                if row != col:
                    a[row, col] = 0.5 * a[row, sire - 1]
                    a[col, row] = a[row, col]
            else:
                if row == col:
                    a[row, row] = a[row, fidx] + 0.5 * a[sire - 1, dam - 1]
                else:
                    a[row, col] = 0.5 * (a[row, sire - 1] + a[row, dam - 1])
                    a[col, row] = a[row, col]

    # Step 3: Compute partial inbreeding coefficients
    partial_f[founder] = {}
    for i in range(len(founderlist), l):
        partial_f[founder][_animals[i]] = a[_sires[i] - 1, _dams[i] - 1]

    if debug:
        np.set_printoptions(precision=4, linewidth=100)
        print(a)
        print(partial_f)

    if pedopts.get('debug_messages') and pedopts['messages'] != 'quiet':
        print(f'\t\t[pyp_nrm/fast_partial_a_matrix()] Finished computing A at {pyp_utils.pyp_nice_time()}')

    return partial_f
