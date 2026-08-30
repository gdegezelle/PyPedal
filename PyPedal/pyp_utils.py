#!/usr/bin/env python3

"""
pyp_utils.py - A module for creating and managing PyPedal pedigrees.

Version: see PyPedal.__version__
Author: John B. Cole (john.b.cole@gmail.com)
License: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later

This module contains utilities for handling pedigree data in PyPedal. It provides 
functions for reordering, renumbering, and modifying pedigrees, along with additional 
utilities for related operations.

Functions:
- set_ancestor_flag()
- set_generation()
- set_age()
- set_species()
- set_sexes()
- assign_sexes()
- set_offspring()
- assign_offspring()
- set_upg()
- assign_upg()
- reorder()
- fast_reorder()
- renumber()
- load_id_map()
- delete_id_map()
- id_map_new_to_old()
- trim_pedigree_to_year()
- pedigree_range()
- sort_dict_by_keys()
- sort_dict_by_values()
- simple_histogram_dictionary()
- reverse_string()
- pyp_nice_time()
- string_to_table_name()
- pyp_datestamp()
- subpedigree()
- founders_from_list()
- founder_allele_dict()
- list_union()
- list_intersection()
- guess_pedformat()
- list_duplicates()
- list_likely_same_animals()
- which()
- remove_missing()
"""

import os
import logging
import copy
import heapq
import math
import time
from typing import Dict, List, Tuple, Union, TYPE_CHECKING

from . import pyp_chronology, pyp_demog, pyp_newclasses
from PyPedal import pyp_errors
from PyPedal.pyp_newclasses import NewPedigree, NewAnimal
from functools import cmp_to_key

if TYPE_CHECKING:
    from PyPedal.pyp_newclasses import NewPedigree

def set_ancestor_flag(pedobj):
    """
    Loops through a pedigree to identify all parents and set ancestor flags.
    Expects a reordered and renumbered pedigree as input.
    
    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    try:
        # Ensure there are enough records to process
        if len(pedobj.pedigree) < 2:
            print('[ERROR]: pedobj.pedigree contains fewer than two records; nothing to do.')
            return False

        parents = set()  # Use a set for faster lookups and to avoid duplicates

        # Process the pedigree from young to old
        for animal in reversed(pedobj.pedigree):
            if pedobj.kw.get('messages') == 'debug':
                print(f"[DEBUG]: animal: {animal.animalID}, sire: {animal.sireID}, dam: {animal.damID}")

            # Add sire to parents if valid
            if animal.sireID != pedobj.kw['missing_parent']:
                sire_id = int(animal.sireID)
                if sire_id not in parents:
                    parents.add(sire_id)
                    pedobj.pedigree[sire_id - 1].ancestor = 1

            # Add dam to parents if valid
            if animal.damID != pedobj.kw['missing_parent']:
                dam_id = int(animal.damID)
                if dam_id not in parents:
                    parents.add(dam_id)
                    pedobj.pedigree[dam_id - 1].ancestor = 1

        # Write ancestors to a file if file I/O is enabled
        if pedobj.kw.get('file_io'):
            a_outputfile = f"{pedobj.kw['filetag']}_ancestors.dat"
            try:
                with open(a_outputfile, 'w') as aout:
                    aout.write(f"# FILE: {a_outputfile}\n")
                    aout.write("# ANCESTOR list produced by PyPedal.\n")
                    aout.writelines(f"{parent_id}\n" for parent_id in sorted(parents))
                logging.info(f"Ancestor list written to {a_outputfile}.")
            except IOError as e:
                logging.error(f"Could not write to {a_outputfile}: {e}")
                return False

        return True

    except Exception as e:
        logging.error(f"An error occurred in set_ancestor_flag: {e}")
        return False


def set_generation(pedobj):
    """
    set_generation() Works through a pedigree to infer the generation to which an animal
    belongs based on founders belonging to generation 1. The igen assigned to an animal
    as the larger of sire.igen+1 and dam.igen+1. This routine assumes that myped is
    reordered and renumbered.
    
    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    try:
        if pedobj.kw.get('messages') == 'debug':
            print(f"[NOTE]: pyp_utils/set_generation() assigning inferred generations in pedigree {pedobj.kw['pedname']}.")

        if pedobj.kw.get('gen_coeff', False):
            raise pyp_errors.PyPedalUsageError(
                "set_generation: kw['gen_coeff'] calculation is outside the "
                "PyPedal 4.0 domain. Pedformat 'p' may store a supplied "
                "generation coefficient; PyPedal does not compute Pattie "
                "(1965) coefficients. Leave gen_coeff at False."
            )

        for i in range(pedobj.metadata.num_records):
            animal = pedobj.pedigree[i]

            if animal.sireID == pedobj.kw['missing_parent'] and animal.damID == pedobj.kw['missing_parent']:
                # Case: Founders
                animal.igen = 1
            elif animal.sireID == pedobj.kw['missing_parent']:
                # Case: Missing sire
                dam = pedobj.pedigree[int(animal.damID) - 1]
                animal.igen = dam.igen + 1
            elif animal.damID == pedobj.kw['missing_parent']:
                # Case: Missing dam
                sire = pedobj.pedigree[int(animal.sireID) - 1]
                animal.igen = sire.igen + 1
            else:
                # Case: Both parents are known
                sire = pedobj.pedigree[int(animal.sireID) - 1]
                dam = pedobj.pedigree[int(animal.damID) - 1]
                animal.igen = max(sire.igen + 1, dam.igen + 1)

        logging.info(f"pyp_utils/set_generation() assigned inferred generations in pedigree {pedobj.kw['pedname']}.")
        return True

    except pyp_errors.PyPedalError:
        raise
    except Exception as e:
        logging.error(f"pyp_utils/set_generation() was unable to assign inferred generations in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def set_age(pedobj):
    """
    Assigns the historical demographic year-offset
    ``animal.age = recorded_year - BASE_DEMOGRAPHIC_YEAR``.

    ``animal.age`` is **not** biological or current age. Unknown recorded
    birth year uses ``kw['missing_age']`` (falling back to ``missing_value``).
    Inferred generation depth (``igen``) is never used as age.

    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    try:
        if pedobj.kw.get('messages') == 'debug':
            print(f"[NOTE]: pyp_utils/set_age() assigning year-offsets in pedigree {pedobj.kw['pedname']}.")

        missing_age = pedobj.kw.get('missing_age', pedobj.kw.get('missing_value', -999))
        for animal in pedobj.pedigree:
            year = pyp_chronology.recorded_year(animal)
            if year is None:
                animal.age = missing_age
            else:
                animal.age = year - pyp_demog.BASE_DEMOGRAPHIC_YEAR

        logging.info(f"pyp_utils/set_age() assigned year-offsets in pedigree {pedobj.kw['pedname']}.")
        return True

    except Exception as e:
        logging.error(f"pyp_utils/set_age() was unable to assign ages in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def set_species(pedobj, species='Unknown'):
    """
    Assigns a species to every animal in the pedigree.

    :param pedobj: A PyPedal NewPedigree object.
    :param species: A PyPedal string specifying the species.
    :return: False for failure and True for success.
    """
    try:
        if pedobj.kw.get('messages') == 'debug':
            print(f"[NOTE]: pyp_utils/set_species() assigning species '{species}' to all animals in pedigree {pedobj.kw['pedname']}.")

        for animal in pedobj.pedigree:
            animal.species = species if species else 'u'

        logging.info(f"pyp_utils/set_species() assigned species in pedigree {pedobj.kw['pedname']}.")
        return True

    except Exception as e:
        logging.error(f"pyp_utils/set_species() was unable to assign species in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def set_sexes(pedobj):
    """
    Assigns a sex to every animal in the pedigree using sire and daughter lists for improved accuracy.

    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    try:
        if pedobj.kw.get('messages') == 'verbose':
            print(f"[NOTE]: pyp_utils/set_sexes() assigning a sex to all animals in pedigree {pedobj.kw['pedname']}.")

        for animal in pedobj.pedigree:
            sire_id = animal.sireID
            dam_id = animal.damID

            if sire_id == pedobj.kw['missing_parent'] and dam_id == pedobj.kw['missing_parent']:
                continue

            if sire_id == pedobj.kw['missing_parent']:
                dam = pedobj.pedigree[int(dam_id) - 1]
                if dam.sex != 'f':
                    if pedobj.kw.get('debug_messages'):
                        print(f"\t\tAnimal {dam_id} sex changed from {dam.sex} to 'f'")
                    dam.sex = 'f'

            elif dam_id == pedobj.kw['missing_parent']:
                sire = pedobj.pedigree[int(sire_id) - 1]
                if sire.sex != 'm':
                    if pedobj.kw.get('debug_messages'):
                        print(f"\t\tAnimal {sire_id} sex changed from {sire.sex} to 'm'")
                    sire.sex = 'm'

            else:
                dam = pedobj.pedigree[int(dam_id) - 1]
                sire = pedobj.pedigree[int(sire_id) - 1]
                
                if dam.sex != 'f':
                    if pedobj.kw.get('debug_messages'):
                        print(f"\t\tAnimal {dam_id} sex changed from {dam.sex} to 'f'")
                    dam.sex = 'f'

                if sire.sex != 'm':
                    if pedobj.kw.get('debug_messages'):
                        print(f"\t\tAnimal {sire_id} sex changed from {sire.sex} to 'm'")
                    sire.sex = 'm'

        logging.info(f"pyp_utils/set_sexes() assigned sexes in pedigree {pedobj.kw['pedname']}.")
        return True

    except Exception as e:
        logging.error(f"pyp_utils/set_sexes() was unable to assign sexes in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def assign_sexes(pedobj):
    """
    assign_sexes() assigns a sex to every animal in the pedigree using sire and daughter
    lists for improved accuracy.

    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    print("[DEPRECATION WARNING]: pyp_utils/assign_sexes() will be replaced by pyp_utils/set_sexes() in version 2.1!")
    try:
        set_sexes(pedobj)
        return True
    except Exception as e:
        logging.error(f"pyp_utils/assign_sexes() was unable to assign sexes in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def set_offspring(pedobj):
    """
    set_offspring() assigns offspring to their parent(s)'s unknown sex offspring list
    (well, dictionary).

    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    try:
        if pedobj.kw.get('messages') == 'debug':
            print(f"[NOTE]: pyp_utils/set_offspring() assigning offspring to all parents in pedigree {pedobj.kw['pedname']}.")
        
        # Initialize offspring dictionaries
        for animal in pedobj.pedigree:
            pedobj.pedigree[int(animal.animalID) - 1].sons = {}
            pedobj.pedigree[int(animal.animalID) - 1].daus = {}
            pedobj.pedigree[int(animal.animalID) - 1].unks = {}

        if 'x' not in pedobj.kw['pedformat']:
            # Assign offspring to unknown sex list if sex information is unavailable
            for animal in pedobj.pedigree:
                if animal.sireID == pedobj.kw['missing_parent'] and animal.damID == pedobj.kw['missing_parent']:
                    continue
                if animal.sireID == pedobj.kw['missing_parent']:
                    pedobj.pedigree[int(animal.damID) - 1].unks[animal.animalID] = animal.animalID
                elif animal.damID == pedobj.kw['missing_parent']:
                    pedobj.pedigree[int(animal.sireID) - 1].unks[animal.animalID] = animal.animalID
                else:
                    pedobj.pedigree[int(animal.damID) - 1].unks[animal.animalID] = animal.animalID
                    pedobj.pedigree[int(animal.sireID) - 1].unks[animal.animalID] = animal.animalID
        else:
            # Assign offspring to specific lists based on known sex
            for animal in pedobj.pedigree:
                if animal.sex in {'m', 'M'}:
                    if animal.sireID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.sireID) - 1].sons[animal.animalID] = animal.animalID
                    if animal.damID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.damID) - 1].sons[animal.animalID] = animal.animalID
                elif animal.sex in {'f', 'F'}:
                    if animal.sireID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.sireID) - 1].daus[animal.animalID] = animal.animalID
                    if animal.damID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.damID) - 1].daus[animal.animalID] = animal.animalID
                else:
                    if animal.sireID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.sireID) - 1].unks[animal.animalID] = animal.animalID
                    if animal.damID != pedobj.kw['missing_parent']:
                        pedobj.pedigree[int(animal.damID) - 1].unks[animal.animalID] = animal.animalID

        logging.info(f"pyp_utils/set_offspring() assigned offspring in pedigree {pedobj.kw['pedname']}")
        return True
    except Exception as e:
        logging.error(f"pyp_utils/set_offspring() was unable to assign offspring in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def assign_offspring(pedobj):
    """
    assign_offspring() assigns offspring to their parent(s)'s specific sex offspring list.

    :param pedobj: A PyPedal NewPedigree object.
    :return: False for failure and True for success.
    """
    print("[DEPRECATION WARNING]: pyp_utils/assign_offspring() will be replaced by pyp_utils/set_offspring() in "
          "version 2.1!")
    try:
        set_offspring(pedobj)
        return True
    except Exception as e:
        logging.error(f"pyp_utils/assign_offspring() was unable to assign sexes in pedigree {pedobj.kw['pedname']}: {e}")
        return False


def set_upg(pedobj, upg_rule='asd'):
    """
    set_upg() assigns pseudo-parents to animals with unknown (missing) parents. The
    pseudo-parents can be used as Westell groups in an animal model.

    :param pedobj: A PyPedal NewPedigree object.
    :param upg_rule: A string indicating how UPG will be formed ['asd'].
    :return: False for failure and True for success.
    """
    try:
        if not pedobj.kw.get('pedigree_is_renumbered', False):
            for p in pedobj:
                if pedobj.kw.get('pedformat', '').lower() == 'asd':
                    if p.sireID == pedobj.kw['missing_parent']:
                        p.sireID = -999999
                        p.sireName = 'Sire_UPG'
                    if p.damID == pedobj.kw['missing_parent']:
                        p.damID = -888888
                        p.damName = 'Dam_UPG'
                else:
                    logging.error(f"No rule for assigning unknown parent group '{upg_rule}'!")
                    if pedobj.kw.get('message') == 'verbose':
                        print(f"[ERROR]: No rule for assigning unknown parent group '{upg_rule}'!")
            return True
        else:
            logging.error("Cannot assign unknown parent groups to a renumbered pedigree!")
            if pedobj.kw.get('message') == 'verbose':
                print("[ERROR]: Cannot assign unknown parent groups to a renumbered pedigree!")
            return False
    except Exception as e:
        logging.error(f"pyp_utils/set_upg() encountered an error: {e}")
        return False


def assign_upg(pedobj, upg_rule='asd'):
    """
    assign_upg() assigns pseudo-parents to animals with unknown (missing) parents. The
    pseudo-parents can be used as Westell groups in an animal model.

    :param pedobj: A PyPedal NewPedigree object.
    :param upg_rule: A string indicating how UPG will be formed ['asd'].
    :return: False for failure and True for success.
    """
    print('[DEPRECATION WARNING]: pyp_utils/assign_upg() will be replaced by pyp_utils/set_upg() in version 2.1!')
    try:
        set_upg(pedobj, upg_rule)
        return True
    except Exception as e:
        logging.error(f"pyp_utils/assign_upg() was unable to assign unknown parent groups in pedigree "
                      f"{pedobj.kw.get('pedname', 'unknown')}: {e}")
        return False


def _order_pedigree(myped, missingparent=0, routine='pyp_utils/reorder()'):
    """
    Order a pedigree so that every known parent precedes its offspring, or
    refuse.

    This is the single ordering engine behind both :func:`reorder` and
    :func:`fast_reorder`. It returns a NEW list holding THE SAME animal objects
    in the declared order; neither the input list nor any animal is mutated,
    which lets the two public wrappers keep their opposite aliasing contracts
    without either of them copying anything.

    The order is:

        founders first, in input order, then a stable topological order using
        original input position as the tie-break.

    Founders-first is a documented PyPedal contract, not an implementation
    detail: the 2.0 methodology notes it, and it exists so partial inbreeding
    code sees founders first. The input-position tie-break is a software
    convention for animals the pedigree graph does not order relative to one
    another; it is not a scientific claim, and it is chosen because it is the
    smallest deterministic rule that makes an already-ordered pedigree a
    fixed point.

    Ordering depends on the parent graph and on nothing else. It does not read
    ``by``, ``bd``, ``gen`` or ``igen``: the pedigree relation already defines
    the ordering constraint. Birth chronology is represented independently
    (unknown recorded years are ``None``). ``tests/test_reorder_correctness.py``
    proves that independence rather than leaving it to be inferred from this
    comment.

    Complexity is O((V + E) log V) from the heap, and a pedigree has at most two
    parental edges per animal, so O(V log V); memory is O(V + E). It replaced an
    insertion sort that was quadratic in pedigree disorder and took 297 s on
    100,000 animals. No claim of linear time is made: the heap costs a
    logarithmic factor and saying otherwise would be wrong.

    :param myped: A list of PyPedal animal objects.
    :param missingparent: The value used to indicate a missing parent.
    :param routine: Name used in refusal messages, so a caller reading the
                    message learns which entry point refused.
    :return: A new list containing the same animal objects, ordered.
    :raises PyPedalPedigreeStructureError: on a duplicated animal ID, a parent
        reference with no matching record, or a pedigree containing a cycle.
    """
    missing = str(missingparent)

    # Stringified comparison throughout, as PyPedal 2.0.4 did
    # (legacy pyp_utils.py:560-565). The Python 3 port dropped the str() calls,
    # which diverges the moment an ID's type differs from the sentinel's -- an
    # int 0 sentinel against string IDs, say. Restored here.
    position = {}
    for index, animal in enumerate(myped):
        key = str(animal.animalID)
        if key in position:
            raise pyp_errors.PyPedalPedigreeStructureError(
                '%s: animal ID %s appears more than once in the pedigree. '
                'Animal IDs must be unique; renumbering cannot be applied to a '
                'pedigree in which one ID names two records, because every '
                'reference to that ID is ambiguous.' % (routine, animal.animalID),
                animals=[animal.animalID])
        position[key] = index

    indegree = [0] * len(myped)
    children = {}
    for index, animal in enumerate(myped):
        for role, parent in (('sire', animal.sireID), ('dam', animal.damID)):
            pkey = str(parent)
            if pkey == missing:
                continue
            if pkey not in position:
                raise pyp_errors.PyPedalPedigreeStructureError(
                    '%s: animal %s names %s as its %s, but no animal with that '
                    'ID has a record in the pedigree. Ordering cannot place a '
                    'parent that is not present, and continuing would silently '
                    'discard the relationship.'
                    % (routine, animal.animalID, parent, role),
                    animals=[animal.animalID, parent])
            indegree[index] += 1
            children.setdefault(position[pkey], []).append(index)

    # Founders first, in input order. They have in-degree zero by definition, so
    # emitting all of them before anything else is always a valid prefix.
    ordered = []
    ready = []
    for index, animal in enumerate(myped):
        if indegree[index] == 0 and (str(animal.sireID) == missing
                                     and str(animal.damID) == missing):
            ordered.append(index)
    emitted = set(ordered)
    for index in ordered:
        for child in children.get(index, ()):
            indegree[child] -= 1
            if indegree[child] == 0:
                heapq.heappush(ready, child)
    # Half-founders and anything else already eligible. Pushed after the founder
    # block so they cannot displace it.
    for index in range(len(myped)):
        if index not in emitted and indegree[index] == 0:
            heapq.heappush(ready, index)

    while ready:
        index = heapq.heappop(ready)
        if index in emitted:
            continue
        ordered.append(index)
        emitted.add(index)
        for child in children.get(index, ()):
            indegree[child] -= 1
            if indegree[child] == 0:
                heapq.heappush(ready, child)

    if len(ordered) != len(myped):
        # Everything still carrying in-degree is either on a cycle or descended
        # from one. Distinguishing the two would need a second pass over the
        # graph to extract strongly connected components, and the requirement
        # here is to refuse rather than to reconstruct the cycle -- so the
        # message says "unorderable", which is true of every animal listed, and
        # does not claim they all lie on the cycle, which is not.
        residue = [myped[i].animalID for i in range(len(myped))
                   if i not in emitted]
        raise pyp_errors.PyPedalPedigreeStructureError(
            '%s: the pedigree contains a cycle or a cyclic dependency, so there '
            'is no order in which every parent precedes its offspring. These '
            '%d animal IDs remain unorderable (some lie on the cycle, others '
            'merely descend from it): %s'
            % (routine, len(residue),
               ', '.join(str(a) for a in residue[:20])
               + (', ...' if len(residue) > 20 else '')),
            animals=residue)

    return [myped[i] for i in ordered]


def reorder(myped, filetag='_reordered_', io='no', missingparent=0, debug=False, max_rounds=100):
    """
    Reorder a pedigree such that parents precede their offspring, in place.

    Founders are grouped at the beginning of the pedigree so that partial
    inbreeding code sees them first. Animals the graph does not order relative
    to one another keep their input order.

    ``myped`` is reordered **in place** and is also returned, and the animal
    objects are the ones the caller passed in. An earlier implementation
    substituted ``copy.deepcopy`` of every animal it moved, so callers holding a
    reference to an animal silently ended up with a stale one -- and a large
    load made more than a million such copies.

    Structurally invalid input is refused rather than resolved. Previously a
    pedigree that could not be ordered was returned unchanged with only a log
    message, and ``renumber()`` then deleted the offending parent links; see
    :class:`PyPedal.pyp_errors.PyPedalPedigreeStructureError`.

    :param myped: A PyPedal pedigree object.
    :param filetag: A descriptor prepended to output file names.
    :param io: Indicates whether or not to write the reordered pedigree to a file ('yes'|'no').
    :param missingparent: Indicates the value used to indicate a missing parent.
    :param debug: Toggles debugging messages on and off.
    :param max_rounds: Retained for backwards compatibility and no longer used.
                       The ordering is computed in a single pass, so there are
                       no rounds to exhaust.
    :return: The reordered PyPedal pedigree, which is ``myped``.
    :raises PyPedalPedigreeStructureError: see :func:`_order_pedigree`.
    """
    ordered = _order_pedigree(myped, missingparent=missingparent,
                              routine='pyp_utils/reorder()')
    if debug:
        print('=' * 70)
        print(f'[DEBUG]: Ordered {len(ordered)} animals; '
              f'parents now precede their offspring.')
    myped[:] = ordered

    if io == 'yes':
        a_outputfile = f'{filetag}_reordered.ped'
        with open(a_outputfile, 'w') as aout:
            aout.write(f'# FILE: {a_outputfile}\n')
            aout.write('# REORDERED pedigree produced by PyPedal.\n')
            aout.write('% asd\n')
            for m in myped:
                aout.write(f'{m.animalID}, {m.sireID}, {m.damID}\n')

    return myped


# def fast_reorder(myped, filetag='_new_reordered_', io='no', debug=False, original_filename=None):
#     """
#     Reorders a pedigree such that parents precede their offspring in the pedigree.

#     :param myped: A list of PyPedal pedigree objects.
#     :param filetag: A descriptor prepended to output file names.
#     :param io: Indicates whether to write the reordered pedigree to a file ('yes'|'no').
#     :param debug: Boolean to toggle debugging messages on and off.
#     :param original_filename: Name of the original file (to derive the output file name).
#     :return: A reordered list of PyPedal pedigree objects.
#     """

#     if not myped:
#         if debug:
#             print("DEBUG: Empty pedigree provided.")
#         return myped

#     if debug:
#         print(f"DEBUG: Starting with {len(myped)} animals in the pedigree.")

#     reordered = []
#     seen = set()
#     processing = set()

#     def process_animal(animal, original, reordered, seen, processing, debug):
#         """Process an individual animal and its parents recursively."""
#         if debug:
#             print(f"[DEBUG]: Processing animal {animal.animalID}, birth year: {animal.by}")

#         if animal.animalID in processing:
#             raise ValueError(f"Circular reference detected for animal {animal.animalID}.")

#         if animal.animalID not in seen:
#             processing.add(animal.animalID)  # Mark as currently being processed
            
#             # Recursively process sire and dam
#             if animal.sireID and animal.sireID != 0:
#                 sire = next((a for a in original if a.animalID == animal.sireID), None)
#                 if sire:
#                     process_animal(sire, original, reordered, seen, processing, debug)
#                 elif debug:
#                     print(f"DEBUG: Sire {animal.sireID} not found for animal {animal.animalID}.")

#             if animal.damID and animal.damID != 0:
#                 dam = next((a for a in original if a.animalID == animal.damID), None)
#                 if dam:
#                     process_animal(dam, original, reordered, seen, processing, debug)
#                 elif debug:
#                     print(f"DEBUG: Dam {animal.damID} not found for animal {animal.animalID}.")

#             reordered.append(animal)  # Add the current animal to the reordered list
#             seen.add(animal.animalID)  # Mark as processed
#             processing.remove(animal.animalID)  # Remove from processing stack

#     for animal in myped:
#         process_animal(animal, myped, reordered, seen, processing, debug)

#     if debug:
#         print(f"DEBUG: Reordering complete. Reordered pedigree contains {len(reordered)} animals.")

#     # Write to output files if `io` is enabled
#     if io.lower() == 'yes':
#         # Determine filenames based on the original filename or fallback to default filetag
#         base_name = os.path.splitext(os.path.basename(original_filename))[0] if original_filename else filetag
#         output_file = f"{base_name}_reordered.ped"
#         id_map_file = f"{base_name}_id_map.map"

#         if debug:
#             print(f"Writing reordered pedigree to {output_file}")
#         try:
#             with open(output_file, 'w', encoding='utf-8') as aout:
#                 aout.write(f"# FILE: {output_file}\n")
#                 aout.write("# REORDERED pedigree produced by PyPedal using fast_reorder().\n")
#                 for animal in reordered:
#                     aout.write(
#                         f"{animal.animalID},{animal.sireID},{animal.damID},{animal.gen},{animal.by}\n"
#                     )
#                     if debug:
#                         print(f"[DEBUG]: Writing animal {animal.animalID}, birth year: {animal.by}")
#             if debug:
#                 print(f"DEBUG: Reordered pedigree written to {output_file}")

#             # Write the ID map to a file
#             with open(id_map_file, 'w', encoding='utf-8') as idmap:
#                 idmap.write(f"# FILE: {id_map_file}\n")
#                 idmap.write("# Renumbered ID to Old ID mapping produced by PyPedal.\n")
#                 for animal in reordered:
#                     idmap.write(f"{animal.originalID},{animal.animalID}\n")

#         except Exception as write_err:
#             print(f"ERROR: Unable to write to file {output_file}. {write_err}")
#             raise

#     return reordered

def fast_reorder(myped, filetag='_new_reordered_', io='no', debug=False,
                 missingparent=0):
    """
    Reorders a pedigree such that parents precede their offspring in the pedigree.

    Returns a **new** list containing **the same** animal objects; ``myped`` is
    not mutated. That is the opposite of :func:`reorder`, which reorders in
    place, and both contracts are preserved deliberately because callers depend
    on each. They now share one ordering engine, so the two names no longer
    disagree about what a valid pedigree is.

    Before the shared ordering engine this routine was a recursive depth-first
    sort whose parent lookup was a linear scan of the whole pedigree per
    parental edge, making it quadratic and bounding pedigree depth by the
    recursion limit. It also hardcoded ``0`` as the missing-parent value,
    silently ignored a parent it could not find, and silently dropped duplicate
    IDs.

    Parameters:
        myped: A list of PyPedal pedigree objects.
        filetag: A descriptor prepended to output file names.
        io: Indicates whether to write the reordered pedigree to a file ('yes'|'no').
        debug: Boolean to toggle debugging messages on and off.
        missingparent: The value used to indicate a missing parent. Added in
            3.0.0 so this routine can honour a configured sentinel; it is the
            last parameter, so existing positional calls are unaffected.

    Returns:
        A reordered list of PyPedal pedigree objects.

    Raises:
        PyPedalPedigreeStructureError: see :func:`_order_pedigree`.
    """

    if not myped:
        if debug:
            print("DEBUG: Empty pedigree provided.")
        return myped

    start_time = time.strftime("%Y-%m-%d %H:%M:%S")

    if debug:
        print(f"DEBUG: Starting with {len(myped)} animals in the pedigree.")

    reordered = _order_pedigree(myped, missingparent=missingparent,
                                routine='pyp_utils/fast_reorder()')

    end_time = time.strftime("%Y-%m-%d %H:%M:%S")

    if debug:
        print(f"DEBUG: Reordering complete. Reordered pedigree contains {len(reordered)} animals.")
    
    # Add a timestamp to file names to prevent overwriting previous runs
    # timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    output_file = f"{filetag}_reordered.ped"
    id_map_file = f"{filetag}_id_map.map"

    # Write to output files if `io` is enabled
    if io.lower() == 'yes':
        if debug:
            print(f"Writing reordered pedigree to {output_file}")

        try:
            with open(output_file, 'w', encoding='utf-8') as aout:
                aout.write(f"# FILE: {output_file}\n")
                aout.write(f"# START TIMESTAMP: {start_time}\n")
                aout.write(f"# END TIMESTAMP: {end_time}\n")
                aout.write("# REORDERED pedigree produced by PyPedal using fast_reorder().\n")
                for animal in reordered:
                    aout.write(
                        f"{animal.animalID},{animal.sireID},{animal.damID},{animal.gen},"
                        f"{pyp_chronology.format_year_token(animal.by)}\n"
                    )
                    if debug:
                        print(f"[DEBUG]: Writing animal {animal.animalID}, birth year: {animal.by}")
            if debug:
                print(f"DEBUG: Reordered pedigree written to {output_file}")

            # Write the ID map to a file
            with open(id_map_file, 'w', encoding='utf-8') as idmap:
                idmap.write(f"# FILE: {id_map_file}\n")
                idmap.write(f"# START TIMESTAMP: {start_time}\n")
                idmap.write(f"# END TIMESTAMP: {end_time}\n")
                # idmap.write(f"# TIMESTAMP: {timestamp}\n")
                idmap.write("# Renumbered ID to Old ID mapping produced by PyPedal.\n")
                for animal in reordered:
                    idmap.write(f"{animal.originalID},{animal.animalID}\n")

        except Exception as write_err:
            print(f"ERROR: Unable to write to file {output_file}. {write_err}")
            raise

    return reordered



def renumber(
    myped,
    filetag="_renumbered_",
    io="no",
    outformat="0",
    debug=False,
    returnmap=False,
    missingparent=0,
    animaltype="new",
    cleanmap=True,
    output_dir=None
):
    """
    Renumber a pedigree for consistent ID mapping. Ensures the oldest animal
    has ID 1, the next oldest has ID 2, etc., while updating parent IDs.

    Parameters:
        myped: PyPedal pedigree object.
        filetag: Prefix for output files.
        io: Whether to write output to file ('yes' or 'no').
        outformat: Output format ('0' for minimal, '1' for detailed).
        debug: Toggle debug messages.
        returnmap: Whether to return an ID map.
        missingparent: Marker for missing parents.
        animaltype: Type of animal ('new' or other).
        cleanmap: Whether to clean up the ID map file.
        original_filename: Name of the original file (to derive the output file name).

    Returns:
        Renumbered pedigree, and optionally an ID map.
    """
    if debug:
        print(f"[DEBUG]: Pedigree of size {len(myped)} passed to renumber()")

    # Check if animals are of type "NewAnimal"
    isnewanimal = animaltype == "new"

    # Initialize ID map
    id_map = {}
    idnum = 1

    # Renumber animals and parents
    for idx, animal in enumerate(myped):
        if debug and idx == 0:
            print("[DEBUG]: Renumbering the pedigree...")
        if debug and idx % 10000 == 0:
            print(f"\t[DEBUG]: Processing animal {idx}...")

        id_map[animal.animalID] = idnum  # Map old ID to new ID
        animal.renumberedID = idnum
        animal.animalID = idnum

        # Update animal name if using "NewAnimal" and name matches original ID
        if isnewanimal and animal.name == animal.originalID:
            animal.name = animal.renumberedID

        # Renumber sire and dam.
        #
        # id_map is built incrementally by this same loop, so a parent resolves
        # if and only if it already appeared earlier in the list. That makes
        # this the point where the topological precondition is enforced. An
        # earlier KeyError handler set the parent to 0, deleting the
        # relationship with no log, no counter and no flag. A pedigree that
        # reorder could not order was therefore silently stripped of exactly the
        # edges that made it unorderable, and every coefficient downstream was
        # computed on a graph the caller never supplied.
        #
        # reorder() now either guarantees the precondition or refuses, so
        # arriving here means a PyPedal invariant was violated rather than that
        # the caller's data is bad -- hence PyPedalInternalError, not an input
        # error. The relationship is never discarded to keep going.
        for role in ('sire', 'dam'):
            parent = getattr(animal, '%sID' % role)
            if str(parent) == str(missingparent):
                continue
            try:
                setattr(animal, '%sID' % role, id_map[parent])
            except KeyError:
                raise pyp_errors.PyPedalInternalError(
                    'pyp_utils/renumber(): animal %s names %s as its %s, but '
                    'that animal has not been renumbered yet, so it does not '
                    'precede its offspring. renumber() requires a pedigree in '
                    'which parents come first; reorder the pedigree before '
                    'renumbering it. The relationship has NOT been discarded.'
                    % (animal.originalID, parent, role))

        idnum += 1

    # Renumber offspring for "NewAnimal"
    if isnewanimal:
        for animal in myped:
            animal.sons = {id_map.get(son, 0): son for son in animal.sons.keys()}
            animal.daus = {id_map.get(dau, 0): dau for dau in animal.daus.keys()}
            animal.unks = {id_map.get(unk, 0): unk for unk in animal.unks.keys()}

    if isinstance(myped, NewPedigree):
        if debug:
            print("[DEBUG] myped is a NewPedigree instance")
        myped.kw["pedigree_is_renumbered"] = True
        myped.updateidmap()

    if debug:
        for animal in myped:
            print(f"[DEBUG] AnimalID: {animal.animalID}, SireID: {animal.sireID}, DamID: {animal.damID}")


    # if output_dir is None:
    #     output_dir = os.getcwd()  # Default to current directory

    # os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
    # ped_outputfile = os.path.join(output_dir, f"{filetag}_renum.ped")
    # map_outputfile = os.path.join(output_dir, f"{filetag}_id_map.map")
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
        ped_outputfile = os.path.join(output_dir, f"{filetag}_renum.ped")
        map_outputfile = os.path.join(output_dir, f"{filetag}_id_map.map")
    else:
        ped_outputfile = f"{filetag}_renum.ped"
        map_outputfile = f"{filetag}_id_map.map"

    if debug:
        print(f"[DEBUG] Writing files to: {output_dir}")
        print(f"[DEBUG] Pedigree file: {ped_outputfile}")
        print(f"[DEBUG] Map file: {map_outputfile}")



    # Write outputs if required
    if io.lower() == 'yes':
        # Write renumbered pedigree
        try:
            with open(ped_outputfile, "w", encoding="utf-8") as pout:
                pout.write(f"# FILE: {ped_outputfile}\n# RENUMBERED pedigree produced by PyPedal.\n")
                for animal in myped:
                    if outformat == "0":
                        pout.write(f"{animal.animalID}, {animal.sireID}, {animal.damID}\n")
                    else:
                        pout.write(
                            f"{animal.animalID}, {animal.sireID}, {animal.damID}, "
                            f"{pyp_chronology.format_year_token(animal.by)}, "
                            f"{animal.sex}, {animal.fa}, {animal.gen}\n"
                        )
        except Exception as write_err:
            print(f"ERROR: Unable to write to file {ped_outputfile}. {write_err}")
            raise

    if io.lower() == 'yes':
        try:
            with open(map_outputfile, "w", encoding="utf-8") as mout:
                mout.write(f"# FILE: {map_outputfile}\n# Old ID to New ID mapping\n")
                for old_id, new_id in id_map.items():
                    mout.write(f"{old_id}, {new_id}\n")
        except Exception as write_err:
            print(f"ERROR: Unable to write to file {map_outputfile}. {write_err}")
            raise

    if cleanmap and map_outputfile and os.path.exists(map_outputfile):
        os.remove(map_outputfile)


    # Return values
    if not returnmap:
        return myped
    else:
        return myped, id_map


def load_id_map(filetag='_renumbered_'):
    """
    load_id_map() reads an ID map from the file generated by pyp_utils/renumber()
    into a dictionary. There is a VERY similar function, pyp_io/id_map_from_file(), that
    is deprecated because it is much more fragile than this procedure.
    
    :param filetag: A descriptor prepended to output file names.
    :return: A dictionary whose keys are renumbered IDs and whose values are original IDs, or an empty dictionary
             on failure.
    """
    try:
        _infile = f'{filetag}_id_map.map'
        idmap = {}

        with open(_infile, 'r') as mapin:
            for line in mapin:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip comments and empty lines
                _line = line.split(',')
                if len(_line) != 2:
                    print(f'[ERROR]: Invalid number of elements in line read from ID map file ({_line})')
                    return {}
                try:
                    renumbered_id = int(_line[1].strip())
                    original_id = int(_line[0].strip())
                    idmap[renumbered_id] = original_id
                except ValueError:
                    print(f'[ERROR]: Failed to parse IDs from line: {_line}')
                    return {}

        return idmap
    except Exception as e:
        print(f'[ERROR]: Failed to load ID map. Exception: {e}')
        return {}


def delete_id_map(filetag='_renumbered_'):
    """
    delete_id_map() checks to see if an ID map for the given filetag exists. If the file
    exists, it is deleted.
    
    :param filetag: A descriptor prepended to output file names.
    :return: True if the file was successfully deleted, False otherwise.
    """
    try:
        _infile = f'{filetag}_id_map.map'
        if os.path.exists(_infile):  # Check if file exists
            os.remove(_infile)  # Delete the file
            return True
        return False  # File does not exist
    except Exception as e:
        print(f'[ERROR]: Unable to delete file {_infile}. Exception: {e}')
        return False


def trim_pedigree_to_year(pedobj, year):
    """
    trim_pedigree_to_year() takes pedigrees and removes all individuals who were not born
    in birthyear 'year'.

    :param pedobj: A PyPedal NewPedigree object.
    :param year: A birthyear.
    :return: A pedigree containing only animals born in the given birthyear or an empty list (on failure).
    """
    try:
        modped = pedobj.pedigree[:]  # Create a copy of the pedigree
        # Filter the pedigree to include only those born in the specified year
        trimmed_pedigree = [
            animal for animal in modped
            if animal.by is not None and int(animal.by) == int(year)
        ]
        return trimmed_pedigree
    except Exception as e:
        print(f"[ERROR]: Failed to trim pedigree. Exception: {e}")
        return []


def pedigree_range(pedobj, n):
    """
    pedigree_range() takes a renumbered pedigree and removes all individuals with a renumbered ID > n.
    The reduced pedigree is returned. Assumes that the input pedigree is sorted on animal key in ascending order.

    :param pedobj: A PyPedal pedigree object.
    :param n: A renumbered animalID.
    :return: A pedigree containing only animals with renumbered IDs <= n or an empty list (on failure).
    """
    try:
        # Use slicing to extract individuals with renumbered IDs <= n
        modped = pedobj.pedigree[:n]
        return modped
    except Exception as e:
        print(f"[ERROR]: Failed to filter pedigree. Exception: {e}")
        return []


def sort_dict_by_keys(mydict):
    """
    sort_dict_by_keys() returns a list of values in the dictionary in the order obtained by sorting the keys.
    Taken from the routine sortedDictValues3 in the "Python Cookbook", p. 39.

    :param mydict: A non-empty Python dictionary.
    :return: A list of values sorted by ascending order of their keys or an empty list (on failure).
    """
    try:
        if not mydict:  # Check if the dictionary is empty
            return []
        else:
            sorted_keys = sorted(mydict.keys())  # Get sorted keys
            return [mydict[key] for key in sorted_keys]  # Return values based on sorted keys
    except Exception as e:
        print(f"[ERROR]: Failed to sort dictionary by keys. Exception: {e}")
        return []


def sort_dict_by_values(dictionary):
    """
    sort_dict_by_values() returns a list of tuples where the keys in the dictionary
    are sorted ascending by value, first on value and then on key within value.

    :param dictionary: A Python dictionary to be sorted by values.
    :return: A list of tuples sorted in ascending order of value, then key.
    """
    try:
        # Sort by values first, then by keys if values are the same
        sorted_items = sorted(dictionary.items(), key=lambda item: (item[1], item[0]))
        return sorted_items
    except Exception as e:
        print(f"[ERROR]: Failed to sort dictionary by values. Exception: {e}")
        return []


def simple_histogram_dictionary(mydict, histchar='*', histstep=5):
    """
    simple_histogram_dictionary() returns a dictionary containing a simple, text
    histogram. The input dictionary is assumed to contain keys which are distinct levels
    and values that are counts.
    
    :param mydict: A dictionary where keys are levels and values are counts.
    :param histchar: The character used to draw the histogram (default is '*').
    :param histstep: Used to determine the number of bins (stars) in the diagram.
    :return: A dictionary containing the histogram by level or an empty dictionary (on failure).
    """
    try:
        if not mydict:
            return {}

        hist_dict = {}
        hist_sum = sum(mydict.values())
        
        if histstep < 1 or histstep > 100:
            histstep = 5

        for key, count in mydict.items():
            # Calculate frequency percentage
            freq_percentage = (float(count) / float(hist_sum)) * 100.0
            # Calculate the number of stars for the histogram
            num_stars = math.ceil(freq_percentage / float(histstep))
            if num_stars > 0:
                hist_dict[key] = f"{histchar * int(num_stars):<20}"
            else:
                hist_dict[key] = " " * 20

        return hist_dict
    except Exception as e:
        print(f"[ERROR]: Exception occurred in simple_histogram_dictionary: {e}")
        return {}


def reverse_string(mystring):
    """
    reverse_string() reverses the input string and returns the reversed version.
    
    :param mystring: A non-empty Python string.
    :return: The input string with the order of its characters reversed, or False (on failure).
    """
    try:
        if not isinstance(mystring, str):
            raise ValueError("Input must be a string.")
        
        return mystring[::-1]
    except Exception as e:
        print(f"[ERROR]: Exception occurred in reverse_string: {e}")
        return False


def pyp_nice_time():
    """
    pyp_nice_time() returns the current date and time formatted as, e.g., Wed Mar 30 10:26:31 2005.
    
    :return: A string containing the formatted date and time, or False (on failure).
    """
    try:
        return time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())
    except Exception as e:
        print(f"[ERROR]: Exception occurred in pyp_nice_time: {e}")
        return False


def string_to_table_name(instring):
    """
    string_to_table_name() takes an arbitrary string and returns a string that
    is safe to use as an SQLite table name.
    
    :param instring: A string that will be converted to an SQLite-safe table name.
    :return: A string that is safe to use as an SQLite table name.
    """
    try:
        # Allow only alphanumeric characters and underscores
        allowed_chars = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')
        outstring = ''.join(c for c in instring if c in allowed_chars)
        return outstring
    except Exception as e:
        print(f"[ERROR]: Exception occurred in string_to_table_name: {e}")
        return instring


def pyp_datestamp():
    """
    pyp_datestamp() returns a datestamp, as a string, of the format YYYYMMDDHHMMSS.
    
    :return: A 14-character string containing the datestamp, or False (on failure).
    """
    try:
        # Format current local time as YYYYMMDDHHMMSS
        return time.strftime('%Y%m%d%H%M%S', time.localtime())
    except Exception as e:
        print(f"[ERROR]: Exception occurred in pyp_datestamp: {e}")
        return False


def subpedigree(pedobj, anlist):
    """
    subpedigree() takes a NewPedigree object and a list of animal IDs and returns 
    a NewPedigree object containing only the animals in the animal list.

    :param pedobj: PyPedal pedigree object.
    :param anlist: A list of animal IDs.
    :return: An instance of a NewPedigree object, or False (on failure).
    """
    try:
        # Create a deep copy of the pedigree object
        NewPed = copy.deepcopy(pedobj)
        order = []
        _tempped = copy.deepcopy(NewPed.pedigree)

        # Create an initial order of animal IDs in the pedigree
        for _p in _tempped:
            order.append(_p.animalID)

        # Remove animals not in the provided list
        for _p in list(NewPed.pedigree):  # Iterate over a copy of the list to allow modification
            if _p.animalID not in anlist:
                _anidx = order.index(_p.animalID)

                # Delete the animal from the maps
                del NewPed.namemap[NewPed.namebackmap[NewPed.backmap[_p.animalID]]]
                del NewPed.namebackmap[NewPed.backmap[_p.animalID]]
                del NewPed.idmap[NewPed.backmap[_p.animalID]]
                del NewPed.backmap[_p.animalID]
                del _tempped[_anidx]
                del order[_anidx]

        # Update the metadata object
        NewPed.pedigree = _tempped
        NewPed.metadata = pyp_newclasses.PedigreeMetadata(NewPed.pedigree, NewPed.kw, NewPed.snp)

        # Renumber the pedigree if required
        if NewPed.kw.get('renumber', False):
            NewPed.renumber()

        return NewPed
    except Exception as e:
        print(f"[ERROR]: Exception occurred in subpedigree: {e}")
        return False


def founders_from_list(anlist, unkid):
    """
    founders_from_list() takes a list of NewAnimal objects and returns a
    list of animalIDs that represent founders in that pedigree (animals
    with an unknown sire and dam).
    
    :param anlist: A list of NewAnimal objects.
    :param unkid: The animalID assigned to unknown parents.
    :return: A list of animalIDs of founders, or an empty list (on failure).
    """
    try:
        # Filter animals with unknown sire and dam IDs
        flist = [animal.animalID for animal in anlist if animal.sireID == unkid and animal.damID == unkid]
        return flist
    except Exception as e:
        print(f"[ERROR]: Exception occurred in founders_from_list: {e}")
        return []


def founder_allele_dict(pedobj):
    """
    founder_allele_dict() takes a pedigree and returns a dictionary containing
    an entry for each unique founder allele.
    
    :param pedobj: PyPedal pedigree.
    :return: A dictionary whose keys are founder alleles and whose values are 0.0, or an empty dictionary (on failure).
    """
    try:
        # Extract founder alleles from the unique founder list
        falist = [pedobj.pedigree[x - 1].alleles for x in pedobj.metadata.unique_founder_list]
        fadict = {}
        
        # Populate the dictionary with founder alleles
        for fa in falist:
            fadict[fa[0]] = 0.0
            fadict[fa[1]] = 0.0
        
        return fadict
    except Exception as e:
        print(f"[ERROR]: Exception occurred in founder_allele_dict: {e}")
        return {}


def founder_allele_map(pedobj):
    """
    founder_allele_map() takes a pedigree and returns a dictionary
    mapping founder alleles to their animal of origin.
    
    :param pedobj: PyPedal pedigree.
    :return: A dictionary whose keys are founder alleles and whose values are the founder of origin, 
             or an empty dictionary (on failure).
    """
    try:
        # Extract founder alleles from the unique founder list
        falist = [pedobj.pedigree[x - 1].alleles for x in pedobj.metadata.unique_founder_list]
        fadict = {}

        # Initialize the dictionary with founder alleles
        for fa in falist:
            fadict[fa[0]] = None
            fadict[fa[1]] = None

        # Assign founder indices to each allele
        for fa in fadict.keys():
            fadict[fa] = list(fadict.keys()).index(fa)

        return fadict
    except Exception as e:
        print(f"[ERROR]: Exception occurred in founder_allele_map: {e}")
        return {}


def founder_allele_freq(pedobj, anlist, allele_map, allele_mat, column):
    """
    founder_allele_freq() takes a pedigree and returns a dictionary containing
    an entry for each unique founder allele and the frequency of that allele in
    a list of provided animals.
    
    :param pedobj: PyPedal pedigree.
    :param anlist: A list of animals in which to track founder allele frequencies.
    :param allele_map: Dictionary mapping founder alleles to their animal of origin.
    :param allele_mat: Dictionary counting founder alleles.
    :param column: Select the correct column to increment.
    :return: A dictionary whose keys are founder alleles and whose values are frequencies, or an empty
             dictionary (on failure).
    """
    try:
        # Get the mapped indices of the animals in `anlist`
        backmap = [pedobj.idmap[an] for an in anlist]
        
        # Extract alleles for each animal in the backmap
        fafreq = [pedobj.pedigree[an - 1].alleles for an in backmap]
        
        # Update allele frequencies in the allele matrix
        for fa in fafreq:
            row0 = allele_map.get(fa[0], None)
            row1 = allele_map.get(fa[1], None)
            if row0 is not None:
                allele_mat[row0, column] = allele_mat.get((row0, column), 0) + 1.0
            if row1 is not None:
                allele_mat[row1, column] = allele_mat.get((row1, column), 0) + 1.0
        
        return allele_mat
    except Exception as e:
        print(f"[ERROR]: Exception occurred in founder_allele_freq: {e}")
        return {}


def list_intersection(pedobjs):
    """
    list_intersection() returns a PyPedal pedigree object which contains the animals that are common to
    both input pedigrees. If there are no animals in common between the two pedigrees, a value
    of False is returned.
    
    :param pedobjs: A list of PyPedal pedigrees.
    :return: A new PyPedal pedigree containing the animals that are common to both input pedigrees, or False.
    """
    if not pedobjs:
        logging.error(
            "An empty list was passed to pyp_utils/list_intersection(). Cannot compute intersection on 0 pedigrees."
        )
        return False

    try:
        if pedobjs[0].kw.get('debug_messages', False):
            print(f"[INFO]: Computing intersection of {len(pedobjs)} pedigrees")
        logging.info("Computing intersection of %s pedigrees", len(pedobjs))

        plen = len(pedobjs)

        # If there's only one pedigree, return it directly
        if plen == 1:
            return pedobjs[0]

        # Initialize with an empty pedigree for intersection computation
        intersected = pyp_newclasses.loadPedigree({'pedfile': 'null'}, pedsource='null')

        # Compute intersection using a nested approach
        for p in range(plen - 1, -1, -1):
            if pedobjs[p].__class__.__name__ != 'NewPedigree':
                if pedobjs[0].kw.get('debug_messages', False):
                    print(
                        f"[ERROR]: Pedigree {p} in pyp_utils.list_intersection() is not a NewPedigree instance! Skipping."
                    )
                logging.error(
                    "Pedigree %s in pyp_utils.list_intersection() is not a NewPedigree instance! Skipping.", p
                )
            else:
                intersected = intersected.intersection(pedobjs[p])

        return intersected

    except Exception as e:
        logging.error("An error occurred in pyp_utils/list_intersection(): %s", str(e))
        return False


def list_union(pedobjs):
    """
    list_union() returns a PyPedal pedigree object which contains all animals included in either or both of the
    input pedigrees.
    :param pedobjs: A list of PyPedal pedigrees.
    :return: A new PyPedal pedigree containing the animals that are included in either or both input pedigrees, or False
             on failure.
    """
    if not pedobjs:
        logging.error("Cannot complete union operation because no pedigrees were provided.")
        return NotImplemented

    try:
        logging.info("Computing union of %s pedigrees", len(pedobjs))
        # Initialize an empty (null) pedigree to build the union
        new_pedigree = pyp_newclasses.loadPedigree({'pedfile': 'null'}, pedsource='null')

        for idx, ped in enumerate(pedobjs):
            if ped.__class__.__name__ != "NewPedigree":
                logging.warning(
                    "Pedigree %s in pyp_utils.list_union() is not a NewPedigree instance! Skipping.", idx
                )
                continue

            logging.info("Using match rule %s to compare pedigrees", ped.kw.get("match_rule", "default"))

            # Renumber the pedigrees if necessary
            if idx > 0 and not new_pedigree.kw.get("pedigree_is_renumbered", False):
                new_pedigree.renumber()
                logging.info("Renumbering pedigree %s", new_pedigree.kw.get("pedname"))

            if ped.kw.get("pedigree_is_renumbered") != 1:
                ped.renumber()
                logging.info("Renumbering pedigree %s", ped.kw.get("pedname"))

            # Compute the union using the NewPedigree::__add__() method
            try:
                new_pedigree += ped
            except Exception as e:
                logging.error("Could not compute union of input pedigrees: %s", str(e))
                return False

        # Return the union of all pedigrees
        return new_pedigree

    except Exception as e:
        logging.error("An error occurred during the union operation: %s", str(e))
        return False


def guess_pedformat(animal, ped_kw):
    """
    guess_pedformat() tries to guess the pedigree format code that best matches a NewAnimal instance provided
    as input based on the animal attributes as compared to the default missing values for those attributes.
    This function is intended primarily for use with functions that take as input NewAnimal objects which may
    come from different pedigrees with different pedformats. Note that the logic of guess_pedformat() depends
    entirely on the single NewAnimal passed as input.

    :param animal: A NewAnimal instance.
    :param ped_kw: The kw dictionary from a NewPedigree object.
    :return: A string containing the best-guess pedformat, or False.
    """
    if isinstance(animal, pyp_newclasses.NewAnimal):
        try:
            # Determine the basic format based on `originalID` type.
            if isinstance(animal.originalID, int):
                pedformat = "asd"
            elif isinstance(animal.originalID, str):
                pedformat = "ASD"
            else:
                logging.error(
                    "pyp_utils/guess_pedformat() cannot process animal.originalID if it is not a string or an integer type!"
                )
                if ped_kw.get("messages") == "verbose":
                    print("[ERROR]: pyp_utils/guess_pedformat() cannot process animal.originalID if it is not a string or an integer type!")
                return False

            # Add additional format codes based on other attributes.
            if animal.bd is not None:
                pedformat += "b"
            if animal.by is not None:
                pedformat += "y"
            if animal.name != ped_kw["missing_name"]:
                pedformat += "n"
            if animal.breed != ped_kw["missing_breed"]:
                pedformat += "r"
            if animal.herd != ped_kw["missing_herd"]:
                pedformat += "H"
            if animal.sex != ped_kw["missing_sex"]:
                pedformat += "x"
            if animal.fa != ped_kw["missing_inbreeding"]:
                pedformat += "f"
            if animal.alive != ped_kw["missing_alive"]:
                pedformat += "l"
            if animal.age != ped_kw["missing_age"]:
                pedformat += "e"
            if animal.gen != ped_kw["missing_gen"]:
                pedformat += "g"
            if animal.gencoeff != ped_kw["missing_gencoeff"]:
                pedformat += "p"
            if animal.alleles != ped_kw["missing_alleles"]:
                pedformat += "L"
            if animal.userField != ped_kw["missing_userfield"]:
                pedformat += "u"

            # Logging the guessed pedformat.
            logging.info("The best-guess pedformat is: '%s'", pedformat)
            if ped_kw.get("messages") == "verbose":
                print(f"[INFO]: The best-guess pedformat is: '{pedformat}'")
            return pedformat
        except KeyError as e:
            logging.error("Missing key in ped_kw: %s", e)
            if ped_kw.get("messages") == "verbose":
                print(f"[ERROR]: Missing key in ped_kw: {e}")
            return False
    else:
        logging.error("You passed an item to pyp_utils/guess_pedformat() that is not a NewAnimal!")
        if ped_kw.get("messages") == "verbose":
            print("[ERROR]: You passed an item to pyp_utils/guess_pedformat() that is not a NewAnimal!")
        return False


def list_duplicates(pedobj, keep_rule=''):
    """
    list_duplicates() identifies animals that are duplicates (have the same animal ID in the input pedigree file) and
    decides which record to retain as the "correct" animal.
    :param pedobj: A PyPedal pedigree.
    :param keep_rule: The rule to use to decide which duplicate record to keep.
    :return: A list of duplicate animals to delete.
    """
    _temp_ped = copy.deepcopy(pedobj.pedigree)
    duplicate_id_count = {}
    duplicates = []
    if pedobj.kw.get('messages') == 'verbose':
        print(f"\t[INFO]: pyp_utils/resolve_duplicates(): There are {len(pedobj.pedigree)} animals in the pedigree.")
    logging.info("pyp_utils/resolve_duplicates(): There are %s animals in the pedigree.", len(pedobj.pedigree))

    _p_idx = 0
    for _p in pedobj.pedigree:
        _p_printed = False
        for _t_idx in range(_p_idx + 1, len(pedobj.pedigree)):
            _t = pedobj.pedigree[_t_idx]
            if _p.originalID == _t.originalID:
                if _p.originalID not in duplicate_id_count:
                    duplicate_id_count[_p.originalID] = 1
                else:
                    duplicate_id_count[_p.originalID] += 1
                if not _p_printed:
                    # Sire
                    if _p.sireID == pedobj.kw['missing_parent']:
                        _orig_sire = pedobj.kw['missing_parent']
                    else:
                        _orig_sire = pedobj.pedigree[_p.sireID - 1].originalID
                    # Dam
                    if _p.damID == pedobj.kw['missing_parent']:
                        _orig_dam = pedobj.kw['missing_parent']
                    else:
                        _orig_dam = pedobj.pedigree[_p.damID - 1].originalID
                    duplicates.append(f"Renumbered ID:\t{_p.animalID}\toriginalID:\t{_p.originalID}\tsire ID:\t{_orig_sire}\tdam ID:\t{_orig_dam}")
                    _p_printed = True
                # Sire
                if _t.sireID == pedobj.kw['missing_parent']:
                    _orig_sire = pedobj.kw['missing_parent']
                else:
                    _orig_sire = pedobj.pedigree[_t.sireID - 1].originalID
                # Dam
                if _t.damID == pedobj.kw['missing_parent']:
                    _orig_dam = pedobj.kw['missing_parent']
                else:
                    _orig_dam = pedobj.pedigree[_t.damID - 1].originalID
                duplicates.append(f"Renumbered ID:\t{_t.animalID}\toriginalID:\t{_t.originalID}\tsire ID:\t{_orig_sire}\tdam ID:\t{_orig_dam}")
        _p_idx += 1

    for k, v in duplicate_id_count.items():
        if pedobj.kw.get('messages') == 'verbose':
            print(f"\t[INFO]: pyp_utils/resolve_duplicates(): originalID {k} occurs {v} times in the pedigree file!")
        logging.info("pyp_utils/resolve_duplicates(): originalID %s occurs %s times in the pedigree file!", k, v)

    return duplicates


def list_likely_same_animals(pedobj, unique_external_field=None):
    """
    Advisory groups of pedigree records that may represent one real animal.

    This does not mutate, merge, or delete. It does not overload
    ``list_duplicates``, which is about duplicate ``originalID`` values.

    A **strong** group requires a matching declared unique external
    identifier (registration, chip/EID, or ``userField`` only when that
    role is declared). A **heuristic** group requires the same call name,
    the same parents, and the same recorded birth date or year, without
    that strong identifier. Same name alone is not enough.

    Parameters
    ----------
    pedobj : NewPedigree
        Pedigree to inspect.
    unique_external_field : str or None
        Attribute that is a unique external identifier. ``'userField'``
        (also ``'u'`` / ``'userfield'``) is accepted. When omitted, the
        pedigree option ``kw['unique_external_field']`` is used.

    Returns
    -------
    list of dict
        Each item has ``animals`` (originalIDs), ``strength``
        (``'strong'`` or ``'heuristic'``), and ``evidence`` (field names).
    """
    if unique_external_field is None:
        unique_external_field = pedobj.kw.get('unique_external_field')

    attr = None
    if unique_external_field in ('userField', 'userfield', 'u'):
        attr = 'userField'

    missing_user = pedobj.kw.get('missing_userfield', 'Unknown')
    missing_name = pedobj.kw.get('missing_name')
    missing_parent = pedobj.kw['missing_parent']

    def _is_missing(value, sentinel):
        return value == sentinel or str(value) == str(sentinel)

    def _original_parent(animal, slot):
        pid = getattr(animal, slot)
        if pid == missing_parent or str(pid) == str(missing_parent):
            return None
        return pedobj.backmap.get(pid)

    def _recorded_chrono(animal):
        if animal.bd is not None:
            return ('bd', str(animal.bd))
        if animal.by is not None:
            return ('by', animal.by)
        return None

    groups = []
    strong_ids = set()

    if attr is not None:
        buckets = {}
        for animal in pedobj.pedigree:
            value = getattr(animal, attr)
            if _is_missing(value, missing_user) or value in ('', None):
                continue
            buckets.setdefault(value, []).append(animal.originalID)
        for value, oids in buckets.items():
            if len(oids) < 2:
                continue
            groups.append({
                'animals': list(oids),
                'strength': 'strong',
                'evidence': [attr],
            })
            strong_ids.update(oids)

    heuristic = {}
    for animal in pedobj.pedigree:
        if _is_missing(animal.name, missing_name) or animal.name in ('', None):
            continue
        chrono = _recorded_chrono(animal)
        if chrono is None:
            continue
        key = (
            str(animal.name),
            _original_parent(animal, 'sireID'),
            _original_parent(animal, 'damID'),
            chrono,
        )
        heuristic.setdefault(key, []).append(animal.originalID)

    for key, oids in heuristic.items():
        if len(oids) < 2:
            continue
        if strong_ids and set(oids) <= strong_ids:
            continue
        groups.append({
            'animals': list(oids),
            'strength': 'heuristic',
            'evidence': ['name', 'parents', key[3][0]],
        })

    if groups and pedobj.kw.get('messages') == 'verbose':
        print(
            "[INFO]: list_likely_same_animals() found %s candidate group(s)."
            % len(groups)
        )
        logging.info(
            "list_likely_same_animals() found %s candidate group(s).",
            len(groups),
        )
    return groups


def which(program):
    """
    which() tries to determine if an executable program exists in the user's path. The code was taken from Stack
    Overflow (http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python).
    :param program: The name of the program to find.
    :return: The name of the program, or False (on failure).
    """
    def is_exe(fpath):
        """
        Checks if the given file path is an executable file.
        :param fpath: File path to check.
        :return: True if the file is executable, False otherwise.
        """
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        # If a direct path to the program is provided
        if is_exe(program):
            return program
    else:
        # Search through system PATH
        for path in os.environ.get("PATH", "").split(os.pathsep):
            path = path.strip('"')  # Remove any extra quotes
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return False


def remove_missing(pedobj):
    """
    remove_missing() takes a NewPedigree object and removes any animals whose ID
    is equal to the missing animal identifier.
    :param pedobj: PyPedal pedigree object.
    :return: An instance of a NewPedigree object with no animals with missing-valued IDs, or False (on failure).
    """
    try:
        n_removed = 0
        for _p in pedobj.pedigree[:]:  # Create a shallow copy to safely iterate while modifying
            if _p.animalID == pedobj.kw['missing_parent']:
                # Delete the animal from the pedigree and the ID maps
                del pedobj.namemap[pedobj.namebackmap[pedobj.backmap[_p.animalID]]]
                del pedobj.namebackmap[pedobj.backmap[_p.animalID]]
                del pedobj.idmap[pedobj.backmap[_p.animalID]]
                del pedobj.backmap[_p.animalID]
                pedobj.pedigree.remove(_p)  # Remove the animal from the pedigree
                n_removed += 1

        # Update pedigree metadata since animal counts, offspring counts, etc., may have changed
        pedobj.metadata = pyp_newclasses.PedigreeMetadata(pedobj.pedigree, pedobj.kw)
        
        # Renumber the pedigree if specified
        if pedobj.kw.get('renumber', False):
            pedobj.renumber()
        
        # Log verbose messages if enabled
        if pedobj.kw.get('messages') == 'verbose':
            print(f"[INFO]: pyp_utils/remove_missing() removed {n_removed} animals with missing IDs from the pedigree!")
        
        return pedobj

    except Exception as e:
        logging.error(f"pyp_utils/remove_missing() encountered an error: {e}")
        return False
