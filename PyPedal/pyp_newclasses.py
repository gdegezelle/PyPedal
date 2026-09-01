#!/usr/bin/env python3

"""
pyp_newclasses: Class structure for PyPedal

This module includes the class structure for PyPedal, forming the backbone of PyPedal.
It includes a master class to which most computational routines will be bound as methods, a `NewAnimal` class, 
and a `PedigreeMetadata` class.


**Version**: see ``PyPedal.__version__``.

**Author**: John B. Cole (john.b.cole@gmail.com)

**License**: LGPL

Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle, 2025-2026. See CHANGELOG.md. SPDX-License-Identifier: LGPL-2.1-or-later.
"""

from __future__ import annotations

import hashlib
import logging

import os
import sys
import warnings
from typing import Any, Dict, List, Literal, Optional, Tuple

from . import _pyp_parse
from networkx import DiGraph
import numpy as np
import math
from configparser import ConfigParser, Error

from . import (
    pyp_chronology,
    pyp_db,
    pyp_errors,
    pyp_io,
    pyp_metrics,
    pyp_nrm,
    pyp_snp,
    pyp_utils,
)

logger = logging.getLogger(__name__)

PYPEDAL_LOGGER_NAME = "PyPedal"
_PYPEDAL_OWNED_HANDLER = "_pypedal_owned"


def filetag_from_pedfile(pedfile) -> str:
    """Return the output prefix for a pedigree path.

    Only the final extension is stripped, so the directory is kept and a
    leading ``./`` does not collapse to ``untitled_pedigree``.
    """
    tag = os.path.splitext(str(pedfile))[0]
    if len(tag) == 0:
        return "untitled_pedigree"
    return tag


def install_pedigree_logfile(logfile: str) -> None:
    """Attach a PyPedal-owned FileHandler for the current pedigree logfile.

    Any previous PyPedal-owned FileHandler is closed and removed first.
    Handlers owned by the host application are not touched. The package
    logger is used, never the root logger.
    """
    package = logging.getLogger(PYPEDAL_LOGGER_NAME)
    package.setLevel(logging.DEBUG)
    for handler in list(package.handlers):
        if getattr(handler, _PYPEDAL_OWNED_HANDLER, False):
            package.removeHandler(handler)
            handler.close()
    handler = logging.FileHandler(logfile, mode="w")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(levelname)-8s %(message)s",
            datefmt="%a, %d %b %Y %H:%M:%S",
        )
    )
    setattr(handler, _PYPEDAL_OWNED_HANDLER, True)
    package.addHandler(handler)


PedigreeSource = Literal[
    "file",
    "db",
    "graph",
    "graphfile",
    "null",
    "animallist",
    "gedcomfile",
    "textstream",
]


class NewPedigree:
    """
    The NewPedigree class is the main data structure for PyPedal,
    handling pedigree data structures, simulation, and metadata.
    """

    def __init__(self, kw: Optional[Dict[str, Any]] = None, kwfile: str = 'pypedal.ini') -> None:
        """
        Initializes a NewPedigree object.

        :param kw: A dictionary of options. If not provided, it's loaded from kwfile.
        :param kwfile: An optional configuration file name.
        """

        # Load the configuration file if `kw` is not provided
        if kw is None:
            kw = self._load_config_file(kwfile)
            if kw.get('debug_messages'):
                print(f"[DEBUG]: Loaded configuration: {kw}")

        # Set defaults for pedigree simulation
        self._set_simulation_defaults(kw)

        if "pedigree_save" not in kw:
            kw["pedigree_save"] = False

        # Set general configuration defaults
        self._set_general_defaults(kw)

        # Default missing values for NewAnimal objects
        self._set_animal_defaults(kw)

        kw.setdefault("file_io", True)
        kw.setdefault("debug_messages", False)
        kw.setdefault("form_nrm", False)
        # Optional declared unique external identifier for same-animal
        # candidate detection (e.g. 'userField'). Generic userField is
        # not assumed to be a registration number.
        kw.setdefault("unique_external_field", None)
        kw.setdefault("nrm_method", "nrm")
        kw.setdefault("nrm_format", "text")
        kw.setdefault("f_computed", False)

        # Computed genomic inbreeding
        kw.setdefault("g_computed", False)
        kw.setdefault("log_ped_lines", 0)
        kw.setdefault("log_long_filenames", False)
        kw.setdefault("pedigree_summary", 1)
        if kw["pedigree_summary"] not in [0, 1, 2]:
            kw["pedigree_summary"] = 1
        kw.setdefault("animal_type", "new")
        if kw["animal_type"] not in ["new", "light"]:
            kw["animal_type"] = "new"

        if "f" in kw["pedformat"]:
            kw["f_computed"] = True
        if "G" in kw["pedformat"]:
            kw["g_computed"] = True

        # Pedigree file handling
        if kw["simulate_pedigree"]:
            kw["pedfile"] = "simulated_pedigree"
        kw["filetag"] = filetag_from_pedfile(kw["pedfile"])

        # Handle database settings
        kw.setdefault("database_name", "pypedal")
        kw.setdefault("database_file", f"{kw['database_name']}.db")
        kw.setdefault("database_table", pyp_utils.string_to_table_name(kw.get("filetag", "untitled_pedigree")))
        kw.setdefault("database_debug", False)
        kw.setdefault("database_type", "sqlite")
        kw.setdefault("database_host", "localhost")
        kw.setdefault("database_user", "anonymous")
        kw.setdefault("database_passwd", "anonymous")
        kw.setdefault("database_port", "")
        kw.setdefault("database_sql", "SELECT * FROM %s")
        kw.setdefault("database_compatibility", "sqlite")

        # This keyword is used by pyp_nrm/fast_a_matrix() to determine if the diagonals
        # of the relationship matrix should be augmented by founder coefficients of
        # inbreeding or not. This is disabled by default.
        kw.setdefault("foundercoi", False)

        # If the user provides a paper size make sure that it is a supported value.
        # Right now, only a4 and letter are supported.  Note that 'a4' is silently
        # changed to 'A4' because ReportLab is case-sensitive.  Also handle setting
        # the default unit of measurement.
        kw.setdefault("paper_size", "letter")
        if kw["paper_size"] == "a4":
            kw["paper_size"] = "A4"
        if kw["paper_size"] not in ["A4", "letter"]:
            kw["paper_size"] = "letter"
        kw.setdefault("default_unit", "inch")
        if kw["default_unit"] not in ["cm", "inch"]:
            kw["default_unit"] = "inch"
        kw.setdefault("default_fontsize", 10)
        try:
            kw["default_fontsize"] = int(kw["default_fontsize"])
        except ValueError:
            kw["default_fontsize"] = 10
        kw.setdefault("default_report", kw["filetag"])

        # This is a hack to fix inbreeding routines to use sparse matrices if necessary
        kw.setdefault("matrix_type", "sparse")
        if kw["matrix_type"] not in ["dense", "sparse"]:
            kw["matrix_type"] = "sparse"

        # SNP genotypes require a separate filename in addition to the 'P' code in
        # the pedformat.
        kw.setdefault("snpfile", False)
        kw.setdefault("snp_sepchar", " ")

        # Internal use only
        kw.setdefault("newanimal_caller", "loader")

        # Options related to traits
        kw.setdefault("trait_names", [])
        kw.setdefault("trait_auto_name", True)
        kw.setdefault("trait_count", 0)

        # Set matching rule to use with the __add__() method
        if "match_rule" not in kw:
            kw["match_rule"] = "ASD" if "A" in kw["pedformat"] and "S" in kw["pedformat"] and "D" in kw["pedformat"] else "asd"

        # After all arguments in the options dictionary are processed, attach it to this object
        self.kw = kw

        # Initialize the Big Main Data Structures to null values
        self.pedigree: list[NewAnimal | LightAnimal] = []
        # Animals PyPedal created itself because they appeared as a sire or dam
        # but had no record of their own. The post-load count is then larger
        # than the number of input rows, which used to be discoverable only by
        # reading thousands of [NOTE] lines. Surfaced in
        # PedigreeMetadata.num_implicit_parents.
        self._implicit_parents: list[int | str] = []
        self.metadata = {}
        self.idmap: dict[int | str, int | str] = {}
        self.backmap: dict[int | str, int | str] = {}
        self.namemap: dict[str, int | str] = {}
        self.namebackmap: dict[int | str, str] = {}
        self.stringmap: dict[Any, Any] = {}
        self.snp = False  # Holds SNP genotype data, indexed by originalID

        # Maybe these will go in a configuration file later
        self.starline = "*" * 80

        # This is the list of valid pedformat codes
        self.pedformat_codes = ["a", "s", "d", "g", "x", "b", "f", "r", "n",
                                "y", "l", "e", "p", "A", "S", "D", "L", "Z",
                                "h", "H", "u", "T", "P", "G", "Y"
                                ]

        # This dictionary maps pedformat codes to NewAnimal attributes
        self.new_animal_attr = {"a": "animalID", "s": "sireID", "d": "damID",
                                "g": "gen", "x": "sex", "b": "bd", "f": "fa",
                                "r": "breed", "n": "name", "y": "by",
                                "l": "alive", "e": "age", "p": "gencoeff",
                                "A": "name", "S": "sireName", "D": "damName",
                                "L": "alleles",  "Z": False, "h": "herd",
                                "H": "originalHerd", "u": "userField",
                                "T": "traits", "P": "SNPgenotype",
                                "G": "genomicInbreeding",
                                "Y": "genomicHomozygosity"
                                }

        # Start logging!
        if "logfile" not in self.kw:
            if self.kw["log_long_filenames"]:
                self.kw["logfile"] = f"{self.kw['filetag']}_{pyp_utils.pyp_datestamp()}.log"
            else:
                self.kw["logfile"] = f"{self.kw['filetag']}.log"
        install_pedigree_logfile(self.kw["logfile"])
        logger.info("Logfile %s instantiated.", self.kw["logfile"])
        if self.kw["messages"] == "verbose" and self.kw["pedigree_summary"]:
            print(f"[INFO]: Logfile {self.kw['logfile']} instantiated.")

        # Deal with aberrant cases of log_ped_lines here.
        _lpl = self.kw["log_ped_lines"]
        try:
            self.kw["log_ped_lines"] = int(self.kw["log_ped_lines"])
        except ValueError:
            self.kw["log_ped_lines"] = 0
            logger.warning(
                "An incorrect value (%s) was provided for the option log_ped_lines, "
                "which must be a number greater than or equal to 0. It has been set to 0.",
                _lpl,
            )
        if self.kw["log_ped_lines"] < 0:
            self.kw["log_ped_lines"] = 0
            logger.warning(
                "A negative value (%s) was provided for the option log_ped_lines, "
                "which must be greater than or equal to 0. It has been set to 0.",
                self.kw["log_ped_lines"],
            )

        if self.kw["messages"] == "verbose" and self.kw["debug_messages"]:
            logger.info("Printing program options for debugging")
            for _k, _v in self.kw.items():
                logger.debug("\t%s\t%s", _k, _v)


    def __len__(self):
        return len(self.pedigree)
    
    
    def __iter__(self):
        """
        Makes the NewPedigree object iterable by returning an iterator
        over the animals in the pedigree.
        """
        return iter(self.pedigree)


    def _load_config_file(self, kwfile: str) -> Dict[str, Any]:
        """
        Load an option file into a flat options dictionary.

        Parameters
        ----------
        kwfile : str
            Path to the INI configuration file.

        Returns
        -------
        dict
            A flat, type-coerced options dict, ready to be used as `kw`.
            Section headers, if present, are flattened away: PyPedal option
            names are global, and every consumer reads kw['renumber'], never
            kw['analysis']['renumber'].

        Raises
        ------
        PyPedalOptionError
            If the configuration file fails to load or is invalid.
        """
        try:
            # Delegated so that there is one definition of what a PyPedal
            # option file means. Returns a FLAT, type-coerced dict -- see
            # pyp_io.read_ini_file() for why all three of those words matter.
            return pyp_io.read_ini_file(kwfile)

        except (Error, ValueError, OSError) as e:
            # Was `raise SystemExit(...)`, which killed the host process over a
            # bad .ini file and could not be caught by `except Exception`.
            raise pyp_errors.PyPedalOptionError(
                f"Failed to load configuration file '{kwfile}': {e}") from e


    def _set_simulation_defaults(self, kw: Dict[str, Any]) -> None:
        """Set the default values for simulation-related keywords."""
        kw.setdefault('simulate_pedigree', False)
        if kw['simulate_pedigree']:
            kw.setdefault('simulate_n', 15)
            kw.setdefault('simulate_g', 3)
            kw.setdefault('simulate_ns', 4)
            kw.setdefault('simulate_nd', 4)
            kw.setdefault('simulate_mp', False)
            kw.setdefault('simulate_po', False)
            kw.setdefault('simulate_fs', False)
            kw.setdefault('simulate_sr', 0.5)
            kw.setdefault('simulate_ir', 0.0)
            kw.setdefault('simulate_pmd', 100)
            kw.setdefault('simulate_save', False)
        else:
            if 'pedfile' not in kw:
                raise ValueError("[ERROR]: 'pedfile' must be provided if simulation is not enabled.")


    def _set_general_defaults(self, kw: Dict[str, Any]) -> None:
        """Set default values for general configuration options."""
        kw.setdefault('pedformat', 'asd')
        kw.setdefault('has_header', False)
        kw.setdefault('pedname', 'Untitled')
        kw.setdefault('animal_type', 'new')
        kw.setdefault('messages', 'verbose')
        kw.setdefault('renumber', True)
        kw.setdefault('reorder', False)
        kw.setdefault('reorder_max_rounds', 100)
        kw.setdefault('pedigree_is_renumbered', False)
        kw.setdefault('set_generations', False)
        kw.setdefault('gen_coeff', False)
        kw.setdefault('set_ancestors', False)
        kw.setdefault('set_alleles', False)
        kw.setdefault('set_offspring', False)
        kw.setdefault('set_sexes', False)
        kw.setdefault('assign_sexes', False)
        kw.setdefault('pedcomp', False)
        kw.setdefault('pedcomp_gens', 3)
        kw.setdefault('sepchar', ' ')
        kw.setdefault('alleles_sepchar', '/')
        kw.setdefault('counter', 1000)
        kw.setdefault('slow_reorder', True)
        kw.setdefault('resolve_duplicates', False)
        # Scientific validation (see PyPedal/pyp_validate.py).
        #   validate    -- cheap runtime postconditions and preconditions on
        #                  values a routine has already computed. On by default:
        #                  raising beats returning a mathematically impossible
        #                  number. Setting this False is unsafe.
        #   diagnostics -- opt-in expensive checking (e.g. full O(n^2) matrix
        #                  scans). Off by default and NOT implied by validate,
        #                  because PyPedal must stay usable at >100,000 animals.
        kw.setdefault('validate', True)
        kw.setdefault('diagnostics', False)


    def _set_animal_defaults(self, kw: Dict[str, Any]) -> None:
        """Default missing values for NewAnimal objects."""
        retired = [
            key for key in ("missing_byear", "missing_bdate") if key in kw
        ]
        for key in retired:
            kw.pop(key, None)
        if retired:
            warnings.warn(
                "missing_byear/missing_bdate no longer define unknown birth "
                "chronology; unknown recorded year and date are None. Old INI "
                "values are ignored and are not copied onto animals. Use "
                "legacy_missing_byear_token or legacy_missing_bdate_token to "
                "import files that stored 1800, 1900, or 01011800 as missing.",
                DeprecationWarning,
                stacklevel=3,
            )
        kw.setdefault("legacy_missing_byear_token", None)
        kw.setdefault("legacy_missing_bdate_token", None)
        kw.setdefault("estimate_birth_dates", False)
        kw.setdefault("vital_rate_profile", None)
        defaults = {
            "missing_parent": 0,
            "missing_name": "Unknown_Name",
            "missing_breed": "Unknown_Breed",
            "missing_herd": "Unknown_Herd",
            "missing_sex": "u",
            "missing_inbreeding": 0.0,
            "missing_genomic_inbreeding": 0.0,
            "missing_homozygosity": -999.0,
            "missing_alive": 0,
            "missing_age": -999,
            "missing_gen": -999.0,
            "missing_gencoeff": -999.0,
            "missing_igen": -999.0,
            "missing_pedcomp": -999.0,
            # The generic missing-DERIVED-result marker, distinct from the
            # input sentinels above. 2.0.4 established it as the float -999.
            # (pyp_newclasses.py:285-286); baa9211 carried it through the first
            # Python 3 conversion; bf0d2b2 dropped it when this helper was
            # extracted from __init__, stranding the three reads in
            # pyp_utils.set_age() that have been there since 2.0.4. Restored so
            # set_age still has a missing-derived-result marker.
            # Float, not the int -999 of missing_age: the two are different keys
            # with different consumers, and the type is part of the contract.
            "missing_value": -999.0,
            "missing_alleles": ["", ""],
            # PyPedal 2.0.4 sets this key twice in NewPedigree.__init__: 'Unknown'
            # at :181-182, which runs first and therefore wins, and '' at :283-284,
            # which is guarded by `not in kw.keys()` and is dead code. 'Unknown' is
            # 2.0.4's effective runtime default, and it is what the differential
            # harness measures on every corpus pedigree without a 'u' column.
            "missing_userfield": "Unknown"
        }
        for key, value in defaults.items():
            kw.setdefault(key, value)


    def __add__(self, other, filename=None, debugLoad=False):
        """
        Method to add two pedigrees and return a new pedigree representing the merged pedigrees.

        Parameters
        ----------
        other : NewPedigree
            Reference to the NewPedigree object of animals to add.
        filename : str, optional
            The file name to be used for saving/loading.
        debugLoad : bool, default: False
            Toggle debugging messages during pedigree loading on (True) or off (False).

        Returns
        -------
        NewPedigree or bool
            A NewPedigree object with the contents of the merged pedigrees or False on failure.
        """
        if isinstance(self, NewPedigree) and isinstance(other, NewPedigree):
            logger.info("Adding pedigrees %s and %s", self.kw['pedname'], other.kw['pedname'])
            if self.kw.get('debug_messages', False):
                print("[DEBUG]: self and other are both NewPedigree objects. Combining them.")
                print(f"[DEBUG]: Using match rule: {self.kw['match_rule']}")
            logger.info("Using match rule %s to merge pedigrees", self.kw['match_rule'])

            # Ensure both pedigrees are renumbered
            if not self.kw['pedigree_is_renumbered']:
                self.renumber()
                logger.info("Renumbering pedigree %s", self.kw['pedname'])
            if not other.kw['pedigree_is_renumbered']:
                other.renumber()
                logger.info("Renumbering pedigree %s", other.kw['pedname'])

            ped_to_write = {"a": {}, "b": {}}
            for a in self.pedigree:
                ped_to_write["a"][a.animalID] = True
                for b in other.pedigree:
                    mismatches = 0  # Count mismatches between animals
                    for match in self.kw["match_rule"]:
                        if match in ["a", "A"]:
                            if a.originalID != b.originalID:
                                mismatches += 1
                        elif match in ["s", "S"]:
                            sire_a = self.pedigree[a.sireID] if a.sireID != self.kw["missing_parent"] else None
                            sire_b = other.pedigree[b.sireID] if b.sireID != other.kw["missing_parent"] else None
                            if sire_a and sire_b and sire_a.originalID != sire_b.originalID:
                                mismatches += 1
                            elif (sire_a and not sire_b) or (sire_b and not sire_a):
                                mismatches += 1
                        elif match in ["d", "D"]:
                            dam_a = self.pedigree[a.damID] if a.damID != self.kw["missing_parent"] else None
                            dam_b = other.pedigree[b.damID] if b.damID != other.kw["missing_parent"] else None
                            if dam_a and dam_b and dam_a.originalID != dam_b.originalID:
                                mismatches += 1
                            elif (dam_a and not dam_b) or (dam_b and not dam_a):
                                mismatches += 1
                        elif getattr(a, self.new_animal_attr[match], None) != getattr(b, self.new_animal_attr[match], None):
                            mismatches += 1

                    if mismatches == 0:
                        ped_to_write["b"][b.animalID] = False
                        if self.kw.get("debug_messages", False):
                            print(f"[DEBUG]: Animals {a.animalID} and {b.animalID} are identical.")
                    else:
                        ped_to_write["b"][b.animalID] = True
                        if self.kw.get("debug_messages", False):
                            print(f"[DEBUG]: Animals {a.animalID} and {b.animalID} are different.")

            # Create the filename for the merged pedigree if not provided
            if not filename:
                filename = f"{self.kw['pedname']}_{other.kw['pedname']}.ped"
                print(f"[INFO]: Using filename = {filename}")

            # Save the unique animals from both pedigrees
            self.save(filename=filename, write_list=ped_to_write['a'], pedformat=self.kw['pedformat'], originalID=True)
            other.save(filename=filename, write_list=ped_to_write['b'], pedformat=self.kw['pedformat'], originalID=True, append=True)

            # Load the new merged pedigree
            merged_pedname = f"Merged Pedigree: {self.kw['pedname']} + {other.kw['pedname']}"
            new_options = {
                "messages": self.kw["messages"],
                "pedname": merged_pedname,
                "renumber": True,
                "pedfile": filename,
                "pedformat": self.kw["pedformat"]
            }

            try:
                new_pedigree = loadPedigree(new_options, debugLoad=debugLoad)
                if self.kw["messages"] == "verbose":
                    print(f"[INFO]: Loaded merged pedigree {merged_pedname} from file {filename}!")
                logger.info("Loaded merged pedigree %s from file %s.", merged_pedname, filename)
                return new_pedigree
            except Exception as e:
                if self.kw['messages'] == 'verbose':
                    print(f"[ERROR]: Could not load merged pedigree {merged_pedname} from file {filename}! Error: {e}")
                logger.error("Could not load merged pedigree %s from file %s! Error: %s", merged_pedname, filename, e)
                return False
        else:
            logger.error("Cannot complete __add__() operation because types do not match.")
            return NotImplemented

    def __sub__(self, other, filename=False, debugLoad=False):
        """
        Method to subtract two pedigrees and return a new pedigree representing the
        first pedigree without any animals shared in common with the second pedigree,
        or: A - B = A - (A ∩ B).

        Parameters
        ----------
        other : NewPedigree
            The pedigree object to subtract from the current pedigree.
        filename : str, optional
            The file name to save the resulting pedigree.
        debugLoad : bool, optional
            Toggle debugging messages during pedigree loading.

        Returns
        -------
        NewPedigree or False
            The resulting pedigree object or False on failure.
        """
        if isinstance(other, NewPedigree):
            logger.info("Subtracting pedigrees %s and %s", self.kw['pedname'], other.kw['pedname'])
            logger.info("Using match rule %s to subtract pedigrees", self.kw['match_rule'])

            # Ensure pedigrees are renumbered
            if not self.kw['pedigree_is_renumbered']:
                self.renumber()
                logger.info("Renumbering pedigree %s", self.kw['pedname'])
            if not other.kw['pedigree_is_renumbered']:
                other.renumber()
                logger.info("Renumbering pedigree %s", other.kw['pedname'])

            # Track animals to write
            ped_to_write = {'a': {}, 'b': {}}
            for a in self.pedigree:
                ped_to_write['a'][a.animalID] = True
            for b in other.pedigree:
                ped_to_write['b'][b.animalID] = False

            # Compare animals based on the match rule
            for a in self.pedigree:
                for b in other.pedigree:
                    mismatches = 0
                    for match in self.kw['match_rule']:
                        if match in ['a', 'A']:
                            if a.originalID != b.originalID:
                                mismatches += 1
                        elif match in ['s', 'S']:
                            if self.pedigree[a.sireID - 1].originalID != other.pedigree[b.sireID - 1].originalID:
                                mismatches += 1
                        elif match in ['d', 'D']:
                            if self.pedigree[a.damID - 1].originalID != other.pedigree[b.damID - 1].originalID:
                                mismatches += 1
                        elif getattr(a, self.new_animal_attr[match]) != getattr(b, self.new_animal_attr[match]):
                            mismatches += 1

                    # If there are no mismatches, animals are identical
                    if mismatches == 0:
                        ped_to_write['a'][a.animalID] = False
                    else:
                        ped_to_write['b'][b.animalID] = True

            # Determine filename
            if not filename:
                filename = f"{self.kw['pedname']}_{other.kw['pedname']}.ped"
                print(f"[INFO]: filename = {filename}")

            # Save and merge unique animals
            self.save(filename=filename, write_list=ped_to_write['a'],
                    pedformat=self.kw['pedformat'], originalID=True)
            other.save(filename=filename, write_list=ped_to_write['b'],
                    pedformat=self.kw['pedformat'], originalID=True, append=True)

            # Load the resulting pedigree
            merged_pedname = f"Subtracted Pedigree: {self.kw['pedname']} - {other.kw['pedname']}"
            new_options = {
                'messages': self.kw['messages'],
                'pedname': merged_pedname,
                'renumber': True,
                'pedfile': filename,
                'pedformat': self.kw['pedformat'],
            }
            try:
                new_pedigree = loadPedigree(new_options, debugLoad=debugLoad)
                if self.kw['messages'] == 'verbose':
                    print(f"[INFO]: Loaded subtracted pedigree {merged_pedname} from file {filename}!")
                logger.info("Loaded subtracted pedigree %s from file %s.", merged_pedname, filename)
                return new_pedigree
            except Exception as e:
                if self.kw['messages'] == 'verbose':
                    print(f"[ERROR]: Could not load subtracted pedigree {merged_pedname} from file {filename}! Error: {e}")
                logger.error("Could not load subtracted pedigree %s from file %s. Error: %s", merged_pedname, filename, e)
                return False
        else:
            logger.error("Cannot complete __sub__() operation because types do not match.")
            return NotImplemented

    def union(self, other, filename=False, debugLoad=False):
        """
        union() is an alias to NewPedigree::__add__().

        Parameters
        ----------
        other : NewPedigree
            The pedigree object to be merged with the current pedigree.
        filename : str, optional
            The file name to save the resulting pedigree.
        debugLoad : bool, optional
            Toggle debugging messages during pedigree loading.

        Returns
        -------
        NewPedigree or NotImplemented
            The resulting merged pedigree or NotImplemented if types do not match.
        """
        return self.__add__(other, filename=filename, debugLoad=debugLoad)

    def intersection(self, other, newpedname='intersected_pedigree'):
        """
        intersection() returns a PyPedal pedigree object which contains the animals that are common to
        both input pedigrees. If there are no animals in common between the two pedigrees, a value
        of False is returned.

        Parameters
        ----------
        other : NewPedigree
            Another PyPedal pedigree object to compare with the current one.
        newpedname : str, optional
            The name of the new pedigree resulting from the intersection.

        Returns
        -------
        NewPedigree or bool
            A new PyPedal pedigree containing the animals common to both input pedigrees,
            or False on failure.
        """
        if isinstance(self, NewPedigree) and isinstance(other, NewPedigree):
            logger.info('Computing intersection of pedigrees %s and %s', self.kw['pedname'], other.kw['pedname'])
            logger.info('Using match rule %s to compare pedigrees', self.kw['match_rule'])

            # Ensure pedigrees are renumbered
            if not self.kw.get('pedigree_is_renumbered', False):
                self.renumber()
                logger.info('Renumbering pedigree %s', self.kw['pedname'])
            if not other.kw.get('pedigree_is_renumbered', False):
                other.renumber()
                logger.info('Renumbering pedigree %s', other.kw['pedname'])

            # Determine common animals based on match rules
            animals_to_write = []
            for a in self.pedigree:
                for b in other.pedigree:
                    matches = 0  # Count the number of matching attributes
                    for match in self.kw['match_rule']:
                        if match == 'a' and a.originalID == b.originalID:
                            matches += 1
                        elif match == 'A' and a.name == b.name:
                            matches += 1
                        elif match == 's' and self.pedigree[a.sireID - 1].originalID == other.pedigree[b.sireID - 1].originalID:
                            matches += 1
                        elif match == 'S' and a.sireName == b.sireName:
                            matches += 1
                        elif match == 'd' and self.pedigree[a.damID - 1].originalID == other.pedigree[b.damID - 1].originalID:
                            matches += 1
                        elif match == 'D' and a.damName == b.damName:
                            matches += 1
                        elif getattr(a, self.new_animal_attr[match]) == getattr(b, other.new_animal_attr[match]):
                            matches += 1

                    if matches == len(self.kw['match_rule']):
                        animals_to_write.append(a)

            # Save intersected pedigree to a file
            filename = f'{newpedname}.ped'
            if self.kw.get('debug_messages', False):
                logger.info('Filename: %s', filename)
                logger.info('There are %s animals in the intersection of pedigrees %s and %s.',
                            len(animals_to_write), self.kw['pedname'], other.kw['pedname'])
            pyp_io.save_newanimals_to_file(animals_to_write, filename, self.kw, self.new_animal_attr)

            # Create new options dictionary for the intersected pedigree
            intersect_pedname = f'Intersected Pedigree: {self.kw["pedname"]} + {other.kw["pedname"]}'
            new_options = {
                'messages': self.kw['messages'],
                'pedname': intersect_pedname,
                'renumber': True,
                'pedfile': filename,
                'pedformat': self.kw['pedformat'],
                'sepchar': self.kw['sepchar'],
            }

            # Load the intersected pedigree
            try:
                new_pedigree = loadPedigree(new_options, debugLoad=True)
                if self.kw['messages'] == 'verbose':
                    logger.info('Loaded intersected pedigree %s from file %s', intersect_pedname, filename)
                return new_pedigree
            except Exception as e:
                if self.kw['messages'] == 'verbose':
                    logger.error('Could not load intersected pedigree %s from file %s: %s', intersect_pedname, filename, e)
                return False
        else:
            logger.error('Cannot complete intersection operation because types do not match.')
            return NotImplemented

    def load(
        self,
        pedsource: PedigreeSource = 'file',
        pedgraph=None,
        pedstream='',
        animallist=None,
    ):
        """
        load() wraps several processes useful for loading and preparing a pedigree for
        use in an analysis, including reading the animals into a list of animal objects,
        forming lists of sires and dams, checking for common errors, setting ancestor
        flags, and renumbering the pedigree.

        Parameters
        ----------
        pedsource : str
            Source of the pedigree ('file', 'db', 'graph', 'graphfile', 'null', 'animallist',
            'gedcomfile', 'textstream').
        pedgraph : networkx.DiGraph, optional
            Directed graph from which to load the pedigree (used when `pedsource='graph'`).
        pedstream : str, optional
            Stream of text from which to load the pedigree (used when `pedsource='textstream'`).
        animallist : list, optional
            List of NewAnimal objects from which to create a pedigree (used when `pedsource='animallist'`).

        Returns
        -------
        None
        """
        valid_sources = [
            'file', 'db', 'graph', 'graphfile', 'null', 'animallist',
            'gedcomfile', 'textstream'
        ]

        if pedsource not in valid_sources:
            logger.error('Invalid pedsource provided: %s', pedsource)
            raise pyp_errors.PyPedalUsageError(
                'Invalid pedsource %r. Valid sources are: %s.'
                % (pedsource, ', '.join(valid_sources)))

        # If simulating a pedigree, generate records and then fall through
        # to the common post-load finalization below (reorder/renumber,
        # PedigreeMetadata, id maps). An early return here was a Python 3
        # port regression: both 2.0.4 and pypedal3 use this if/elif shape
        # so simulated pedigrees reach the same lifecycle as file loads.
        if self.kw.get('simulate_pedigree', False):
            self.simulate()

        elif pedsource == 'db':
            self.kw['pedformat'] = 'ASDx'
            self.kw['sepchar'] = ','
            logger.info('Loading from database %s.%s at %s.',
                        self.kw['database_name'], self.kw['database_table'], pyp_utils.pyp_nice_time())

            try:
                if pyp_db.doesTableExist(self):
                    conn = pyp_db.connectToDatabase(self)
                    if conn:
                        sql = self.kw['database_sql'] % self.kw['database_table']
                        cursor = conn.cursor()
                        cursor.execute(sql)
                        dbstream = cursor.fetchall()
                        conn.close()
                        self.preprocess(dbstream=dbstream)
                    else:
                        logger.error('Failed to connect to database %s', self.kw['database_name'])
                        raise pyp_errors.PyPedalConfigurationError(
                            'Could not connect to database %r.'
                            % self.kw['database_name'])
            except pyp_errors.PyPedalError:
                raise
            except Exception as e:
                logger.error('Database load failed: %s', e)
                raise pyp_errors.PyPedalConfigurationError(
                    'Loading from database %r failed: %s'
                    % (self.kw.get('database_name'), e)) from e

        elif pedsource == 'graph':
            if pedgraph is not None:
                self.fromgraph(pedgraph)
            else:
                logger.error('No pedgraph provided for graph source')
                raise pyp_errors.PyPedalUsageError(
                    "pedsource='graph' requires a pedgraph argument, but none "
                    'was given.')

        elif pedsource == 'graphfile':
            try:
                import networkx as nx
                if self.kw['pedfile']:
                    pedgraph = nx.read_adjlist(self.kw['pedfile'])
                    self.fromgraph(pedgraph)
                else:
                    logger.error('No filename provided for graphfile source')
                    raise pyp_errors.PyPedalConfigurationError(
                        "pedsource='graphfile' requires kw['pedfile'] to be "
                        'set, but it is empty.')
            except ImportError as e:
                logger.error('NetworkX module not available')
                raise pyp_errors.PyPedalDependencyError(
                    "pedsource='graphfile' needs NetworkX, which is not "
                    "installed. Install it with: pip install networkx") from e
            except pyp_errors.PyPedalError:
                raise
            except Exception as e:
                logger.error('Failed to load graph from file: %s', e)
                raise pyp_errors.PyPedalPedigreeSourceError(
                    'Could not read an adjacency list from %r: %s'
                    % (self.kw['pedfile'], e)) from e

        elif pedsource == 'null':
            try:
                self.fromnull()
            except Exception as e:
                logger.error('Failed to create null pedigree: %s', e)
                raise pyp_errors.PyPedalInternalError(
                    'Constructing a null pedigree failed: %s' % e) from e

        elif pedsource == 'animallist':
            if animallist and len(animallist) > 0:
                try:
                    self.fromanimallist(animallist)
                except Exception as e:
                    logger.error('Failed to create pedigree from animallist: %s', e)
                    raise pyp_errors.PyPedalInputError(
                        'Building a pedigree from the supplied animallist '
                        'failed: %s' % e) from e
            else:
                logger.error('Empty or missing animallist provided')
                raise pyp_errors.PyPedalUsageError(
                    "pedsource='animallist' requires a non-empty animallist, "
                    'but none was given.')

        elif pedsource == 'gedcomfile':
            if self.kw['pedfile']:
                pedformat = pyp_io.load_from_gedcom(
                    infilename=self.kw['pedfile'],
                    standalone=False,
                    messages=self.kw['messages'],
                    missing_sex=self.kw['missing_sex'],
                    missing_parent=self.kw['missing_parent'],
                    missing_name=self.kw['missing_name'],
                )
                if pedformat != 'xxxx':
                    self.kw.update({
                        'pedformat': pedformat,
                        'sepchar': ',',
                        'pedfile': f"{self.kw['pedfile']}.tmp"
                    })
                    self.preprocess()
                else:
                    logger.error('Invalid pedigree format from GEDCOM file')
                    raise pyp_errors.PyPedalPedigreeSourceError(
                        'Could not determine a pedigree format from the GEDCOM '
                        'file %r.' % self.kw['pedfile'])
            else:
                logger.error('No filename provided for GEDCOM file')
                raise pyp_errors.PyPedalConfigurationError(
                    "pedsource='gedcomfile' requires kw['pedfile'] to be set, "
                    'but it is empty.')

        elif pedsource == 'textstream':
            self.kw.update({'pedformat': 'ASD', 'sepchar': ','})
            logger.info('Preprocessing a textstream')
            self.preprocess(textstream=pedstream)

        else:  # Default is loading from a file
            if not self.kw['pedfile']:
                print(f"pedfile is None or not set: {self.kw.get('pedfile')}")
                raise ValueError("The key 'pedfile' must be set in self.kw before calling preprocess().")

            logger.info('Preprocessing %s', self.kw['pedfile'])
            self.preprocess()

        # Post-processing: renumbering, assigning generations, etc.
        if self.kw.get('reorder') and not self.kw.get('renumber'):
            logger.info('Reordering pedigree')
            self.pedigree = (
                pyp_utils.fast_reorder(
                    self.pedigree,
                    missingparent=self.kw['missing_parent'],
                )
                if not self.kw.get('slow_reorder')
                # Keywords, not positions. reorder()'s signature is
                # (myped, filetag, io, missingparent, debug, max_rounds), so the
                # positional form bound missing_parent to filetag and
                # reorder_max_rounds to io, leaving missingparent and max_rounds
                # at their defaults. Latent at the defaults,
                # silently wrong as soon as missing_parent is not 0.
                else pyp_utils.reorder(
                    self.pedigree,
                    missingparent=self.kw['missing_parent'],
                    max_rounds=self.kw['reorder_max_rounds'],
                )
            )

        # Read SNP genotypes before renumbering, so renumber_snp_ids() below has
        # something to renumber. Until genotypes were loaded here, nothing
        # read kw['snpfile']: it had a default and `self.snp = False` was the only
        # assignment to the attribute, so every pyp_snp entry point guarded on
        # `pedobj.snp is False` always took the failure path and
        # form_grm_from_snp() was unreachable in production.
        if self.kw.get('snpfile'):
            logger.info('Reading SNP genotypes from %s', self.kw['snpfile'])
            pyp_snp.load_snp_file(self)

        if self.kw.get('renumber'):
            self.renumber()
            pyp_snp.renumber_snp_ids(self)

        if self.kw.get('set_ancestors'):
            logger.info('Setting ancestor flags')
            pyp_utils.set_ancestor_flag(self)

        if self.kw.get('set_sexes') or self.kw.get('assign_sexes'):
            logger.info('Assigning sexes')
            pyp_utils.set_sexes(self)

        if self.kw.get('set_alleles'):
            logger.info('Gene dropping to compute founder genome equivalents')
            pyp_metrics.effective_founder_genomes(self)

        if self.kw.get('form_nrm'):
            logger.info('Forming numerator relationship matrix')
            self.nrm = NewAMatrix(self.kw)
            self.nrm.form_a_matrix(self.pedigree)

        if self.kw.get('set_offspring') and not self.kw.get('renumber'):
            logger.info('Assigning offspring')
            pyp_utils.set_offspring(self)

        # Metadata and summary
        logger.info('Creating pedigree metadata')
        self.metadata = PedigreeMetadata(self.pedigree, self.kw, self.snp)
        # Carry the implicit-parent record onto the metadata so callers can
        # reconcile input rows against post-load counts.
        self.metadata.implicit_parent_ids = list(self._implicit_parents)
        self.metadata.num_implicit_parents = len(self._implicit_parents)

        # set_generation reads pedobj.metadata.num_records. Until this
        # point metadata is the empty dict __init__ installed, so a
        # load-time request used to raise AttributeError, get swallowed,
        # and leave igen at the initialisation sentinel.
        if self.kw.get('set_generations'):
            logger.info('Assigning generations')
            pyp_utils.set_generation(self)

        pyp_chronology.validate_recorded_chronology(self)
        if self.kw.get('estimate_birth_dates') and self.kw.get('vital_rate_profile'):
            pyp_chronology.estimate_birth_date_ranges(
                self, self.kw.get('vital_rate_profile')
            )

        if self.kw.get('messages') != 'quiet' and self.kw.get('pedigree_summary'):
            self.metadata.printme()

        if self.kw.get('pedcomp'):
            logger.info('Calculating pedigree completeness for %s generations', self.kw['pedcomp_gens'])
            pyp_metrics.pedigree_completeness(self, self.kw['pedcomp_gens'])

    def oldsave(self, filename='', outformat='o', idformat='o'):
        """
        oldsave() writes a PyPedal pedigree to a user-specified file. The saved pedigree
        includes all fields recognized by PyPedal, not just the original fields read
        from the input pedigree file.

        Parameters
        ----------
        filename : str, optional
            The file to which the pedigree should be written (default is '').
        outformat : str, optional
            The format in which the pedigree should be written: 'o' for original (as read),
            'l' for long version (all available variables) (default is 'o').
        idformat : str, optional
            Write 'o' (original) or 'r' (renumbered) animal, sire, and dam IDs (default is 'o').

        Returns
        -------
        bool
            True on success, False on failure.
        """
        # Default filename to avoid overwriting user data
        if not filename:
            filename = f"{self.kw['filetag']}_saved.ped"
            if self.kw.get('messages') == 'verbose':
                print(f"[WARNING]: Saving pedigree to file {filename} to avoid overwriting {self.kw['pedfile']}.")
            logger.warning("Saving pedigree to file %s to avoid overwriting %s.", filename, self.kw['pedfile'])

        try:
            with open(filename, 'w', encoding='utf-8') as ofh:
                if self.kw.get('messages') == 'verbose':
                    print(f"[INFO]: Opened file {filename} for pedigree save at {pyp_utils.pyp_nice_time()}.")
                logger.info("Opened file %s for pedigree save at %s.", filename, pyp_utils.pyp_nice_time())

                # Determine the new pedigree format
                if outformat == 'l':
                    _newpedformat = 'asdgx'
                    _newpedformat += 'y' if 'y' in self.kw['pedformat'] else 'b'
                    _newpedformat += 'frnleh'
                else:
                    _newpedformat = f"{self.kw['pedformat']}f" if self.kw.get('f_computed') else self.kw['pedformat']

                # Write file header
                ofh.write(f"# {filename} created by PyPedal at {pyp_utils.pyp_nice_time()}\n")
                ofh.write("# Current pedigree metadata:\n")
                ofh.write(f"#\tpedigree file: {filename}\n")
                ofh.write(f"#\tpedigree name: {self.kw['pedname']}\n")
                ofh.write(f"#\tpedigree format: '{_newpedformat}'\n")
                if idformat == 'o':
                    ofh.write("#\tNOTE: Animal, sire, and dam IDs are RENUMBERED IDs, not original IDs!\n")
                ofh.write("# Original pedigree metadata:\n")
                ofh.write(f"#\tpedigree file: {self.kw['pedfile']}\n")
                ofh.write(f"#\tpedigree name: {self.kw['pedname']}\n")
                ofh.write(f"#\tpedigree format: {self.kw['pedformat']}\n")

                # Write pedigree records.
                # idformat='o' looks parents up by 1-based animalID. A missing
                # parent is the configured sentinel (default 0); indexing it
                # is pedigree[-1] and writes the last real animal as a parent.
                # Do not fall back to missing on a
                # failed real-parent lookup -- that would invent the inverse
                # lie. PyPedalPedigreeStructureError is re-raised below.
                missing = self.kw['missing_parent']

                def original_parent(animal, parent_id, role):
                    if parent_id == missing:
                        return missing
                    try:
                        idx = int(parent_id) - 1
                    except (TypeError, ValueError):
                        idx = -1
                    if idx < 0 or idx >= len(self.pedigree):
                        raise pyp_errors.PyPedalPedigreeStructureError(
                            "oldsave(): animal %s names %s as its %s, but no "
                            "animal with that ID has a record in the pedigree."
                            % (animal.originalID, parent_id, role),
                            animals=[animal.originalID, parent_id],
                        )
                    return self.pedigree[idx].originalID

                for _a in self.pedigree:
                    if idformat == 'o':
                        _outstring = (
                            f"{_a.originalID} "
                            f"{original_parent(_a, _a.sireID, 'sire')} "
                            f"{original_parent(_a, _a.damID, 'dam')}"
                        )
                    else:
                        _outstring = f"{_a.animalID} {_a.sireID} {_a.damID}"

                    if 'g' in _newpedformat:
                        _outstring += f" {_a.gen}"
                    if 'p' in _newpedformat:
                        _outstring += f" {_a.gencoeff}"
                    if 'x' in _newpedformat:
                        _outstring += f" {_a.sex}"
                    if 'y' in _newpedformat:
                        _outstring += f" {pyp_chronology.format_year_token(_a.by)}"
                    if 'b' in _newpedformat:
                        _outstring += (
                            f" {pyp_chronology.format_date_token(_a.bd, _a.by, allow_year_only='y' not in _newpedformat)}"
                        )
                    if 'y' not in _newpedformat and 'b' not in _newpedformat:
                        _outstring += f" {pyp_chronology.format_year_token(_a.by)}"
                    if 'f' in _newpedformat:
                        _outstring += f" {_a.fa}"
                    if 'r' in _newpedformat:
                        _outstring += f" {_a.breed}"
                    if 'n' in _newpedformat:
                        _outstring += f" {_a.name}"
                    if 'l' in _newpedformat:
                        _outstring += f" {_a.alive}"
                    if 'e' in _newpedformat:
                        _outstring += f" {_a.age}"
                    if 'h' in _newpedformat or 'H' in _newpedformat:
                        _outstring += f" {_a.herd} {_a.originalHerd}"
                    if 'u' in _newpedformat:
                        _outstring += f" {_a.userField}"

                    ofh.write(f"{_outstring}\n")

            if self.kw.get('messages') == 'verbose':
                print(f"[INFO]: Closed file {filename} after pedigree save at {pyp_utils.pyp_nice_time()}.")
            logger.info("Closed file %s after pedigree save at %s.", filename, pyp_utils.pyp_nice_time())
            return True

        except pyp_errors.PyPedalError:
            raise
        except Exception as e:
            if self.kw.get('messages') == 'verbose':
                print(f"[ERROR]: Unable to open file {filename} for pedigree save!")
            logger.error("Unable to open file %s for pedigree save. Error: %s", filename, str(e))
            return False

    def save(self, filename='', pedformat='asd', sepchar=' ', append=False, write_list=None, originalID=False):
        """
        save() writes a PyPedal pedigree to a user-specified file. The saved pedigree
        includes all fields recognized by PyPedal, not just the original fields read
        from the input pedigree file.

        Parameters
        ----------
        filename : str, optional
            The file to which the pedigree should be written (default is '').
        pedformat : str, optional
            Pedigree format string for the pedigree to be written (default is 'asd').
        sepchar : str, optional
            Character used to separate columns in the output pedigree file (default is ' ').
        append : bool, optional
            Add animal records to an existing file instead of creating a new one (default is False).
        write_list : dict, optional
            Optional dictionary of animal records to save. Defaults to saving all animals.
        originalID : bool, optional
            Save original IDs or renumbered IDs (default is False).

        Returns
        -------
        bool
            True on success, False on failure.
        """
        # Validate and process the pedformat string
        pedformat_in = str(pedformat)
        pedformat = ''.join(pf for pf in pedformat_in if pf in self.pedformat_codes)
        if any(pf not in self.pedformat_codes for pf in pedformat_in):
            invalid_codes = [pf for pf in pedformat_in if pf not in self.pedformat_codes]
            if self.kw.get('messages') == 'verbose':
                print(f"[WARNING]: Invalid pedigree format codes {invalid_codes} in NewPedigree::save().")
            logger.warning("Invalid pedigree format codes %s in NewPedigree::save().", invalid_codes)

        # Warn about incomplete information in the pedformat
        if 'asd' not in pedformat.lower():
            if self.kw.get('messages') == 'verbose':
                print(f"[WARNING]: Pedigree format {pedformat} may lead to incomplete parentage information.")
            logger.warning("Pedigree format %s may lead to incomplete parentage information.", pedformat)

        # Validate sepchar
        if not sepchar:
            sepchar = self.kw['sepchar']
            if self.kw.get('messages') == 'verbose':
                print(f"[WARNING]: Invalid sepchar ''. Changed to '{sepchar}'.")
            logger.warning("Invalid sepchar ''. Changed to '%s'.", sepchar)

        # Ensure filename is valid
        if not filename:
            filename = f"{self.kw['filetag']}_saved.ped"
            if self.kw.get('messages') == 'verbose':
                print(f"[WARNING]: Saving pedigree to file {filename} to avoid overwriting {self.kw['pedfile']}.")
            logger.warning("Saving pedigree to file %s to avoid overwriting %s.", filename, self.kw['pedfile'])

        mode = 'a' if append else 'w'
        try:
            with open(filename, mode, encoding='utf-8') as ofh:
                if self.kw.get('messages') == 'verbose':
                    print(f"[INFO]: Opened file {filename} for saving at {pyp_utils.pyp_nice_time()}.")
                logger.info("Opened file %s for saving at %s.", filename, pyp_utils.pyp_nice_time())

                if append:
                    ofh.write(f"# {filename} created by PyPedal at {pyp_utils.pyp_nice_time()}\n")
                    ofh.write("# Current pedigree metadata:\n")
                    ofh.write(f"#\tpedigree file: {filename}\n")
                    ofh.write(f"#\tpedigree name: {self.kw['pedname']}\n")
                    ofh.write(f"#\tpedigree format: {pedformat}\n")
                    if self.kw['pedigree_is_renumbered']:
                        note = "original IDs" if originalID else "renumbered IDs"
                        ofh.write(f"#\tNOTE: Animal, sire, and dam IDs are {note}!\n")
                    ofh.write("# Original pedigree metadata:\n")
                    ofh.write(f"#\tpedigree file: {self.kw['pedfile']}\n")
                    ofh.write(f"#\tpedigree name: {self.kw['pedname']}\n")
                    ofh.write(f"#\tpedigree format: {self.kw['pedformat']}\n")

                for _a in self.pedigree:
                    if not write_list or write_list.get(_a.animalID, True):
                        _outstring = []
                        for pf in pedformat:
                            formatted = pyp_chronology.format_pedigree_field(pf, _a, pedformat)
                            if formatted is not None:
                                value = formatted
                            elif not originalID:
                                value = getattr(_a, self.new_animal_attr[pf], '')
                            else:
                                if pf in ['a', 'A']:
                                    value = _a.originalID
                                elif pf in ['s', 'S']:
                                    value = (
                                        self.pedigree[_a.sireID - 1].originalID
                                        if _a.sireID != self.kw['missing_parent'] else 0
                                    )
                                elif pf in ['d', 'D']:
                                    value = (
                                        self.pedigree[_a.damID - 1].originalID
                                        if _a.damID != self.kw['missing_parent'] else 0
                                    )
                                else:
                                    value = getattr(_a, self.new_animal_attr[pf], '')
                            _outstring.append(str(value))
                        ofh.write(sepchar.join(_outstring) + '\n')

                if self.kw.get('messages') == 'verbose':
                    print(f"[INFO]: Closed file {filename} after saving at {pyp_utils.pyp_nice_time()}.")
                logger.info("Closed file %s after saving at %s.", filename, pyp_utils.pyp_nice_time())
            return True

        except Exception as e:
            if self.kw.get('messages') == 'verbose':
                print(f"[ERROR]: Unable to save pedigree to file {filename}. Error: {e}")
            logger.error("Unable to save pedigree to file %s. Error: %s", filename, e)
            return False

    def estimate_birth_date_ranges(self, profile=None):
        """Fill ``birth_date_estimate`` from an explicit vital-rate profile."""
        return pyp_chronology.estimate_birth_date_ranges(self, profile=profile)

    def savegraph(self, pedoutfile=None, pedgraph=None):
        """
        Save a pedigree to a file as an adjacency list.

        Parameters
        ----------
        pedoutfile : str, optional
            The name of the file to which the graph is written. If not provided,
            a default name is derived from the pedigree file name.
        pedgraph : networkx.DiGraph, optional
            A NetworkX directed graph object representing the pedigree. If not provided,
            the method will attempt to generate one.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        _retval = False

        try:
            # Use a default output file name if none is provided
            if not pedoutfile:
                pedoutfile = f"{self.kw['pedfile']}.adjlist"

            # If no graph is provided, attempt to create one
            if not pedgraph:
                try:
                    from . import pyp_network
                    pedgraph = pyp_network.ped_to_graph(self)
                    logger.info("[savegraph]: Converted pedigree to a directed graph.")
                except Exception as e:
                    logger.error("[savegraph]: Unable to convert pedigree to a directed graph. Error: %s", e)
                    return False

            # Save the graph to a file as an adjacency list
            try:
                import networkx as nx
                nx.write_adjlist(pedgraph, pedoutfile)
                logger.info("[savegraph]: Saved directed graph to file %s.", pedoutfile)
                _retval = True
            except Exception as e:
                logger.error("[savegraph]: Unable to save directed graph to file %s. Error: %s", pedoutfile, e)
                _retval = False

        except Exception as e:
            logger.error("[savegraph]: An unexpected error occurred. Error: %s", e)
            _retval = False

        return _retval

    def savegedcom(self, pedoutfile=None):
        """
        Save a pedigree to a file in GEDCOM 5.5 format.

        Parameters
        ----------
        pedoutfile : str, optional
            The name of the file to which the pedigree is written. If not provided,
            a default name is derived from the pedigree file name.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        _retval = False

        try:
            # Use a default output file name if none is provided
            if not pedoutfile:
                pedoutfile = f"{self.kw['pedfile']}.ged"

            try:
                # Save the pedigree to a GEDCOM file
                from . import pyp_io
                pyp_io.save_to_gedcom(self, pedoutfile)
                logger.info("[savegedcom]: Saved GEDCOM pedigree to the file %s.", pedoutfile)
                _retval = True
            except Exception as e:
                logger.error("[savegedcom]: Unable to save GEDCOM pedigree to the file %s. Error: %s", pedoutfile, e)
                _retval = False

        except Exception as e:
            logger.error("[savegedcom]: An unexpected error occurred. Error: %s", e)
            _retval = False

        return _retval

    def savedb(self, drop=False):
        """
        Saves a pedigree to a database table in ASDx format for NewAnimals and LightAnimals.

        Parameters
        ----------
        drop : bool
            Indicates if existing data should be deleted (True) or kept (False). Default is False.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        _retval = False

        try:
            _savedb_status = False
            _table_loaded = False
            _table_created = False

            if self.kw['database_compatibility'] == 'sqlite':
                # Handle SQLite database operations
                if pyp_db.doesTableExist(self):
                    if drop:
                        pyp_db.deleteTable(self)

                conn = pyp_db.connectToDatabase(self)
                if conn:
                    cursor = conn.cursor()
                    if not pyp_db.doesTableExist(self):
                        try:
                            sql = f"""
                            CREATE TABLE {self.kw['database_table']} (
                                animalName VARCHAR(128) PRIMARY KEY,
                                sireName   VARCHAR(128),
                                damName    VARCHAR(128),
                                sex        CHAR(1)
                            );
                            """
                            cursor.execute(sql)
                            conn.commit()
                            _table_created = True
                        except Exception as e:
                            logger.error('[savedb]: Error creating table. Error: %s', e)

                    else:
                        if self.kw['messages'] == 'verbose':
                            print(f"[WARNING]: The table {self.kw['database_table']} already exists in database "
                                f"{self.kw['database_name']} and you chose to keep existing data. This may lead to "
                                "duplicate data or multiple pedigrees stored in the same table!")
                        logger.warning("The table %s already exists in database %s and you chose to keep existing data. "
                                        "This may lead to duplicate data or multiple pedigrees stored in the same table!",
                                        self.kw['database_table'], self.kw['database_name'])

                    # Load the pedigree data into the table
                    if pyp_db.doesTableExist(self):
                        try:
                            for p in self.pedigree:
                                an = p.name
                                si = p.sireName if p.sireName != self.kw['missing_name'] else self.kw['missing_parent']
                                da = p.damName if p.damName != self.kw['missing_name'] else self.kw['missing_parent']

                                sql = f"""
                                INSERT INTO {self.kw['database_table']} (animalName, sireName, damName, sex)
                                VALUES (?, ?, ?, ?)
                                """
                                cursor.execute(sql, (an, si, da, p.sex))
                                conn.commit()
                            _table_loaded = True
                        except Exception as e:
                            logger.error('[savedb]: Error loading pedigree data into table. Error: %s', e)

                    conn.close()

            else:
                # Only SQLite is supported in Python 3 version
                logger.error('[savedb]: Only SQLite database is supported. database_compatibility must be "sqlite".')
                if self.kw['messages'] == 'verbose':
                    print('[ERROR]: Only SQLite database is supported. Set database_compatibility to "sqlite".')

            if _table_loaded:
                if self.kw['messages'] == 'verbose':
                    print(f"[INFO]: Saved pedigree to {self.kw['database_name']}.{self.kw['database_table']} at "
                        f"{pyp_utils.pyp_nice_time()}.")
                logger.info('Saved pedigree to %s.%s at %s.', self.kw['database_name'],
                            self.kw['database_table'], pyp_utils.pyp_nice_time())
                _retval = True
            else:
                if self.kw['messages'] == 'verbose':
                    print(f"[ERROR]: Could not save pedigree to {self.kw['database_name']}.{self.kw['database_table']} "
                        f"at {pyp_utils.pyp_nice_time()}.")
                logger.error('Could not save pedigree to %s.%s at %s.', self.kw['database_name'],
                            self.kw['database_table'], pyp_utils.pyp_nice_time())

        except Exception as e:
            logger.error('[savedb]: An unexpected error occurred. Error: %s', e)
            _retval = False

        return _retval


    def preprocess(self, textstream: str = '', dbstream = '') -> bool:
        """Read animal records into ``self.pedigree``.

        This method owns pedformat interpretation, record iteration,
        NewAnimal construction, and implicit-parent materialization.
        ``load()`` still owns later lifecycle stages: reorder, SNP,
        renumbering, metadata, and chronology.

        :param textstream: Optional string containing animal records.
        :param dbstream: Optional list of tuples of animal records.
        :return: True on success, False on failure.
        :rtype: bool
        """
        _retval = False  # Initialize return value for the procedure
        line_counter = 0  # Count the number of lines in the pedigree file
        animal_counter = 0  # Count the number of animal records in the pedigree file
        critical_count = 0  # Number of critical errors encountered
        pedformat_locations = {}  # Stores columns numbers for input data

        # We need to track the sires and dams read from the pedigree
        # file in order to insert records for any parents that do not
        # have their own records in the pedigree file.
        # A variable, 'pedformat, is passed as a parameter that indicates the format of the
        # pedigree in the input file.  Note that A PEDIGREE FORMAT STRING IS NO LONGER
        # REQUIRED in the input file, and any found will be ignored.  The index of the single-
        # digit code in the format string indicates the column in which the corresponding
        # variable is found.  Duplicate values in the pedformat atring are ignored.
        _sires = {}  # Track sires in the pedigree
        _dams = {}  # Track dams in the pedigree

        try:
            if not self.kw.get('pedfile'):
                raise ValueError("pedfile is not set in self.kw")
            logger.info(f"Processing file: {self.kw['pedfile']}")
        except Exception as e:
            logger.error(f"Error during preprocessing: {e}")
            raise

        try:
            # Ensure pedformat is defined
            if not self.kw.get('pedformat'):
                self.kw['pedformat'] = 'asd'
                logger.error('Null pedigree format string assigned a default value of %s.', self.kw['pedformat'])
                if self.kw['messages'] == 'verbose':
                    print(f"[ERROR]: Null pedigree format string assigned a default value of {self.kw['pedformat']}.")
                raise ValueError("pedformat is not set in self.kw")

            source = None

            # pedformat interpretation is a one-shot map. load() still owns
            # reorder, SNP, renumbering, metadata, and chronology.
            _pedformat, _format_events = _pyp_parse.canonicalize_pedformat(
                self.kw['pedformat'], self.pedformat_codes
            )
            for _kind, _char in _format_events:
                if _kind == "Z":
                    if self.kw['messages'] == 'verbose':
                        print("[INFO]: Skipping one or more columns in the input file.")
                    logger.info(
                        'Skipping one or more columns in the input file as requested '
                        'by the pedigree format string %s',
                        self.kw['pedformat'],
                    )
                else:
                    if self.kw['messages'] == 'verbose':
                        print(f"[DEBUG]: Invalid format code, {_char}, encountered!")
                    logger.error(
                        'Invalid column format code %s found while reading pedigree '
                        'format string %s',
                        _char,
                        self.kw['pedformat'],
                    )

            (
                pedformat_locations,
                critical_count,
                _alleles_collision,
                _optional_debug,
            ) = _pyp_parse.build_pedformat_locations(
                _pedformat,
                alleles_sepchar=self.kw['alleles_sepchar'],
                sepchar=self.kw['sepchar'],
                pedformat=self.kw['pedformat'],
            )
            if 'animal' not in pedformat_locations:
                print(
                    f"[CRITICAL]: No animal identification code was specified in the "
                    f"pedigree format string {_pedformat}! This is a critical error "
                    f"and the program will halt."
                )
            if 'sire' not in pedformat_locations:
                print(
                    f"[CRITICAL]: No sire identification code was specified in the "
                    f"pedigree format string {_pedformat}! This is a critical error "
                    f"and the program will halt."
                )
            if 'dam' not in pedformat_locations:
                print(
                    f"[CRITICAL]: No dam identification code was specified in the "
                    f"pedigree format string {_pedformat}! This is a critical error "
                    f"and the program will halt."
                )
            if self.kw['messages'] == 'all':
                for _msg in _optional_debug:
                    print(_msg)
            if _alleles_collision:
                if self.kw['messages'] == 'all':
                    print(
                        "[DEBUG]: The same separating character was specified for both "
                        "columns of input "
                        "(option sepchar) and alleles (option alleles_sepchar) in an "
                        "animal's "
                        "allelotype. The allelotypes will not be used in this pedigree."
                    )
                logger.warning(
                    'The same separating character was specified for both columns of '
                    'input (option sepchar) and alleles (option alleles_sepchar) in '
                    'an animal\'s allelotype. The allelotypes will not be used in '
                    'this pedigree.'
                )
            # Pedformat 'p' stores an input gencoeff. Pattie calculation is
            # outside the 4.0 domain; do not treat 'p' as "compute gen_coeff".
            if self.kw.get('gen_coeff'):
                raise pyp_errors.PyPedalUsageError(
                    "kw['gen_coeff'] calculation is outside the PyPedal 4.0 "
                    "domain. Pedformat 'p' may store a supplied generation "
                    "coefficient; PyPedal does not compute Pattie (1965) "
                    "coefficients. Leave gen_coeff at False."
                )

            # If the pedigree file includes coefficients of inbreeding flag the pedigree
            if 'f' in self.kw['pedformat']:
                self.kw['f_computed'] = True
            if 'G' in self.kw['pedformat']:
                self.kw['g_computed'] = True

            if critical_count > 0:
                # The pedformat is missing an animal, sire or dam code. Each
                # such omission already printed a [CRITICAL] diagnostic above.
                raise pyp_errors.PyPedalConfigurationError(
                    'The pedigree format string %r is missing %d required '
                    'field code(s); see the [CRITICAL] messages above. A '
                    'pedigree cannot be read without animal, sire and dam.'
                    % (self.kw['pedformat'], critical_count))
            else:
                if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                    print(f"[INFO]: Opening pedigree file {self.kw['pedfile']}")
                    logger.info('Opening pedigree file %s', self.kw['pedfile'])

                # File, textstream, and db records are normalized to raw lines.
                # A textstream without a trailing newline still drops the last
                # record. Database tuples are joined with a comma.
                source = _pyp_parse.PedigreeRecordSource(
                    self.kw['pedfile'],
                    textstream=textstream,
                    dbstream=dbstream,
                )

                while True:
                    line = source.readline(line_counter, logger)
                    if line is False:
                        break

                    # Log the raw line
                    logger.debug(f"Raw line read: {repr(line)}")
                    
                    if not line or line.strip() == '':
                        logger.info('Skipping empty or null line. Reached end-of-line in %s after reading %s lines.', self.kw['pedfile'], line_counter)
                        # logger.info('Skipping empty or null line.')
                        # sys.exit(0)
                        # logger.info('Reached end-of-line in %s after reading %s lines.', self.kw['pedfile'], line_counter)
                        break

                    else:
                        logger.debug(f"Processing line: {line.strip()}")
                        # Handling and processing each line
                        line_counter += 1
                        if line_counter <= self.kw['log_ped_lines']:
                            logger.info('Pedigree (line %s): %s', line_counter, line.strip())

                        # Handle header lines
                        if line_counter == 1 and self.kw['has_header']:
                            logger.info('Converted the first line in the input file into a comment because the pedigree file has a header row.')
                            if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                                print(f"[INFO]: Converted the first line in the input file into a comment because the pedigree file has a header row.")
                            line = f"# {line}"

                        # Handle comment lines
                        if line[0] == '#':
                            logger.info('Pedigree comment (line %s): %s', line_counter, line.strip())
                            continue

                        # Handle deprecated pedigree format strings
                        elif line[0] == '%':
                            self.kw['old_pedformat'] = line[1:].strip()  # Store the format string
                            logger.warning('Encountered deprecated pedigree format string (%s) on line %s of the pedigree file.', line.strip(), line_counter)
                            continue

                        # Handle empty or blank lines
                        elif len(line.strip()) == 0:
                            logger.warning('Encountered an empty (blank) record on line %s of the pedigree file.', line_counter)
                            continue
                        else:
                            animal_counter += 1
                            if animal_counter % self.kw['counter'] == 0:
                                logger.info('Records read: %s', animal_counter)

                            # Split the line based on the separator and validate fields
                            lfields = line.strip().split(self.kw['sepchar'])
                            if (len(self.kw['pedformat']) == len(lfields)) or (len(self.kw['pedformat']) == len(lfields) + 1 and 'P' in self.kw['pedformat']):
                                self.namemap = {}
                                self.namebackmap = {}
                                lfields = [field.strip() for field in lfields]
                                if len(lfields) < 3:
                                    error_msg = (f"The record on line {line_counter} of file {self.kw['pedfile']} is too short - all records "
                                                f"must contain animal, sire, and dam ID numbers ({len(lfields)} fields detected).")
                                    print(f"[ERROR]: {error_msg}")
                                    print(f"[ERROR]: {line}")
                                    raise pyp_errors.PyPedalPedigreeFormatError(
                                        error_msg)
                                else:
                                    # Process and validate animal records
                                    if lfields[0] != self.kw['missing_parent']:
                                        if self.kw['animal_type'] == 'light':
                                            an = LightAnimal(pedformat_locations, lfields, self.kw) 
                                        else:
                                             an = NewAnimal(pedformat_locations, lfields, self.kw)
                                    else:
                                        error_msg = (f"The record on line {line_counter} of file {self.kw['pedfile']} has an animal ID that "
                                                     f"is the same as the missing value code specified for the pedigree. This animal is being "
                                                     f"skipped and will not have an entry in the pedigree.")
                                        print(f"[ERROR]: {error_msg}")
                                        logger.error(error_msg)
                                        continue

                                    # Track sires and dams
                                    if 'S' in self.kw['pedformat']:
                                        if an.sireName != self.kw['missing_name']:
                                            _sires[an.sireName] = an.sireName
                                    else:
                                        if str(an.sireID) != str(self.kw['missing_parent']):
                                            _sires[an.sireID] = an.sireID

                                    if 'D' in self.kw['pedformat']:
                                        if an.damName != self.kw['missing_name']:
                                            _dams[an.damName] = an.damName
                                    else:
                                        if str(an.damID) != str(self.kw['missing_parent']):
                                            _dams[an.damID] = an.damID

                                    self.pedigree.append(an)
                                    # Map IDs and names
                                    # If strings are used for animals names we need to put those names in idmap
                                    # and backmap so that the missing sire and dam assignment code will
                                    # work correctly.
                                    #
                                    # !!! Note that this is broken if you renumber the pedigree. In that case, use the
                                    # namemap to map from name to original ID, and the backmap to go from the original
                                    # ID to the renumbered ID.
                                    if 'A' in self.kw['pedformat']:
                                        self.idmap[an.name] = an.name
                                        self.backmap[an.name] = an.name
                                        if self.kw['animal_type'] == 'new':
                                            self.namemap[an.name] = an.name
                                            self.namebackmap[an.name] = an.name
                                    else:
                                        self.idmap[an.animalID] = an.animalID
                                        self.backmap[an.animalID] = an.animalID
                                        if self.kw['animal_type'] == 'new':
                                            self.namemap[an.name] = an.animalID
                                            self.namebackmap[an.animalID] = an.name
                            else:
                                error_msg = (
                                    f"The record on line {line_counter} of file {self.kw['pedfile']} has {len(lfields)} columns, "
                                    f"but the pedigree format string ({self.kw['pedformat']}) says that it should have "
                                    f"{len(self.kw['pedformat'])} columns "
                                    f"(separator {self.kw['sepchar']!r}). Please check your pedigree file, "
                                    f"the pedigree format string, and the column separator."
                                )
                                print(f"[ERROR]: {error_msg}")
                                # kw['debug'] has no default anywhere; this is a
                                # typo for debug_messages and raised KeyError
                                # while reporting a column-count mismatch, so
                                # the user got an unrelated KeyError instead of
                                # the diagnostic above.
                                if self.kw.get('debug_messages'):
                                    print(f"[DEBUG]: {self.kw['pedformat']}")
                                    print(f"[DEBUG]: {lfields}")
                                logger.error(error_msg)
                                raise pyp_errors.PyPedalPedigreeFormatError(
                                    error_msg)

                # This is where we deal with parents with no pedigree file entry.
                # Things are kind of tricky when we are working with the S and D codes.
                _null_locations = _pyp_parse.implicit_parent_locations(
                    pedformat_locations
                )

                # Source tokens already given a record by the loops below.
                #
                # The `self.idmap[...]` guards cannot do this on their own: the
                # sire loop registers the animal it creates under
                # `an.animalID`, which for an 'A'/'S'/'D' pedigree is the HASHED
                # integer, while both loops probe with `_s`/`_d`, the raw source
                # token. The key domains differ, so a parent that is both a sire
                # and a dam -- the usual case for a shared missing-parent
                # placeholder -- missed the guard in the dam loop and got a
                # second record with the same animalID. A string-ID pedigree
                # can load with more records than distinct IDs because of it.
                #
                # Deliberately narrow: this prevents a duplicate the loader
                # itself creates. It is not general record merging, and it does
                # not make duplicate animal IDs acceptable -- a genuine
                # duplicate still refuses in pyp_utils/_order_pedigree().
                for _role, _token in _pyp_parse.iter_implicit_parent_tokens(
                    _sires,
                    _dams,
                    self.idmap,
                    self.kw['pedformat'],
                    self.kw['missing_parent'],
                    self.kw['missing_name'],
                ):
                    if self.kw['messages'] == 'verbose':
                        print(f'[NOTE]: Adding {_role} {_token} to the pedigree')
                    an = NewAnimal(
                        _null_locations,
                        [_token, self.kw['missing_parent'], self.kw['missing_parent']],
                        self.kw
                    )
                    self.pedigree.append(an)
                    self.idmap[an.animalID] = an.animalID
                    self.backmap[an.animalID] = an.animalID
                    self.namemap[an.name] = an.animalID
                    self.namebackmap[an.animalID] = an.name
                    logger.info(f'Added pedigree entry for {_role} {_token}')
                    if self.kw['messages'] == 'verbose':
                        print(f'[NOTE]: Added pedigree entry for {_role} {_token}')
                    self._implicit_parents.append(an.animalID)

            # Finish up
            #
            logger.info('Closing pedigree file')
            if source is not None:
                source.close()
            _retval = True

        except PyPedalError:
            raise
        except Exception as e:
            # This handler used to set
            # _retval = False and return, which produced a pedigree object
            # containing ZERO animals and no indication that anything had gone
            # wrong -- and every caller of preprocess() ignores the return value
            # anyway, so the flag was never even consulted. A caller who asks
            # for a pedigree and receives an empty one silently is worse off
            # than one who receives an exception.
            logger.error("Error during preprocessing: %s", e)
            raise PyPedalError(
                "Failed to read pedigree %r after %d record(s): %s: %s. The "
                "pedigree has NOT been loaded. Check the pedigree format string "
                "(kw['pedformat']) against the file's columns."
                % (self.kw.get('pedfile'), animal_counter,
                   type(e).__name__, e)) from e

        return _retval


    def fromgraph(self, pedgraph):
        """
        Loads the animals to populate the pedigree from a DiGraph object.

        Parameters
        ----------
        pedgraph : networkx.DiGraph
            A directed graph containing the pedigree.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        _retval = False
        missing = ['sex', 'generation', 'gencoeff', 'birthyear', 'inbreeding', 'breed', 'name',
                'birthdate', 'alive', 'age', 'alleles', 'herd', 'userfield']
        pedformat_locations = {
            'animal': 0,
            'sire': 1,
            'dam': 2,
        }

        try:
            for _m in missing:
                pedformat_locations[_m] = -999

            for _n in pedgraph.nodes():
                _s = pedgraph.nodes[_n].get('sire', self.kw['missing_parent'])
                _d = pedgraph.nodes[_n].get('dam', self.kw['missing_parent'])

                an = NewAnimal(pedformat_locations, [_n, _s, _d], self.kw)
                self.pedigree.append(an)
                self.idmap[an.animalID] = an.animalID
                self.backmap[an.animalID] = an.animalID
                self.namemap[an.name] = an.animalID
                self.namebackmap[an.animalID] = an.name

                if self.kw.get('debug_messages', False):
                    logger.info('Added pedigree entry for animal %s', _n)

            if self.kw.get('debug_messages', False):
                logger.info('Added %s animals to the pedigree in NewPedigree.fromgraph()', len(pedgraph.nodes))

            if self.kw.get('messages', False) == 'verbose':
                print(f"[INFO]: Added {len(pedgraph.nodes)} animals to the pedigree in NewPedigree.fromgraph()!")

            _retval = True

        except Exception as e:
            if self.kw.get('debug_messages', False):
                logger.error('Unable to add animals to the pedigree in NewPedigree.fromgraph(). Error: %s', e)
            if self.kw.get('messages', False) == 'verbose':
                print(f"[ERROR]: Unable to add animals to the pedigree in NewPedigree.fromgraph()! Error: {e}")
            _retval = False

        return _retval

    def fromnull(self):
        """
        Creates a new pedigree with no animal records in it.

        Returns
        -------
        bool
            Always True, indicating success.
        """
        # Log the creation of an empty pedigree
        logger.info('Created a null (empty) pedigree.')
        return True

    def fromanimallist(self, animallist):
        """
        Populates a NewPedigree with instances of NewAnimal objects.

        Parameters
        ----------
        animallist : list
            A list of instances of NewAnimal.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        _retval = False

        if len(animallist) > 0:
            # Guess the pedigree format from the first animal in the list
            self.kw['pedformat'] = pyp_utils.guess_pedformat(animallist[0], self.kw)

            for an in animallist:
                if an.__class__.__name__ == 'NewAnimal':
                    self.pedigree.append(an)
                    self.idmap[an.animalID] = an.animalID
                    self.backmap[an.animalID] = an.animalID
                    self.namemap[an.name] = an.animalID
                    self.namebackmap[an.animalID] = an.name
                    if self.kw['debug_messages']:
                        logger.info('Added pedigree entry for animal %s', an.originalID)
                else:
                    logger.error('An entry in the animallist was not a NewAnimal object, skipping!')

            _retval = True
        else:
            if self.kw['messages']:
                print('[ERROR]: Could not create a pedigree from an empty animal list!')
            logger.error('Could not create a pedigree from an empty animal list!')
            _retval = False

        return _retval

    def tostream(self):
        """
        Creates a text stream from a pedigree.

        Returns
        -------
        str
            A text stream representing the pedigree.
        """
        streamout = ''
        try:
            for p in self.pedigree:
                an = p.name
                si = p.sireName
                da = p.damName

                # Replace missing names with the missing parent indicator
                if si == self.kw['missing_name']:
                    si = self.kw['missing_parent']
                if da == self.kw['missing_name']:
                    da = self.kw['missing_parent']

                # Append the current animal record to the output stream
                streamout += f"{an},{si},{da}\n"

            if self.kw['messages'] == 'verbose':
                print('[INFO]: Created text stream from pedigree.')
            logger.info('Created text stream from pedigree.')
        except Exception as e:
            if self.kw['messages'] == 'verbose':
                print(f'[ERROR]: Could not create text stream from pedigree! Exception: {e}')
            logger.error('Could not create text stream from pedigree!', exc_info=True)

        return streamout

    def renumber(self):
        """
        Updates the ID map after a pedigree has been renumbered so that all references
        are to renumbered rather than original IDs.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        try:
            if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                print(f'\t[INFO]: Renumbering pedigree at {pyp_utils.pyp_nice_time()}')
                print(f'\t\t[INFO]: Reordering pedigree at {pyp_utils.pyp_nice_time()}')
            logger.info('Reordering pedigree')

            # Determine whether to use fast or slow reorder based on format and configuration
            if ('b' in self.kw['pedformat'] or 'y' in self.kw['pedformat']) and not self.kw['slow_reorder']:
                self.pedigree = pyp_utils.fast_reorder(
                    self.pedigree,
                    missingparent=self.kw['missing_parent'],
                )
            else:
                self.pedigree = pyp_utils.reorder(
                    self.pedigree,
                    missingparent=self.kw['missing_parent'],
                    max_rounds=self.kw['reorder_max_rounds']
                )

            if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                print(f'\t\t[INFO]: Renumbering at {pyp_utils.pyp_nice_time()}')
            logger.info('Renumbering pedigree')

            # Renumber the pedigree
            self.pedigree = pyp_utils.renumber(
                self.pedigree,
                missingparent=self.kw['missing_parent'],
                animaltype=self.kw.get('animal_type', 'new')  # Provide a fallback value
            )

            if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                print(f'\t\t[INFO]: Updating ID map at {pyp_utils.pyp_nice_time()}')
            logger.info('Updating ID map')

            # Update the ID map
            self.updateidmap()

            if self.kw['messages'] == 'verbose' and self.kw['pedigree_summary']:
                print(f'\t[INFO]: Assigning offspring at {pyp_utils.pyp_nice_time()}')
            logger.info('Assigning offspring')

            # Assign offspring
            pyp_utils.set_offspring(self)

            # Update flags
            self.kw['pedigree_is_renumbered'] = True
            self.kw['assign_offspring'] = 1
            if getattr(self, 'nrm', None) is not None:
                self._invalidate_relationship_cache()
            return True
        except pyp_errors.PyPedalError:
            # Structural refusals from reorder()/renumber() must reach the
            # caller. Collapsing them into `return False` hid structural
            # refusals: every call site here discards the return value,
            # so a pedigree left half-renumbered by a failure was indistinguish-
            # able from one that renumbered cleanly, and the analysis layer went
            # on to compute coefficients from it.
            logger.error('Renumbering refused; the pedigree was not modified '
                          'past the point of failure', exc_info=True)
            raise
        except Exception as e:
            logger.error('Error during renumbering', exc_info=True)
            if self.kw['messages'] == 'verbose':
                print(f'[ERROR]: An error occurred during renumbering: {e}')
            return False

    def _require_renumbered_pedigree(self, operation):
        """Refuse mutation APIs that assume 1-based current IDs and idmap keys."""
        if not self.kw.get('pedigree_is_renumbered'):
            raise pyp_errors.PyPedalUsageError(
                f"{operation} requires a supported renumbered pedigree "
                "(pedigree_is_renumbered is false). Renumber first; this "
                "call did not mutate."
            )

    def _invalidate_relationship_cache(self):
        """Drop NRM / relationship cache. Never reform a matrix here."""
        self.nrm = None

    def _is_missing_parent(self, parent_id):
        missing = self.kw['missing_parent']
        return parent_id == missing or str(parent_id) == str(missing)

    def _rebuild_offspring_maps(self):
        """Rebuild sons/daus/unks by current animalID. Safe with ID gaps."""
        missing = self.kw['missing_parent']
        by_aid = {animal.animalID: animal for animal in self.pedigree}
        for animal in self.pedigree:
            animal.sons = {}
            animal.daus = {}
            animal.unks = {}
        pedformat = self.kw.get('pedformat', '')
        for animal in self.pedigree:
            parents = []
            if not self._is_missing_parent(animal.sireID):
                parents.append(by_aid.get(animal.sireID))
            if not self._is_missing_parent(animal.damID):
                parents.append(by_aid.get(animal.damID))
            for parent in parents:
                if parent is None:
                    continue
                if 'x' in pedformat and animal.sex in {'m', 'M'}:
                    parent.sons[animal.animalID] = animal.animalID
                elif 'x' in pedformat and animal.sex in {'f', 'F'}:
                    parent.daus[animal.animalID] = animal.animalID
                else:
                    parent.unks[animal.animalID] = animal.animalID

    def _refresh_parent_names(self):
        """Rebuild sireName/damName from surviving parent records."""
        missing_name = self.kw['missing_name']
        by_aid = {animal.animalID: animal for animal in self.pedigree}
        for animal in self.pedigree:
            if self._is_missing_parent(animal.sireID):
                animal.sireName = missing_name
            else:
                sire = by_aid.get(animal.sireID)
                animal.sireName = sire.name if sire is not None else missing_name
            if self._is_missing_parent(animal.damID):
                animal.damName = missing_name
            else:
                dam = by_aid.get(animal.damID)
                animal.damName = dam.name if dam is not None else missing_name

    def _rebuild_metadata(self):
        self.metadata = PedigreeMetadata(self.pedigree, self.kw, self.snp)
        self.metadata.implicit_parent_ids = list(self._implicit_parents)
        self.metadata.num_implicit_parents = len(self._implicit_parents)

    def _rebuild_structural_state(self):
        """Rebuild maps, metadata, and offspring after a structural mutation.

        Runs the flag-gated post-path used after delete, except it never
        forms an NRM from ``kw['form_nrm']``. Always invalidates cache.
        """
        if self.kw.get('reorder', 0) == 1 and not self.kw.get('renumber', 0):
            if not self.kw.get('slow_reorder', False):
                self.pedigree = pyp_utils.fast_reorder(
                    self.pedigree,
                    missingparent=self.kw['missing_parent'],
                )
            else:
                self.pedigree = pyp_utils.reorder(
                    self.pedigree,
                    missingparent=self.kw.get('missing_parent', ''),
                    max_rounds=self.kw.get('reorder_max_rounds', 100),
                )

        if self.kw.get('renumber', False):
            self.renumber()
        else:
            self.updateidmap()
            self._rebuild_offspring_maps()

        self._refresh_parent_names()
        self._rebuild_metadata()

        if self.kw.get('set_generations', False):
            pyp_utils.set_generation(self)

        if self.kw.get('set_ancestors', False):
            pyp_utils.set_ancestor_flag(self)

        if self.kw.get('set_sexes', False) or self.kw.get('assign_sexes', False):
            pyp_utils.set_sexes(self)

        if self.kw.get('set_alleles', False):
            pyp_metrics.effective_founder_genomes(self)

        if self.kw.get('set_offspring', False) and not self.kw.get('renumber', False):
            self._rebuild_offspring_maps()

        if self.kw.get('resolve_duplicates', False):
            self.duplicates = pyp_utils.list_duplicates(self)
        else:
            self.duplicates = False

        if self.duplicates:
            logger.warning(
                "There are duplicate animals in the pedigree file, requiring manual intervention to correct!"
            )
            if self.kw.get('messages', '') != 'quiet':
                print(
                    "[WARNING]: There are duplicate animals in the pedigree file, requiring manual intervention to correct!"
                )

        self._invalidate_relationship_cache()

    def addanimal(self, animalID, sireID, damID):
        """
        Adds a new animal of class NewAnimal to the pedigree.

        Parameters
        ----------
        animalID : int or str
            ID of the new animal to be added to the pedigree.
        sireID : int or str
            Sire ID of the new animal to be added to the pedigree.
        damID : int or str
            Dam ID of the new animal to be added to the pedigree.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        added = False
        appended = None
        prior_caller = self.kw.get('newanimal_caller')

        if not self.kw['pedigree_is_renumbered']:
            logger.warning(
                "Adding an animal to an unrenumbered pedigree using NewPedigree::addanimal() is unsafe!"
            )

        try:
            # Define missing attributes and pedigree format locations
            missing = [
                'sex', 'generation', 'gencoeff', 'birthyear', 'inbreeding', 'breed',
                'name', 'birthdate', 'alive', 'age', 'alleles', 'herd', 'userfield'
            ]
            pedformat_locations = {'animal': 0, 'sire': 1, 'dam': 2}
            for m in missing:
                pedformat_locations[m] = -999

            # Construct the list of animal data
            animal_data = [animalID]
            if sireID == 0:
                animal_data.append(self.kw['missing_parent'])
            elif 'S' in self.kw['pedformat']:
                animal_data.append(self.namebackmap[self.backmap[sireID]])
            else:
                animal_data.append(str(sireID))

            if damID == 0:
                animal_data.append(self.kw['missing_parent'])
            elif 'D' in self.kw['pedformat']:
                animal_data.append(self.namebackmap[self.backmap[damID]])
            else:
                animal_data.append(str(damID))

            # Temporarily override the caller setting for NewAnimal initialization
            self.kw['newanimal_caller'] = 'addanimal'
            animal = NewAnimal(pedformat_locations, animal_data, self.kw)
            self.kw['newanimal_caller'] = 'loader'

            # Add the new animal to the pedigree
            self.pedigree.append(animal)
            appended = animal

            # Renumber the new animal if the pedigree is renumbered
            animal.animalID = max(self.idmap.values()) + 1
            animal.renumberedID = animal.animalID
            if 'n' in self.kw['pedformat']:
                animal.name = self.kw['missing_name']
            if animal.name == animal.originalID:
                animal.name = animal.renumberedID
            animal.sireID = self.idmap[animal.sireID]
            animal.damID = self.idmap[animal.damID]

            # Update ID and name maps
            self.idmap[animal.originalID] = animal.animalID
            self.backmap[animal.animalID] = animal.originalID
            self.namemap[animal.name] = animal.originalID
            self.namebackmap[animal.originalID] = animal.name

            added = True
            self._invalidate_relationship_cache()
        except Exception as e:
            logger.error("Error adding animal to the pedigree: %s", e, exc_info=True)
            added = False
            if appended is not None:
                for i in range(len(self.pedigree) - 1, -1, -1):
                    if self.pedigree[i] is appended:
                        del self.pedigree[i]
                        break
                oid = getattr(appended, 'originalID', None)
                aid = getattr(appended, 'animalID', None)
                name = getattr(appended, 'name', None)
                if oid is not None and self.idmap.get(oid) == aid:
                    del self.idmap[oid]
                if aid is not None and self.backmap.get(aid) == oid:
                    del self.backmap[aid]
                if name is not None and self.namemap.get(name) == oid:
                    del self.namemap[name]
                if oid is not None and self.namebackmap.get(oid) == name:
                    del self.namebackmap[oid]
            if prior_caller is not None:
                self.kw['newanimal_caller'] = prior_caller

        return added

    def delanimal(self, animalID):
        """
        Deletes an animal from the pedigree. Note that this method does not update the
        metadata attached to the pedigree and should only be used if that is not important.
        As of 04/10/2006, delanimal() is intended for use by pyp_metrics/mating_coi()
        rather than directly by users.

        Parameters
        ----------
        animalID : int or str
            ID of the animal to delete.

        Returns
        -------
        bool
            True on success, False on failure.
        """
        deleted = False

        if not self.kw['pedigree_is_renumbered']:
            logger.warning(
                "Deleting an animal from an unrenumbered pedigree using NewPedigree::delanimal() is unsafe!"
            )

        try:
            # Locate the animal's index in the pedigree
            anidx = self.idmap[animalID] - 1

            # Remove animal references from the various maps
            del self.namebackmap[self.pedigree[anidx].originalID]
            del self.namemap[self.pedigree[anidx].name]
            del self.backmap[self.pedigree[anidx].renumberedID]
            del self.idmap[animalID]

            # Remove the animal from the pedigree list
            del self.pedigree[anidx]

            deleted = True
            self._invalidate_relationship_cache()
            logger.info("Successfully deleted animal with ID %s from the pedigree.", animalID)

        except KeyError as e:
            logger.error(
                "Failed to delete animal with ID %s: Animal not found in the pedigree. Error: %s",
                animalID,
                e,
            )
        except Exception as e:
            logger.error(
                "An unexpected error occurred while deleting animal with ID %s: %s",
                animalID,
                e,
                exc_info=True,
            )

        return deleted

    def delete_animals(self, animal_list):
        """
        Atomically delete one or more animals from a renumbered pedigree.

        Request IDs are ``idmap`` keys (``originalID`` after a normal load).
        Duplicate entries have set semantics: ``[A, A]`` deletes ``A`` once.
        An empty list returns ``True`` with no mutation. The call succeeds
        only when every requested ID exists and no survivor retains a known
        ``sireID`` / ``damID`` pointing at a deleted animal. There is no
        cascade, no orphan rewrite, and no automatic NRM rebuild.

        Parameters
        ----------
        animal_list : list
            List of animal IDs to delete.

        Returns
        -------
        bool
            True on success.

        Raises
        ------
        PyPedalUsageError
            Unrenumbered pedigree, or a requested ID is not in ``idmap``.
        PyPedalPedigreeStructureError
            A survivor would retain a known parent in the deletion set.
        """
        self._require_renumbered_pedigree('delete_animals')
        if not animal_list:
            return True

        unique_ids = list(dict.fromkeys(animal_list))
        missing_ids = [animal_id for animal_id in unique_ids if animal_id not in self.idmap]
        if missing_ids:
            raise pyp_errors.PyPedalUsageError(
                "delete_animals refused: requested identit(y/ies) %s do not "
                "exist in idmap. No animals were deleted." % (missing_ids,)
            )

        by_original = {animal.originalID: animal for animal in self.pedigree}
        delete_original_ids = set()
        delete_current_ids = set()
        for animal_id in unique_ids:
            animal = by_original.get(animal_id)
            if animal is None:
                raise pyp_errors.PyPedalUsageError(
                    "delete_animals refused: requested identity %s is in "
                    "idmap but has no pedigree record. No animals were "
                    "deleted." % (animal_id,)
                )
            delete_original_ids.add(animal.originalID)
            delete_current_ids.add(animal.animalID)

        dangling = []
        for animal in self.pedigree:
            if animal.originalID in delete_original_ids:
                continue
            if (
                animal.sireID in delete_current_ids
                or animal.damID in delete_current_ids
            ):
                dangling.append(animal.originalID)
        if dangling:
            raise pyp_errors.PyPedalPedigreeStructureError(
                "delete_animals refused: surviving animal(s) %s still "
                "reference a requested parent. No animals were deleted."
                % (dangling,),
                animals=dangling,
            )

        self.pedigree = [
            animal for animal in self.pedigree
            if animal.originalID not in delete_original_ids
        ]
        self._rebuild_structural_state()
        return True

    _FACTUAL_MERGE_ATTRS = (
        'name', 'displayName', 'sex', 'bd', 'by', 'breed',
        'herd', 'originalHerd', 'alive', 'userField',
    )

    def _is_ancestor(self, ancestor, descendant):
        """True if ancestor is in descendant's sire/dam chain (current IDs)."""
        missing = self.kw['missing_parent']
        by_aid = {animal.animalID: animal for animal in self.pedigree}
        target = ancestor.animalID
        seen = set()
        stack = [descendant.sireID, descendant.damID]
        while stack:
            pid = stack.pop()
            if pid == missing or str(pid) == str(missing) or pid in seen:
                continue
            if pid == target:
                return True
            seen.add(pid)
            parent = by_aid.get(pid)
            if parent is None:
                continue
            stack.append(parent.sireID)
            stack.append(parent.damID)
        return False

    def _factual_missing(self, attr, value):
        sentinels = {
            'name': self.kw.get('missing_name'),
            'displayName': self.kw.get('missing_name'),
            'sex': self.kw.get('missing_sex', 'u'),
            'breed': self.kw.get('missing_breed'),
            'herd': self.kw.get('missing_herd'),
            'originalHerd': self.kw.get('missing_herd'),
            'alive': self.kw.get('missing_alive', 0),
            'userField': self.kw.get('missing_userfield', 'Unknown'),
        }
        if attr in ('bd', 'by'):
            return value is None
        sentinel = sentinels[attr]
        return value == sentinel or str(value) == str(sentinel)

    def _reconcile_merge_value(self, attr, keep_val, drop_val):
        keep_missing = self._factual_missing(attr, keep_val)
        drop_missing = self._factual_missing(attr, drop_val)
        if keep_missing and drop_missing:
            return None
        if keep_missing and not drop_missing:
            return drop_val
        if not keep_missing and drop_missing:
            return None
        if keep_val == drop_val or str(keep_val) == str(drop_val):
            return None
        raise pyp_errors.PyPedalPedigreeStructureError(
            "merge_animals refused: conflicting known %s values (%r vs %r). "
            "The pedigree was not modified." % (attr, keep_val, drop_val)
        )

    def _reconcile_merge_parent(self, keep_pid, drop_pid, slot):
        keep_missing = self._is_missing_parent(keep_pid)
        drop_missing = self._is_missing_parent(drop_pid)
        if keep_missing and drop_missing:
            return None
        if keep_missing and not drop_missing:
            return drop_pid
        if not keep_missing and drop_missing:
            return None
        if keep_pid == drop_pid or str(keep_pid) == str(drop_pid):
            return None
        raise pyp_errors.PyPedalPedigreeStructureError(
            "merge_animals refused: conflicting known %s (%r vs %r). "
            "The pedigree was not modified." % (slot, keep_pid, drop_pid)
        )

    def merge_animals(self, keep, drop):
        """
        Merge two records that represent one real animal.

        Request IDs are ``idmap`` keys. Offspring of ``drop`` are redirected
        to ``keep``. Factual conflicts and ancestor relationships refuse with
        no mutation. Does not auto-form an NRM.

        Parameters
        ----------
        keep : int or str
            Identity of the surviving record.
        drop : int or str
            Identity of the record to merge into ``keep`` and remove.

        Returns
        -------
        bool
            True on success.

        Raises
        ------
        PyPedalUsageError
            Unrenumbered pedigree, missing ID, or ``keep == drop``.
        PyPedalPedigreeStructureError
            Ancestor relationship, conflicting parents, or conflicting facts.
        """
        self._require_renumbered_pedigree('merge_animals')
        if keep == drop:
            raise pyp_errors.PyPedalUsageError(
                "merge_animals refused: keep and drop are the same identity "
                "(%r). The pedigree was not modified." % (keep,)
            )
        missing = [ident for ident in (keep, drop) if ident not in self.idmap]
        if missing:
            raise pyp_errors.PyPedalUsageError(
                "merge_animals refused: requested identit(y/ies) %s do not "
                "exist in idmap. The pedigree was not modified." % (missing,)
            )

        by_original = {animal.originalID: animal for animal in self.pedigree}
        keep_animal = by_original.get(keep)
        drop_animal = by_original.get(drop)
        if keep_animal is None or drop_animal is None:
            raise pyp_errors.PyPedalUsageError(
                "merge_animals refused: keep %r or drop %r is in idmap but "
                "has no pedigree record. The pedigree was not modified."
                % (keep, drop)
            )

        if self._is_ancestor(keep_animal, drop_animal) or self._is_ancestor(
            drop_animal, keep_animal
        ):
            raise pyp_errors.PyPedalPedigreeStructureError(
                "merge_animals refused: keep %r and drop %r are in an "
                "ancestor relationship. The pedigree was not modified."
                % (keep, drop),
                animals=[keep, drop],
            )

        parent_updates = {}
        sire_copy = self._reconcile_merge_parent(
            keep_animal.sireID, drop_animal.sireID, 'sireID')
        if sire_copy is not None:
            parent_updates['sireID'] = sire_copy
        dam_copy = self._reconcile_merge_parent(
            keep_animal.damID, drop_animal.damID, 'damID')
        if dam_copy is not None:
            parent_updates['damID'] = dam_copy

        field_updates = {}
        for attr in self._FACTUAL_MERGE_ATTRS:
            copied = self._reconcile_merge_value(
                attr, getattr(keep_animal, attr), getattr(drop_animal, attr))
            if copied is not None:
                field_updates[attr] = copied

        for attr, value in field_updates.items():
            setattr(keep_animal, attr, value)
        for attr, value in parent_updates.items():
            setattr(keep_animal, attr, value)
        if parent_updates:
            keep_animal.founder = (
                'y'
                if self._is_missing_parent(keep_animal.sireID)
                and self._is_missing_parent(keep_animal.damID)
                else 'n'
            )

        drop_aid = drop_animal.animalID
        keep_aid = keep_animal.animalID
        for animal in self.pedigree:
            if animal is drop_animal:
                continue
            if animal.sireID == drop_aid:
                animal.sireID = keep_aid
            if animal.damID == drop_aid:
                animal.damID = keep_aid

        self.pedigree = [
            animal for animal in self.pedigree if animal is not drop_animal
        ]
        self._rebuild_structural_state()
        return True

    def updateidmap(self):
        """
        Updates the ID map after a pedigree has been renumbered so that
        all references are to renumbered rather than original IDs.

        :param self: NewPedigree object.
        :return: None
        """
        logger.info("Updating ID map...")

        # Initialize or reset the ID maps
        self.idmap = {}
        self.backmap = {}
        self.namemap = {}
        self.namebackmap = {}

        # Iterate through each animal in the pedigree and update mappings
        for animal in self.pedigree:
            try:
                self.idmap[animal.originalID] = animal.animalID
                self.backmap[animal.renumberedID] = animal.originalID

                if self.kw.get('animal_type') == 'new':
                    self.namemap[animal.name] = animal.originalID
                    self.namebackmap[animal.originalID] = animal.name
            except KeyError as e:
                logger.warning(
                    "A KeyError occurred while updating ID map for animal: %s", animal
                )
            except Exception as e:
                logger.error(
                    "An unexpected error occurred while updating ID map: %s", e, exc_info=True
                )

        logger.info("ID map update complete.")

    def printoptions(self):
        """
        printoptions() prints the contents of the options dictionary.
        :param self: NewPedigree object.
        :return: None
        """
        print(f"{self.kw['pedname']} OPTIONS")
        for key, value in self.kw.items():
            if len(key) <= 14:
                print(f"\t{key}:\t\t{value}")
            else:
                print(f"\t{key}:\t{value}")

    def simulate(self):
        """
        Simulate an arbitrary pedigree of size n with g generations
        starting from n_s base sires and n_d base dams. This method is based
        on the concepts and algorithms in the Pedigree::sample method from
        Matvec 1.1a.
        :param self: NewPedigree object.
        :return: None
        """
        import numpy as np
        from random import randint

        if self.kw.get('messages') == 'verbose':
            print('[SIMULATE]: Preparing to simulate a pedigree')
            logger.info('Preparing to simulate a pedigree')

        # Ensure the pedigree is empty before simulation
        if len(self.pedigree) > 0:
            logger.error(
                'The simulate() method did not create a new randomly-generated pedigree because the pedigree %s '
                'has already been populated with animals.',
                self.kw['pedname']
            )
            if self.kw.get('messages') == 'verbose':
                print(f'[ERROR]: The pedigree {self.kw["pedname"]} already has animals.')
            return

        # Validate simulation parameters
        self.kw['simulate_n'] = max(1, int(self.kw.get('simulate_n', 15)))
        self.kw['simulate_g'] = max(1, int(self.kw.get('simulate_g', 3)))
        self.kw['simulate_ns'] = max(1, int(self.kw.get('simulate_ns', 4)))
        self.kw['simulate_nd'] = max(1, int(self.kw.get('simulate_nd', 4)))
        self.kw['simulate_sr'] = max(0.0, min(1.0, float(self.kw.get('simulate_sr', 0.5))))
        self.kw['simulate_ir'] = max(0.0, min(1.0, float(self.kw.get('simulate_ir', 0.5))))
        self.kw['simulate_pmd'] = max(1, int(self.kw.get('simulate_pmd', 100)))

        # Derived parameters
        _snt = self.kw['simulate_n'] + 2
        _sng = self.kw['simulate_g']
        _sns = self.kw['simulate_ns']
        _snd = self.kw['simulate_nd']
        _ssr = self.kw['simulate_sr']
        _sir = self.kw['simulate_ir']
        _spmd = self.kw['simulate_pmd']

        # Initialize variables
        _smd = int(round((_snt - _sns - _snd) * (1 - _ssr)))  # Max daughters
        _sms = int(round((_snt - _sns - _snd) * _ssr))        # Max sons
        _smf = _snd + _smd                                    # Max females
        _smm = _sns + _sms                                    # Max males
        _totalna = _sns + _snd
        _npg = ((_snt - _snd - _sns) // _sng) + 1
        males = list(range(_snd, _snd + _sns))
        females = list(range(_snd))
        _pedholder = [None] * (_snt + 1)

        # Initialize founders
        for i in females:
            _pedholder[i] = SimAnimal(i, 0, 0, 'f', 0)
        for i in males:
            _pedholder[i] = SimAnimal(i, 0, 0, 'm', 0)

        # Generate animals generation by generation
        for g in range(_sng):
            if _totalna >= _snt:
                break
            for _ in range(_npg):
                if _totalna >= _snt:
                    break

                sire = 0
                dam = 0
                if np.random.rand() > _sir:  # Not an immigrant
                    for _ in range(_spmd):  # Retry loop for valid parents
                        sire = np.random.choice(males)
                        dam = np.random.choice(females)

                        sire_ok = sire != dam
                        if sire_ok:
                            break
                    else:
                        sire = 0
                        dam = 0

                _totalna += 1
                if np.random.rand() < _ssr:
                    females.append(_totalna)
                    _pedholder[_totalna] = SimAnimal(_totalna, sire, dam, 'f', g + 1)
                else:
                    males.append(_totalna)
                    _pedholder[_totalna] = SimAnimal(_totalna, sire, dam, 'm', g + 1)

        # Convert to NewAnimal objects and populate pedigree
        pedformat_locations = {
            'animal': 0, 'sire': 1, 'dam': 2, 'sex': 3, 'generation': 4,
            'gencoeff': -999, 'birthyear': -999, 'inbreeding': -999,
            'genomic_inbreeding': -999, 'homozygosity': -999, 'breed': -999,
            'name': -999, 'birthdate': -999, 'alive': -999, 'age': -999,
            'alleles': -999, 'herd': -999, 'userfield': -999
        }

        for entry in _pedholder:
            if entry and entry.animalID != 0:
                lightans = [
                    str(entry.animalID),
                    str(entry.sireID) if entry.sireID else self.kw['missing_parent'],
                    str(entry.damID) if entry.damID else self.kw['missing_parent'],
                    entry.sex,
                    entry.gen,
                ]
                if self.kw['animal_type'] == 'light':
                    animal = LightAnimal(pedformat_locations, lightans, self.kw)
                else:
                    animal = NewAnimal(pedformat_locations, lightans, self.kw)

                self.pedigree.append(animal)
                self.idmap[animal.animalID] = animal.animalID
                self.backmap[animal.animalID] = animal.animalID
                self.namemap[animal.name] = animal.animalID
                self.namebackmap[animal.animalID] = animal.name

        if self.kw.get('pedigree_save'):
            try:
                filename = f"{self.kw['filetag']}.ped"
                with open(filename, 'w') as f:
                    for entry in _pedholder:
                        if entry and entry.animalID != 0:
                            f.write(f"{entry.stringme()}\n")
            except Exception as e:
                logger.error(f"Failed to save pedigree to file: {e}")


class SimAnimal:
    """The SimAnimal class is a placeholder used for simulating animals."""

    def __init__(self, animalID, sireID=0, damID=0, sex='u', gen=0):
        """
        Initialize a SimAnimal object.
        :param animalID: Animal's ID.
        :param sireID: Sire's ID.
        :param damID: Dam's ID.
        :param sex: Sex of the animal.
        :param gen: Generation to which an animal belongs.
        """
        self.animalID = animalID
        self.sireID = sireID
        self.damID = damID
        self.sex = sex
        self.gen = gen

    def printme(self):
        """
        Print a summary of the data stored in the SimAnimal object.
        """
        try:
            print(f'\t{self.animalID}\t{self.sireID}\t{self.damID}\t{self.sex}\t{self.gen}')
        except Exception as e:
            print(f"Error while printing SimAnimal object: {e}")

    def stringme(self):
        """
        Return the data stored in the SimAnimal object as a string.
        """
        try:
            return f'{self.animalID} {self.sireID} {self.damID} {self.sex} {self.gen}'
        except Exception as e:
            print(f"Error while converting SimAnimal object to string: {e}")
            return False


class NewAnimal:
    """A simple class to hold animal records read from a pedigree file."""

    def __init__(self, locations, data, mykw):
        """
        Initialize a NewAnimal object.

        Parameters
        ----------
        locations : dict
            A dictionary mapping attribute names to their positions in the input data.
        data : list
            The input data for the current animal.
        mykw : dict
            A dictionary containing options and default values.
        """
        def safe_strip(value):
            """Safely strip strings, converting non-strings to strings first."""
            return str(value).strip() if isinstance(value, str) else str(value)

        # Required fields: animal, sire, dam
        #
        # The pedformat codes are case-sensitive: lowercase 'a'/'s'/'d' mean the
        # column holds an integer ID, uppercase 'A'/'S'/'D' mean it holds a
        # string that must be hashed to an integer. The uppercase branch was
        # missing in an earlier Python 3 port, so int() was applied
        # unconditionally, every string pedigree raised ValueError, and
        # preprocess() swallowed it and handed back an EMPTY pedigree.
        pedformat = mykw.get('pedformat', '')

        if locations.get('animal', -999) == -999:
            raise ValueError("The 'animal' field is missing or not defined in the pedigree format.")
        if 'A' in pedformat:
            if mykw.get('newanimal_caller') == 'addanimal':
                # Adding an animal programmatically: the caller supplies the ID
                # it wants and it is kept verbatim, as in PyPedal 2.0.4.
                self.animalID = data[locations['animal']]
                self.originalID = data[locations['animal']]
            else:
                self.animalID = hashed_string_id(
                    data[locations['animal']], mykw, 'animal')
                self.originalID = self.animalID
            self.name = safe_strip(data[locations['animal']])
        else:
            self.animalID = int(safe_strip(data[locations['animal']]))
            self.originalID = self.animalID
            self.name = (
                safe_strip(data[locations['name']]) if locations.get('name', -999) != -999
                else safe_strip(data[locations['animal']])
            )

        if locations.get('sire', -999) == -999:
            raise ValueError("The 'sire' field is missing or not defined in the pedigree format.")
        if str(data[locations['sire']]) != str(mykw['missing_parent']):
            if 'S' in pedformat:
                self.sireID = hashed_string_id(data[locations['sire']], mykw, 'sire')
            else:
                self.sireID = int(safe_strip(data[locations['sire']]))
            self.sireName = safe_strip(data[locations['sire']])
        else:
            self.sireID = mykw['missing_parent']
            self.sireName = mykw['missing_name']

        if locations.get('dam', -999) == -999:
            raise ValueError("The 'dam' field is missing or not defined in the pedigree format.")
        if str(data[locations['dam']]) != str(mykw['missing_parent']):
            if 'D' in pedformat:
                self.damID = hashed_string_id(data[locations['dam']], mykw, 'dam')
            else:
                self.damID = int(safe_strip(data[locations['dam']]))
            self.damName = safe_strip(data[locations['dam']])
        else:
            self.damID = mykw['missing_parent']
            self.damName = mykw['missing_name']

        # Optional fields
        self.gen = (
            data[locations['generation']] if 'generation' in locations and locations['generation'] != -999
            else mykw.get('missing_gen', -1)
        )
        self.gencoeff = (
            data[locations['gencoeff']] if 'gencoeff' in locations and locations['gencoeff'] != -999
            else mykw.get('missing_gencoeff', -1.0)
        )
        self.sex = (
            safe_strip(data[locations['sex']]).lower() if 'sex' in locations and locations['sex'] != -999
            else mykw.get('missing_sex', 'unknown')
        )
        has_date = locations.get('birthdate', -999) != -999
        has_year = locations.get('birthyear', -999) != -999
        date_token = data[locations['birthdate']] if has_date else None
        year_token = data[locations['birthyear']] if has_year else None
        parsed_bd, parsed_by = pyp_chronology.parse_animal_chronology(
            year_token,
            date_token,
            has_year_column=has_year,
            has_date_column=has_date,
            legacy_missing_byear_token=mykw.get('legacy_missing_byear_token'),
            legacy_missing_bdate_token=mykw.get('legacy_missing_bdate_token'),
            animal_id=self.animalID,
        )
        self.bd = parsed_bd
        self.by = parsed_by
        self.birth_date_estimate = pyp_chronology.BirthDateEstimate()

        if mykw.get('debug_messages'):
            print(f"[DEBUG]: Animal {self.animalID} initialized with birth year {self.by}")

        self.fa = (
            float(safe_strip(data[locations['inbreeding']]))
            if 'inbreeding' in locations and locations['inbreeding'] != -999
            else mykw.get('missing_inbreeding', 0.0)
        )
        self.fg = (
            float(safe_strip(data[locations['genomic_inbreeding']]))
            if 'genomic_inbreeding' in locations and locations['genomic_inbreeding'] != -999
            else mykw.get('missing_genomic_inbreeding', 0.0)
        )
        self.homozygosity = (
            float(data[locations['homozygosity']])
            if locations.get('homozygosity', -999) != -999
            else mykw.get('missing_homozygosity', 0.0)
        )
        self.displayName = (
            safe_strip(data[locations['name']]) if 'name' in locations and locations['name'] != -999
            else mykw.get('missing_name', str(self.animalID))
        )
        self.breed = (
            safe_strip(data[locations['breed']]) if 'breed' in locations and locations['breed'] != -999
            else mykw.get('missing_breed', 'unknown')
        )
        self.age = (
            int(safe_strip(data[locations['age']]))
            if 'age' in locations and locations['age'] != -999
            else mykw.get('missing_age', -1)
        )
        self.alive = (
            int(safe_strip(data[locations['alive']]))
            if 'alive' in locations and locations['alive'] != -999
            else mykw.get('missing_alive', -1)
        )
        self.herd = (
            int(safe_strip(data[locations['herd']])) if 'herd' in locations and locations['herd'] != -999
            else mykw.get('missing_herd', -1)
        )
        self.originalHerd = (
            safe_strip(data[locations['herd']]) if 'herd' in locations and locations['herd'] != -999
            else mykw.get('missing_herd', 'unknown')
        )
        self.renumberedID = -999
        self.igen = mykw.get('missing_igen', -1)

        # Founder determination
        self.founder = (
            'y' if str(self.sireID) == str(mykw['missing_parent']) and str(self.damID) == str(mykw['missing_parent'])
            else 'n'
        )
        self.paddedID = self.pad_id()

        # Allele assignment
        if 'alleles' in locations and locations['alleles'] != -999:
            self.alleles = [
                safe_strip(data[locations['alleles']]).split(mykw['alleles_sepchar'])[0],
                safe_strip(data[locations['alleles']]).split(mykw['alleles_sepchar'])[1],
            ]
        else:
            if self.founder == 'y':
                self.alleles = [f"{self.paddedID}__1", f"{self.paddedID}__2"]
            elif self.sireID == mykw['missing_parent']:
                self.alleles = [f"{self.paddedID}__1", ""]
            elif self.damID == mykw['missing_parent']:
                self.alleles = ["", f"{self.paddedID}__2"]
            else:
                self.alleles = mykw.get('missing_alleles', ['unknown', 'unknown'])

        # Final attributes
        self.pedcomp = mykw.get('missing_pedcomp', 'unknown')
        # pyp_utils.set_ancestor_flag() only ever writes 1, and only onto animals
        # that are a parent, so without this initialiser every non-parent had no
        # 'ancestor' attribute at all. PyPedal 2.0.4 sets it here,
        # unconditionally, at pyp_newclasses.py:3206.
        self.ancestor = 0
        # The 'u' pedformat code was accepted, its column position computed and
        # stored, and its default declared -- but never used to set the attribute,
        # while ten direct and five indirect consumers read it.
        # PyPedal 2.0.4 assigns it here, at pyp_newclasses.py:3234-3238.
        self.userField = (
            safe_strip(data[locations['userfield']])
            if 'userfield' in locations and locations['userfield'] != -999
            else mykw.get('missing_userfield', 'Unknown')
        )
        self.sons = {}
        self.daus = {}
        self.unks = {}

    def pad_id(self):
        """
        Generate a padded ID for the animal.

        :return: A padded string representing the animal's ID.
        """
        return pyp_chronology.padded_identity(self.animalID, self.by)

    def string_to_int(self, idstring, mymaxint=9223372036854775807):
        """
        Convert a string ID to an integer using a hash function.

        PyPedal 2.0.4 provided this method on NewAnimal as well as on
        LightAnimal; an earlier Python 3 port kept it only on LightAnimal.
        Restored, delegating to the
        shared module-level implementation.
        """
        return string_to_int(idstring, mymaxint)


class LightAnimal:
    """
    The LightAnimal class holds animal records read from a pedigree file. It
    is a much simpler object than the NewAnimal class and is intended for use
    with the graph-theoretic routines in pyp_network. The only attributes of these
    objects are: animal ID, sire ID, dam ID, original ID, birth year, and sex.
    """

    def __init__(self, locations, data, mykw):
        """
        Initialize a LightAnimal object.
        :param locations: A dictionary containing the locations of variables in the input line.
        :param data: The line of input read from the pedigree file.
        :param mykw: A dictionary of keyword arguments.
        """
        if locations['animal'] != -999:
            if 'A' in mykw['pedformat']:
                self.animalID = self.string_to_int(data[locations['animal']])
                self.originalID = self.string_to_int(data[locations['animal']])
            else:
                self.animalID = safe_strip(data[locations['animal']])
                self.originalID = safe_strip(data[locations['animal']])

        if locations['sire'] != -999 and safe_strip(data[locations['sire']]) != mykw['missing_parent']:
            if 'S' in mykw['pedformat']:
                self.sireID = (
                    self.string_to_int(data[locations['sire']])
                    if data[locations['sire']] != mykw['missing_parent']
                    else data[locations['sire']]
                )
            else:
                self.sireID = safe_strip(data[locations['sire']])
        else:
            self.sireID = mykw['missing_parent']

        if locations['dam'] != -999 and safe_strip(data[locations['dam']]) != mykw['missing_parent']:
            if 'D' in mykw['pedformat']:
                self.damID = (
                    self.string_to_int(data[locations['dam']])
                    if data[locations['dam']] != mykw['missing_parent']
                    else data[locations['dam']]
                )
            else:
                self.damID = safe_strip(data[locations['dam']])
        else:
            self.damID = mykw['missing_parent']

        self.sex = safe_strip(data[locations['sex']]) if locations['sex'] != -999 else 'u'

        has_year = locations.get('birthyear', -999) != -999
        has_date = locations.get('birthdate', -999) != -999
        year_token = data[locations['birthyear']] if has_year else None
        date_token = data[locations['birthdate']] if has_date else None
        _bd, self.by = pyp_chronology.parse_animal_chronology(
            year_token,
            date_token,
            has_year_column=has_year,
            has_date_column=has_date,
            legacy_missing_byear_token=mykw.get('legacy_missing_byear_token'),
            legacy_missing_bdate_token=mykw.get('legacy_missing_bdate_token'),
            animal_id=getattr(self, 'animalID', None),
        )

        self.paddedID = pad_id(self.by, self.animalID)

    def printme(self):
        """
        Print the contents of an animal record - used for debugging.
        """
        try:
            print(f'ANIMAL {self.animalID} RECORD')
            print(f'\tAnimal ID:\t{self.animalID}')
            print(f'\tSire ID:\t{self.sireID}')
            print(f'\tDam ID:\t\t{self.damID}')
            print(f'\tBirth Year:\t{self.by}')
            print(f'\tSex:\t\t{self.sex}')
            print(f'\tOriginal ID:\t{self.originalID}')
        except Exception as e:
            print(f"Error printing LightAnimal object: {e}")

    def stringme(self):
        """
        Return the contents of an animal record as a string.
        """
        try:
            return (
                f'ANIMAL {self.animalID} RECORD\n'
                f'\tAnimal ID:\t{self.animalID}\n'
                f'\tSire ID:\t{self.sireID}\n'
                f'\tDam ID:\t\t{self.damID}\n'
                f'\tBirth Year:\t{self.by}\n'
                f'\tSex:\t\t{self.sex}\n'
                f'\tOriginal ID:\t{self.originalID}\n'
            )
        except Exception as e:
            print(f"Error converting LightAnimal object to string: {e}")
            return ""

    def dictme(self):
        """
        Return the contents of an animal record as a dictionary.
        """
        try:
            return {
                'animalID': self.animalID,
                'sireID': self.sireID,
                'damID': self.damID,
                'by': self.by,
                'sex': self.sex,
                'originalID': self.originalID,
            }
        except Exception as e:
            print(f"Error converting LightAnimal object to dictionary: {e}")
            return {}

    def trap(self):
        """
        Trap common errors in pedigree file entries.
        """
        try:
            if int(self.animalID) == int(self.sireID):
                print(f"[ERROR]: Animal {self.animalID} has an ID equal to its sire's ID (sire ID {self.sireID}).")
            if int(self.animalID) == int(self.damID):
                print(f"[ERROR]: Animal {self.animalID} has an ID equal to its dam's ID (dam ID {self.damID}).")
            if int(self.animalID) < int(self.sireID):
                print(f"[ERROR]: Animal {self.animalID} is older than its sire (sire ID {self.sireID}).")
            if int(self.animalID) < int(self.damID):
                print(f"[ERROR]: Animal {self.animalID} is older than its dam (dam ID {self.damID}).")
        except ValueError as e:
            print(f"[ERROR]: Invalid ID data encountered: {e}")


    def string_to_int(self, idstring, mymaxint=9223372036854775807):
        """
        Convert a string ID to an integer using a hash function.

        Retained as a method because PyPedal 2.0.4 exposed it as one and
        callers may use it; it now delegates to the module-level
        ``string_to_int`` so that NewAnimal and LightAnimal cannot drift apart.

        The previous body caught every exception and returned 0 on failure,
        which mapped every unhashable ID to the same animal. Hashing failures
        now propagate.
        """
        return string_to_int(idstring, mymaxint)


class PedigreeMetadata:
    """
    A class to hold metadata about pedigrees. Hopefully this will help improve performance
    in some procedures, as well as provide some useful summary data.
    """

    def __init__(self, myped, kw, snp):
        """
        Initialize a PedigreeMetadata object.
        """
        self.kw = kw
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Instantiating a new PedigreeMetadata() object...')
            print('\t[INFO]: Naming the Pedigree()...')
        self.name = kw.get('pedname', '')
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Assigning a filename...')
        self.filename = kw.get('pedfile', '')
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Attaching a pedigree...')
        self.myped = myped
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Setting the pedcode...')
        self.pedcode = kw.get('pedformat', '')
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting the number of animals in the pedigree...')
        self.num_records = len(self.myped)
        # Animals PyPedal created because they were named as a parent but had no
        # record of their own. Without this, num_records
        # exceeding the input row count looks like a parsing error and there is
        # no programmatic way to reconcile the two. Set by NewPedigree after
        # construction; defaults to 0 for metadata built from other sources.
        self.num_implicit_parents = 0
        self.implicit_parent_ids = []
        self._calculate_metadata(snp)

    def _calculate_metadata(self, snp):
        """
        Helper method to calculate all metadata attributes.
        """
        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique sires...')
        self.num_unique_sires, self.unique_sire_list = self.nus()

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique dams...')
        self.num_unique_dams, self.unique_dam_list = self.nud()

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Setting renumbered flag...')
        self.renumbered = self.kw.get('renumber', False)

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique generations...')
        self.num_unique_gens, self.unique_gen_list = self.nug()

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique birthyears...')
        self.num_unique_years, self.unique_year_list = self.nuy()

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique founders...')
        self.num_unique_founders, self.unique_founder_list = self.nuf()

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Counting and finding unique herds...')
        self.num_unique_herds, self.unique_herd_list = self.nuherds()

        if 'u' in self.kw.get('pedformat', ''):
            if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
                print('\t[INFO]: Counting and finding unique userField items...')
            self.num_unique_fields, self.unique_field_list = self.nufields()

        if 'P' in self.pedcode and snp is not None:
            if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
                print('\t[INFO]: Counting SNP...')
            self.snp_count = snp['n_snps'].iloc[0]
        else:
            self.snp_count = None

        if self.kw.get('messages') == 'verbose' and self.kw.get('pedigree_summary'):
            print('\t[INFO]: Detaching pedigree...')
        self.myped = []

    def printme(self):
        """
        Print the pedigree metadata.
        """
        print(f'Metadata for {self.name} ({self.filename})')
        print(f'\tRecords:\t\t{self.num_records}')
        print(f'\tUnique Sires:\t\t{self.num_unique_sires}')
        print(f'\tUnique Dams:\t\t{self.num_unique_dams}')
        print(f'\tUnique Gens:\t\t{self.num_unique_gens}')
        print(f'\tUnique Years:\t\t{self.num_unique_years}')
        print(f'\tUnique Founders:\t{self.num_unique_founders}')
        print(f'\tUnique Herds:\t\t{self.num_unique_herds}')
        if 'u' in self.pedcode:
            print(f'\tUnique userFields:\t{self.num_unique_fields}')
        if 'P' in self.pedcode and self.snp_count is not None:
            print(f'\tSNP Count:\t\t{self.snp_count}')
        print(f'\tPedigree Code:\t\t{self.pedcode}')

    def stringme(self):
        """Return pedigree metadata as a string."""
        lines = [
            f'Metadata for {self.name} ({self.filename})',
            f'\tRecords:\t\t{self.num_records}',
            f'\tUnique Sires:\t\t{self.num_unique_sires}',
            f'\tUnique Dams:\t\t{self.num_unique_dams}',
            f'\tUnique Gens:\t\t{self.num_unique_gens}',
            f'\tUnique Years:\t\t{self.num_unique_years}',
            f'\tUnique Founders:\t{self.num_unique_founders}',
            f'\tUnique Herds:\t\t{self.num_unique_herds}',
        ]
        if 'u' in self.pedcode:
            lines.append(f'\tUnique userFields:\t{self.num_unique_fields}')
        if 'P' in self.pedcode and self.snp_count is not None:
            lines.append(f'\tSNP Count:\t\t{self.snp_count}')
        lines.append(f'\tPedigree Code:\t\t{self.pedcode}')
        return '\n'.join(lines)

    def nus(self):
        """
        Count unique sires.
        """
        sirelist = {x.sireID for x in self.myped if x.sireID != self.kw.get('missing_parent')}
        return len(sirelist), sirelist

    def nud(self):
        """
        Count unique dams.
        """
        damlist = {x.damID for x in self.myped if x.damID != self.kw.get('missing_parent')}
        return len(damlist), damlist

    def nug(self):
        """
        Count unique generations.
        """
        genlist = {x.gen for x in self.myped} if self.kw.get('animal_type') != 'light' else []
        return len(genlist), genlist

    def nuy(self):
        """
        Count unique recorded birth years. Unknown chronology is omitted.
        """
        known = [x.by for x in self.myped if x.by is not None]
        bylist = set(known)
        self.num_unknown_birth_years = sum(1 for x in self.myped if x.by is None)
        return len(bylist), bylist

    def nuf(self):
        """
        Count unique founders.
        """
        if self.kw.get('animal_type') == 'light':
            flist = {x.animalID for x in self.myped if
                     x.sireID == self.kw.get('missing_parent') and x.damID == self.kw.get('missing_parent')}
        else:
            flist = {x.animalID for x in self.myped if x.founder == 'y'}
        return len(flist), flist

    def nuherds(self):
        """
        Count unique herds.
        """
        herdlist = {x.originalHerd for x in self.myped} if self.kw.get('animal_type') != 'light' else []
        return len(herdlist), herdlist

    def nufields(self):
        """
        Count unique user fields.
        """
        fieldlist = {x.userField for x in self.myped} if self.kw.get('animal_type') != 'light' else []
        return len(fieldlist), fieldlist


class NewAMatrix:
    """
    NewAMatrix provides an instance of a numerator relationship matrix (NRM) as a NumPy array
    of floats with some convenience methods. The idea here is to provide a wrapper around a NRM
    so that it is easier to work with. For large pedigrees, it can take a long time to compute
    the elements of A, so there is real value in providing an easy way to save and retrieve a
    NRM once it has been formed.
    """

    def __init__(self, kw):
        """
        Initialize a new numerator relationship matrix.
        """
        self.kw = kw
        self.kw.setdefault('messages', 'verbose')
        self.kw.setdefault('nrm_method', 'nrm')
        self.kw.setdefault('nrm_format', 'text')
        self.nrm = None

    def form_a_matrix(self, pedigree):
        """
        form_a_matrix() calls pyp_nrm/fast_a_matrix() or pyp_nrm/fast_a_matrix_r()
        to form a NRM from a pedigree.
        """
        if self.kw['nrm_method'] not in ['nrm', 'frm']:
            self.kw['nrm_method'] = 'nrm'
        if self.kw['messages'] == 'verbose':
            print(f"[INFO]: Forming A-matrix from pedigree at {pyp_utils.pyp_nice_time()}.")
        logger.info("Forming A-matrix from pedigree")

        try:
            if self.kw['nrm_method'] == 'nrm':
                self.nrm = pyp_nrm.fast_a_matrix(pedigree, self.kw)
                if self.kw['messages'] == 'verbose':
                    print(f"[INFO]: Formed A-matrix from pedigree using pyp_nrm.fast_a_matrix() at {pyp_utils.pyp_nice_time()}.")
                logger.info("Formed A-matrix from pedigree using pyp_nrm.fast_a_matrix()")
            else:
                self.nrm = pyp_nrm.fast_a_matrix_r(pedigree, self.kw)
                if self.kw['messages'] == 'verbose':
                    print(f"[INFO]: Formed A-matrix from pedigree using pyp_nrm.fast_a_matrix_r() at {pyp_utils.pyp_nice_time()}.")
                logger.info("Formed A-matrix from pedigree using pyp_nrm.fast_a_matrix_r()")
            return True
        except Exception as e:
            if self.kw['messages'] == 'verbose':
                print(f"[ERROR]: Unable to form A-matrix from pedigree using {self.kw['nrm_method']} at {pyp_utils.pyp_nice_time()}. Error: {e}")
            logger.error(f"Unable to form A-matrix from pedigree using {self.kw['nrm_method']}")
            return False

    def load(self, nrm_filename):
        """
        load() uses NumPy's `fromfile()` to load an array from a binary file.
        If the load is successful, self.nrm contains the matrix.
        """
        if self.kw['messages'] == 'verbose':
            print(f"[INFO]: Loading A-matrix from file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")
        logger.info(f"Loading A-matrix from file {nrm_filename}")

        try:
            self.nrm = np.fromfile(nrm_filename, dtype='float64', sep=self.kw.get('sepchar', ' '))
            size = int(math.sqrt(self.nrm.shape[0]))
            self.nrm = np.reshape(self.nrm, (size, size))
            if self.kw['messages'] == 'verbose':
                print(f"[INFO]: A-matrix successfully loaded from file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")
            logger.info(f"A-matrix successfully loaded from file {nrm_filename}")
            return True
        except Exception as e:
            if self.kw['messages'] == 'verbose':
                print(f"[ERROR]: Unable to load A-matrix from file {nrm_filename} at {pyp_utils.pyp_nice_time()}. Error: {e}")
            logger.error(f"Unable to load A-matrix from file {nrm_filename}")
            return False

    def save(self, nrm_filename, nrm_format=''):
        """
        save() uses NumPy's `tofile()` to save an array to a binary or text file.
        """
        if self.kw['messages'] == 'verbose':
            print(f"[INFO]: Saving A-matrix to file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")
        logger.info(f"Saving A-matrix to file {nrm_filename}")

        try:
            if not nrm_format:
                nrm_format = self.kw.get('nrm_format', 'text')
            if nrm_format == 'binary':
                self.nrm.tofile(nrm_filename)
            else:
                self.nrm.tofile(nrm_filename, sep=self.kw.get('sepchar', ' '))
            if self.kw['messages'] == 'verbose':
                print(f"[INFO]: A-matrix successfully saved to file {nrm_filename} at {pyp_utils.pyp_nice_time()}.")
            logger.info(f"A-matrix successfully saved to file {nrm_filename}")
            return True
        except Exception as e:
            if self.kw['messages'] == 'verbose':
                print(f"[ERROR]: Unable to save A-matrix to file {nrm_filename} at {pyp_utils.pyp_nice_time()}. Error: {e}")
            logger.error(f"Unable to save A-matrix to file {nrm_filename}")
            return False

    def printme(self):
        """
        printme() prints the NRM to the screen.
        """
        try:
            print(self.nrm)
        except Exception as e:
            print(f"[ERROR]: Unable to print A-matrix. Error: {e}")


def loadPedigree(
    options=None,
    optionsfile="pypedal.ini",
    pedsource="file",
    pedgraph=None,
    pedstream="",
    debugLoad=False,
):
    """
    loadPedigree() wraps pedigree creation and loading into a one-step process. 
    If the user passes both a dictionary and a filename, the dictionary will be 
    used instead of the filename unless the dictionary is empty.

    Parameters
    ----------
    options : dict, optional
        Dictionary of pedigree options.
    optionsfile : str, optional
        File from which pedigree options should be read.
    pedsource : str, optional
        Source of the pedigree ('file', 'graph', 'graphfile', 'db', 'gedcomfile',
        'textstream'). Default is 'file'.
    pedgraph : DiGraph, optional
        DiGraph from which to load the pedigree.
    pedstream : str, optional
        String of tuples to unpack into a pedigree.
    debugLoad : bool, optional
        When True, print debugging messages while loading.

    Returns
    -------
    NewPedigree or bool
        An instance of a NewPedigree object on success, False on failure.
    """
    if debugLoad:
        print("[DEBUG]: Debugging pyp_newclasses/loadPedigree()...")

    try:
        # Pass None, not {}, when no options dict was supplied. NewPedigree
        # reads the option file only when kw IS None, so normalising None to an
        # empty dict here meant loadPedigree(optionsfile=...) silently ignored
        # the file and then failed for want of a pedfile.
        # An explicitly supplied options dict still takes precedence.
        _pedigree = NewPedigree(kw=options if options else None,
                                kwfile=optionsfile)
        
        if debugLoad:
            print(f"[DEBUG]: Loading pedigree from source: {pedsource}")

        _pedigree.load(pedsource=pedsource, pedgraph=pedgraph, pedstream=pedstream)
        return _pedigree

    except PyPedalError:
        # A PyPedal-level failure is already specific and already explains
        # itself. Downgrading it to `return False` here hid string-ID load
        # failures: preprocess() raising was pointless while this handler turned
        # the raise back into a falsy return that callers then treated as a
        # pedigree.
        raise

    except TypeError as e:
        error_msg = f"[ERROR]: loadPedigree() failed due to a type mismatch. Error: {e}"
        print(error_msg)
        logger.error(error_msg)
        return False

    except Exception as e:
        error_msg = f"[ERROR]: loadPedigree() encountered an unexpected error. Error: {e}"
        print(error_msg)
        logger.error(error_msg)
        return False


def load_pedigree(*args, **kwargs):
    """Snake-case wrapper around loadPedigree() used by examples and the GUI."""
    if "options_file" in kwargs:
        kwargs["optionsfile"] = kwargs.pop("options_file")
    if "debug_load" in kwargs:
        kwargs["debugLoad"] = kwargs.pop("debug_load")
    return loadPedigree(*args, **kwargs)


def pad_id(by, animalID):
    """
    Generate a padded ID by combining birth year and animal ID.

    Unknown birth year uses identity padding, never a fake year.
    """
    return pyp_chronology.light_padded_id(by, animalID)


def safe_strip(value):
    return str(value).strip() if isinstance(value, str) else str(value)


# Key under which the per-pedigree string-ID collision registry is kept. It
# lives in the pedigree's own `kw` dict, so it is scoped to a single pedigree
# and a single load: two pedigrees loaded in one process cannot interfere, and
# reloading resets it.
_STRING_ID_REGISTRY = '_string_id_registry'


def string_to_int(idstring, mymaxint=9223372036854775807):
    """
    Map a string ID to an integer with md5 modulo 2**63 - 1.

    Shared by NewAnimal and LightAnimal. It previously existed only as a method
    on LightAnimal, which is why NewAnimal -- the default animal class -- had no
    way to handle the A/S/D pedformat codes at all.

    The hash is taken over the field exactly as it appears in the record,
    including any surrounding quotes, which is what PyPedal 2.0.4 does. Keeping
    that identical means originalID values still match the reference
    implementation and the differential harness continues to work on string
    pedigrees.
    """
    return int(hashlib.md5(str(idstring).encode()).hexdigest(), 16) % mymaxint


def hashed_string_id(idstring, mykw, field):
    """
    Hash a string ID and record it, raising if a distinct string already
    claimed the same integer.

    md5 modulo 2**63 makes a collision vanishingly unlikely, but the
    consequence of one is severe and entirely silent: two animals would merge
    into a single record and every relationship computed afterwards would be
    wrong, with nothing to indicate it. The registry is scoped to the pedigree
    being loaded and covers animal, sire and dam IDs together, because a
    collision between an animal ID and a sire ID merges records just as
    effectively as one between two animal IDs.
    """
    value = str(idstring)
    hashed = string_to_int(value)
    registry = mykw.setdefault(_STRING_ID_REGISTRY, {})
    previous = registry.get(hashed)
    if previous is not None and previous != value:
        raise PyPedalStringIDCollisionError(previous, value, hashed, field)
    registry[hashed] = value
    return hashed


# The exception hierarchy moved to PyPedal.pyp_errors, which imports nothing
# from PyPedal so that pyp_validate can raise PyPedal exceptions without a
# circular import back through this module. Re-exported here so that existing
# `from PyPedal.pyp_newclasses import PyPedalError` imports keep working.
from .pyp_errors import (  # noqa: E402,F401
    PyPedalError,
    PyPedalConfigurationError,
    PyPedalDependencyError,
    PyPedalInputError,
    PyPedalInternalError,
    PyPedalNotImplementedError,
    PyPedalOptionError,
    PyPedalPedigreeFormatError,
    PyPedalPedigreeInputFileNameError,
    PyPedalPedigreeSourceError,
    PyPedalStringIDCollisionError,
    PyPedalUsageError,
    PyPedalValidationError,
)
