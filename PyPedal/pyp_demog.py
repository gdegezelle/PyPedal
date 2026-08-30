#!/usr/bin/env python3

"""
pyp_demog.py - A module for demographic calculations on the population described in a pedigree.

Version: see PyPedal.__version__
Author: John B. Cole (john.b.cole@gmail.com)
License: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later

This module contains utilities for demographic calculations on pedigree data.
"""

import logging
from . import pyp_utils

# Define some globals in case the user forgets to call set_base_year() and
# set_age_units().
BASE_DEMOGRAPHIC_YEAR = 1800
BASE_DEMOGRAPHIC_UNIT = 'year'
SEX_CODE_MAP = {'m': 'Male', 'f': 'Female', 'u': 'Unk'}


def set_base_year(year=1800):
    """
    set_base_year() defines a global variable, BASE_DEMOGRAPHIC_YEAR.
    """
    global BASE_DEMOGRAPHIC_YEAR
    BASE_DEMOGRAPHIC_YEAR = year


def set_age_units(units='year'):
    """
    set_age_units() defines a global variable, BASE_DEMOGRAPHIC_UNIT.
    """
    global BASE_DEMOGRAPHIC_UNIT
    valid_units = ['year', 'month', 'day']
    BASE_DEMOGRAPHIC_UNIT = units if units in valid_units else 'year'


def age_distribution(pedobj, sex=1):
    """
    Print histograms of the demographic year-offset (``animal.age``), not
    biological age. Animals whose year-offset is the configured missing-age
    marker are omitted.
    """
    age_dict = {}
    age_freq_total = 0.0
    missing_age = pedobj.kw.get('missing_age', -999)
    missing_value = pedobj.kw.get('missing_value', -999.0)

    def _is_missing_age(value):
        return value is None or value == missing_age or value == missing_value

    if not pedobj.pedigree:
        return
    if _is_missing_age(pedobj.pedigree[0].age):
        pyp_utils.set_age(pedobj)

    def _known_ages(animals):
        out = {}
        for animal in animals:
            if _is_missing_age(animal.age):
                continue
            out[animal.age] = out.get(animal.age, 0) + 1
        return out

    if not sex:
        age_dict = _known_ages(pedobj.pedigree)
        age_hist = pyp_utils.simple_histogram_dictionary(age_dict) if age_dict else {}
        known_n = sum(age_dict.values())

        if pedobj.kw.get('debug_messages'):
            print("-" * 80)
            print("Population year-offset from base year")
            print("-" * 80)
            print("\tYear-offset\tCount\tFrequency\tHistogram")
            for key, count in age_dict.items():
                freq = count / known_n if known_n else 0.0
                age_freq_total += freq
                print(f"\t{key}\t{count}\t{freq:.2f}\t{age_hist[key]}")
            print(f"\tTOTAL\t{known_n}\t{age_freq_total:.2f}")
            print("-" * 80)
    else:
        males, females, unknowns = [], [], []
        for animal in pedobj.pedigree:
            if animal.sex == 'm':
                males.append(animal)
            elif animal.sex == 'f':
                females.append(animal)
            else:
                unknowns.append(animal)

        male_dict = _known_ages(males)
        female_dict = _known_ages(females)
        unknown_dict = _known_ages(unknowns)

        male_hist = pyp_utils.simple_histogram_dictionary(male_dict) if male_dict else {}
        female_hist = pyp_utils.simple_histogram_dictionary(female_dict) if female_dict else {}
        unknown_hist = pyp_utils.simple_histogram_dictionary(unknown_dict) if unknown_dict else {}

        if pedobj.kw.get('messages') == 'verbose':
            print("-" * 80)
            print("Population year-offset from base year by sex")
            print("-" * 80)
            print("Males")
            print("\tYear-offset\tCount\tFrequency\tHistogram")
            _print_histogram(male_dict, sum(male_dict.values()) or 1, male_hist)

            print("Females")
            print("\tYear-offset\tCount\tFrequency\tHistogram")
            _print_histogram(female_dict, sum(female_dict.values()) or 1, female_hist)

            print("Unknowns")
            print("\tYear-offset\tCount\tFrequency\tHistogram")
            _print_histogram(unknown_dict, sum(unknown_dict.values()) or 1, unknown_hist)


def sex_ratio(pedobj):
    """
    Returns a dictionary containing the proportion of males and females in the population.

    :param pedobj: An instance of a PyPedal NewPedigree object.
    :return: A dictionary containing entries for each sex/gender code.
    """
    sexratiodict = {s: 0 for s in SEX_CODE_MAP.keys()}

    for animal in pedobj.pedigree:
        sexratiodict[animal.sex] = sexratiodict.get(animal.sex, 0) + 1

    if pedobj.kw.get('messages') == 'verbose':
        print("-" * 80)
        print("Overall Sex Ratio")
        print("-" * 80)
        print(f"(n = {len(pedobj.pedigree)})")
        print("Sex\tCount\tFrequency")
        for sex, count in sexratiodict.items():
            freq = count / len(pedobj.pedigree)
            print(f"{SEX_CODE_MAP.get(sex, 'Unknown')}:\t{count}\t{freq:.2f}")

        if sexratiodict['u'] > 0:
            marginal = sexratiodict['m'] + sexratiodict['f']
            print("Conditional Sex Ratio")
            print("-" * 80)
            print(f"(n = {marginal})")
            print("Sex\tCount\tFrequency")
            print(f"{SEX_CODE_MAP['m']}:\t{sexratiodict['m']}\t{sexratiodict['m'] / marginal:.2f}")
            print(f"{SEX_CODE_MAP['f']}:\t{sexratiodict['f']}\t{sexratiodict['f'] / marginal:.2f}")

    return sexratiodict


def founders_by_year(pedobj):
    """
    Count founders by known recorded birth year.

    Unknown-year founders are omitted from the year-keyed dict. Gap filling
    runs only between actual known integer years.
    """
    founderbyyeardict = {}
    for founder in pedobj.metadata.unique_founder_list:
        mapped = pedobj.idmap.get(founder, founder)
        birth_year = pedobj.pedigree[int(mapped) - 1].by
        if birth_year is None:
            continue
        founderbyyeardict[birth_year] = founderbyyeardict.get(birth_year, 0) + 1

    years = sorted(founderbyyeardict.keys())
    if years:
        for year in range(years[0], years[-1] + 1):
            founderbyyeardict.setdefault(year, 0)

    return founderbyyeardict


def _print_histogram(age_dict, total_count, hist_dict):
    """
    Helper function to print a histogram.

    :param age_dict: Dictionary with year-offset as key and count as value.
    :param total_count: Total count for frequency calculation.
    :param hist_dict: Histogram dictionary.
    """
    age_freq_total = 0.0
    if not age_dict or not total_count:
        print("\tTOTAL\t0\t0.00")
        print("-" * 80)
        return
    for age, count in age_dict.items():
        freq = count / total_count
        age_freq_total += freq
        print(f"\t{age}\t{count}\t{freq:.2f}\t{hist_dict[age]}")
    print(f"\tTOTAL\t{total_count}\t{age_freq_total:.2f}")
    print("-" * 80)
