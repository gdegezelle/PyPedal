#!/usr/bin/env python3

###############################################################################
# NAME: __init__.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################

"""
Package initialization for the PyPedal library.
"""

import logging

logging.getLogger("PyPedal").addHandler(logging.NullHandler())

__all__ = [
    "application",
    "pyp_chronology",
    "pyp_db",
    "pyp_demog",
    "pyp_graphics",
    "pyp_io",
    "pyp_metrics",
    "pyp_network",
    "pyp_newclasses",
    "pyp_nrm",
    "pyp_reports",
    "pyp_tests",
    "pyp_utils",
    "pyp_jbc",
    "pyp_reports_templates",
    "pyp_app",
]

from .__version__ import version as __version__

from . import (
    pyp_db,
    pyp_demog,
    pyp_graphics,
    pyp_io,
    pyp_metrics,
    pyp_network,
    pyp_newclasses,
    pyp_nrm,
    pyp_reports,
    pyp_tests,
    pyp_utils,
    pyp_jbc,
    pyp_reports_templates,
)