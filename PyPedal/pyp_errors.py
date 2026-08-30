###############################################################################
# NAME: pyp_errors.py
# VERSION: see PyPedal.__version__
# LICENSE: LGPL
# Written for PyPedal 4.0; this module is not part of the original 2.0.4 sources.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
"""
PyPedal's exception hierarchy.

This module imports nothing from PyPedal, by design. The exception base class
was previously defined in ``pyp_newclasses``, which every other module depends
on; anything wanting to raise a PyPedal exception therefore had to import the
package's largest module and could not itself be imported by it.
``pyp_validate`` needs exactly that, so the hierarchy lives here instead.

``pyp_newclasses`` re-exports ``PyPedalError`` and
``PyPedalPedigreeInputFileNameError``, so existing imports of those names from
``pyp_newclasses`` continue to work unchanged.

The hierarchy
-------------

::

    PyPedalError
    |-- PyPedalInputError                    the data is malformed
    |   |-- PyPedalPedigreeFormatError       record does not match pedformat
    |   |-- PyPedalPedigreeSourceError       GEDCOM/adjlist parse failure
    |   +-- PyPedalPedigreeStructureError    the graph cannot be ordered
    |-- PyPedalConfigurationError            the run is misconfigured
    |   |-- PyPedalOptionError               bad .ini, missing required option
    |   +-- PyPedalPedigreeInputFileNameError    no pedfile and no simulation
    |-- PyPedalDependencyError               an optional dependency is absent
    |-- PyPedalUsageError                    the calling code is wrong (also ValueError)
    |-- PyPedalInternalError                 a PyPedal invariant was violated
    |-- PyPedalNotImplementedError           deliberately unavailable
    |-- PyPedalStringIDCollisionError        two string IDs hashed alike
    +-- PyPedalValidationError               a computed quantity is impossible

Why more than one class
-----------------------

These conditions call for different responses, and a single generic exception
would erase distinctions the call sites already make. Malformed data means fix
the file; a missing option means fix the configuration; an absent dependency
means install something; a bad argument means fix the calling code; a violated
invariant means report a PyPedal bug. The CLI maps each to its own exit status
so a wrapping script can branch on the difference.

**Library code must never terminate the process.** Nothing under ``PyPedal/``
raises ``SystemExit`` or calls ``sys.exit`` except the three top-level entry
points (``pyp_app.py`` twice, ``__main__.py`` once). Deciding that the program
is over belongs to the application, not to a library, and
``tests/test_error_contract.py`` enforces this with an AST scan so it cannot
regress quietly.
"""


class PyPedalError(Exception):
    """PyPedalError is the base class for exceptions in PyPedal."""
    pass


class PyPedalConfigurationError(PyPedalError):
    """
    The run is misconfigured: an option is missing, unreadable or inconsistent.

    Distinct from :class:`PyPedalInputError` because the fix is different --
    change the options, not the data.
    """


class PyPedalPedigreeInputFileNameError(PyPedalConfigurationError):
    """PyPedalPedigreeInputFileNameError is raised when a simulated pedigree
    is not requested and a pedigree file name is not provided.
    """

    def __init__(self):
        super().__init__()
        self.message = (
            'You did not request that a pedigree be simulated, and you did not provide the name of a '
            'pedigree file to be read.'
        )

    def __str__(self):
        return self.message


class PyPedalStringIDCollisionError(PyPedalError):
    """
    Raised when two distinct string IDs in one pedigree hash to the same
    integer.

    String IDs (pedformat codes ``A``/``S``/``D``) are mapped to integers with
    md5 modulo 2**63 - 1. The collision probability is negligible in practice,
    but a collision that did occur would silently merge two animals into one,
    corrupting every relationship computed afterwards with nothing to indicate
    it. Detecting it costs one dictionary lookup per ID.
    """

    def __init__(self, first, second, hashed, field):
        message = (
            'Two distinct %s IDs in this pedigree hash to the same integer '
            '(%d): %r and %r. Continuing would silently merge them into one '
            'animal and corrupt every relationship computed from this '
            'pedigree. Rename one of the IDs.' % (field, hashed, first, second)
        )
        super().__init__(message)
        self.message = message
        self.first = first
        self.second = second
        self.hashed = hashed
        self.field = field

    def __str__(self):
        return self.message


class PyPedalInputError(PyPedalError):
    """
    The data PyPedal was given is malformed.

    The user's *input* is at fault, not the way PyPedal was configured and not
    the calling code. A caller that catches this should expect to report a
    problem with the file, not to retry with different options.
    """


class PyPedalPedigreeFormatError(PyPedalInputError):
    """A pedigree record does not match the declared ``pedformat``."""


class PyPedalPedigreeSourceError(PyPedalInputError):
    """A pedigree source could not be parsed -- GEDCOM or an adjacency list."""


class PyPedalPedigreeStructureError(PyPedalInputError, ValueError):
    """
    The pedigree graph cannot be unambiguously ordered.

    Also a ``ValueError``, for the same reason ``PyPedalUsageError`` is:
    ``fast_reorder()`` already raised a bare ``ValueError`` on a circular
    reference, so callers that wrap an ordering call in ``except ValueError``
    keep working. A cyclic pedigree is also what ``ValueError`` means in Python
    -- the argument has the right type and an unusable value.

    A cycle, a self-parent, a duplicated animal ID, or a parent reference to an
    animal that has no record. PyPedal 2.0.4's manual already treated this as an
    error rather than something reordering resolves: PyPedal can reorder any
    pedigree unless there is an error in it that would prevent unambiguously
    placing parents before offspring.

    ``PyPedalInputError`` because the user's *data* is at fault -- not their
    configuration, and not the calling code. Retrying with different options
    will not help; the pedigree has to be corrected.

    Raised **before** anything is reordered or renumbered. Silently dropping
    the offending parent links and continuing would strip the edges that made
    the pedigree unorderable, and every coefficient downstream would then be
    computed on a graph the user never supplied.

    ``animals`` carries the IDs the refusal is about, so a caller can report
    them without parsing the message.
    """

    def __init__(self, message, animals=None):
        super().__init__(message)
        self.message = message
        self.animals = list(animals) if animals else []

    def __str__(self):
        return self.message


class PyPedalOptionError(PyPedalConfigurationError):
    """An option file could not be read, or a required option is absent."""


class PyPedalDependencyError(PyPedalError):
    """
    An optional dependency needed for the requested operation is not installed.

    Recoverable by installing something, which is why it is not lumped in with
    configuration errors: the message can name the extra to install.
    """


class PyPedalUsageError(PyPedalError, ValueError):
    """
    The calling code is wrong -- an argument outside its documented domain.

    Also a ``ValueError``, deliberately. Passing an unrecognised ``pedsource``
    is exactly what ``ValueError`` means in Python, and existing callers that
    write ``except ValueError`` around a PyPedal call keep working. Catching
    ``PyPedalError`` catches it too, so neither audience is surprised.
    """


class PyPedalInternalError(PyPedalError):
    """
    An invariant inside PyPedal was violated.

    Always a PyPedal defect. A caller cannot fix it by changing data, options or
    arguments, and should report it rather than handle it.
    """


class PyPedalNotImplementedError(PyPedalError):
    """
    Raised when a documented routine is deliberately unavailable because the
    reference defining it is not in hand.

    This is not "unfinished work". It marks a routine whose published algorithm
    PyPedal does not have, where the code that used to stand in its place was
    demonstrably wrong. Returning a plausible-looking wrong number is worse than
    returning nothing, because a caller cannot tell the difference; refusing is
    the only honest option until the reference arrives.

    Deliberately **not** a subclass of the builtin ``NotImplementedError``. That
    one conventionally means "an abstract method a subclass must override", and
    code that catches it is looking for a very different condition.

    ``blocked_item`` is an optional internal label for the refusal.
    """

    def __init__(self, routine, reason, blocked_item=None):
        message = '%s is unavailable: %s' % (routine, reason)
        if blocked_item:
            message += ' This domain is not supported in PyPedal 4.0 (%s).' % blocked_item
        super().__init__(message)
        self.message = message
        self.routine = routine
        self.reason = reason
        self.blocked_item = blocked_item

    def __str__(self):
        return self.message


class PyPedalExportFormatError(PyPedalError):
    """
    An in-memory value cannot be encoded into a declared fixed-width output
    format.

    Deliberately its own class rather than a reuse of an existing one. It is not
    ``PyPedalInputError``: the reproduced case is PyPedal's own ``missing_age``
    default colliding with PyPedal's own ``AGE`` field declaration, so telling
    the caller their pedigree is malformed would be wrong. It is not
    ``PyPedalValidationError`` either: nothing scientific has been computed, and
    diluting that class would cost it its meaning.

    Raised **before** anything is written. A fixed-width record whose fields do
    not match the declared layout desynchronises every reader after it, and the
    writer returns success -- so refusing is the only honest outcome. The value
    is never truncated or coerced to fit.

    ``declared_width`` and ``encoded_width`` are **byte** counts, not character
    counts: the target is a byte-oriented format, and a value can satisfy a
    character-count check while overflowing the field once encoded.
    """

    def __init__(self, fmt, field, declared_width, value, encoded_width,
                 note=None):
        message = (
            '%s export: field %s is declared %d bytes wide, but the value %r '
            'encodes to %d bytes. Refusing to write a record that does not '
            'match the declared layout; the value is not truncated to fit.'
            % (fmt, field, declared_width, value, encoded_width)
        )
        if note:
            message += ' ' + note
        super().__init__(message)
        self.message = message
        self.fmt = fmt
        self.field = field
        self.declared_width = declared_width
        self.value = value
        self.encoded_width = encoded_width
        self.note = note

    def __str__(self):
        return self.message


class PyPedalValidationError(PyPedalError):
    """
    Raised when a computed scientific quantity violates a mathematical property
    it must satisfy by definition.

    This always indicates a defect -- in PyPedal, or in input data reaching a
    routine whose documented preconditions it does not meet. It is never a
    normal outcome, and must not be handled by substituting a default value.
    """

    def __init__(self, message):
        super().__init__(message)
        self.message = message

    def __str__(self):
        return self.message
