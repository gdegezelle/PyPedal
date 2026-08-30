"""
The matplotlib-missing error paths in ``pyp_graphics``.

Found by static analysis (Ruff ``F821``) during an earlier revision stage E3, not by any
test -- these branches only execute when matplotlib is absent, which it never is
in a development environment.

The migration defect: PyPedal 2.0.4 gave ``pcolor_matrix_pylab``,
``spy_matrix_pylab`` and ``plot_line_xy`` a ``debug=False`` parameter and gated
their diagnostics on ``if debug:``. The Python 3 port dropped the parameter and
replaced the guard with ``if pedobj.kw["messages"] == "verbose":`` -- the
project-wide convention, but wrong here, because all three take a MATRIX or a
DICT and no pedigree. The result was ``NameError: name 'pedobj' is not defined``
raised from inside the handler that was supposed to report the missing module:
the error path destroyed the error message.

These tests force the ``ImportError`` branch by making ``import matplotlib``
fail, and assert that each function returns rather than raising.
"""
import builtins
import unittest

import numpy as np

from PyPedal import pyp_graphics


class _NoMatplotlib:
    """Context manager that makes ``import matplotlib`` raise ImportError."""

    def __enter__(self):
        self._real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "matplotlib" or name.startswith("matplotlib."):
                raise ImportError("matplotlib is not installed (simulated)")
            return self._real_import(name, *args, **kwargs)

        builtins.__import__ = fake_import
        return self

    def __exit__(self, *exc):
        builtins.__import__ = self._real_import
        return False


class TestMatplotlibMissingIsReportedNotCrashed(unittest.TestCase):

    def test_pcolor_matrix_pylab_returns_instead_of_raising(self):
        A = np.eye(3)
        with _NoMatplotlib():
            result = pyp_graphics.pcolor_matrix_pylab(A, fname="unused")
        self.assertFalse(result)

    def test_spy_matrix_pylab_returns_instead_of_raising(self):
        A = np.eye(3)
        with _NoMatplotlib():
            result = pyp_graphics.spy_matrix_pylab(A, fname="unused")
        self.assertFalse(result)

    def test_plot_line_xy_returns_instead_of_raising(self):
        with _NoMatplotlib():
            result = pyp_graphics.plot_line_xy({1: 2, 3: 4}, gfilename="unused")
        self.assertFalse(result)

    def test_the_debug_flag_is_accepted_again(self):
        """
        PyPedal 2.0.4 had ``debug=False`` on all three. The port dropped it,
        which is what left the guard referring to a name that was never in
        scope. Passing it must work, and must not change the return.
        """
        A = np.eye(3)
        with _NoMatplotlib():
            self.assertFalse(pyp_graphics.pcolor_matrix_pylab(A, "unused", debug=True))
            self.assertFalse(pyp_graphics.spy_matrix_pylab(A, "unused", debug=True))
            self.assertFalse(
                pyp_graphics.plot_line_xy({1: 2}, gfilename="unused", debug=True))

    def test_no_undefined_name_survives_in_these_handlers(self):
        """
        Structural backstop. The three functions must not reference ``pedobj``
        at all -- none of them receives one, so any occurrence is the defect
        returning.
        """
        import ast
        import inspect

        for func in (pyp_graphics.pcolor_matrix_pylab,
                     pyp_graphics.spy_matrix_pylab,
                     pyp_graphics.plot_line_xy):
            with self.subTest(function=func.__name__):
                tree = ast.parse(inspect.getsource(func))
                names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
                self.assertNotIn(
                    "pedobj", names,
                    f"{func.__name__} refers to pedobj, which it never receives")


if __name__ == "__main__":
    unittest.main()
