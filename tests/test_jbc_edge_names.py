"""
``new_draw_colored_pedigree`` edge construction (Ruff F821).

``_anim_node`` was read twice and assigned nowhere -- in this port **and** in
PyPedal 2.0.4, so this is an inherited defect rather than a migration one. Both
reads sit in the ``garrow=0`` branch, the documented way to ask for edges
without arrowheads:

    g.add_edge(_sire_edge, _node_name)
    if not _tf[garrow]:
        e = g.get_edge(_sire_edge, _anim_node)    # NameError

The correction is forced rather than chosen: ``AGraph.get_edge(u, v)`` raises
``KeyError`` unless both endpoints are existing nodes, and the line immediately
above added the edge ``(_sire_edge, _node_name)``. ``_node_name`` is the only
value that retrieves it.

THE END-TO-END TEST
-------------------
This branch was long unreachable. ``new_draw_colored_pedigree`` reads
``_m.userField`` at pyp_jbc.py:594, before it ever gets to the edges, and
``NewAnimal.userField`` did not exist in the Python 3 port -- it is set in
PyPedal 2.0.4 at pyp_newclasses.py:3235 and was dropped in the migration, while
ten direct and five indirect consumer sites still read it. Any call raised
``AttributeError`` first, so the fix could only be covered structurally.

 restored the attribute, so ``TestGarrowZeroDrawsUndirectedEdges``
below now executes the branch for real: it draws a pedigree with ``garrow=0``
and asserts the undirected edges are produced. The structural tests are kept --
they check the defect *class* (a name read but never bound) rather than this one
instance of it, and they run without Graphviz installed.
"""
import ast
import inspect
import os
import unittest

import pytest

from PyPedal import pyp_jbc

from _pedhelpers import chdir_tmp, load_corpus


class TestNoUnboundNamesInDrawColoredPedigree(unittest.TestCase):

    FUNC = pyp_jbc.new_draw_colored_pedigree

    def _unbound_locals(self, func):
        """
        Names read inside ``func`` that are neither assigned in it, nor
        parameters, nor module globals, nor builtins.
        """
        source = inspect.getdedent(inspect.getsource(func)) \
            if hasattr(inspect, "getdedent") else inspect.getsource(func)
        tree = ast.parse(source)
        fn = tree.body[0]

        assigned = {a.arg for a in fn.args.args}
        assigned |= {a.arg for a in fn.args.kwonlyargs}
        for node in ast.walk(fn):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                assigned.add(node.id)
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                for alias in node.names:
                    assigned.add((alias.asname or alias.name).split(".")[0])
            elif isinstance(node, ast.ExceptHandler) and node.name:
                assigned.add(node.name)
            elif isinstance(node, (ast.For, ast.comprehension)):
                target = getattr(node, "target", None)
                for sub in ast.walk(target) if target is not None else []:
                    if isinstance(sub, ast.Name):
                        assigned.add(sub.id)

        module_globals = set(vars(pyp_jbc)) | set(dir(__builtins__)) | set(dir(int))
        import builtins
        module_globals |= set(dir(builtins))

        read = {n.id for n in ast.walk(fn)
                if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
        return sorted(read - assigned - module_globals)

    def test_no_name_is_read_without_ever_being_bound(self):
        unbound = self._unbound_locals(self.FUNC)
        self.assertEqual(
            [], unbound,
            "new_draw_colored_pedigree reads names it never binds: %r" % unbound)

    def test_anim_node_specifically_is_gone(self):
        """
        Checked against the AST, not the source text: the fix carries an
        explanatory comment that names ``_anim_node``, and a substring match
        would fire on the comment explaining the defect rather than on the
        defect. What matters is that no *identifier* by that name is read.
        """
        tree = ast.parse(inspect.getsource(self.FUNC))
        names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
        self.assertNotIn(
            "_anim_node", names,
            "_anim_node is read as an identifier, and it is never assigned "
            "anywhere in PyPedal or in PyPedal 2.0.4")

    def test_each_get_edge_matches_the_add_edge_above_it(self):
        """
        The property that makes the correction forced: get_edge() must ask for
        the same endpoints add_edge() just created, or it raises KeyError.
        """
        source = inspect.getsource(self.FUNC)
        added, fetched = [], []
        for line in source.splitlines():
            line = line.strip()
            if line.startswith("g.add_edge("):
                added.append(line[len("g.add_edge("):].rstrip(")"))
            elif "g.get_edge(" in line:
                fetched.append(line.split("g.get_edge(", 1)[1].rstrip(")"))
        self.assertTrue(added and fetched)
        self.assertEqual(
            len(added), len(fetched),
            "every add_edge in this function is followed by a get_edge of the "
            "same edge; a mismatch means one of them changed")
        for a, f in zip(added, fetched):
            self.assertEqual(a, f, f"get_edge({f}) does not match add_edge({a})")


class TestGarrowZeroDrawsUndirectedEdges(unittest.TestCase):
    """
    The end-to-end cover for the ``_anim_node`` repair, unblocked by .

    ``garrow=0`` maps to ``_tf[0] is False`` (pyp_jbc.py:505), so ``not _tf[garrow]``
    is true and the ``get_edge`` branch executes -- the branch that used to raise
    ``NameError`` on ``_anim_node``. Reaching it at all is most of the assertion;
    the ``dir=none`` attribute it sets is the observable result.
    """

    def _draw(self, tmp, **kwargs):
        ped = load_corpus("mrode.ped")
        shading = {a.animalID: len(a.sons) + len(a.daus) for a in ped.pedigree}
        dot = os.path.join(tmp, "garrow")
        pyp_jbc.new_draw_colored_pedigree(
            ped, shading, gfilename=dot, gtitle="garrow",
            gformat="dot", gdot=1, **kwargs)
        return ped, dot + ".dot"

    def test_undirected_edges_are_produced(self):
        pytest.importorskip("pygraphviz")
        with chdir_tmp() as tmp:
            ped, dotfile = self._draw(tmp, garrow=0)
            self.assertTrue(os.path.exists(dotfile), "no .dot file was written")
            text = open(dotfile).read()

        # One edge per known parent in mrode.ped, and every one undirected.
        expected_edges = sum(
            (a.sireID != ped.kw["missing_parent"])
            + (a.damID != ped.kw["missing_parent"])
            for a in ped.pedigree)
        self.assertTrue(expected_edges)
        self.assertEqual(
            expected_edges, text.count("dir=none"),
            "garrow=0 must mark every parent edge dir=none; got:\n" + text)

    def test_garrow_one_leaves_edges_directed(self):
        """
        The control. Without it, a build that silently wrote dir=none on every
        edge regardless of garrow would pass the test above.
        """
        pytest.importorskip("pygraphviz")
        with chdir_tmp() as tmp:
            _, dotfile = self._draw(tmp, garrow=1)
            text = open(dotfile).read()
        self.assertNotIn("dir=none", text)


if __name__ == "__main__":
    unittest.main()
