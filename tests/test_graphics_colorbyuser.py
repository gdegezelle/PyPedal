"""
``colorByUser`` emits the colormap's own RGB, and needs no undeclared dependency.

WHAT CHANGED
------------
``new_draw_pedigree(colorByUser=True)`` used to pass each colormap RGB through
``pyp_graphics.get_colour_name()``, which asked ``webcolors`` for the nearest CSS
colour *name*. That was two problems at once:

* **it could not run.** ``closest_colour()`` called ``webcolors.css3_hex_to_names``
  -- lowercase -- an API predating webcolors 1.11.1, whose raw mappings were
  withdrawn entirely in 24.6.0. And no extra in ``pyproject.toml`` declared
  ``webcolors`` at all, so the path raised in every supported installation,
  including the full optional profile.
* **it was lossy.** An exact CSS name was returned only when ``rgb_to_name()``
  happened to match; otherwise it returned the *nearest* name by Euclidean RGB
  distance, which is a different colour from the one the colormap chose.

Graphviz accepts ``#RRGGBB`` directly, so the name lookup is gone.

WHAT THESE TESTS ASSERT -- AND WHAT THEY MUST NOT
-------------------------------------------------
The expected colour is recomputed **from the colormap**, the same source the
drawing code reads. It is deliberately *not* compared against the old
nearest-CSS-name output: that approximation is the thing being removed, and
pinning it would re-enshrine the defect as the specification.
"""
import os
import re
import unittest

import pytest

from _pedhelpers import chdir_tmp, load_corpus


def _expected_hex(pedobj):
    """
    Recompute the fill colours the way ``new_draw_pedigree`` does: an "Accent"
    colormap with one entry per unique userField, indexed by that field's
    position in ``metadata.unique_field_list``.
    """
    import matplotlib

    # Same call the drawing code makes. matplotlib.cm.get_cmap was removed in
    # 3.11; .resampled() is the documented replacement and is bit-identical.
    cmap = matplotlib.colormaps["Accent"].resampled(
        pedobj.metadata.num_unique_fields)
    fields = list(pedobj.metadata.unique_field_list)
    out = {}
    for animal in pedobj.pedigree:
        rgb = cmap(fields.index(animal.userField))
        out[animal.animalID] = "#%02x%02x%02x" % (
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    return out


class TestColorByUserEmitsColormapRGB(unittest.TestCase):

    def _draw(self, tmp):
        from PyPedal import pyp_graphics

        ped = load_corpus("userfield.ped", "asdu")
        base = os.path.join(tmp, "cbu")
        pyp_graphics.new_draw_pedigree(
            ped, gfilename=base, gtitle="colorByUser", colorByUser=True,
            gdot=1)
        return ped, base + ".dot"

    def test_every_node_is_filled_with_the_exact_colormap_colour(self):
        pytest.importorskip("matplotlib")
        pytest.importorskip("pygraphviz")

        with chdir_tmp() as tmp:
            ped, dotfile = self._draw(tmp)
            self.assertTrue(os.path.exists(dotfile), "no .dot file was written")
            text = open(dotfile).read()
            expected = _expected_hex(ped)

        emitted = set(re.findall(r'fillcolor="(#[0-9a-f]{6})"', text))
        self.assertTrue(emitted, "no hex fillcolor was emitted:\n" + text)
        # Every colour in the drawing is one the colormap actually produced...
        self.assertTrue(emitted <= set(expected.values()),
                        "emitted %r, colormap gives %r"
                        % (sorted(emitted), sorted(set(expected.values()))))
        # ...and every distinct userField got its own colour.
        self.assertEqual(len(set(expected.values())), len(emitted))

    def test_no_css_colour_names_are_emitted(self):
        """
        The regression guard for the lossy step. A CSS name in the output means
        the nearest-colour approximation came back.
        """
        pytest.importorskip("matplotlib")
        pytest.importorskip("pygraphviz")

        with chdir_tmp() as tmp:
            _ped, dotfile = self._draw(tmp)
            text = open(dotfile).read()

        for match in re.findall(r'fillcolor="([^"]+)"', text):
            self.assertRegex(
                match, r"^#[0-9a-f]{6}$",
                "fillcolor %r is not a hex triplet" % match)


class TestWebcolorsIsGone(unittest.TestCase):
    """
    The dependency must not creep back. It is not declared in any extra, and no
    installable version provides the API the old code called, so a re-introduced
    import would fail only in the environments that actually draw -- which is
    how it went unnoticed in the first place.
    """

    def test_pyp_graphics_does_not_reference_webcolors(self):
        """
        Checked against the AST, not the source text. The removal carries an
        explanatory comment that necessarily names ``webcolors.css3_hex_to_names``,
        and a substring match would fire on the comment explaining the defect
        rather than on the defect. What matters is that no *code* imports the
        module or reads an attribute from it.
        """
        import ast

        from PyPedal import pyp_graphics

        tree = ast.parse(open(pyp_graphics.__file__).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    self.assertNotEqual("webcolors", alias.name.split(".")[0])
            elif isinstance(node, ast.ImportFrom):
                self.assertNotEqual("webcolors", (node.module or "").split(".")[0])
            elif isinstance(node, ast.Attribute):
                if isinstance(node.value, ast.Name):
                    self.assertNotEqual(
                        "webcolors", node.value.id,
                        "pyp_graphics reads webcolors.%s" % node.attr)

    def test_the_removed_helpers_are_gone(self):
        from PyPedal import pyp_graphics

        for name in ("closest_colour", "get_colour_name"):
            self.assertFalse(
                hasattr(pyp_graphics, name),
                "%s is back; it cannot run against any installable webcolors"
                % name)


if __name__ == "__main__":
    unittest.main()
