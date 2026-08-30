"""
The published worked examples, checked against the values their papers print.

These fixtures exist so that later work is validated against **published**
numbers rather than against our own reimplementation of the same idea. That
only holds if the transcription is right, so each fixture is checked here, at
the point it enters the repository, using machinery that is ALREADY trusted:

* ``pyp_nrm.inbreeding`` -- cross-validated against Mrode and an independent
  oracle since an earlier revision;
* ``pyp_metrics.ballou_ancestral_inbreeding`` -- confirmed to implement
  Ballou (1997) eq. 1 exactly (a).

Nothing here validates the routines under adjudication. The Boichard fixtures
carry expected values for Appendix A and Appendix B, and those are checked when
the paper-derived oracle exists -- not by the code they are meant to judge.

A transcription is only as good as its checks, so the checks are deliberately
over-determined: the Suwanlee fixture alone is pinned by 13 published values.
"""
import unittest

from _pedhelpers import load_corpus

from PyPedal import pyp_metrics, pyp_nrm

# Published values are printed rounded to three decimals, so a computed value
# agrees with them exactly when it falls inside the rounding interval
# [published - 0.0005, published + 0.0005]. The bound is INCLUSIVE: several of
# these fixtures land precisely on it (0.3125 against a published 0.313), which
# is agreement, not disagreement. A tiny epsilon absorbs binary representation
# error at the endpoint.
THREE_DP = 5e-4 + 1e-9


def _renumbered(ped, paper_label):
    """
    Map a paper's animal label onto the ID the loaded pedigree actually uses.

    Loading reorders and renumbers, and it does NOT preserve the paper's
    numbering: on these fixtures it reverses the founders among themselves
    (suwanlee_fig1 founders 1,2,3,4 become 4,3,2,1). Non-founders happen to
    keep their IDs, but nothing guarantees that, and checking a published value
    against the wrong animal is exactly the failure these fixtures exist to
    prevent. So every lookup goes through ``idmap`` -- the canonical
    originalID -> renumberedID map -- rather than assuming identity.
    """
    idmap = ped.idmap
    for key in (paper_label, str(paper_label)):
        if key in idmap:
            return idmap[key]
    raise KeyError(f"animal {paper_label!r} is not in the pedigree's idmap")


def _f_and_fa(name, pedformat):
    ped = load_corpus(name, pedformat)
    fx = pyp_nrm.inbreeding(ped, output=False)["fx"]
    fa = pyp_metrics.ballou_ancestral_inbreeding(ped)
    # Re-key both results by the paper's own labels so the expectations below
    # can be written exactly as the figures print them.
    by_label = {}
    for animal in ped.pedigree:
        label = int(animal.originalID)
        rid = _renumbered(ped, label)
        by_label[label] = (fx.get(rid, 0.0), fa[rid])
    return ped, {k: v[0] for k, v in by_label.items()}, {k: v[1] for k, v in by_label.items()}


class TestFixtureLabelMapping(unittest.TestCase):
    """
    The published expectations are keyed to each paper's numbering, so the
    mapping from that numbering onto loaded IDs has to be total and injective.
    If it were not, a value could be silently checked against the wrong animal.
    """

    FIXTURES = (("boichard_fig1.ped", "asdg"), ("boichard_fig2.ped", "asdg"),
                ("ballou_fig1.ped", "asd"), ("suwanlee_fig1.ped", "asd"))

    def test_every_paper_label_maps_to_exactly_one_animal(self):
        for name, pedformat in self.FIXTURES:
            with self.subTest(pedigree=name):
                ped = load_corpus(name, pedformat)
                labels = [int(a.originalID) for a in ped.pedigree]
                mapped = [_renumbered(ped, label) for label in labels]
                self.assertEqual(sorted(labels), sorted(int(a.originalID) for a in ped.pedigree))
                self.assertEqual(len(set(mapped)), len(mapped),
                                 f"{name}: idmap is not injective")
                self.assertEqual(sorted(int(m) for m in mapped),
                                 sorted(int(a.animalID) for a in ped.pedigree),
                                 f"{name}: idmap does not cover the pedigree")


class TestBallouFigure1(unittest.TestCase):
    """
    Ballou (1997) Figure 1, article p.170. Values are printed beside each node.

    ``LL`` is also printed there and is deliberately NOT tested: it comes from
    Ballou's separate lethal-recessive Monte-Carlo model, is not ancestral
    inbreeding, and PyPedal does not compute it.
    """

    PUBLISHED_F = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.25, 6: 0.0, 7: 0.375}
    PUBLISHED_FA = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.125, 7: 0.125}

    def test_published_inbreeding(self):
        _, fx, _ = _f_and_fa("ballou_fig1.ped", "asd")
        for animal, want in self.PUBLISHED_F.items():
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, fx[animal], delta=THREE_DP)

    def test_published_ancestral_inbreeding(self):
        _, _, fa = _f_and_fa("ballou_fig1.ped", "asd")
        for animal, want in self.PUBLISHED_FA.items():
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, fa[animal], delta=THREE_DP)


class TestSuwanleeFigure1(unittest.TestCase):
    """
    Suwanlee et al. (2007) Figure 1, article p.491.

    The figure prints one row of values per generation. ``f`` and ``f_a_b``
    are deterministic and are tested; ``f_a_g`` is a stochastic gene-drop
    estimate published rounded to three decimals and is NOT tested here -- it
    is the subject of the predeclared denominator adjudication (b).

    Published values are given to three decimals, so they are compared to
    three decimals. The exact rational values are recorded in the fixture
    header.
    """

    PUBLISHED_F = {11: 0.125, 12: 0.125,
                   13: 0.313, 14: 0.313,
                   15: 0.438, 16: 0.438,
                   17: 0.547, 18: 0.547, 19: 0.547, 20: 0.547,
                   21: 0.633, 22: 0.633,
                   23: 0.633}

    PUBLISHED_FA_BALLOU = {11: 0.000, 12: 0.000,
                           13: 0.125, 14: 0.125,
                           15: 0.399, 16: 0.399,
                           17: 0.662, 18: 0.662, 19: 0.662, 20: 0.662,
                           21: 0.847, 22: 0.847,
                           23: 0.944}

    # Animals 15 and 16 are the ONE row where the published f_a does not lie in
    # the rounding interval of the exact recursion. See
    # TestSuwanleePublicationRounding below: the paper reproduces its own 0.399
    # only if it carries its printed three-decimal inbreeding coefficients back
    # into eq. 1 instead of the exact values. That is a publication artefact,
    # not a defect here, and it is measured rather than tolerated away.
    ROUNDING_ARTEFACT = (15, 16)

    EXACT_FA_BALLOU = {15: 0.3984375, 16: 0.3984375}

    def test_published_inbreeding(self):
        _, fx, _ = _f_and_fa("suwanlee_fig1.ped", "asd")
        for animal, want in self.PUBLISHED_F.items():
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, fx[animal], delta=THREE_DP)

    def test_published_ballou_ancestral_inbreeding(self):
        _, _, fa = _f_and_fa("suwanlee_fig1.ped", "asd")
        for animal, want in self.PUBLISHED_FA_BALLOU.items():
            if animal in self.ROUNDING_ARTEFACT:
                continue
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, fa[animal], delta=THREE_DP)

    def test_the_one_divergent_row_matches_the_exact_recursion(self):
        """
        Animals 15 and 16 must equal the exact value of Ballou eq. 1, not the
        paper's printed 0.399. Asserting the exact value keeps the discrepancy
        visible: if our implementation ever drifted TO 0.399 that would be a
        regression, not an improvement.
        """
        _, _, fa = _f_and_fa("suwanlee_fig1.ped", "asd")
        for animal, want in self.EXACT_FA_BALLOU.items():
            with self.subTest(animal=animal):
                self.assertAlmostEqual(want, fa[animal], places=9)
                self.assertNotAlmostEqual(0.399, fa[animal], delta=THREE_DP)

    def test_founders_and_first_generations_are_not_inbred(self):
        """
        The figure prints no values for animals 1-10, but the topology fixes
        them: four founders, two full-sib pairs from unrelated parents, and one
        half-sib pair. All must be non-inbred, and 9 and 10 being half sibs is
        what produces the published 0.125 one generation later.
        """
        _, fx, fa = _f_and_fa("suwanlee_fig1.ped", "asd")
        for animal in range(1, 11):
            with self.subTest(animal=animal):
                self.assertAlmostEqual(0.0, fx[animal], places=9)
                self.assertAlmostEqual(0.0, fa[animal], places=9)


class TestSuwanleePublicationRounding(unittest.TestCase):
    """
    Why one published Ballou value cannot be reproduced exactly.

    Suwanlee Figure 1 prints f_a = 0.399 for animals 15 and 16. Ballou eq. 1
    applied to exact inbreeding coefficients gives 0.3984375, which displays as
    0.398 -- outside the rounding interval of 0.399. Every other published row
    agrees.

    The explanation is that the paper fed its own THREE-DECIMAL inbreeding
    coefficients back into the recursion rather than the exact ones. This test
    demonstrates that, because an unexplained mismatch in a transcription is
    indistinguishable from a bad transcription, and the difference matters: one
    means the fixture is wrong, the other means the paper rounded.

    Running eq. 1 on the paper's printed f values reproduces all six published
    f_a values including 0.399; running it on exact f values reproduces five.
    That asymmetry is the evidence.
    """

    PUBLISHED_F_ROUNDED = (0.125, 0.313, 0.438, 0.547, 0.633)
    EXACT_F = (0.125, 0.3125, 0.4375, 0.546875, 0.6328125)
    PUBLISHED_FA = (0.125, 0.399, 0.662, 0.847, 0.944)

    @staticmethod
    def _ballou_chain(inbreeding_by_generation, round_to=None):
        """
        Ballou eq. 1 down a chain in which both parents are full sibs with the
        same f and f_a, which is what Figure 1's pedigree is from animal 11 on.
        The equation then collapses to f_a' = f_a + (1 - f_a) * f.
        """
        f_a, out = 0.0, []
        for f in inbreeding_by_generation:
            f_a = f_a + (1.0 - f_a) * f
            if round_to is not None:
                f_a = round(f_a, round_to)
            out.append(f_a)
        return out

    def test_rounded_inputs_reproduce_every_published_value(self):
        got = self._ballou_chain(self.PUBLISHED_F_ROUNDED, round_to=3)
        for want, have in zip(self.PUBLISHED_FA, got):
            self.assertAlmostEqual(want, have, places=9)

    def test_exact_inputs_reproduce_all_but_the_one_row(self):
        got = [round(v, 3) for v in self._ballou_chain(self.EXACT_F)]
        disagreements = [i for i, (w, h) in enumerate(zip(self.PUBLISHED_FA, got))
                         if abs(w - h) > 1e-9]
        self.assertEqual([1], disagreements,
                         "exactly one published row (animals 15/16) should "
                         "disagree with the exact recursion; if this changes, "
                         "the transcription or the recursion has changed")
        self.assertAlmostEqual(0.398, got[1], places=9)


class TestBallouIsNotGeneDropping(unittest.TestCase):
    """
    Suwanlee Figure 1 publishes Ballou's coefficient and the gene-drop estimate
    side by side ON THE SAME PEDIGREE, and they differ from the third
    generation onward. That is the paper's result, not sampling noise.

    This test guards the conclusion in both directions, because both errors are
    easy to make later:

    * asserting the two estimators AGREE is refuted by the published table;
    * asserting a universal ordering between them is not something Suwanlee
      proves -- it reports an observation in the scenarios it simulated, and
      separately reports gene dropping overestimating under lethal models.

    So what is pinned is the specific published divergence, not a law.
    """

    # (animal, Ballou, gene drop) rows where the paper's two columns differ.
    # Rows where the paper's two columns differ. Animals 15/16 are omitted:
    # their Ballou column carries the publication-rounding artefact examined in
    # TestSuwanleePublicationRounding, so they are not a clean comparison.
    PUBLISHED_DIVERGENCES = ((17, 0.662, 0.585), (21, 0.847, 0.728),
                             (23, 0.944, 0.822))

    def test_published_columns_are_not_equal(self):
        for animal, ballou, genedrop in self.PUBLISHED_DIVERGENCES:
            with self.subTest(animal=animal):
                self.assertNotAlmostEqual(
                    ballou, genedrop, places=2,
                    msg="Suwanlee Figure 1 publishes these as different estimators")

    def test_our_ballou_matches_the_ballou_column_not_the_gene_drop_column(self):
        _, _, fa = _f_and_fa("suwanlee_fig1.ped", "asd")
        for animal, ballou, genedrop in self.PUBLISHED_DIVERGENCES:
            with self.subTest(animal=animal):
                self.assertAlmostEqual(ballou, fa[animal], delta=THREE_DP)
                self.assertNotAlmostEqual(genedrop, fa[animal], delta=THREE_DP)


if __name__ == "__main__":
    unittest.main()
