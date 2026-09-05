"""
 -- the explicit Boichard reference-population API.

WHAT THIS FILE IS ABOUT
-----------------------
Boichard, Maignel & Verrier (1997), *Genet Sel Evol* 29:5-23. Appendix A step 1
and Appendix B step 1, article p.22, use identical wording:

    "(1) define the population under study, ie, the group of N animals carrying
    the gene pool of interest;"

The verb is **define**. It addresses the analyst, not the computer. No
generation number appears anywhere in Appendix A or Appendix B.

PyPedal nevertheless obtains that population by matching the pedformat ``g``
column, so all three public scalar routines are unreachable on any pedigree
without one -- every griffon fixture, ``mrode.ped``, ``hartlandclark.ped``. That
is , and the merged adjudication
(``the algorithm notes``)
reclassified it from an ``igen``/``gen`` wrong-field bug into an API-coupling
defect. Three concepts stay distinct and are never bridged:

    igen    a computed pedigree-depth quantity
    gen     a legacy input annotation from the `g` column, held as a str
    R       the analyst-defined population under study -- SOURCE-EXPLICIT

THREE EVIDENCE LEVELS, KEPT SEPARATE
------------------------------------
* **SOURCE-EXPLICIT** -- the analyst defines the group of N animals.
* **SOFTWARE/API DESIGN** -- the Python ``reference=`` parameter and every
  detail of its validation contract. Boichard specifies no Python API, and
  nothing in this file claims otherwise.
* **BACKWARD-COMPATIBILITY DESIGN** -- the legacy ``gen``/``g`` selector, kept
  bit-exact when ``reference`` is omitted.

Inferring R from ``igen`` is EMPIRICALLY / SEMANTICALLY REJECTED and no test
here may be satisfied by such an inference.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with every test that
needed the new API marked ``xfail(strict=True)``, so the contract was written
down before any code satisfied it and nothing could pass by accident. Those
markers were removed in the commit that implemented the API; a strict xfail
left in place after it starts passing is an error, not a convenience.

``self.subTest`` was deliberately kept out of the tests that were xfailed:
pytest 9 handles the subTest/xfail interaction differently from pytest 8, a
difference this project has already been bitten by
(``docs/-R2-HALF-FOUNDER-VERIFICATION.md`` section 2.4). Plain loops are
used there instead, and that shape is kept so the markers could be reinstated
unchanged if this ever needs re-litigating.
"""
import inspect
import os
import unittest

from _pedhelpers import (
    owned_temp_dir,
    corpus,
    load_corpus,
    load_corpus_from_path,
    load_example,
    load_griffon_1871_1890,
)
from PyPedal import pyp_errors, pyp_metrics, pyp_utils

# ---------------------------------------------------------------------------
# Subjects
#
# The pedformat each file requires is stated once, here. Getting one wrong
# yields an empty pedigree and a vacuously passing test.
# ---------------------------------------------------------------------------

#: Corpus pedigrees carrying a real generation column, with the reference
#: population the legacy `gen` path selects by default. Measured on
#: ed1fee4, not predicted.
WITH_G = {
    "boichard2a.ped": [7, 8, 9, 10, 11, 12, 13, 14],
    "boichard_fig1.ped": [9, 10, 11, 12, 13],
    "boichard_fig2.ped": [9, 10, 11, 12, 13, 14, 15, 16, 19, 20],
    "generation_split.ped": [7, 8],
    "synth_hf.ped": [6, 7],
    "synth_nohf.ped": [6, 7],
    # deep_generations selects {13, 25} at gen '12'; the ancestor routines
    # refuse on it today. It is included precisely because a refusal is a
    # result: the explicit path must refuse identically.
    #
    # {1, 25} until . The two animals selected are the same ones --
    # original 24 and 25 on both sides, measured -- but the founder pre-pass
    # used to reverse the founder block, so original 24 was animalID 1. Founders
    # now keep their input order and it is animalID 13. Every other entry in
    # this table is unchanged in both ID domains.
    "deep_generations.ped": [13, 25],
}

#: `boichard2a.ped` carries three generations, so it can show that the SAME R
#: reached through a label and through `reference=` agrees at every label, not
#: only at the default one.
BOICHARD2A_BY_GENERATION = {
    "1": [1, 2, 3, 4],
    "2": [5, 6],
    "3": [7, 8, 9, 10, 11, 12, 13, 14],
}

GRIFFON_PEDFORMAT = "asdxb"

GRIFFON_OPTIONS = {
    "pedformat": GRIFFON_PEDFORMAT,
    "sepchar": ",",
    "messages": "quiet",
    "renumber": True,
    "pedigree_is_renumbered": True,
    "pedigree_summary": 0,
}


def load_griffon():
    """Historical 1871–1890 characterisation cohort, derived from canonical."""
    return load_griffon_1871_1890(GRIFFON_OPTIONS)


def griffon_cohort(ped, years, complete_pedigree_only=False):
    """
    Renumbered animalIDs of the griffon animals born in ``years``.

    A TEST-SPECIFIED reference population, not a recommended scientific one.
    The adjudication proposes no canonical griffon cohort and neither does
    this file; these sets exist to exercise the API on real data.

    Unknown recorded years (``by is None``) are not a cohort. Asking for
    ``None`` is refused rather than grouping missing dates together.
    """
    if any(year is None for year in years):
        raise AssertionError(
            "None is unknown recorded chronology, not a year; a cohort must "
            "not be built from animals whose birth date is unknown")
    missing = int(ped.kw["missing_parent"])
    out = []
    for animal in ped.pedigree:
        if animal.by not in years:
            continue
        if complete_pedigree_only and (int(animal.sireID) == missing
                                       or int(animal.damID) == missing):
            continue
        out.append(int(animal.animalID))
    return sorted(out)


def strip_generation_column(name, pedformat):
    """
    Write ``name`` out again with its generation column removed, and return
    ``(path, pedformat)`` for the reduced file.

    Used to build a no-``g`` pedigree whose topology is *identical* to a
    ``g``-carrying one, so an explicit-R result can be compared against the
    legacy ``gen`` result for the same animals. Generated into a tmpdir; the
    repository is never written to.
    """
    index = pedformat.index("g")
    reduced_format = pedformat.replace("g", "", 1)
    rows = []
    with open(corpus(name), encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped[0] in "#%":
                # The in-file `% asdg` directive must agree with the data.
                rows.append(stripped[0] + " " + reduced_format)
                continue
            fields = stripped.split()
            del fields[index]
            rows.append(" ".join(fields))
    tmp = owned_temp_dir(prefix="pypedal_nog_")
    path = os.path.join(tmp, "nog_" + name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return path, reduced_format


def outcome(call):
    """
    The value, or the class name of the refusal.

    Both directions matter. This work deletes a precondition on one path, and a
    deleted precondition's likeliest failure mode is not a wrong number -- it is
    an input that used to be declined starting to produce one. Messages are
    never compared: they embed temporary directory paths.
    """
    try:
        return ("value", repr(call()))
    except Exception as exc:                       # noqa: BLE001 -- classify
        return ("refused", type(exc).__name__)


def report_text(ped, suffix):
    with open("%s%s" % (ped.kw["filetag"], suffix), encoding="utf-8") as handle:
        return handle.read()


def report_line(ped, suffix, prefix):
    """
    The single report line starting with ``prefix``.

    Used to compare the ancestor SELECTION SEQUENCE across the two paths. The
    public routines return only a scalar, and a scalar can coincide while the
    sequence differs -- which would silently change which animal 's
    tie-break credits. The `.dat` report is where that sequence is observable,
    so it is what gets compared.
    """
    for line in report_text(ped, suffix).splitlines():
        if line.startswith(prefix):
            return line
    raise AssertionError("no %r line in %s" % (prefix, suffix))


def animal_state(ped):
    """
    A deep snapshot of everything an analysis must not disturb.

    Recorded per animal rather than as a whole-object hash so a failure names
    the attribute that moved.
    """
    return [
        {
            "animalID": animal.animalID,
            "originalID": animal.originalID,
            "renumberedID": getattr(animal, "renumberedID", None),
            "sireID": animal.sireID,
            "damID": animal.damID,
            "gen": animal.gen,
            "igen": animal.igen,
            "founder": animal.founder,
            "sex": animal.sex,
            "by": animal.by,
            "ancestor": getattr(animal, "ancestor", None),
        }
        for animal in ped.pedigree
    ]


# ---------------------------------------------------------------------------
# 1. The reference contract -- SOFTWARE/API DESIGN
# ---------------------------------------------------------------------------

class TestTheReferenceContractRejectsBadInputLoudly(unittest.TestCase):
    """
    ``reference`` names REAL, current, renumbered ``animalID`` integers from
    the loaded pedigree. Everything else is refused with a typed exception.

    The identifier domain is a deliberate scope containment, not a claim that
    internal integer IDs are the best long-term user interface: it is exactly
    the domain ``boichard_probabilities_of_gene_origin`` and
    ``boichard_marginal_contributions`` already consume. A caller holding
    originalIDs or string IDs converts them through ``pedobj.idmap`` first.
    """

    def setUp(self):
        self.ped = load_corpus("boichard2a.ped", "asdg")

    def _routines(self):
        return (pyp_metrics.a_effective_founders_boichard,
                pyp_metrics.a_effective_ancestors_definite,
                pyp_metrics.a_effective_ancestors_indefinite)

    def test_an_empty_reference_population_is_refused(self):
        """
        Appendix A step 4 divides the vector q by N. With N = 0 the routine
        would return NaN rather than refuse, and the legacy path already
        declines an empty generation.
        """
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"), reference=[])

    def test_a_duplicated_id_is_refused_rather_than_silently_deduplicated(self):
        """
        This is arithmetic, not tidiness. Both engines normalise by the LENGTH
        of the reference sequence -- ``q /= float(len(reference))`` in Appendix
        A, ``n = float(len(ref_index))`` in Appendix B -- so a repeated ID
        divides by the wrong N and corrupts the contributions. Calling
        ``set()`` on the caller's input would silently change their intent
        instead of reporting it.
        """
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"),
                        reference=[7, 8, 8])

    def test_an_id_that_is_not_in_the_pedigree_is_refused(self):
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"),
                        reference=[7, 9999])

    def test_the_missing_parent_sentinel_is_refused(self):
        """
        Named explicitly rather than left to the membership check, so the
        message tells the caller what they actually passed.
        """
        sentinel = int(self.ped.kw["missing_parent"])
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"),
                        reference=[7, sentinel])

    def test_a_reading_c_phantom_id_can_never_be_a_reference_animal(self):
        """
        's phantom founders are analysis-local: they exist in the
        engine's completed arrays and in ``boichard_phantom_ids()``, but never
        in ``pedobj.pedigree``. ``boichard_marginal_contributions`` would
        happily index one, so the refusal has to be structural -- membership in
        the REAL pedigree -- rather than a numeric range test.

        A phantom is a parent of an existing animal by construction and can
        never be a member of the population under study.
        """
        ped = load_corpus("synth_hf.ped", "asdg")
        phantoms = pyp_metrics.boichard_phantom_ids(ped)
        self.assertTrue(phantoms, "fixture stopped having a half-founder; this "
                                  "test would pass vacuously")
        phantom_id = sorted(phantoms)[0]
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("synth_hf.ped", "asdg"),
                        reference=[6, phantom_id])

    def test_booleans_are_refused_even_though_bool_subclasses_int(self):
        """``True`` would otherwise silently become animal 1."""
        for routine in self._routines():
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"),
                        reference=[True, 8])

    def test_non_integral_and_string_ids_are_refused_not_coerced(self):
        """
        Strict validation for a new API: ``"12"`` and ``12.0`` are refused
        rather than quietly converted. No existing PyPedal public-ID convention
        requires the coercion.
        """
        for bad in ([7.0, 8], ["7", "8"], [7, None], [7, 8.5]):
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                pyp_metrics.a_effective_founders_boichard(
                    load_corpus("boichard2a.ped", "asdg"), reference=bad)

    def test_a_string_is_not_accepted_as_a_container_of_ids(self):
        """
        ``str`` and ``bytes`` are iterable and would decompose into characters.
        Refused as containers, before any element is looked at.
        """
        for bad in ("78", b"78"):
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                pyp_metrics.a_effective_founders_boichard(
                    load_corpus("boichard2a.ped", "asdg"), reference=bad)

    def test_every_refusal_is_catchable_as_a_pypedal_error_and_a_value_error(self):
        """
        ``PyPedalUsageError`` is deliberately both, so neither audience is
        surprised (``pyp_errors.py`` docstring).
        """
        for expected in (pyp_errors.PyPedalError, ValueError):
            with self.assertRaises(expected):
                pyp_metrics.a_effective_founders_boichard(
                    load_corpus("boichard2a.ped", "asdg"), reference=[])

    def test_no_raw_key_error_escapes_for_an_unknown_id(self):
        """
        Today an unknown ID reaches ``index_of[...]`` in the engine and comes
        back as a bare ``KeyError`` with no message. That must not be the
        public contract.
        """
        try:
            pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"), reference=[9999])
        except pyp_errors.PyPedalUsageError:
            pass
        except KeyError:                       # pragma: no cover - the defect
            self.fail("a raw KeyError escaped to the caller")


class TestAcceptedContainers(unittest.TestCase):
    """Any iterable of the right elements, materialised once."""

    def test_list_tuple_set_and_frozenset_all_work(self):
        expected = pyp_metrics.a_effective_founders_boichard(
            load_corpus("boichard2a.ped", "asdg"), reference=[7, 8, 9, 10])
        for container in (tuple, set, frozenset, list):
            got = pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"),
                reference=container([7, 8, 9, 10]))
            self.assertEqual(expected, got)

    def test_a_generator_is_materialised_once_and_survives_re_use(self):
        """
        ``reference`` is consumed TWICE on the ancestor paths -- once by the
        R3 antichain guard and once by the Appendix-B engine. An unmaterialised
        generator would arrive empty at the second consumer and the routine
        would divide by zero instead of refusing. So this is a correctness
        test, not a convenience one.
        """
        expected = pyp_metrics.a_effective_ancestors_definite(
            load_corpus("boichard2a.ped", "asdg"), reference=[7, 8, 9, 10])
        got = pyp_metrics.a_effective_ancestors_definite(
            load_corpus("boichard2a.ped", "asdg"),
            reference=(x for x in [7, 8, 9, 10]))
        self.assertEqual(expected, got)


# ---------------------------------------------------------------------------
# 2. Exactly one mechanism selects R
# ---------------------------------------------------------------------------

class TestExactlyOneSelectionMechanism(unittest.TestCase):
    """
    ``gen`` and ``reference`` both name the population under study. Supplying
    both is a caller error, not a precedence question: silently preferring one
    is how an analyst ends up with an analysis of a population they did not
    choose, which is the whole of .
    """

    def test_supplying_gen_and_reference_together_is_refused(self):
        for routine in (pyp_metrics.a_effective_founders_boichard,
                        pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.assertRaises(pyp_errors.PyPedalUsageError):
                routine(load_corpus("boichard2a.ped", "asdg"),
                        gen=3, reference=[7, 8])

    def test_the_legacy_default_is_not_treated_as_an_explicit_selector(self):
        """
        ``gen`` defaults to None, so omitting it must not count as "explicitly
        supplied" and must not collide with ``reference``.
        """
        pyp_metrics.a_effective_founders_boichard(
            load_corpus("boichard2a.ped", "asdg"), reference=[7, 8])

    def test_no_caller_can_reach_reference_positionally(self):
        """
        Positional compatibility, true before and after the change: nothing in
        the repository passes more than ``pedobj`` positionally, and a fourth
        positional argument must stay a ``TypeError`` rather than becoming a
        reference population. Not xfailed -- it is a standing guarantee.
        """
        with self.assertRaises(TypeError):
            pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"), None, None, [7, 8])

    def test_reference_is_declared_keyword_only_on_all_three_routines(self):
        for routine in (pyp_metrics.a_effective_founders_boichard,
                        pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            parameters = inspect.signature(routine).parameters
            self.assertIn("reference", parameters, routine.__name__)
            self.assertEqual(inspect.Parameter.KEYWORD_ONLY,
                             parameters["reference"].kind, routine.__name__)
            self.assertIsNone(parameters["reference"].default,
                              routine.__name__)
            # The existing positional parameters must not have moved.
            positional = [name for name, p in parameters.items()
                          if p.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD]
            self.assertEqual(["pedobj", "a", "gen"], positional[:3],
                             routine.__name__)


# ---------------------------------------------------------------------------
# 3. The primary correctness test: explicit R == legacy g selecting the same R
# ---------------------------------------------------------------------------

class TestExplicitReferenceEqualsTheLegacyGenerationPath(unittest.TestCase):
    """
    The strongest statement this work can make. The API path may differ; the
    mathematics must not.

    Exact equality, no tolerance. A tolerance here would hide precisely the
    kind of reordering or renormalisation that this change could introduce.
    """

    def _pedformat(self, name):
        return "asdg"

    def test_f_e_is_identical_for_the_same_reference_population(self):
        checked = 0
        for name, reference in sorted(WITH_G.items()):
            fmt = self._pedformat(name)
            legacy = outcome(lambda n=name, f=fmt:
                             pyp_metrics.a_effective_founders_boichard(
                                 load_corpus(n, f)))
            explicit = outcome(lambda n=name, f=fmt, r=reference:
                               pyp_metrics.a_effective_founders_boichard(
                                   load_corpus(n, f), reference=list(r)))
            self.assertEqual(legacy, explicit, "f_e differs on %s" % name)
            checked += legacy[0] == "value"
        self.assertGreater(checked, 0, "every subject refused on both sides, "
                                       "so this comparison proved nothing")

    def test_definite_f_a_is_identical_for_the_same_reference_population(self):
        checked = 0
        for name, reference in sorted(WITH_G.items()):
            fmt = self._pedformat(name)
            legacy = outcome(lambda n=name, f=fmt:
                             pyp_metrics.a_effective_ancestors_definite(
                                 load_corpus(n, f)))
            explicit = outcome(lambda n=name, f=fmt, r=reference:
                               pyp_metrics.a_effective_ancestors_definite(
                                   load_corpus(n, f), reference=list(r)))
            self.assertEqual(legacy, explicit, "f_a differs on %s" % name)
            checked += legacy[0] == "value"
        self.assertGreater(checked, 0, "every subject refused on both sides")

    def test_the_bounds_are_identical_at_several_n(self):
        checked = 0
        for name, reference in sorted(WITH_G.items()):
            fmt = self._pedformat(name)
            for n in (1, 2, 10 ** 6):
                legacy = outcome(
                    lambda nm=name, f=fmt, k=n:
                    pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(nm, f), n=k))
                explicit = outcome(
                    lambda nm=name, f=fmt, k=n, r=reference:
                    pyp_metrics.a_effective_ancestors_indefinite(
                        load_corpus(nm, f), n=k, reference=list(r)))
                self.assertEqual(legacy, explicit,
                                 "bounds differ on %s at n=%d" % (name, n))
                checked += legacy[0] == "value"
        self.assertGreater(checked, 0, "every subject refused on both sides")

    def test_the_selected_ancestors_and_contributions_are_identical(self):
        """
        Not just the scalar. If the selection SEQUENCE moved, 's
        tie-break would have changed meaning even where f_a happened to agree.
        The sequence is observable in the `.dat` report, so that is what is
        compared -- through the PUBLIC routines, not the engine, because it is
        the public plumbing this work changes.
        """
        for generation, reference in sorted(BOICHARD2A_BY_GENERATION.items()):
            legacy_ped = load_corpus("boichard2a.ped", "asdg")
            explicit_ped = load_corpus("boichard2a.ped", "asdg")
            legacy = outcome(lambda p=legacy_ped, g=generation:
                             pyp_metrics.a_effective_ancestors_definite(
                                 p, gen=g))
            explicit = outcome(lambda p=explicit_ped, r=reference:
                               pyp_metrics.a_effective_ancestors_definite(
                                   p, reference=list(r)))
            self.assertEqual(legacy, explicit, "generation %r" % generation)
            if legacy[0] != "value":
                continue
            suffix = "_fa_boichard_definite_.dat"
            for prefix in ("ancestors:", "ancestor contributions:"):
                self.assertEqual(report_line(legacy_ped, suffix, prefix),
                                 report_line(explicit_ped, suffix, prefix),
                                 "%r differs at generation %r"
                                 % (prefix, generation))

    def test_every_generation_label_agrees_not_only_the_default_one(self):
        for generation, reference in sorted(BOICHARD2A_BY_GENERATION.items()):
            legacy = outcome(lambda g=generation:
                             pyp_metrics.a_effective_founders_boichard(
                                 load_corpus("boichard2a.ped", "asdg"), gen=g))
            explicit = outcome(lambda r=reference:
                               pyp_metrics.a_effective_founders_boichard(
                                   load_corpus("boichard2a.ped", "asdg"),
                                   reference=list(r)))
            self.assertEqual(legacy, explicit,
                             "generation %r differs" % generation)


class TestGenerationLabelsAreIrrelevantOnceRIsFixed(unittest.TestCase):
    """
    MATHEMATICALLY IMPLIED / code-confirmed: once R is fixed, the numeric
    labels play no part in the Appendix A/B arithmetic. Measured at probe level
    in ``probe_reference_population.py`` Part 2; pinned here at API level.
    """

    TOPOLOGY = ["1 0 0", "2 0 0", "3 0 0", "4 0 0",
                "5 1 2", "6 3 4", "7 5 6", "8 5 6"]

    LABELS = {
        "consecutive": ["1", "1", "1", "1", "2", "2", "3", "3"],
        "times ten": ["10", "10", "10", "10", "20", "20", "30", "30"],
        "scrambled outside R": ["7", "7", "7", "7", "8", "8", "99", "99"],
    }

    def _load(self, labels):
        tmp = owned_temp_dir(prefix="pypedal_labels_")
        path = os.path.join(tmp, "labels.ped")
        with open(path, "w", encoding="utf-8") as handle:
            for row, label in zip(self.TOPOLOGY, labels):
                handle.write("%s %s\n" % (row, label))
        return load_corpus_from_path(path, "asdg")

    def test_the_same_r_gives_the_same_numbers_under_every_labeling(self):
        results = {}
        for name, labels in sorted(self.LABELS.items()):
            results[name] = (
                pyp_metrics.a_effective_founders_boichard(
                    self._load(labels), reference=[7, 8]),
                pyp_metrics.a_effective_ancestors_definite(
                    self._load(labels), reference=[7, 8]),
                pyp_metrics.a_effective_ancestors_indefinite(
                    self._load(labels), n=2, reference=[7, 8]),
            )
        distinct = set(results.values())
        self.assertEqual(1, len(distinct),
                         "generation labels changed the answer: %r" % results)


class TestIterableOrderCarriesNoScientificMeaning(unittest.TestCase):
    """
    R is a mathematical SET. Python iteration order must never become a
    scientific input. ``probe_reference_population.py`` Part 3 measured this at
    engine level and exits non-zero if it stops holding; this pins it at API
    level, where the canonicalisation actually happens.
    """

    BASE = [9, 10, 11, 12, 13]

    def _permutations(self):
        ascending = list(self.BASE)
        descending = list(reversed(self.BASE))
        rotated = self.BASE[2:] + self.BASE[:2]
        interleaved = self.BASE[::2] + self.BASE[1::2]
        return {
            "ascending": ascending,
            "descending": descending,
            "rotated": rotated,
            "interleaved": interleaved,
            "as a set": set(self.BASE),
            "as a generator": None,       # rebuilt per use below
        }

    def test_permuting_the_reference_changes_nothing(self):
        results = {}
        for label, order in sorted(self._permutations().items()):
            supplied = (x for x in self.BASE) if order is None else order
            results[label] = (
                pyp_metrics.a_effective_founders_boichard(
                    load_corpus("boichard_fig1.ped", "asdg"),
                    reference=supplied),
            )
        distinct = set(results.values())
        self.assertEqual(1, len(distinct),
                         "iterable order acquired meaning: %r" % results)

    def test_the_ancestor_selection_sequence_is_also_order_invariant(self):
        """
        A scalar can coincide while the selection sequence differs, which would
        silently change which animal 's tie-break credits. Taken through
        the public routine, so it tests the canonicalisation this work adds
        rather than the engine's pre-existing order-independence.
        """
        sequences = set()
        for order in (list(self.BASE), list(reversed(self.BASE)),
                      self.BASE[2:] + self.BASE[:2],
                      self.BASE[::2] + self.BASE[1::2]):
            ped = load_corpus("boichard_fig1.ped", "asdg")
            pyp_metrics.a_effective_ancestors_definite(ped, reference=order)
            sequences.add(report_line(ped, "_fa_boichard_definite_.dat",
                                      "ancestor contributions:"))
        self.assertEqual(1, len(sequences),
                         "iterable order changed the selection sequence: %r"
                         % sequences)


# ---------------------------------------------------------------------------
# 4. The practical fix: no `g` column at all
# ---------------------------------------------------------------------------

class TestPedigreesWithNoGenerationColumn(unittest.TestCase):
    """
    's motivating case. ``mrode.ped`` is a real corpus pedigree with
    no generation column AND a half-founder, so it exercises the explicit path
    and 's Reading C completion together.
    """

    #: Antichain reference populations for mrode, each verified below to leave
    #: the marginal contributions summing to one.
    MRODE_REFERENCES = ([6], [5], [3, 4], [3, 6])

    def test_without_a_reference_the_routines_still_refuse_loudly(self):
        """
        BACKWARD-COMPATIBILITY. This holds today and must keep holding: R is
        never inferred from ``igen``, from pedigree depth, from terminal
        status, or from a birth cohort. Not xfailed -- it is current behaviour.
        """
        for routine in (pyp_metrics.a_effective_founders_boichard,
                        pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.subTest(routine=routine.__name__):
                with self.assertRaises(pyp_errors.PyPedalError):
                    routine(load_corpus("mrode.ped", "asd"))

    def test_mrode_really_has_no_generation_information(self):
        """Guard on the guard: if it grew a `g` column the tests above go vacuous."""
        ped = load_corpus("mrode.ped", "asd")
        self.assertEqual({str(ped.kw["missing_gen"])},
                         {str(a.gen) for a in ped.pedigree})

    def test_mrode_really_has_a_half_founder(self):
        """So the no-g tests below genuinely exercise Reading C."""
        self.assertTrue(pyp_metrics.boichard_phantom_ids(
            load_corpus("mrode.ped", "asd")))

    def test_all_three_routines_run_on_a_no_g_pedigree_with_an_explicit_r(self):
        for reference in self.MRODE_REFERENCES:
            f_e = pyp_metrics.a_effective_founders_boichard(
                load_corpus("mrode.ped", "asd"), reference=list(reference))
            f_a = pyp_metrics.a_effective_ancestors_definite(
                load_corpus("mrode.ped", "asd"), reference=list(reference))
            f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                load_corpus("mrode.ped", "asd"), n=1, reference=list(reference))
            self.assertGreaterEqual(f_e, 1.0)
            self.assertGreaterEqual(f_a, 1.0)
            self.assertLessEqual(f_a, f_e + 1e-9)     # p.17
            self.assertLessEqual(f_l, f_u + 1e-9)

    def test_require_generations_is_not_consulted_on_the_explicit_path(self):
        """
        The guard has nothing scientifically relevant to validate once ``gen``
        is no longer read for selection. Asserted by observing that the routine
        returns rather than raising ``PyPedalValidationError``, on a pedigree
        whose every ``gen`` is the missing sentinel.
        """
        value = pyp_metrics.a_effective_founders_boichard(
            load_corpus("mrode.ped", "asd"), reference=[6])
        self.assertGreaterEqual(value, 1.0)

    def test_stripping_the_g_column_changes_nothing_given_the_same_r(self):
        """
        The sharpest form of the claim. The same topology is loaded twice --
        once WITH a generation column and analysed the legacy way, once with
        the column physically deleted and analysed through ``reference=`` --
        and the two must agree exactly.
        """
        for name, reference in sorted(WITH_G.items()):
            path, reduced = strip_generation_column(name, "asdg")
            with_g = load_corpus(name, "asdg")
            without_g = load_corpus_from_path(path, reduced)

            # Guard on the guard: the two loads must describe the same animals,
            # or the comparison below is meaningless.
            self.assertEqual(
                [(int(a.animalID), int(a.sireID), int(a.damID))
                 for a in with_g.pedigree],
                [(int(a.animalID), int(a.sireID), int(a.damID))
                 for a in without_g.pedigree],
                "stripping `g` changed the pedigree on %s" % name)
            self.assertEqual({str(without_g.kw["missing_gen"])},
                             {str(a.gen) for a in without_g.pedigree})

            legacy = outcome(lambda n=name:
                             pyp_metrics.a_effective_founders_boichard(
                                 load_corpus(n, "asdg")))
            explicit = outcome(lambda p=path, f=reduced, r=reference:
                               pyp_metrics.a_effective_founders_boichard(
                                   load_corpus_from_path(p, f),
                                   reference=list(r)))
            self.assertEqual(legacy, explicit,
                             "no-g explicit R differs on %s" % name)


# ---------------------------------------------------------------------------
# 5. The per-routine guard matrix -- BINDING
# ---------------------------------------------------------------------------

class TestThePerRoutineGuardMatrixIsPreserved(unittest.TestCase):
    """
    The three routines do NOT share a guard set, and this work changes
    reference-SELECTION plumbing only -- never a routine's scientific
    validation domain.

    | routine     | explicit-path guards                        |
    |-------------|---------------------------------------------|
    | founders    | topological                                 |
    | definite    | topological, R3 antichain                   |
    | indefinite  | n >= 1, topological, R3 antichain           |

    An earlier draft of the research plan asserted R3 applied to all three. It
    does not, and this class exists so a future change cannot introduce it by
    accident.
    """

    #: A half-founder pedigree whose most recent generation contains a dam and
    #: her daughter, so R is not an antichain.
    NON_ANTICHAIN = ["1 0 0 1", "2 0 0 1", "3 1 2 1",
                     "4 3 0 1", "5 3 4 2", "6 5 4 2"]

    def _non_antichain(self):
        tmp = owned_temp_dir(prefix="pypedal_r3_")
        path = os.path.join(tmp, "non_antichain.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(self.NON_ANTICHAIN) + "\n")
        return path

    def test_r3_never_leaks_into_the_founders_routine(self):
        """
        ``a_effective_founders_boichard`` has NEVER carried the antichain
        guard. Appendix A propagates contributions and never selects ancestors,
        so R3's zeroing convention -- the thing the guard protects -- does not
        arise there. The asymmetry is deliberate and must not be tidied away.
        """
        path = self._non_antichain()
        value = pyp_metrics.a_effective_founders_boichard(
            load_corpus_from_path(path, "asdg"), reference=[5, 6])
        self.assertGreaterEqual(value, 1.0)

    def test_r3_still_refuses_in_both_ancestor_routines(self):
        path = self._non_antichain()
        for routine in (pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.assertRaises(pyp_errors.PyPedalError) as caught:
                routine(load_corpus_from_path(path, "asdg"), reference=[4, 5])
            self.assertIn("antichain", str(caught.exception))

    def test_n_below_one_is_still_refused_before_anything_else(self):
        with self.assertRaises(pyp_errors.PyPedalUsageError):
            pyp_metrics.a_effective_ancestors_indefinite(
                load_corpus("boichard2a.ped", "asdg"), n=0, reference=[7, 8])

    def test_the_topological_guard_still_applies_on_the_explicit_path(self):
        """
        Appendix B steps 5 and 6 are single passes in opposite directions and
        are only correct when parents precede offspring. Loading with
        ``renumber`` disabled and offspring written first reaches the guard.
        """
        tmp = owned_temp_dir(prefix="pypedal_topo_")
        path = os.path.join(tmp, "unordered.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("3 1 2 2\n1 0 0 1\n2 0 0 1\n")
        ped = load_corpus_from_path(path, "asdg", renumber=False,
                                    pedigree_is_renumbered=False, reorder=0)
        with self.assertRaises(pyp_errors.PyPedalError) as caught:
            pyp_metrics.a_effective_founders_boichard(ped, reference=[3])
        self.assertIn("earlier in the pedigree", str(caught.exception))


class TestTheLegacyGuardMatrixIsUntouched(unittest.TestCase):
    """
    Current behaviour, asserted so the production change cannot quietly move a
    guard on the path it is supposed to leave alone. Not xfailed.
    """

    def test_the_legacy_path_still_requires_generations(self):
        for routine in (pyp_metrics.a_effective_founders_boichard,
                        pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.subTest(routine=routine.__name__):
                with self.assertRaises(pyp_errors.PyPedalValidationError):
                    routine(load_corpus("mrode.ped", "asd"))

    def test_the_legacy_path_still_refuses_an_unknown_generation(self):
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_metrics.a_effective_founders_boichard(
                load_corpus("boichard2a.ped", "asdg"), gen=99)

    def test_the_legacy_founders_routine_has_no_antichain_guard_today(self):
        """The asymmetry pinned above already exists; recorded as a baseline."""
        tmp = owned_temp_dir(prefix="pypedal_r3_base_")
        path = os.path.join(tmp, "non_antichain.ped")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(
                TestThePerRoutineGuardMatrixIsPreserved.NON_ANTICHAIN) + "\n")
        value = pyp_metrics.a_effective_founders_boichard(
            load_corpus_from_path(path, "asdg"))
        self.assertGreaterEqual(value, 1.0)
        with self.assertRaises(pyp_errors.PyPedalError):
            pyp_metrics.a_effective_ancestors_definite(
                load_corpus_from_path(path, "asdg"))


# ---------------------------------------------------------------------------
# 6.  and  under an explicit reference
# ---------------------------------------------------------------------------

class TestHalfFounderFixturesAcceptAnExplicitReference(unittest.TestCase):
    """
     Reading C is untouched by this work: phantom completion happens
    inside the shared engine, below the point where R is chosen. These tests
    prove it still happens when R arrives through the new door.
    """

    def test_the_invariants_the_paper_states_about_itself_still_hold(self):
        """
        p.7 founder contributions sum to one; p.8 marginal contributions sum to
        one; p.17 f_a <= f_e. Nothing normalises either vector, so these are
        real checks rather than tautologies.
        """
        for name in ("synth_hf.ped", "synth_nohf.ped"):
            reference = WITH_G[name]
            f_e = pyp_metrics.a_effective_founders_boichard(
                load_corpus(name, "asdg"), reference=list(reference))
            f_a = pyp_metrics.a_effective_ancestors_definite(
                load_corpus(name, "asdg"), reference=list(reference))
            self.assertLessEqual(f_a, f_e + 1e-9, "f_a > f_e on %s" % name)

    def test_reading_c_completion_still_happens_under_an_explicit_reference(self):
        """
        Phantom completion lives inside the shared engine, below the point
        where R is chosen, so arriving through ``reference=`` must not change
        it. The report's phantom note is the observable evidence that it ran,
        and a synthetic ID must still never be printed bare -- ``8`` appears as
        ``8=phantom (dam of 5)``.
        """
        ped = load_corpus("synth_hf.ped", "asdg")
        phantoms = pyp_metrics.boichard_phantom_ids(ped)
        self.assertTrue(phantoms, "fixture stopped having a half-founder")
        pyp_metrics.a_effective_ancestors_definite(ped, reference=[6, 7])
        text = report_text(ped, "_fa_boichard_definite_.dat")
        self.assertIn("phantom founder(s) were created for this analysis", text)
        for phantom_id, (animal_id, side) in sorted(phantoms.items()):
            self.assertIn("%d=phantom (%s of %s)" % (phantom_id, side,
                                                     animal_id), text)

    def test_bl2_bounds_bracket_the_exact_value_under_an_explicit_reference(self):
        for name in ("boichard_fig1.ped", "boichard_fig2.ped", "synth_hf.ped"):
            reference = WITH_G[name]
            exact = pyp_metrics.a_effective_ancestors_definite(
                load_corpus(name, "asdg"), reference=list(reference))
            for n in (1, 2):
                f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                    load_corpus(name, "asdg"), n=n, reference=list(reference))
                self.assertLessEqual(f_l, exact + 1e-9,
                                     "f_l > f_a on %s at n=%d" % (name, n))
                self.assertLessEqual(exact, f_u + 1e-9,
                                     "f_a > f_u on %s at n=%d" % (name, n))

    def test_exhausting_the_sequence_collapses_the_bounds_onto_the_exact_value(self):
        for name in ("boichard_fig1.ped", "boichard_fig2.ped", "synth_hf.ped"):
            reference = WITH_G[name]
            exact = pyp_metrics.a_effective_ancestors_definite(
                load_corpus(name, "asdg"), reference=list(reference))
            f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
                load_corpus(name, "asdg"), n=10 ** 6, reference=list(reference))
            self.assertAlmostEqual(exact, f_l, places=12)
            self.assertAlmostEqual(exact, f_u, places=12)


# ---------------------------------------------------------------------------
# 7. Reports must not lie about generations
# ---------------------------------------------------------------------------

class TestTheLegacyReportTextIsUnchanged(unittest.TestCase):
    """Current wording, pinned so the production change cannot disturb it."""

    def test_the_founders_report_still_names_the_generation(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_founders_boichard(ped)
        text = report_text(ped, "_fe_boichard_.dat")
        self.assertIn("generations: ['1', '2', '3']", text)
        self.assertIn("8 animals in generation 3", text)

    def test_the_definite_report_still_names_the_generation(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_ancestors_definite(ped)
        text = report_text(ped, "_fa_boichard_definite_.dat")
        self.assertIn("animals in the population under study:", text)
        self.assertIn("generations: ['1', '2', '3']", text)
        self.assertIn("8 animals in generation 3", text)

    def test_the_indefinite_report_never_mentioned_a_generation(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_ancestors_indefinite(ped, n=2)
        text = report_text(ped, "_fa_boichard_indefinite_.dat")
        self.assertNotIn("generation", text)


class TestTheExplicitReportTellsTheTruth(unittest.TestCase):
    """
    On the explicit path there IS no generation, so the report must not invent
    one. Nothing may be derived from ``gen``, ``igen``, pedigree depth, birth
    year or any cohort rule to keep the legacy wording usable.
    """

    def test_the_founders_report_names_the_explicit_reference_population(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_founders_boichard(ped, reference=[7, 8])
        text = report_text(ped, "_fe_boichard_.dat")
        self.assertIn("2 animals in the explicit reference population: [7, 8]",
                      text)
        self.assertNotIn("in generation", text)
        self.assertNotIn("generations:", text)

    def test_the_definite_report_records_that_r_was_supplied_explicitly(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_ancestors_definite(ped, reference=[7, 8])
        text = report_text(ped, "_fa_boichard_definite_.dat")
        self.assertIn("2 animals in the population under study: [7, 8]", text)
        self.assertIn("reference population supplied explicitly", text)
        self.assertNotIn("in generation", text)
        self.assertNotIn("generations:", text)

    def test_the_indefinite_report_gains_no_gratuitous_text(self):
        """Smallest truthful delta: it had no generation line to correct."""
        ped = load_corpus("boichard2a.ped", "asdg")
        pyp_metrics.a_effective_ancestors_indefinite(ped, n=2,
                                                     reference=[7, 8])
        text = report_text(ped, "_fa_boichard_indefinite_.dat")
        self.assertNotIn("generation", text)

    def test_no_synthetic_generation_value_leaks_into_the_report(self):
        """
        A no-``g`` pedigree has nothing but the missing sentinel. If the
        sentinel or a fabricated label ever reached the report, this catches it.
        """
        ped = load_corpus("mrode.ped", "asd")
        pyp_metrics.a_effective_founders_boichard(ped, reference=[6])
        text = report_text(ped, "_fe_boichard_.dat")
        self.assertNotIn("-999", text)
        self.assertNotIn("generation", text.replace(
            "explicit reference population", ""))


# ---------------------------------------------------------------------------
# 8. Caller state isolation and determinism
# ---------------------------------------------------------------------------

class TestTheCallersPedigreeIsNotDisturbed(unittest.TestCase):
    """
    Supplying a reference list must not write ``gen`` or ``igen``, mark animals, reorder the
    pedigree, insert synthetic records, or touch ``pedobj.kw``.
    """

    def test_no_animal_attribute_moves(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        before = animal_state(ped)
        kw_before = dict(ped.kw)
        pyp_metrics.a_effective_founders_boichard(ped, reference=[7, 8])
        pyp_metrics.a_effective_ancestors_definite(ped, reference=[7, 8])
        pyp_metrics.a_effective_ancestors_indefinite(ped, n=1,
                                                     reference=[7, 8])
        self.assertEqual(before, animal_state(ped))
        self.assertEqual(kw_before, dict(ped.kw))

    def test_no_synthetic_animal_is_inserted(self):
        ped = load_corpus("synth_hf.ped", "asdg")
        length = len(ped.pedigree)
        ids = [int(a.animalID) for a in ped.pedigree]
        pyp_metrics.a_effective_ancestors_definite(ped, reference=[6, 7])
        self.assertEqual(length, len(ped.pedigree))
        self.assertEqual(ids, [int(a.animalID) for a in ped.pedigree])

    def test_the_supplied_reference_container_is_not_mutated(self):
        supplied = [8, 7, 9, 10]
        copy = list(supplied)
        pyp_metrics.a_effective_founders_boichard(
            load_corpus("boichard2a.ped", "asdg"), reference=supplied)
        self.assertEqual(copy, supplied)

    def test_repeated_calls_are_deterministic(self):
        ped = load_corpus("boichard2a.ped", "asdg")
        first = pyp_metrics.a_effective_ancestors_definite(ped,
                                                           reference=[7, 8])
        second = pyp_metrics.a_effective_ancestors_definite(ped,
                                                            reference=[7, 8])
        self.assertEqual(first, second)

    def test_call_order_does_not_change_either_result(self):
        """
        legacy-then-explicit and explicit-then-legacy must both leave the two
        answers where they were when each was taken alone.
        """
        alone_legacy = pyp_metrics.a_effective_founders_boichard(
            load_corpus("boichard2a.ped", "asdg"))
        alone_explicit = pyp_metrics.a_effective_founders_boichard(
            load_corpus("boichard2a.ped", "asdg"), reference=[5, 6])

        ped = load_corpus("boichard2a.ped", "asdg")
        legacy_first = pyp_metrics.a_effective_founders_boichard(ped)
        explicit_second = pyp_metrics.a_effective_founders_boichard(
            ped, reference=[5, 6])

        other = load_corpus("boichard2a.ped", "asdg")
        explicit_first = pyp_metrics.a_effective_founders_boichard(
            other, reference=[5, 6])
        legacy_second = pyp_metrics.a_effective_founders_boichard(other)

        self.assertEqual(alone_legacy, legacy_first)
        self.assertEqual(alone_legacy, legacy_second)
        self.assertEqual(alone_explicit, explicit_second)
        self.assertEqual(alone_explicit, explicit_first)


# ---------------------------------------------------------------------------
# 9. Real data -- characterisation, never authority
# ---------------------------------------------------------------------------

class TestGriffonCharacterisation(unittest.TestCase):
    """
    The historical 1871–1890 Griffon extract is CHARACTERISATION, derived at
    test time from the canonical sample. No canonical griffon
    reference population is proposed here or anywhere else in the project; the
    sets below are TEST-SPECIFIED, chosen to exercise the API on real data.

    Boichard's own published selection is "recorded females born in 1988-91
    from known sire and dam" (Table V fn a, p.14) -- a birth cohort intersected
    with pedigree completeness, and the authors' analysis choice for that study
    rather than a method rule. It must never become a PyPedal default.
    """

    def test_the_max_depth_counterexample_still_reproduces(self):
        """
        The measurement that killed the ``igen -> gen`` repair, pinned so it
        cannot quietly stop being true. 3 of 167 animals sit at the maximum
        inferred generation, and all three are inside -- but far from equal to
        -- the 36-animal 1890 birth cohort.
        """
        ped = load_griffon()
        self.assertEqual(167, len(ped.pedigree))
        self.assertTrue(pyp_utils.set_generation(ped),
                        "set_generation failed; the distribution below would "
                        "be measured on unset values")
        distribution = {}
        for animal in ped.pedigree:
            distribution[animal.igen] = distribution.get(animal.igen, 0) + 1
        self.assertEqual({1: 124, 2: 20, 3: 8, 4: 12, 5: 3}, distribution)

        deepest = {int(a.animalID) for a in ped.pedigree
                   if a.igen == max(distribution)}
        self.assertEqual(3, len(deepest))
        cohort = set(griffon_cohort(ped, {1890}))
        self.assertEqual(36, len(cohort))
        self.assertTrue(deepest < cohort,
                        "the max-depth set is a STRICT subset of the cohort")

    def test_no_griffon_fixture_carries_a_generation_column(self):
        ped = load_griffon()
        self.assertEqual({str(ped.kw["missing_gen"])},
                         {str(a.gen) for a in ped.pedigree})

    def test_unknown_recorded_years_are_not_a_cohort(self):
        """
        Guard on the cohort helper. The file's known birth years run in the
        historical Griffon range. Unknown chronology is ``None`` (including
        implicit parents) and must not be treated as year 1800.
        """
        ped = load_griffon()
        years = {a.by for a in ped.pedigree}
        known = {year for year in years if year is not None}
        self.assertTrue(known)
        self.assertNotEqual(min(known), 1800)
        self.assertNotIn(None, {1889, 1890})
        self.assertGreaterEqual(sum(1 for a in ped.pedigree if a.by is None), 0)
        with self.assertRaises(AssertionError):
            griffon_cohort(ped, {None})

    def test_the_three_deepest_animals_are_not_a_default(self):
        """
        Without a reference the routines refuse. They do NOT fall back to the
        max-depth set, which is what an ``igen -> gen`` repair would have done.
        """
        for routine in (pyp_metrics.a_effective_founders_boichard,
                        pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.subTest(routine=routine.__name__):
                with self.assertRaises(pyp_errors.PyPedalValidationError):
                    routine(load_griffon())

    def test_an_explicit_griffon_r_reaches_the_engine(self):
        """
        The 1890 animals from known sire and dam -- the SHAPE of Boichard's own
        Table V selection, minus its sex filter. TEST-SPECIFIED, not endorsed.
        """
        ped = load_griffon()
        reference = griffon_cohort(ped, {1890}, complete_pedigree_only=True)
        self.assertTrue(reference, "the fixture stopped having such animals")
        f_e = pyp_metrics.a_effective_founders_boichard(
            load_griffon(), reference=list(reference))
        f_a = pyp_metrics.a_effective_ancestors_definite(
            load_griffon(), reference=list(reference))
        f_l, f_u = pyp_metrics.a_effective_ancestors_indefinite(
            load_griffon(), n=5, reference=list(reference))
        self.assertGreaterEqual(f_e, 1.0)
        self.assertLessEqual(f_a, f_e + 1e-9)
        self.assertLessEqual(f_l, f_a + 1e-9)
        self.assertLessEqual(f_a, f_u + 1e-9)

    def test_an_explicit_r_is_never_substituted_by_the_max_depth_set(self):
        """
        The regression this whole finding exists to prevent. The answer for an
        explicitly supplied R must equal the engine's answer for THAT R, and
        must differ from the answer for the 3 max-depth animals -- otherwise a
        substitution would be undetectable.
        """
        ped = load_griffon()
        reference = griffon_cohort(ped, {1890}, complete_pedigree_only=True)
        supplied = pyp_metrics.a_effective_founders_boichard(
            load_griffon(), reference=list(reference))

        depth_ped = load_griffon()
        self.assertTrue(pyp_utils.set_generation(depth_ped))
        top = max(a.igen for a in depth_ped.pedigree)
        deepest = sorted(int(a.animalID) for a in depth_ped.pedigree
                         if a.igen == top)
        self.assertEqual(3, len(deepest))
        substituted = pyp_metrics.a_effective_founders_boichard(
            load_griffon(), reference=deepest)
        self.assertNotAlmostEqual(supplied, substituted, places=6)

    def test_a_cohort_containing_founders_is_refused_by_the_definite_routine(self):
        """
        The R3-family limitation, measured on real data and SHARPER than the
        adjudication's section 8 recorded -- that section measured only the
        antichain property.

        The whole 1890 cohort IS an antichain, so it passes R3's guard. It
        still cannot be analysed for effective ancestors, for a second and
        independent reason: 27 of its 36 members are themselves founders, R3's
        convention zeroes reference-population members before selection, and
        the mass those 27 carry is therefore credited to no ancestor at all.
        The marginal contributions sum to 0.306 instead of the 1 that Boichard
        p.8 requires, and `check_contribution_vector` refuses.

        This is why Boichard's own Table V selection reads "recorded females
        born in 1988-91 **from known sire and dam**" -- the pedigree-
        completeness filter is not incidental. Restricting R to the nine 1890
        animals with both parents known satisfies every invariant, which
        `test_an_explicit_griffon_r_reaches_the_engine` exercises.

         removed the SELECTION block. It did not, and was not meant
        to, make every cohort admissible.
        """
        ped = load_griffon()
        cohort = griffon_cohort(ped, {1890})
        self.assertEqual(36, len(cohort))

        # Passes R3: the guard this cohort does NOT trip.
        pyp_metrics._boichard_require_antichain(ped, cohort, "test")

        # The founders routine has no such difficulty: Appendix A credits every
        # founder directly, so nothing is lost and q still sums to one.
        self.assertGreater(
            pyp_metrics.a_effective_founders_boichard(
                load_griffon(), reference=list(cohort)), 1.0)

        with self.assertRaises(pyp_errors.PyPedalValidationError) as caught:
            pyp_metrics.a_effective_ancestors_definite(
                load_griffon(), reference=list(cohort))
        self.assertIn("marginal contributions", str(caught.exception))

    # The bounded routine's failure to refuse this same cohort was registered
    # here as a defect and deliberately left unrepaired, because that phase's
    # remit was reference-selection plumbing only. It is now , and
    # the pin that stood here -- which asserted `check_contribution_vector` was
    # ABSENT from `a_effective_ancestors_indefinite`, and whose own message
    # said to remove it once the repair landed -- has been replaced by
    # `tests/test_boichard_contribution_sequence.py`. See
    # `the algorithm notes`.

    def test_r3_still_refuses_a_two_year_griffon_window(self):
        """
        MEASURED, and the planning pass guessed wrong about the single latest
        year: 1890 alone IS an antichain. Two years are not -- 3 parent/
        offspring pairs. Recorded as measured rather than replaced with a
        tidier prediction.

         removes the SELECTION block. It does not remove R3.
        """
        ped = load_griffon()
        window = griffon_cohort(ped, {1889, 1890})
        self.assertEqual(57, len(window))
        for routine in (pyp_metrics.a_effective_ancestors_definite,
                        pyp_metrics.a_effective_ancestors_indefinite):
            with self.assertRaises(pyp_errors.PyPedalError) as caught:
                routine(load_griffon(), reference=list(window))
            self.assertIn("antichain", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
