"""
 -- the contract of ``pyp_utils.set_age`` and the ``missing_value`` option.

WHAT THIS FILE IS ABOUT
-----------------------
``pyp_utils.set_age`` reads ``pedobj.kw['missing_value']`` three times
(``pyp_utils.py:195-199``). PyPedal 2.0.4 established that option with
setdefault semantics; the first Python 3 conversion preserved it; commit
``bf0d2b2`` dropped it when the flat ``kw.setdefault`` run was extracted into
``NewPedigree._set_animal_defaults``. Fourteen of the sixteen deleted keys were
re-established in the new helper and one moved to ``__init__``. Two were not:
``update_sexes``, which has no readers, and ``missing_value``, which has three.

So every call raised ``KeyError('missing_value')`` on the first animal, the
routine's broad handler swallowed it, and the caller got ``False`` into a return
value nobody reads. Ages stayed at their initialisation sentinel and
``pyp_demog.age_distribution`` histogrammed that sentinel.

The repair restores the one dropped default. ``set_age``'s executable logic, its
broad handler and ``pyp_demog`` are unchanged.

HOW THIS FILE WAS BUILT
-----------------------
It landed one commit *ahead* of the production change, with every test asserting
the repaired behaviour marked ``xfail(strict=True)``, so the contract was written
down before any code satisfied it and nothing could pass by accident. Those
markers were removed in the commit that implemented the repair.

The characterisation tests that pinned the *broken* behaviour were kept rather
than deleted, restated as assertions that the defect no longer reproduces
(``TestTheReproducerNoLongerReproduces``). A reproducer thrown away once it stops
reproducing leaves nothing to stop the defect returning.

``self.subTest`` was deliberately kept out of the tests that were xfailed: a
strict xfail must fail as ONE outcome, and pytest 8 and pytest 9 handle the
subTest/xfail interaction differently -- a difference this project has already
been bitten by (``docs/-R2-HALF-FOUNDER-VERIFICATION.md`` section 2.4). Those
tests still use plain loops.

THE ORACLES
-----------
Two of them, and neither is PyPedal 2.0.4's ``age_distribution``.

*For ``set_age`` itself*, the oracle is **PyPedal 2.0.4**, measured in the
pinned ``pypedal-py27:audit`` container against the read-only legacy checkout.
With the option supplied and the base year normalised, the two implementations
agree bit-exactly on the return value, the computed ages, both sentinel
branches, the custom-marker path and idempotence. See
``the independent oracle``.

*For ``age_distribution``*, 2.0.4 **cannot** serve as the oracle: its own
``age_distribution`` passes ``pedobj.pedigree`` -- a list -- to
``pyp_utils.set_generation``, whose bare ``except:`` handler then dereferences
``.kw`` on that list and raises an uncaught ``AttributeError``. 2.0.4 never
reaches ``set_age`` through that path at all. The Python 3 port already fixed
the argument. The oracle for AGE-7 is therefore the port's own arithmetic
identity -- ``age == by - pyp_demog.BASE_DEMOGRAPHIC_YEAR`` -- and the claim
being made is that restoring the option lets the already-repaired demographic
call path consume populated ages. It is NOT a claim that Python 3
``age_distribution`` equals Python 2 ``age_distribution``.

WHAT IS *NOT* CLAIMED HERE
--------------------------
 pins that a pedigree with no birth-year column ends up with every age
equal to 0. That is a **legacy-compatibility pin, not a scientific invariant**
-- see the class docstring.  owns the birth-year sentinel model and
may deliberately supersede it.
"""
import contextlib
import io
import os
import tempfile
import unittest
import warnings

from _pedhelpers import chdir_tmp, corpus, load_corpus, load_corpus_from_path

from PyPedal import pyp_demog, pyp_newclasses, pyp_utils


# ---------------------------------------------------------------------------
# Fixtures and oracles
#
# hartlandclark.ped carries a real birth-year column (pedformat 'asdb', years
# 1900-1970) and is the known-birth-year fixture. new_lacy.ped has no birth
# column at all and is the unknown-birth-year fixture. Both are already in the
# differential corpus, so the same bytes drive the Python 2 comparison.
# ---------------------------------------------------------------------------

# Measured on hartlandclark.ped: by - BASE_DEMOGRAPHIC_YEAR, base year 1800.
HARTL_BY = [1900, 1900, 1900, 1910, 1910, 1920, 1920,
            1930, 1930, 1930, 1930, 1940, 1950, 1960, 1970]
HARTL_AGES = [100, 100, 100, 110, 110, 120, 120,
              130, 130, 130, 130, 140, 150, 160, 170]

# The eight buckets the restored histogram must show, as {age: count}.
HARTL_BUCKETS = {100: 3, 110: 2, 120: 2, 130: 4, 140: 1, 150: 1, 160: 1, 170: 1}

# The single bucket the BROKEN histogram shows: every animal at missing_age.
BROKEN_BUCKETS = {-999: 15}

MIXED_ROWS = ["1 0 0 1950", "2 0 0 0", "3 1 2 1970", "4 1 2 0", "5 3 4 1990"]
MIXED_BY = [1950, None, 1970, None, 1990]
MIXED_AGES = [150, -999, 170, -999, 190]

# 2.0.4's default, and the exact value the repair must restore.
LEGACY_MISSING_VALUE = -999.0


def hartl(**overrides):
    return load_corpus("hartlandclark.ped", **overrides)


def lacy(**overrides):
    return load_corpus("new_lacy.ped", **overrides)


def rows_to_ped(rows, pedformat="asdb", **overrides):
    """
    Load an inline pedigree with every generated file confined to a tmpdir.

    Written through ``load_corpus_from_path`` rather than by hand so the
    repository-delta guard in ``conftest.py`` stays satisfied.
    """
    tmp = tempfile.mkdtemp(prefix="age33_")
    path = os.path.join(tmp, "age33.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    return load_corpus_from_path(path, pedformat, **overrides)


def histogram(text, section):
    """
    Parse one ``age_distribution`` section into ``{age: count}``.

    The routine prints rather than returns, so the printed table is the only
    observable it has. Rows look like ``\\tAGE\\tCOUNT\\tFREQ\\tHIST`` and the
    section ends at the ``TOTAL`` row.
    """
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() != section:
            continue
        out = {}
        for row in lines[i + 2:]:
            fields = row.split("\t")
            if len(fields) < 3 or not fields[0] == "":
                break
            if fields[1] == "TOTAL":
                break
            out[int(float(fields[1]))] = int(fields[2])
        return out
    raise AssertionError("section %r not found in age_distribution output" % section)


def run_age_distribution(pedobj, sex=1):
    """Call age_distribution and return everything it printed."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pyp_demog.age_distribution(pedobj, sex=sex)
    return buf.getvalue()


def mutable_state(pedobj, skip=("age",)):
    """A comparable snapshot of every animal attribute except the named ones."""
    return [
        {k: repr(v) for k, v in sorted(a.__dict__.items()) if k not in skip}
        for a in pedobj.pedigree
    ]


# ===========================================================================
# THE REPRODUCER, INVERTED.
#
# On the phase baseline 2f8ceb4 each of these asserted the BROKEN behaviour and
# passed. They are kept rather than deleted -- a reproducer thrown away once it
# stops reproducing leaves nothing to stop the defect returning -- and restated
# as the assertion that it no longer reproduces.
# ===========================================================================
class TestTheReproducerNoLongerReproduces(unittest.TestCase):
    """The Finding-33 reproducer, now inverted. Baseline: 2f8ceb4."""

    def test_the_option_is_present_on_a_freshly_loaded_pedigree(self):
        """Was: assertNotIn -- the option was absent from every pedigree."""
        self.assertIn("missing_value", hartl().kw)

    def test_set_age_reports_success_on_a_pedigree_with_real_birth_years(self):
        """Was: assertFalse -- set_age failed on every call."""
        ped = hartl()
        self.assertTrue(pyp_utils.set_age(ped))

    def test_ages_are_no_longer_left_at_the_initialisation_sentinel(self):
        """Was: every age stayed at kw['missing_age']."""
        ped = hartl()
        pyp_utils.set_age(ped)
        ages = [a.age for a in ped.pedigree]
        self.assertNotEqual(ages, [ped.kw["missing_age"]] * 15)
        self.assertEqual(ages, HARTL_AGES)

    def test_the_comparison_that_used_to_raise_now_succeeds(self):
        """
        The cause, isolated from the handler that hid it: the same comparison
        set_age performs, run outside the try block. Was: assertRaises(KeyError)
        with args[0] == 'missing_value'.
        """
        ped = hartl()
        self.assertFalse(ped.pedigree[0].by == ped.kw["missing_value"])

    def test_age_distribution_no_longer_histograms_the_sentinel(self):
        """Was: eight real buckets collapsed into the single BROKEN_BUCKETS."""
        ped = hartl(messages="verbose")
        buckets = histogram(run_age_distribution(ped), "Unknowns")
        self.assertNotEqual(buckets, BROKEN_BUCKETS)
        self.assertEqual(buckets, HARTL_BUCKETS)


# ===========================================================================
# PERMANENT -- true now, and required to stay true after the repair.
#
# These are the guard against a repair that assigns the default
# unconditionally instead of through setdefault.
# ===========================================================================
class TestAGE4CallerSuppliedValueIsPreserved(unittest.TestCase):
    """AGE-4 -- a caller-supplied missing_value survives into kw and is used."""

    def test_a_custom_value_reaches_kw_untouched(self):
        ped = hartl(missing_value=-1234.5)
        self.assertEqual(ped.kw["missing_value"], -1234.5)

    def test_a_custom_value_makes_set_age_succeed_today(self):
        """
        The positive evidence that only the DEFAULT is missing. The option is
        read, not computed, so supplying it is already a working workaround --
        which is why the repair belongs in the defaults block and not in
        set_age.
        """
        ped = hartl(missing_value=-1234.5)
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], HARTL_AGES)

    def test_a_custom_value_is_not_overwritten_by_the_defaults_pass(self):
        """
        Stated separately from the two above because this is the specific
        property a `kw['missing_value'] = -999.0` repair would break while
        still making the symptom disappear.
        """
        for supplied in (-1234.5, -1.0, 0.0):
            ped = hartl(missing_value=supplied)
            self.assertEqual(ped.kw["missing_value"], supplied)


class TestTheFixturesStillDiscriminate(unittest.TestCase):
    """
    Fixture guards. True on the baseline and after the repair, so they carry no
    xfail marker -- a fixture that quietly stopped discriminating would make
    the /AGE-5 assertions vacuous without failing anything.
    """

    def test_new_lacy_really_has_no_birth_information(self):
        ped = lacy()
        self.assertNotIn("b", ped.kw["pedformat"])
        self.assertNotIn("y", ped.kw["pedformat"])
        self.assertEqual([a.by for a in ped.pedigree], [None] * 7)

    def test_the_mixed_fixture_actually_mixes_the_two_cases(self):
        ped = rows_to_ped(MIXED_ROWS)
        known = [a for a in ped.pedigree if a.by is not None]
        unknown = [a for a in ped.pedigree if a.by is None]
        self.assertEqual(len(known), 3)
        self.assertEqual(len(unknown), 2)

    def test_the_base_year_is_not_a_missing_birth_year(self):
        self.assertEqual(pyp_demog.BASE_DEMOGRAPHIC_YEAR, 1800)
        self.assertNotIn("missing_byear", lacy().kw)
        self.assertTrue(all(a.by is None for a in lacy().pedigree))


class TestAGE8OptionFileRoundTrip(unittest.TestCase):
    """AGE-8 -- missing_value supplied through an .ini file."""

    def _ini_pedigree(self, value):
        with chdir_tmp() as tmp:
            with open(os.path.join(tmp, "hartlandclark.ped"), "wb") as out:
                with open(corpus("hartlandclark.ped"), "rb") as src:
                    out.write(src.read())
            with open(os.path.join(tmp, "age33.ini"), "w", encoding="utf-8") as handle:
                handle.write(
                    "pedfile = hartlandclark.ped\n"
                    "pedformat = asdb\n"
                    "pedname = age33\n"
                    "messages = quiet\n"
                    "missing_value = %s\n" % value)
            return pyp_newclasses.load_pedigree(optionsfile="age33.ini")

    def test_an_ini_value_is_coerced_to_float_and_survives(self):
        """
        ``pyp_io.read_ini_file`` has no key allowlist and runs at
        ``NewPedigree.__init__:45``, before ``_set_animal_defaults`` at :59,
        so an option-file value precedes every setdefault.
        """
        ped = self._ini_pedigree("-1234.5")
        self.assertEqual(ped.kw["missing_value"], -1234.5)
        self.assertIsInstance(ped.kw["missing_value"], float)

    def test_an_ini_value_is_consumed_by_set_age(self):
        ped = self._ini_pedigree("-1234.5")
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], HARTL_AGES)


# ===========================================================================
# STRICT XFAIL -- the contract the repair must deliver.
# Remove the markers when it lands. Do not delete the tests.
# ===========================================================================
class TestAGE1TheDefaultExists(unittest.TestCase):
    """AGE-1 -- the option must be established by the canonical defaults pass."""

    def test_the_default_is_established_with_the_legacy_value_and_type(self):
        ped = hartl()
        self.assertIn("missing_value", ped.kw)
        self.assertEqual(ped.kw["missing_value"], LEGACY_MISSING_VALUE)
        # 2.0.4 wrote `-999.`, a float. missing_age is the int -999; the two
        # keys are not interchangeable and the type is part of the contract.
        self.assertIsInstance(ped.kw["missing_value"], float)

    def test_the_default_is_established_for_a_pedigree_with_no_birth_column(self):
        self.assertEqual(lacy().kw["missing_value"], LEGACY_MISSING_VALUE)


class TestAGE2KnownBirthYears(unittest.TestCase):
    """AGE-2 -- known birth years produce ages, matching PyPedal 2.0.4."""

    def test_set_age_reports_success(self):
        self.assertTrue(pyp_utils.set_age(hartl()))

    def test_ages_match_the_measured_legacy_result(self):
        ped = hartl()
        pyp_utils.set_age(ped)
        self.assertEqual([a.by for a in ped.pedigree], HARTL_BY)
        self.assertEqual([a.age for a in ped.pedigree], HARTL_AGES)

    def test_ages_follow_the_base_year_identity(self):
        """
        Stated as an identity as well as a literal list, so the test still
        means something if pyp_demog.BASE_DEMOGRAPHIC_YEAR is ever changed
        deliberately rather than silently.
        """
        ped = hartl()
        pyp_utils.set_age(ped)
        for animal in ped.pedigree:
            self.assertEqual(animal.age,
                             animal.by - pyp_demog.BASE_DEMOGRAPHIC_YEAR)


class TestAGE3aUnknownBirthYearCompatibility(unittest.TestCase):
    """
     -- unknown recorded birth year uses the missing-age marker.

     / RC4 an earlier revision: unknown is ``None``, not 1800, so ``set_age``
    no longer produces a fake year-offset of 0.
    """

    def test_unknown_birth_year_uses_missing_age_not_zero(self):
        ped = lacy()
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], [ped.kw["missing_age"]] * 7)
        self.assertNotEqual([a.age for a in ped.pedigree], [0] * 7)


class TestAGE3bTheSentinelBranches(unittest.TestCase):
    """
     -- unknown recorded year uses missing age, never ``igen``.
    """

    def test_unknown_birth_year_gives_the_missing_age_marker(self):
        ped = lacy()
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], [ped.kw["missing_age"]] * 7)

    def test_unknown_birth_year_does_not_fall_back_to_generation(self):
        ped = lacy()
        pyp_utils.set_generation(ped)
        expected_igen = [a.igen for a in ped.pedigree]
        self.assertNotIn(ped.kw["missing_age"], expected_igen)
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], [ped.kw["missing_age"]] * 7)
        self.assertNotEqual([a.age for a in ped.pedigree], expected_igen)


class TestAGE5MixedPedigree(unittest.TestCase):
    """AGE-5 -- known and unknown birth years in one pedigree, one pass."""

    def test_each_animal_gets_its_own_answer(self):
        ped = rows_to_ped(MIXED_ROWS)
        self.assertEqual([a.by for a in ped.pedigree], MIXED_BY)
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], MIXED_AGES)


class TestAGE6StateAndMutationContract(unittest.TestCase):
    """AGE-6 -- what set_age is allowed to touch, and how it repeats."""

    def test_repeated_calls_are_idempotent(self):
        ped = hartl()
        self.assertTrue(pyp_utils.set_age(ped))
        first = [a.age for a in ped.pedigree]
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], first)

    def test_no_attribute_other_than_age_is_changed(self):
        ped = hartl()
        before = mutable_state(ped)
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual(mutable_state(ped), before)

    def test_two_pedigree_objects_do_not_interfere(self):
        """No process-global state: a custom marker on one must not leak."""
        a = hartl()
        b = hartl(missing_value=-4242.0)
        self.assertTrue(pyp_utils.set_age(b))
        self.assertTrue(pyp_utils.set_age(a))
        self.assertEqual([x.age for x in a.pedigree], HARTL_AGES)
        self.assertEqual(a.kw["missing_value"], LEGACY_MISSING_VALUE)
        self.assertEqual(b.kw["missing_value"], -4242.0)


class TestAGE7AgeDistributionIntegration(unittest.TestCase):
    """
    AGE-7 -- the demographic consumer, against the PORT's contract.

    NOT a Python 2 differential. PyPedal 2.0.4's own ``age_distribution`` is
    fatally broken before it reaches ``set_age``: it passes ``pedobj.pedigree``
    -- a list -- to ``set_generation``, whose bare handler then dereferences
    ``.kw`` on that list and raises an uncaught ``AttributeError``. The port
    already fixed the argument. What is claimed here is only that restoring the
    option lets the already-repaired call path consume populated ages.
    """

    def test_the_histogram_shows_real_ages(self):
        ped = hartl(messages="verbose")
        out = run_age_distribution(ped)
        self.assertEqual(histogram(out, "Unknowns"), HARTL_BUCKETS)

    def test_the_histogram_agrees_with_the_ages_on_the_animals(self):
        """
        The port's own oracle: whatever the buckets are, they must be the
        tally of the age attribute the same call assigned.

        The non-vacuity assertion is load-bearing. Without it this test passes
        on the BROKEN baseline, where the one-bucket histogram faithfully
        tallies fifteen unassigned sentinels -- agreement with nothing.
        """
        ped = hartl(messages="verbose")
        out = run_age_distribution(ped)
        tally = {}
        for animal in ped.pedigree:
            tally[animal.age] = tally.get(animal.age, 0) + 1
        self.assertNotEqual(tally, {ped.kw["missing_age"]: 15})
        self.assertEqual(histogram(out, "Unknowns"), tally)

    def test_age_distribution_populates_ages_through_its_own_call_path(self):
        """age_distribution discards set_age's return value; check the effect."""
        ped = hartl(messages="verbose")
        self.assertEqual([a.age for a in ped.pedigree], [ped.kw["missing_age"]] * 15)
        run_age_distribution(ped)
        self.assertEqual([a.age for a in ped.pedigree], HARTL_AGES)


class TestAGE9IndependenceFromTheBirthYearSentinel(unittest.TestCase):
    """
    AGE-9 -- the generic missing-result marker is not the birth-year sentinel.

    Deliberately a UNIT test on ``set_age``, driven by direct attribute
    assignment rather than through the loader. The parser path cannot express
    this comparison honestly: for a pedigree with no ``b``/``y`` column the port
    derives ``by`` from ``missing_bdate`` (``'01011800'`` -> ``int('1800')``),
    so a supplied ``missing_byear`` has no effect. That is a Finding-36
    property, and routing AGE-9 through it would test , not
    .
    """

    def test_the_restored_default_marker_is_what_lands_in_age(self):
        ped = lacy()
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], [ped.kw["missing_age"]] * 7)

    def test_the_two_options_are_distinct_keys_with_distinct_types(self):
        ped = lacy()
        self.assertIsInstance(ped.kw["missing_value"], float)
        self.assertIsInstance(ped.kw["missing_age"], int)
        self.assertNotIn("missing_byear", ped.kw)


class TestAGE9IndependenceAtTheSuppliedOptionLevel(unittest.TestCase):
    """
    AGE-9, the half that is already provable on the baseline, because
    supplying the option bypasses the missing default.

    Permanent: it states the Finding-33/Finding-36 boundary itself -- the
    generic missing-result marker follows ``missing_value`` and is unmoved by
    ``missing_byear``. That must remain true whatever  later decides
    about birth-year sentinels.
    """

    def test_the_marker_written_into_age_follows_missing_age(self):
        for marker in (-999, -7777):
            ped = lacy(missing_age=marker)
            self.assertTrue(pyp_utils.set_age(ped))
            self.assertEqual([a.age for a in ped.pedigree], [marker] * 7)

    def test_retired_missing_byear_does_not_move_the_marker(self):
        marker = -7777
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            ped = lacy(missing_age=marker, missing_byear=1234)
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual([a.age for a in ped.pedigree], [marker] * 7)
        self.assertNotIn("missing_byear", ped.kw)


if __name__ == "__main__":
    unittest.main()
