"""RC4 an earlier revision: true unknown birth chronology."""
import datetime
import os
import pickle
import tempfile
import unittest
import warnings

from PyPedal import pyp_chronology, pyp_db, pyp_demog, pyp_io, pyp_metrics, pyp_nrm, pyp_utils
from PyPedal.pyp_chronology import VitalRateProfile
from PyPedal.pyp_errors import PyPedalUsageError, PyPedalValidationError
from PyPedal.pyp_newclasses import NewAnimal

from _pedhelpers import (
    chdir_tmp,
    load_canonical_griffon,
    load_corpus,
    load_corpus_from_path,
    load_example,
    load_griffon_1871_1890,
    load_griffon_test_small,
)


def rows_to_ped(rows, pedformat="asdy", **overrides):
    tmp = tempfile.mkdtemp(prefix="chrono_")
    path = os.path.join(tmp, "chrono.ped")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    options = {"messages": "quiet", "pedigree_summary": 0}
    options.update(overrides)
    return load_corpus_from_path(path, pedformat, **options)


class TestParseMissingAndMalformed(unittest.TestCase):
    def test_blank_and_dot_are_none(self):
        for token in ("", "   ", ".", None):
            self.assertIsNone(pyp_chronology.parse_recorded_year(token))
            bd, by = pyp_chronology.parse_recorded_date(token)
            self.assertIsNone(bd)
            self.assertIsNone(by)

    def test_year_zero_is_missing_on_input(self):
        self.assertIsNone(pyp_chronology.parse_recorded_year(0))
        self.assertIsNone(pyp_chronology.parse_recorded_year("0"))
        bd, by = pyp_chronology.parse_recorded_date("0")
        self.assertIsNone(bd)
        self.assertIsNone(by)

    def test_malformed_year_raises(self):
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.parse_recorded_year("20XX")
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.parse_recorded_year("foobar")

    def test_malformed_and_impossible_dates_raise(self):
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.parse_recorded_date("foobar")
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.parse_recorded_date("2026-99-88")
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.parse_recorded_date("02312026")

    def test_year_only_in_b_does_not_invent_a_day(self):
        bd, by = pyp_chronology.parse_recorded_date("1991")
        self.assertIsNone(bd)
        self.assertEqual(by, 1991)

    def test_mmddyyyy_becomes_datetime_date(self):
        bd, by = pyp_chronology.parse_recorded_date("04251991")
        self.assertEqual(bd, datetime.date(1991, 4, 25))
        self.assertEqual(by, 1991)

    def test_seven_digit_mmdyyyy_is_unpadded_day(self):
        bd, by = pyp_chronology.parse_recorded_date("0111957")
        self.assertEqual(bd, datetime.date(1957, 1, 1))
        self.assertEqual(by, 1957)


class TestRealYearsSurvive(unittest.TestCase):
    def test_1800_and_1900_are_real_years_by_default(self):
        ped = rows_to_ped(["1 0 0 1800", "2 0 0 1900"], "asdy")
        years = {animal.originalID: animal.by for animal in ped.pedigree}
        self.assertEqual(years[1], 1800)
        self.assertEqual(years[2], 1900)
        self.assertTrue(all(animal.bd is None for animal in ped.pedigree))

    def test_legacy_token_1800_does_not_consume_1900(self):
        ped = rows_to_ped(
            ["1 0 0 1800", "2 0 0 1900"],
            "asdy",
            legacy_missing_byear_token=1800,
        )
        by_oid = {animal.originalID: animal.by for animal in ped.pedigree}
        self.assertIsNone(by_oid[1])
        self.assertEqual(by_oid[2], 1900)

    def test_legacy_token_1900_does_not_consume_1800(self):
        ped = rows_to_ped(
            ["1 0 0 1800", "2 0 0 1900"],
            "asdy",
            legacy_missing_byear_token=1900,
        )
        by_oid = {animal.originalID: animal.by for animal in ped.pedigree}
        self.assertEqual(by_oid[1], 1800)
        self.assertIsNone(by_oid[2])

    def test_legacy_bdate_token_maps_only_that_string(self):
        ped = rows_to_ped(
            ["1 0 0 01011800", "2 0 0 01011900"],
            "asdb",
            legacy_missing_bdate_token="01011800",
        )
        by_oid = {a.originalID: (a.bd, a.by) for a in ped.pedigree}
        self.assertEqual(by_oid[1], (None, None))
        self.assertEqual(by_oid[2][0], datetime.date(1900, 1, 1))
        self.assertEqual(by_oid[2][1], 1900)


class TestMrodeAndScience(unittest.TestCase):
    def test_mrode_asd_loads_with_unknown_chronology(self):
        ped = load_corpus("mrode.ped")
        self.assertEqual(len(ped.pedigree), 6)
        self.assertTrue(all(animal.by is None for animal in ped.pedigree))
        self.assertTrue(all(animal.bd is None for animal in ped.pedigree))
        self.assertEqual(ped.metadata.num_unique_years, 0)
        self.assertEqual(ped.metadata.num_unknown_birth_years, 6)

    def test_mrode_inbreeding_unchanged(self):
        ped = load_corpus("mrode.ped")
        result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
        self.assertAlmostEqual(float(result["fx"][5]), 0.125, places=6)

    def test_mating_coi_unchanged(self):
        ped = load_corpus("mrode.ped")
        from PyPedal.pyp_metrics import mating_coi
        self.assertAlmostEqual(mating_coi(1, 2, ped), 0.0, places=6)

    def test_lacy_control_unchanged(self):
        ped = load_corpus("new_lacy.ped")
        out = pyp_metrics.a_effective_founders_lacy(ped)
        self.assertAlmostEqual(float(out["fa_effective_founders"]), 2.909090909, places=6)

    def test_retired_missing_byear_is_ignored(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ped = load_corpus("mrode.ped", missing_byear=1800)
        self.assertNotIn("missing_byear", ped.kw)
        self.assertTrue(all(animal.by is None for animal in ped.pedigree))
        self.assertTrue(any(item.category is DeprecationWarning for item in caught))


class TestParentOrdering(unittest.TestCase):
    def test_exact_dates_parent_earlier_is_valid(self):
        rows_to_ped(
            ["1 0 0 01011990", "2 1 0 02011991"],
            "asdb",
        )

    def test_exact_same_day_is_invalid(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(
                ["1 0 0 01011990", "2 1 0 01011990"],
                "asdb",
            )

    def test_exact_parent_later_is_invalid(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(
                ["1 0 0 01011992", "2 1 0 01011991"],
                "asdb",
            )

    def test_year_only_same_year_is_valid(self):
        rows_to_ped(["1 0 0 1990", "2 1 0 1990"], "asdy")

    def test_year_only_parent_later_is_invalid(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(["1 0 0 1992", "2 1 0 1991"], "asdy")

    def test_mixed_resolution_same_year_is_valid(self):
        rows_to_ped(
            ["1 0 0 1990 .", "2 1 0 1990 06151990"],
            "asdyb",
        )

    def test_unknown_parent_chronology_is_valid(self):
        rows_to_ped(["1 0 0 .", "2 1 0 1991"], "asdy")

    def test_year_and_date_conflict_raises(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(["1 0 0 1990 01011991"], "asdyb")

    def test_year_only_b_conflicts_with_y(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(["1 0 0 1990 1991"], "asdyb")

    def test_mixed_parent_later_year_is_invalid(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(
                ["1 0 0 1992 .", "2 1 0 1991 06151991"],
                "asdyb",
            )

    def test_malformed_year_on_load_raises(self):
        with self.assertRaises(PyPedalValidationError):
            rows_to_ped(["1 0 0 20XX"], "asdy")

    def test_blank_year_on_load_is_none(self):
        ped = rows_to_ped(["1 0 0 ."], "asdy")
        self.assertIsNone(ped.pedigree[0].by)
        self.assertIsNone(ped.pedigree[0].bd)


class TestSetAgeAndDemography(unittest.TestCase):
    def test_set_age_uses_year_offset_not_igen(self):
        ped = load_corpus("new_lacy.ped")
        pyp_utils.set_generation(ped)
        self.assertTrue(pyp_utils.set_age(ped))
        missing = ped.kw["missing_age"]
        for animal in ped.pedigree:
            self.assertEqual(animal.age, missing)
            self.assertNotEqual(animal.age, animal.igen)

    def test_known_year_offset(self):
        ped = rows_to_ped(["1 0 0 1900"], "asdy")
        self.assertTrue(pyp_utils.set_age(ped))
        self.assertEqual(ped.pedigree[0].age, 100)

    def test_generation_intervals_unknown_are_none(self):
        ped = load_corpus("mrode.ped")
        out = pyp_metrics.generation_intervals(ped)
        self.assertEqual(set(out), {"ss", "sd", "ds", "dd", "mean"})
        self.assertTrue(all(value is None for value in out.values()))

    def test_generation_intervals_known_years_are_numeric(self):
        ped = load_example(
            "generations.ped",
            {
                "pedformat": "asdbx",
                "messages": "quiet",
                "pedigree_summary": 0,
            },
        )
        out = pyp_metrics.generation_intervals(ped)
        self.assertTrue(all(value is not None for value in out.values()))
        self.assertGreater(out["mean"], 0)

    def test_founders_by_year_omits_unknown(self):
        ped = load_corpus("mrode.ped")
        self.assertEqual(pyp_demog.founders_by_year(ped), {})

    def test_founders_by_year_does_not_bucket_at_1800(self):
        ped = rows_to_ped(["1 0 0 1990", "2 0 0 .", "3 1 2 1991"], "asdy")
        by_year = pyp_demog.founders_by_year(ped)
        self.assertNotIn(1800, by_year)
        self.assertNotIn(None, by_year)
        self.assertEqual(by_year[1990], 1)


class TestPadId(unittest.TestCase):
    def test_unknown_year_does_not_embed_none_or_1800(self):
        ped = load_corpus("mrode.ped")
        for animal in ped.pedigree:
            self.assertNotIn("None", animal.paddedID)
            self.assertFalse(animal.paddedID.startswith("1800"))
            self.assertFalse(animal.paddedID.startswith("1900"))

    def test_known_year_keeps_year_prefix(self):
        locations = {"animal": 0, "sire": 1, "dam": 2, "birthyear": 3}
        animal = NewAnimal(
            locations,
            [123, 0, 0, 2005],
            {"missing_parent": 0, "missing_name": "unknown", "pedformat": "asdy"},
        )
        self.assertEqual(animal.pad_id(), "200500000000123")


class TestRoundTrip(unittest.TestCase):
    def test_text_save_load_preserves_unknown_and_real_years(self):
        ped = rows_to_ped(["1 0 0 1800", "2 0 0 .", "3 0 0 1900"], "asdy")
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "rt.ped")
            self.assertTrue(ped.save(filename=out, pedformat="asdy"))
            with open(out, encoding="utf-8") as handle:
                body = [
                    line.strip()
                    for line in handle
                    if line.strip() and not line.startswith("#")
                ]
            self.assertTrue(any(line.endswith(" .") or line.split()[-1] == "." for line in body))
            self.assertFalse(any("None" in line for line in body))
            reloaded = load_corpus_from_path(out, "asdy", messages="quiet")
        by_oid = {a.originalID: a.by for a in reloaded.pedigree}
        self.assertEqual(by_oid[1], 1800)
        self.assertIsNone(by_oid[2])
        self.assertEqual(by_oid[3], 1900)

    def test_pickle_preserves_none_and_date(self):
        ped = rows_to_ped(["1 0 0 04251991", "2 0 0 ."], "asdb")
        payload = pickle.dumps(ped)
        restored = pickle.loads(payload)
        by_oid = {a.originalID: a for a in restored.pedigree}
        self.assertEqual(by_oid[1].bd, datetime.date(1991, 4, 25))
        self.assertEqual(by_oid[1].by, 1991)
        self.assertIsNone(by_oid[2].bd)
        self.assertIsNone(by_oid[2].by)

    def test_sqlite_stores_null_for_unknown_year(self):
        ped = rows_to_ped(["1 0 0 1800", "2 0 0 ."], "asdy")
        with chdir_tmp() as tmp:
            dbfile = os.path.join(tmp, "chrono.db")
            ped.kw["database_name"] = dbfile
            ped.kw["database_file"] = dbfile
            ped.kw["database_table"] = "chrono"
            conn = pyp_db.connect_to_database(ped)
            self.assertIsNotNone(conn)
            try:
                self.assertTrue(pyp_db.create_pedigree_table(ped, conn=conn, drop=True))
                self.assertTrue(pyp_db.populate_pedigree_table(ped, conn=conn))
                rows = conn.execute(
                    "SELECT originalID, birthyear FROM chrono ORDER BY originalID"
                ).fetchall()
            finally:
                conn.close()
        values = {str(oid): year for oid, year in rows}
        self.assertEqual(values["1"], 1800)
        self.assertIsNone(values["2"])

    def test_gedcom_omits_missing_date(self):
        ped = rows_to_ped(["1 0 0 1991", "2 0 0 ."], "asdy")
        with chdir_tmp() as tmp:
            ged = os.path.join(tmp, "chrono.ged")
            self.assertTrue(pyp_io.save_to_gedcom(ped, ged))
            text = open(ged, encoding="utf-8").read()
        self.assertIn("DATE 1991", text)
        self.assertNotIn("0001", text)
        self.assertNotIn("DATE 1800", text)

    def test_gedcom_tmp_keeps_eight_digit_dates(self):
        assembled = {
            "A": {
                "indi": "A",
                "sire": 0,
                "dam": 0,
                "sex": "M",
                "birth": "06151991",
                "name": "A",
            },
            "B": {
                "indi": "B",
                "sire": 0,
                "dam": 0,
                "sex": "F",
                "birth": None,
                "name": "B",
            },
        }
        with chdir_tmp() as tmp:
            out = os.path.join(tmp, "fromged.tmp")
            self.assertEqual(pyp_io.save_from_gedcom(out, assembled), "ASDxbu")
            body = open(out, encoding="utf-8").read()
        self.assertIn("06151991", body)
        self.assertIn(".", body)
        self.assertNotIn("None", body)


class TestEstimation(unittest.TestCase):
    def _profile(self):
        return VitalRateProfile(
            name="unit",
            gestation_days=60,
            sire_min_age_at_conception_days=365,
            sire_max_age_at_conception_days=3650,
            dam_min_age_at_conception_days=365,
            dam_max_age_at_conception_days=3650,
            founder_typical_age_at_progeny_days=730,
        )

    def test_default_load_does_not_estimate(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 06152000"],
            "asdb",
        )
        self.assertTrue(ped.pedigree[0].birth_date_estimate.is_empty())

    def test_profile_without_call_does_not_estimate(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 06152000"],
            "asdb",
            vital_rate_profile=self._profile(),
        )
        self.assertTrue(ped.pedigree[0].birth_date_estimate.is_empty())

    def test_estimate_on_load_requires_both_flags(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 06152000"],
            "asdb",
            estimate_birth_dates=True,
            vital_rate_profile=self._profile(),
        )
        parent = next(a for a in ped.pedigree if a.originalID == 1)
        self.assertIsNone(parent.bd)
        self.assertIsNone(parent.by)
        self.assertIsNotNone(parent.birth_date_estimate.earliest)

    def test_explicit_estimate_fills_range_not_facts(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 06152000"],
            "asdb",
        )
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        parent = next(a for a in ped.pedigree if a.originalID == 1)
        self.assertIsNone(parent.by)
        self.assertIsNone(parent.bd)
        self.assertIsNotNone(parent.birth_date_estimate.earliest)
        self.assertIsNotNone(parent.birth_date_estimate.latest)
        self.assertEqual(parent.birth_date_estimate.source, "vital-rate-profile")
        self.assertEqual(parent.birth_date_estimate.profile, "unit")

    def test_inconsistent_bounds_raise(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 01011990", "3 1 0 01012020"],
            "asdb",
        )
        tight = VitalRateProfile(
            name="tight",
            gestation_days=1,
            sire_min_age_at_conception_days=0,
            sire_max_age_at_conception_days=10,
        )
        with self.assertRaises(PyPedalValidationError):
            pyp_chronology.estimate_birth_date_ranges(ped, tight)

    def test_founder_typical_is_separate(self):
        ped = rows_to_ped(
            ["1 0 0 .", "2 1 0 06152000"],
            "asdb",
        )
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        parent = next(a for a in ped.pedigree if a.originalID == 1)
        self.assertIsNotNone(parent.birth_date_estimate.typical)
        self.assertIsNone(parent.bd)
        self.assertIsNone(parent.by)

    def test_half_founder_has_no_founder_typical(self):
        ped = rows_to_ped(
            ["1 0 0 01011980", "2 1 0 .", "3 2 0 06152000"],
            "asdb",
        )
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        half = next(a for a in ped.pedigree if a.originalID == 2)
        self.assertIsNone(half.birth_date_estimate.typical)

    def test_known_facts_are_not_overwritten(self):
        ped = rows_to_ped(
            ["1 0 0 01011980", "2 1 0 06152000"],
            "asdb",
        )
        before = [(a.bd, a.by) for a in ped.pedigree]
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        self.assertEqual([(a.bd, a.by) for a in ped.pedigree], before)

    def test_no_builtin_dog_profile(self):
        self.assertFalse(hasattr(pyp_chronology, "DOG_PROFILE"))
        self.assertFalse(hasattr(pyp_chronology, "dog_profile"))

    def test_year_only_offspring_does_not_invent_a_date(self):
        ped = rows_to_ped(["1 0 0 .", "2 1 0 2000"], "asdy")
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        parent = next(a for a in ped.pedigree if a.originalID == 1)
        self.assertTrue(parent.birth_date_estimate.is_empty())

    def test_estimate_does_not_change_inbreeding(self):
        ped = load_corpus("mrode.ped")
        before = pyp_nrm.inbreeding(ped, method="tabular", output=False)["fx"][5]
        pyp_chronology.estimate_birth_date_ranges(ped, self._profile())
        after = pyp_nrm.inbreeding(ped, method="tabular", output=False)["fx"][5]
        self.assertEqual(before, after)

    def test_invalid_profile_raises_usage_error(self):
        with self.assertRaises(PyPedalUsageError):
            VitalRateProfile(gestation_days=-1).validate()


class TestGeneDropLabels(unittest.TestCase):
    def test_founder_genomes_still_numeric(self):
        ped = load_corpus("new_lacy.ped")
        value = pyp_metrics.effective_founder_genomes(ped, rounds=2, quiet=True)
        self.assertGreater(float(value), 0)


class TestGriffonUnknownNotSentinel(unittest.TestCase):
    def test_implicit_parent_is_none_not_1800(self):
        ped = load_griffon_test_small()
        years = {animal.by for animal in ped.pedigree}
        self.assertIn(None, years)
        known = {year for year in years if year is not None}
        self.assertTrue(any(year >= 1900 for year in known) or min(known) >= 1870)
        self.assertNotEqual(min(known), 1800)
        by_year = pyp_demog.founders_by_year(ped)
        self.assertNotIn(1800, by_year)
        self.assertNotIn(None, by_year)


class TestHistoricalImpossibleChronologyIsRefused(unittest.TestCase):
    """Dated historical files keep their chronology columns.

    RC4 refuses files whose recorded parent dates are on or after offspring
    (``new_amatrix.ped``, same-day registry rows). The corrected canonical
    Griffon sample satisfies parent/child order and loads as ``asdxb``.
    """

    def test_new_amatrix_asdgb_refuses_dam_later_than_offspring(self):
        with self.assertRaises(PyPedalValidationError) as raised:
            load_example(
                "new_amatrix.ped",
                {
                    "pedformat": "asdgb",
                    "sepchar": " ",
                    "messages": "quiet",
                    "pedigree_summary": 0,
                },
            )
        message = str(raised.exception)
        self.assertIn("2047", message)
        self.assertIn("1997", message)
        self.assertIn("2048", message)
        self.assertIn("1992", message)

    def test_griffon_characterisation_cohort_asdxb_loads(self):
        ped = load_griffon_1871_1890()
        self.assertEqual(len(ped.pedigree), 167)
        self.assertTrue(any(isinstance(animal.bd, datetime.date) for animal in ped.pedigree))

    def test_griffonbruxellois_2026_pyp_asdxb_loads(self):
        ped = load_canonical_griffon()
        self.assertEqual(len(ped.pedigree), 98001)
        self.assertEqual(ped.metadata.num_records, 98001)
        self.assertEqual(ped.metadata.num_implicit_parents, 0)
        original_ids = [str(animal.originalID) for animal in ped.pedigree]
        self.assertEqual(len(original_ids), len(set(original_ids)))
        self.assertNotIn("51627", set(original_ids))
        self.assertNotIn("45936", set(original_ids))
        self.assertIn("57922", set(original_ids))
        self.assertEqual(ped.metadata.num_unknown_birth_years, 3997)
        years = {animal.by for animal in ped.pedigree if animal.by is not None}
        self.assertNotIn(1800, years)
        self.assertEqual(sum(1 for animal in ped.pedigree if animal.by == 1900), 152)
        self.assertEqual(min(years), 1870)
        self.assertEqual(max(years), 2025)
        self.assertTrue(any(isinstance(animal.bd, datetime.date) for animal in ped.pedigree))
        by_oid = {str(a.originalID): a for a in ped.pedigree}
        animal = by_oid["51614"]
        sire = by_oid["97427"]
        dam = by_oid["24636"]
        self.assertEqual(animal.bd, datetime.date(1898, 1, 1))
        self.assertEqual(sire.bd, datetime.date(1896, 8, 25))
        self.assertEqual(dam.bd, datetime.date(1896, 1, 1))
        self.assertLess(sire.bd, animal.bd)
        self.assertLess(dam.bd, animal.bd)
        self.assertEqual(ped.metadata.num_unique_founders, 6689)
        missing = ped.kw["missing_parent"]
        half = sum(
            1
            for a in ped.pedigree
            if a.founder == "n"
            and (
                (a.sireID == missing and a.damID != missing)
                or (a.sireID != missing and a.damID == missing)
            )
        )
        self.assertEqual(half, 915)
        id_to_index = {a.animalID: idx for idx, a in enumerate(ped.pedigree)}
        parent_after_child = 0
        for a in ped.pedigree:
            idx = id_to_index[a.animalID]
            for parent in (a.sireID, a.damID):
                if parent != missing and parent in id_to_index and id_to_index[parent] > idx:
                    parent_after_child += 1
        self.assertEqual(parent_after_child, 0)

    def test_registry_same_day_seven_digit_date_is_refused(self):
        with self.assertRaises(PyPedalValidationError) as raised:
            rows_to_ped(
                ["88634 0 0 m 0111915", "57922 88634 0 f 0111915"],
                "asdxb",
            )
        message = str(raised.exception)
        self.assertIn("1915-01-01", message)
        self.assertIn("88634", message)
        self.assertIn("57922", message)
