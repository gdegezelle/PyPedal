"""
GENES import/export is gone from the supported product.
"""
import ast
import inspect
import os
import unittest

from PyPedal import pyp_app, pyp_errors, pyp_io, pyp_metrics, pyp_newclasses, pyp_nrm
from PyPedal.pyp_newclasses import NewPedigree

from _pedhelpers import EXAMPLES, REPO, chdir_tmp, load_corpus, load_example

GENES_IO_SYMBOLS = (
    "load_from_genes",
    "save_from_genes",
    "save_to_genes",
    "_GENES_FIELD_MAP",
    "_format_genes_field",
    "_encode_genes_record",
)


def _source(obj):
    return inspect.getsource(obj)


class TestGenesProductionSurfaceRemoved(unittest.TestCase):
    """GEN-RM-1, GEN-RM-3, GEN-RM-7: the executable GENES surface is gone."""

    def test_gen_rm_3_pyp_io_has_no_genes_symbols(self):
        for name in GENES_IO_SYMBOLS:
            self.assertFalse(
                hasattr(pyp_io, name),
                "removed GENES symbol still present: pyp_io.%s" % name)

    def test_gen_rm_3_newpedigree_has_no_savegenes(self):
        self.assertFalse(hasattr(NewPedigree, "savegenes"))

    def test_gen_rm_1_and_7_load_does_not_advertise_or_dispatch_genesfile(self):
        src = _source(NewPedigree.load)
        self.assertNotIn("genesfile", src)
        self.assertNotIn("load_from_genes", src)
        self.assertNotIn("GENES", src)

    def test_gen_rm_1_loadpedigree_does_not_document_genesfile(self):
        src = _source(pyp_newclasses.loadPedigree)
        self.assertNotIn("genesfile", src)

    def test_gen_rm_1_genesfile_is_an_invalid_pedsource(self):
        with chdir_tmp() as tmp:
            dummy = os.path.join(tmp, "dummy.ped")
            open(dummy, "w").write("1 0 0\n")
            ped = NewPedigree({
                "pedfile": dummy,
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
            })
            with self.assertRaises(pyp_errors.PyPedalUsageError) as caught:
                ped.load(pedsource="genesfile")
            self.assertNotIn("genesfile", str(caught.exception).split(
                "Valid sources are:", 1)[-1])

    def test_gen_rm_7_obsolete_genes_export_tests_are_gone(self):
        path = os.path.join(REPO, "tests", "test_genes_export.py")
        self.assertFalse(os.path.isfile(path))

    def test_gen_rm_7_pyp_io_no_longer_imports_struct_or_decimal(self):
        path = os.path.join(REPO, "PyPedal", "pyp_io.py")
        tree = ast.parse(open(path, encoding="utf-8").read())
        imported = set()
        for node in tree.body:
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imported.add(alias.name)
        self.assertNotIn("struct", imported)
        self.assertNotIn("decimal", imported)

    def test_gen_rm_7_no_genes_dispatch_in_production_python(self):
        root = os.path.join(REPO, "PyPedal")
        hits = []
        for dirpath, _dirs, files in os.walk(root):
            for name in files:
                if not name.endswith(".py"):
                    continue
                path = os.path.join(dirpath, name)
                text = open(path, encoding="utf-8", errors="replace").read()
                for token in ("genesfile", "savegenes", "save_to_genes",
                              "load_from_genes", "save_from_genes",
                              "_GENES_FIELD_MAP"):
                    if token in text:
                        hits.append("%s: %s" % (path, token))
        self.assertEqual([], hits)


class TestGenesRemovalAlreadyTrue(unittest.TestCase):
    """GEN-RM checks that already hold on rc1 and must remain true."""

    def test_gen_rm_2_gui_does_not_advertise_genes(self):
        src = _source(pyp_app)
        lowered = src.lower()
        for token in ("genesfile", "savegenes", "save_to_genes",
                      "load_from_genes", "genes 1.20"):
            self.assertNotIn(token, lowered)

    def test_gen_rm_2_no_config_or_package_metadata_selector(self):
        ini = os.path.join(REPO, "PyPedal", "PyPedal.ini")
        if os.path.isfile(ini):
            text = open(ini, encoding="utf-8").read().lower()
            self.assertNotIn("genes", text)
        toml = open(os.path.join(REPO, "pyproject.toml"), encoding="utf-8").read()
        self.assertNotIn("genesfile", toml.lower())
        self.assertNotIn("savegenes", toml.lower())

    def test_gen_rm_3_package_imports(self):
        import PyPedal
        self.assertTrue(PyPedal.__version__)
        self.assertEqual("4.0.0", PyPedal.__version__)

    def test_gen_rm_8_export_format_error_remains(self):
        self.assertTrue(issubclass(
            pyp_errors.PyPedalExportFormatError, pyp_errors.PyPedalError))
        exc = pyp_errors.PyPedalExportFormatError(
            "FORMAT", "FIELD", 3, -999, 4)
        self.assertEqual(73, pyp_app.exit_status_for(exc))

    def test_missing_age_is_not_genes_only(self):
        ped = load_corpus("mrode.ped", "asd")
        self.assertIn("missing_age", ped.kw)
        self.assertEqual(-999, ped.kw["missing_age"])

    def test_gedcom_symbols_are_independent_of_genes(self):
        for name in ("load_from_gedcom", "save_to_gedcom", "save_from_gedcom"):
            self.assertTrue(hasattr(pyp_io, name), name)
        self.assertTrue(hasattr(NewPedigree, "savegedcom"))

    def test_no_genes_working_example_or_dbf_fixture(self):
        for root, _dirs, files in os.walk(EXAMPLES):
            for name in files:
                lowered = name.lower()
                self.assertFalse(
                    lowered.endswith(".dbf"),
                    "GENES fixture packaged as an example: %s" %
                    os.path.join(root, name))
                self.assertNotIn("genes", lowered)


class TestSupportedIOStillWorks(unittest.TestCase):
    """GEN-RM-4 / GEN-RM-5: remaining I/O must survive GENES deletion."""

    def test_gen_rm_4_file_load_mrode_and_lacy(self):
        mrode = load_corpus("mrode.ped", "asd")
        self.assertEqual(6, len(mrode.pedigree))
        lacy = load_example("new_lacy.ped", {
            "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "renumber": True,
        })
        self.assertEqual(7, len(lacy.pedigree))

    def test_gen_rm_4_textstream_load(self):
        with chdir_tmp() as tmp:
            dummy = os.path.join(tmp, "stream.ped")
            open(dummy, "w").write("")
            ped = NewPedigree({
                "pedfile": dummy,
                "pedformat": "ASD",
                "messages": "quiet",
                "pedigree_summary": 0,
                "renumber": True,
            })
            ped.load(pedsource="textstream",
                     pedstream="a1,s1,d1\na2,s2,d2\na3,a1,a2\n")
            self.assertGreaterEqual(len(ped.pedigree), 3)

    def test_gen_rm_4_null_load(self):
        with chdir_tmp() as tmp:
            dummy = os.path.join(tmp, "null.ped")
            open(dummy, "w").write("")
            ped = NewPedigree({
                "pedfile": dummy,
                "messages": "quiet",
                "pedigree_summary": 0,
            })
            ped.load(pedsource="null")
            self.assertEqual(0, len(ped.pedigree))

    def test_gen_rm_5_save_and_oldsave(self):
        ped = load_corpus("mrode.ped", "asd")
        with chdir_tmp() as tmp:
            saved = os.path.join(tmp, "mrode_saved.ped")
            self.assertTrue(ped.save(filename=saved, pedformat="asd"))
            self.assertTrue(os.path.isfile(saved))
            old = os.path.join(tmp, "mrode_oldsave.ped")
            self.assertTrue(ped.oldsave(filename=old))
            self.assertTrue(os.path.isfile(old))

    def test_gen_rm_5_gedcom_export(self):
        # Write-only: the GEDCOM *reader* is a remaining supported path, but
        # round-tripping savegedcom output through load(pedsource='gedcomfile')
        # is not a GENES-removal obligation and the reader can hang on some
        # files. Export must still succeed after GENES is deleted.
        ped = load_corpus("mrode.ped", "asd")
        with chdir_tmp() as tmp:
            ged = os.path.join(tmp, "mrode.ged")
            self.assertTrue(ped.savegedcom(pedoutfile=ged))
            self.assertTrue(os.path.isfile(ged))
            body = open(ged, encoding="utf-8").read()
            self.assertIn("INDI", body)
            self.assertGreater(os.path.getsize(ged), 0)


class TestCurrentDocsDoNotClaimGenesSupport(unittest.TestCase):
    """Current-facing docs must not advertise GENES as supported."""

    def test_user_manual_does_not_advertise_genes(self):
        manual = os.path.join(REPO, "docs", "manual")
        for name in os.listdir(manual):
            if not name.endswith(".md"):
                continue
            body = open(os.path.join(manual, name), encoding="utf-8").read()
            if "GENES" not in body:
                continue
            lowered = body.replace("*", "").lower()
            self.assertTrue(
                "not supported" in lowered
                or "not a supported" in lowered
                or "removed" in lowered
                or "gone" in lowered,
                name,
            )

    def test_readme_does_not_advertise_genes(self):
        body = open(os.path.join(REPO, "README.md"), encoding="utf-8").read()
        if "GENES" in body:
            lowered = body.lower()
            self.assertTrue("removed" in lowered or "not" in lowered)


class TestScientificOraclesUnchanged(unittest.TestCase):
    """GEN-RM-9: published scientific pins must not move."""

    def test_mrode_animal_5_coi_is_one_eighth(self):
        ped = load_example("mrode.ped", {
            "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "renumber": True,
            "pedigree_is_renumbered": True,
        })
        result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
        fx = {int(k): float(v) for k, v in result["fx"].items()}
        self.assertAlmostEqual(fx[5], 0.125, places=3)

    def test_lacy_effective_founders_is_2_91(self):
        ped = load_example("new_lacy.ped", {
            "pedformat": "asd", "sepchar": " ",
            "messages": "quiet", "renumber": True,
            "pedigree_is_renumbered": True,
        })
        result = pyp_metrics.effective_founders_lacy(ped)
        fe = result.get("fa_effective_founders", -999.9)
        self.assertEqual(2.91, round(fe, 2))
