"""RC4 an earlier revision: current product-surface contract.

High-level truths after the RC4 feature phases: advertised capabilities
work; bounded option values use domain validation; examples stay out of
the wheel; the Introduction workflow does not need an installed examples
tree.
"""
from __future__ import annotations

import os
import shutil
import site
import subprocess
import sys
import tarfile
import textwrap
import venv
import zipfile
from pathlib import Path

import pytest

from PyPedal import pyp_metrics, pyp_network, pyp_snp, pyp_utils
from PyPedal.pyp_errors import PyPedalUsageError
from PyPedal.pyp_newclasses import load_pedigree

from _pedhelpers import REPO, chdir_tmp, load_corpus
from test_manual_pages import USER

MRODE_PEDIGREE = """\
# Pedigree from Mrode (2005) Table 2.1
1 0 0
2 0 0
3 1 2
4 1 0
5 4 3
6 5 2
"""


def _plain(name):
    return (USER / name).read_text(encoding="utf-8").replace("*", "").lower()


def test_chrometype_autosome_is_the_supported_domain():
    ped = load_corpus("new_lacy.ped")
    ng = pyp_metrics.effective_founder_genomes(
        ped, rounds=2, seed=31, chrometype="autosome", output=False, quiet=True
    )
    assert ng > 0
    with pytest.raises(PyPedalUsageError):
        pyp_metrics.effective_founder_genomes(
            ped, rounds=2, chrometype="sex", output=False, quiet=True
        )
    with pytest.raises(PyPedalUsageError):
        pyp_metrics.effective_founder_genomes(
            ped, rounds=2, chrometype="x", output=False, quiet=True
        )


def test_gen_coeff_true_is_outside_the_supported_domain(tmp_path):
    pedfile = tmp_path / "mrode.ped"
    pedfile.write_text(MRODE_PEDIGREE, encoding="utf-8")
    with chdir_tmp():
        with pytest.raises(PyPedalUsageError):
            load_pedigree(
                options={
                    "pedfile": str(pedfile),
                    "pedformat": "asd",
                    "messages": "quiet",
                    "pedigree_summary": 0,
                    "gen_coeff": True,
                }
            )
        ped = load_corpus("mrode.ped")
        ped.kw["gen_coeff"] = True
        with pytest.raises(PyPedalUsageError):
            pyp_utils.set_generation(ped)


def test_pedformat_p_stores_without_enabling_calculation(tmp_path):
    rows = "1 0 0 1.5\n2 1 0 2.0\n"
    pedfile = tmp_path / "gencoeff.ped"
    pedfile.write_text(rows, encoding="utf-8")
    with chdir_tmp():
        ped = load_pedigree(
            options={
                "pedfile": str(pedfile),
                "pedformat": "asdp",
                "messages": "quiet",
                "pedigree_summary": 0,
                "gen_coeff": False,
            }
        )
    assert ped.kw["gen_coeff"] is False
    values = [float(a.gencoeff) for a in ped.pedigree]
    assert values == pytest.approx([1.5, 2.0])


def test_dyad_census_is_not_a_current_api():
    with pytest.raises(PyPedalUsageError):
        pyp_network.dyad_census(None)
    assert not hasattr(pyp_snp, "read_agil_pedigree_file")
    header = Path(pyp_snp.__file__).read_text(encoding="utf-8").split("FUNCTIONS:")[1]
    header = header.split("## @package", 1)[0]
    assert "read_agil_pedigree_file" not in header


def test_restored_capabilities_are_not_called_unsupported_in_current_docs():
    current = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(USER.glob("*.md"))
    ).lower()
    readme = (Path(REPO) / "README.md").read_text(encoding="utf-8").lower()
    blob = current + "\n" + readme
    stale = (
        "pdf pedigree reports are **not**",
        "pdf pedigree reports are not part",
        "mating_coi` unsupported",
        "mating coefficient of inbreeding is unsupported",
        "missing_byear = 1800",
        "not generally transactional",
        "can partly mutate then return false",
        "pdf reports removed",
    )
    for needle in stale:
        assert needle not in blob, needle
    intro = _plain("index.md")
    assert "pdf_pedigree_metadata" in _plain("pdf-reports.md")
    assert "mating_coi" in intro or "test mating" in intro or "prospective" in intro
    assert "1 0 0" in (USER / "first-pedigree.md").read_text(encoding="utf-8")
    assert "Path(PyPedal.__file__)" not in (USER / "index.md").read_text(
        encoding="utf-8"
    )
    readme_text = (Path(REPO) / "README.md").read_text(encoding="utf-8")
    assert "4.0.0-rc8" not in readme_text
    assert readme_text.splitlines()[0] == "# PyPedal 4.1.0"
    assert "not been published" in readme_text.lower() or "not published" in readme_text.lower()
    prod = Path(REPO) / "PyPedal"
    hits = []
    for path in prod.glob("*.py"):
        text = path.read_text(encoding="utf-8")
        if "raise PyPedalNotImplementedError" in text or "raise NotImplementedError" in text:
            hits.append(path.name)
    assert hits == [], hits


def test_introduction_inline_mrode_workflow(tmp_path):
    pedfile = tmp_path / "mrode.ped"
    pedfile.write_text(MRODE_PEDIGREE, encoding="utf-8")
    with chdir_tmp():
        from PyPedal import pyp_nrm

        ped = load_pedigree(
            options={
                "pedfile": str(pedfile),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
            }
        )
        result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
    assert len(ped.pedigree) == 6
    assert result["fx"][5] == 0.125


def _build_setuptools_artifact(directory, kind):
    """Build a wheel or sdist via the setuptools PEP 517 backend.

    A stale gitignored ``build/`` tree previously packed historical modules and
    ``PyPedal/examples`` into the wheel even when pyproject excluded them.
    Wipe that cache before building so the artifact matches current config.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(Path(REPO) / "build", ignore_errors=True)
    cwd = os.getcwd()
    os.chdir(REPO)
    try:
        from setuptools.build_meta import build_sdist, build_wheel

        name = build_wheel(str(directory)) if kind == "wheel" else build_sdist(
            str(directory)
        )
    finally:
        os.chdir(cwd)
    return directory / name


def test_wheel_contains_no_examples(tmp_path):
    wheel = _build_setuptools_artifact(tmp_path / "dist", "wheel")
    names = zipfile.ZipFile(wheel).namelist()
    example_members = [n for n in names if "/examples/" in n or n.endswith("/examples")]
    assert example_members == [], example_members
    assert not any("griffon" in n.lower() for n in names)
    assert any(
        n.endswith("PyPedal/pyp_newclasses.py") or n.endswith("pyp_newclasses.py")
        for n in names
    )
    assert any(
        n.endswith("PyPedal/application/__init__.py")
        or n.endswith("application/__init__.py")
        for n in names
    )
    assert any(
        n.endswith("PyPedal/desktop/__init__.py")
        or n.endswith("desktop/__init__.py")
        for n in names
    )
    assert any("LICENSE" in n for n in names)
    assert wheel.stat().st_size > 1024
    assert any(
        n.endswith("PyPedal.ini") or n.endswith("PEDIGREE_FORMAT_CODES.txt")
        for n in names
    )


def test_sdist_may_include_examples(tmp_path):
    manifest = (Path(REPO) / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include PyPedal/examples" in manifest
    sdist = _build_setuptools_artifact(tmp_path / "dist", "sdist")
    with tarfile.open(sdist) as archive:
        names = archive.getnames()
    example_members = [n for n in names if "/examples/" in n]
    assert example_members, "sdist should retain repository examples"
    griffon_peds = [
        n for n in names if "griffon" in n.lower() and n.endswith(".ped")
    ]
    assert len(griffon_peds) == 1, griffon_peds
    assert griffon_peds[0].endswith("griffonbruxellois_2026_pyp.ped")


def test_wheel_install_runs_inline_mrode(tmp_path):
    wheel = _build_setuptools_artifact(tmp_path / "dist", "wheel")

    venv_dir = tmp_path / "venv"
    venv.create(venv_dir, with_pip=True, system_site_packages=False)
    py = venv_dir / "bin" / "python"
    if not py.exists():
        py = venv_dir / "Scripts" / "python.exe"
    install = subprocess.run(
        [str(py), "-m", "pip", "install", "--no-deps", str(wheel)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert install.returncode == 0, install.stdout + install.stderr

    venv_sites = []
    probe = subprocess.run(
        [
            str(py),
            "-c",
            "import site; print('\\n'.join(site.getsitepackages()))",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    for line in probe.stdout.splitlines():
        line = line.strip()
        if line:
            venv_sites.append(line)
    host_sites = site.getsitepackages()
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(venv_sites + host_sites)
    env["PYTHONNOUSERSITE"] = "1"

    work = tmp_path / "work"
    work.mkdir()
    script = work / "intro.py"
    script.write_text(
        textwrap.dedent(
            f"""
            from pathlib import Path
            from PyPedal.pyp_newclasses import load_pedigree
            from PyPedal import pyp_nrm
            import PyPedal
            pkg = Path(PyPedal.__file__).resolve().parent
            assert not (pkg / "examples").exists(), pkg
            pedfile = Path({str(work / "mrode.ped")!r})
            pedfile.write_text({MRODE_PEDIGREE!r})
            ped = load_pedigree(options={{
                "pedfile": str(pedfile),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
            }})
            result = pyp_nrm.inbreeding(ped, method="tabular", output=False)
            print(len(ped.pedigree), result["fx"][5])
            """
        ),
        encoding="utf-8",
    )
    run = subprocess.run(
        [str(py), str(script)],
        cwd=str(work),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert "6 0.125" in run.stdout
