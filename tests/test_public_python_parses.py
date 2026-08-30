"""Public Python sources in the sdist must parse on Python 3.12+."""
import py_compile
from pathlib import Path

from _pedhelpers import REPO

ROOT = Path(REPO)
TREES = ("PyPedal", "tests", "tools")


def _public_python_files():
    files = []
    for tree in TREES:
        base = ROOT / tree
        if not base.is_dir():
            continue
        files.extend(sorted(base.rglob("*.py")))
    return files


def test_public_python_sources_parse():
    files = _public_python_files()
    assert files, "no public Python files found"
    errors = []
    for path in files:
        try:
            py_compile.compile(str(path), doraise=True)
        except (SyntaxError, TabError, IndentationError, py_compile.PyCompileError) as exc:
            errors.append("%s: %s" % (path.relative_to(ROOT), exc))
    assert errors == [], "\n".join(errors)


def test_newfoundland_pedigree_is_data_not_python():
    ped = ROOT / "PyPedal" / "examples" / "newfoundland.ped"
    assert ped.is_file()
    assert not (ROOT / "PyPedal" / "examples" / "newfoundland.ped.py").exists()
    text = ped.read_text(encoding="utf-8")
    assert text.lstrip().startswith("# Pedigree")
    assert "newfoundlanddog-database.net" in text
