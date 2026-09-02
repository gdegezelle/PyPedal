#!/usr/bin/env python3
"""
Quick verification script to test PyPedal setup.
Run this to verify your installation is working correctly.

Usage:
    python verify_setup.py
"""

import sys
import traceback

from pathlib import Path

REPO = Path(__file__).resolve().parent

REQUIRED_PUBLIC_FILES = (
    "README.md",
    "CHANGELOG.md",
    "MIGRATION.md",
    "THIRD-PARTY.md",
    "CONTRIBUTING.md",
    "LICENSE",
    "pyproject.toml",
    "docs/manual/index.md",
)

FORBIDDEN_PUBLIC_PATHS = (
    "docs/user",
    "docs/archive",
    "docs/manual-source",
    "docs/algorithms",
    "CLAUDE.md",
    "tools/difftest",
    "ADODB-LICENSE.txt",
)


def check_public_tree():
    """Confirm the publication tree has current docs and no archaeology."""
    print("=" * 60)
    print("Testing public tree...")
    print("=" * 60)
    ok = True
    for rel in REQUIRED_PUBLIC_FILES:
        path = REPO / rel
        if path.is_file() and path.stat().st_size > 0:
            print(f"    [OK] {rel}")
        else:
            print(f"    [FAIL] missing {rel}")
            ok = False
    manual = REPO / "docs" / "manual"
    if manual.is_dir() and any(manual.glob("*.md")):
        print("    [OK] docs/manual/")
    else:
        print("    [FAIL] docs/manual/ is missing")
        ok = False
    for rel in FORBIDDEN_PUBLIC_PATHS:
        path = REPO / rel
        if path.exists():
            print(f"    [FAIL] leftover {rel}")
            ok = False
        else:
            print(f"    [OK] no {rel}")
    if ok:
        print("\n[OK] Public tree looks current.\n")
    else:
        print("\n[FAIL] Public tree check failed.\n")
    return ok


def test_imports():
    """Test that all core modules can be imported."""
    print("=" * 60)
    print("Testing imports...")
    print("=" * 60)
    
    try:
        print("  - Importing PyPedal core modules...")
        from PyPedal import pyp_newclasses, pyp_nrm, pyp_metrics, pyp_utils
        print("    [OK] Core modules imported")

        print("  - Importing Qt-free application layer...")
        from PyPedal import application
        assert application.PedigreeSession is not None
        print("    [OK] PyPedal.application imported")
        
        print("  - Testing NewPedigree class...")
        from PyPedal.pyp_newclasses import NewPedigree
        print("    [OK] NewPedigree class available")
        
        print("  - Testing NRM functions...")
        from PyPedal.pyp_nrm import fast_a_matrix, inbreeding
        print("    [OK] NRM functions available")
        
        print("\n[OK] All imports successful!\n")
        return True
    except ImportError as e:
        print(f"\n[FAIL] Import error: {e}\n")
        traceback.print_exc()
        return False

def test_basic_functionality():
    """Test basic functionality with a minimal example."""
    print("=" * 60)
    print("Testing basic functionality...")
    print("=" * 60)
    
    try:
        from PyPedal.pyp_newclasses import NewPedigree
        from PyPedal import pyp_utils
        import shutil
        import tempfile

        print("  - Creating a pedigree object...")
        src = REPO / "PyPedal" / "examples" / "boichard2.ped"
        work = Path(tempfile.mkdtemp())
        pedfile = work / "boichard2.ped"
        shutil.copy(src, pedfile)
        options = {
            'pedfile': str(pedfile),
            'pedformat': 'asdx',
            'messages': 'quiet'
        }
        
        print(f"    Attempting to load: {options['pedfile']}")
        ped = NewPedigree(options)
        ped.load()
        
        print(f"    [OK] Loaded {len(ped.pedigree)} animals")
        
        # Test that we can access basic attributes
        if len(ped.pedigree) > 0:
            first_animal = ped.pedigree[0]
            print(f"    [OK] First animal ID: {first_animal.animalID}")
            print(f"    [OK] First animal name: {first_animal.name}")
        
        print("\n[OK] Basic functionality works!\n")
        return True
        
    except FileNotFoundError as e:
        print(f"\n[WARN] File not found: {e}")
        print("   This is okay if test data files are missing.")
        print("   Try running examples from PyPedal/examples/ instead.\n")
        return False
    except Exception as e:
        print(f"\n[FAIL] Error: {e}\n")
        traceback.print_exc()
        return False

# The core runtime dependencies -- the `dependencies` list in pyproject.toml.
# PyPedal cannot load or analyse a pedigree without these, so a missing one is
# a failed installation.
CORE_DEPENDENCIES = {
    'numpy': 'NumPy',
    'pandas': 'Pandas',
    'scipy': 'SciPy',
    'networkx': 'NetworkX',
}

# Optional feature dependencies, grouped by the extra that provides them.
# A missing one disables a FEATURE; it does not break the installation, and it
# must not fail setup verification.
#
# matplotlib belongs here, not above. It moved to the `graphics` extra because
# pyp_graphics imports it lazily inside each plotting function, so the library
# is fully usable without it -- this script previously listed it as required and
# so declared a perfectly good `pip install -e ".[dev]"` broken.
OPTIONAL_DEPENDENCIES = {
    'graphics': {
        'matplotlib': 'Matplotlib',
        'pydot': 'PyDot',
        'graphviz': 'Graphviz',
        'PIL': 'Pillow',
    },
    'reports': {'reportlab': 'ReportLab'},
    'gui': {'customtkinter': 'CustomTkinter'},
    'desktop': {'PySide6': 'PySide6'},
    'graphviz-extra': {'pygraphviz': 'PyGraphviz'},
}


def _importable(module):
    try:
        __import__(module)
        return True
    except ImportError:
        return False


def check_dependencies():
    """
    Verify the CORE dependencies and report on the optional ones.

    Returns True when the core install is sound, regardless of which optional
    features happen to be available.
    """
    print("=" * 60)
    print("Testing dependencies...")
    print("=" * 60)

    all_ok = True
    print("\n  Core (required):")
    for module, name in CORE_DEPENDENCIES.items():
        if _importable(module):
            print(f"    [OK] {name} installed")
        else:
            print(f"    [FAIL] {name} NOT installed")
            all_ok = False

    for extra, modules in OPTIONAL_DEPENDENCIES.items():
        available = [n for m, n in modules.items() if _importable(m)]
        missing = [n for m, n in modules.items() if not _importable(m)]
        if not missing:
            status = "available"
        elif available:
            status = "partly available"
        else:
            status = "not installed"
        print(f"\n  Optional feature '{extra}' -- {status}:")
        for name in available:
            print(f"    [OK] {name}")
        for name in missing:
            print(f"    [WARN] {name} not installed"
                  f"  (install with: pip install 'PyPedal[{extra}]')")

    if all_ok:
        print("\n[OK] Core dependencies installed. PyPedal is usable.")
        print("  Any feature listed above as not installed is optional; it")
        print("  disables that feature only.\n")
    else:
        print("\n[FAIL] A CORE dependency is missing -- PyPedal will not work.")
        print("   Run: pip install -e .\n")

    return all_ok


# Kept so that anything calling the old name still works.
test_dependencies = check_dependencies

def main():
    """Run all verification tests."""
    print("\n" + "=" * 60)
    print("PyPedal Setup Verification")
    print("=" * 60)
    print(f"Python version: {sys.version}\n")
    
    results = {
        'dependencies': check_dependencies(),
        'imports': test_imports(),
        'functionality': test_basic_functionality(),
        'public_tree': check_public_tree(),
    }
    
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    
    for test, passed in results.items():
        status = "[OK] PASS" if passed else "[FAIL] FAIL"
        print(f"  {status}: {test}")
    
    all_passed = all(results.values())
    
    if all_passed:
        print("\n[OK] All tests passed! PyPedal is ready to use.")
        print("\nNext steps:")
        print("  1. Try running: python -m pytest tests/ -m \"not integration\" -q")
        print("  2. Try examples: cd PyPedal/examples && python new_methods.py")
        print("  3. Read: docs/manual/index.md")
    else:
        print("\n[WARN] Some tests failed. See errors above.")
        print("\nTroubleshooting:")
        print("  1. Make sure conda environment is activated: conda activate pypedal3_env")
        print("  2. Reinstall: pip install -e .")
        print("  3. See CONTRIBUTING.md")
    
    print("=" * 60 + "\n")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())

