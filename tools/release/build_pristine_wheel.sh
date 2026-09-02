#!/usr/bin/env bash
# Build a local wheel from a PRISTINE tracked tree and inspect it.
#
# Final artifacts MUST be built this way, never from a dirty
# working copy. A working-tree `pip wheel .` will pack ignored generated
# docs that are not part of the release.
#
# Usage (from the repository root):
#
#   ./tools/release/build_pristine_wheel.sh
#   ./tools/release/build_pristine_wheel.sh /tmp/pypedal-pristine-wheel
#   EXPECTED_VERSION=<pep440> ./tools/release/build_pristine_wheel.sh
#
# The expected wheel version defaults to [project].version in the archived
# pyproject.toml (Python 3.12+ tomllib). Override with EXPECTED_VERSION.
# The built METADATA Version must match that value.
#
# Exit status is 0 only if:
#   - exactly one wheel is produced
#   - its version matches the derived (or supplied) project version
#   - LICENSE is present
#   - production modules are present
#   - examples, Griffon data, and GENES HTML are absent
#
# The wheel is written outside the repository. Nothing is published.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT="${1:-/tmp/pypedal-pristine-wheel}"
SRC="$(mktemp -d /tmp/pypedal-pristine-src.XXXXXX)"

mkdir -p "$OUT"
git -C "$REPO" archive HEAD | tar -x -C "$SRC"

if [ -z "${EXPECTED_VERSION:-}" ]; then
  EXPECTED_VERSION="$(python3 - "$SRC" <<'PY'
import sys
import tomllib
from pathlib import Path

with (Path(sys.argv[1]) / "pyproject.toml").open("rb") as fh:
    print(tomllib.load(fh)["project"]["version"])
PY
)"
fi

python3 -m pip wheel "$SRC" -w "$OUT" --no-deps --no-build-isolation

python3 - "$OUT" "$EXPECTED_VERSION" <<'PY'
import glob
import os
import sys
import zipfile

out, expected = sys.argv[1], sys.argv[2]
wheels = sorted(glob.glob(os.path.join(out, "*.whl")))
if len(wheels) != 1:
    sys.stderr.write("expected exactly one wheel in %s, got %s\n" % (out, wheels))
    sys.exit(1)
path = wheels[0]
print("wheel:", path)
print("size:", os.path.getsize(path))
z = zipfile.ZipFile(path)
names = z.namelist()
print("members:", len(names))
meta = next(n for n in names if n.endswith("METADATA"))
version = None
for line in z.read(meta).decode().splitlines():
    if line.startswith("Version:"):
        version = line.split(":", 1)[1].strip()
        break
print("version:", version)
print("expected:", expected)
if version != expected:
    sys.stderr.write("unexpected wheel version %r (expected %r)\n" % (version, expected))
    sys.exit(1)
examples = [n for n in names if "/examples/" in n or n.endswith("/examples")]
if examples:
    sys.stderr.write("examples present in wheel:\n")
    for n in examples:
        sys.stderr.write("  %s\n" % n)
    sys.exit(1)
griffon = [n for n in names if "griffon" in n.lower()]
if griffon:
    sys.stderr.write("Griffon data present in wheel:\n")
    for n in griffon:
        sys.stderr.write("  %s\n" % n)
    sys.exit(1)
if not any("LICENSE" in n for n in names):
    sys.stderr.write("LICENSE missing from wheel\n")
    sys.exit(1)
if not any(n.endswith("pyp_newclasses.py") for n in names):
    sys.stderr.write("production module pyp_newclasses.py missing from wheel\n")
    sys.exit(1)
if not any(n.endswith("application/__init__.py") for n in names):
    sys.stderr.write("application package missing from wheel\n")
    sys.exit(1)
if not any(n.endswith("desktop/__init__.py") for n in names):
    sys.stderr.write("desktop package missing from wheel\n")
    sys.exit(1)
genes = [n for n in names if "GENES.html" in n or n.endswith("/GENES.html") or "/GENES" in n]
if genes:
    sys.stderr.write("stale GENES material present:\n")
    for n in genes:
        sys.stderr.write("  %s\n" % n)
    sys.exit(1)
print("examples: 0")
print("Griffon: 0")
print("LICENSE: present")
print("GENES: absent")
print("PRISTINE WHEEL OK")
PY
