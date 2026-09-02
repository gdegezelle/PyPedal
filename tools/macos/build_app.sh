#!/usr/bin/env bash
# Build the 4.2-B macOS .app spike outside the repository.
#
# Usage (from the repository root, with the macos-app extra installed):
#
#   ./tools/macos/build_app.sh
#   ./tools/macos/build_app.sh /tmp/pypedal-macos-app
#
# Output is a PyPedal.app under <out>/dist/. Nothing is signed, notarized,
# or published. The canonical Griffon dataset is not included.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT="${1:-${TMPDIR:-/tmp}/pypedal-macos-app}"

mkdir -p "$OUT"

python -m PyInstaller \
  --noconfirm \
  --clean \
  --distpath "$OUT/dist" \
  --workpath "$OUT/work" \
  "$HERE/PyPedal.spec"

APP="$OUT/dist/PyPedal.app"
if [[ ! -d "$APP" ]]; then
  echo "expected $APP" >&2
  exit 1
fi

python - "$APP" <<'PY'
import os
import sys

app = sys.argv[1]
names = []
for root, dirs, files in os.walk(app):
    for name in files:
        names.append(os.path.join(root, name))
blob = "\n".join(names).lower()
if "griffon" in blob:
    sys.stderr.write("Griffon dataset leaked into the app bundle\n")
    sys.exit(1)
plugins = [n for n in names if "plugins" in n.lower() and "platforms" in n.lower()]
if not plugins:
    sys.stderr.write("Qt platform plugins not found in bundle\n")
    sys.exit(1)
print("app:", app)
print("files:", len(names))
print("qt platform plugins:", len(plugins))
PY

echo "MACOS APP SPIKE OK: $APP"
