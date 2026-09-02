#!/usr/bin/env bash
# Compatibility wrapper around tools/macos/build_app.py.
#
# Usage (from the repository root, with the macos-app extra installed):
#
#   ./tools/macos/build_app.sh
#   ./tools/macos/build_app.sh /tmp/pypedal-macos-app
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python "$HERE/build_app.py" "${1:-${TMPDIR:-/tmp}/pypedal-macos-app}"
