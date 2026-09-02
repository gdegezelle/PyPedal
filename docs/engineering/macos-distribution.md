# macOS distribution (PyPedal 4.2-D)

Engineering notes for building `PyPedal.app` and `PyPedal.dmg`. This file
is not part of the user manual and is not in MkDocs navigation.

The Python package version remains **4.1.0** until 4.2-E. CustomTkinter
launchers (`pypedal`, `pypedal-gui`, `python -m PyPedal`) are unchanged.

## Automated now

From a checkout with `pip install -e ".[macos-app]"`:

```bash
python tools/macos/build_app.py /tmp/pypedal-macos-app
python tools/macos/build_dmg.py \
  --app /tmp/pypedal-macos-app/dist/PyPedal.app \
  --out /tmp/pypedal-macos-dmg
```

`./tools/macos/build_app.sh` is a compatibility wrapper around the same
Python entry.

The spec reads `[project].version` from `pyproject.toml`. Rebuilding at
4.2.0 does not require a separate spec edit if that field (and
`PyPedal/__version__.py`) are bumped together.

Default output is outside the repository. Do not commit `.app` or `.dmg`
artifacts.

### Bundle identity

| Key | Value |
|---|---|
| `CFBundleName` | PyPedal |
| `CFBundleDisplayName` | PyPedal |
| `CFBundleIdentifier` | `org.pypedal.PyPedal` |
| `CFBundleExecutable` | PyPedal |
| `CFBundlePackageType` | APPL |
| Version strings | project version |

There is no `.ped` document-type registration.

### Icon

Place an authoritative `tools/macos/PyPedal.icns` when one exists. The
build wires it in automatically. No project-owned logo is in the tree
today; do not invent one. See the D report: **FINAL APPLICATION ICON
ASSET REQUIRED**.

### Smoke

From a directory that is not the checkout, with `PYTHONPATH` unset:

```bash
cd /tmp
/tmp/pypedal-macos-app/dist/PyPedal.app/Contents/MacOS/PyPedal --version
/tmp/pypedal-macos-app/dist/PyPedal.app/Contents/MacOS/PyPedal \
  --self-test /absolute/path/to/tiny.ped
```

`--self-test` loads the pedigree and runs Meuwissen–Luo, Lacy founders,
relationship, mating CoI, and theoretical Ne with `output=False`. It does
not start the GUI.

A frozen app changes its working directory to
`~/Library/Application Support/PyPedal` so library log files are not
written relative to Finder’s cwd. Pedigree paths remain absolute.

## Maintainer credentials (not in this tree)

Developer ID signing and Apple notarization are **not** performed by the
build scripts. They require a local signing identity and Apple notary
credentials that must never be committed.

Expected release-time sequence (placeholders only):

1. Sign `PyPedal.app` with **Developer ID Application**
   (`Developer ID Application: Example Org (TEAMID)`).
2. Enable hardened runtime where Apple’s current guidance requires it.
3. Build `PyPedal.dmg` from the signed app.
4. Submit the disk image with `notarytool submit` using an App Store
   Connect API key or an equivalent notary profile.
5. Staple the ticket (`xcrun stapler staple`).
6. Confirm Gatekeeper with `spctl --assess` on a clean Mac.

Do not store Apple IDs, app-specific passwords, API private keys, or
keychain passwords in the repository.

## Deferred to release time (4.2-E)

- Version bump to 4.2.0 and rebuild of app + DMG
- Application icon asset, if still missing
- Developer ID signing and notarization, if a distribution identity exists
- Finder double-click human confirmation (checklist in the D report)
- Public launcher cutover away from CustomTkinter

An unsigned or ad-hoc-signed development app is expected to fail
Gatekeeper. That is a signature-state result, not an invitation to disable
macOS security settings.
