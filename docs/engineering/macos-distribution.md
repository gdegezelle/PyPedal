# macOS distribution (PyPedal 4.2)

Engineering notes for building `PyPedal.app` and `PyPedal.dmg`. This file
is not part of the user manual and is not in MkDocs navigation.

Official launchers (`pypedal`, `pypedal-gui`, `python -m PyPedal`) start
the PySide6 desktop. CustomTkinter is not a runtime dependency.

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

The spec reads `[project].version` from `pyproject.toml`.

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
today; do not invent one. **FINAL APPLICATION ICON ASSET STILL REQUIRED
FOR PROFESSIONAL BINARY RELEASE.**

### Smoke

From a directory that is not the checkout, with `PYTHONPATH` unset:

```bash
cd /tmp
/tmp/pypedal-macos-app/dist/PyPedal.app/Contents/MacOS/PyPedal --version
/tmp/pypedal-macos-app/dist/PyPedal.app/Contents/MacOS/PyPedal \
  --self-test /absolute/path/to/tiny.ped
```

A frozen app changes its working directory to
`~/Library/Application Support/PyPedal`. Pedigree paths remain absolute.

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

## Deferred to a signed binary release

- Application icon asset, if still missing
- Developer ID signing and notarization, if a distribution identity exists
- Finder double-click human confirmation
- Public macOS download publication

An unsigned or ad-hoc-signed development app is expected to fail
Gatekeeper. That is a signature-state result, not an invitation to disable
macOS security settings.
