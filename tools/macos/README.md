# macOS application packaging (PyPedal 4.2-D)

Canonical commands:

```bash
python tools/macos/build_app.py /tmp/pypedal-macos-app
python tools/macos/build_dmg.py \
  --app /tmp/pypedal-macos-app/dist/PyPedal.app \
  --out /tmp/pypedal-macos-dmg
```

- Bundle name: `PyPedal.app`
- `CFBundleName` / `CFBundleDisplayName`: PyPedal
- Bundle identifier: `org.pypedal.PyPedal`
- Console: false (GUI bundle)
- Version strings: `[project].version` in `pyproject.toml`
- Canonical Griffon dataset: not included
- CustomTkinter: excluded from the `.app` (source extra unchanged)
- Icon: `tools/macos/PyPedal.icns` when a maintainer supplies it

See `docs/engineering/macos-distribution.md` for signing/notarization
status and the human Finder checklist.
