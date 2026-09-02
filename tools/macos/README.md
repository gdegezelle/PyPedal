# macOS application spike (PyPedal 4.2-B)

Engineering spike only. This is not the 4.2-D distribution package.

- Bundle name: `PyPedal.app`
- `CFBundleName`: PyPedal
- Bundle identifier: `org.pypedal.PyPedal`
- Console: false (GUI bundle)
- Signing/notarization: not performed in 4.2-B
- Canonical Griffon dataset: not included

Build (requires `pip install -e ".[macos-app]"`):

```bash
./tools/macos/build_app.sh /tmp/pypedal-macos-app
```

The `.app` is written outside the repository. Smoke the embedded
executable with `PYTHONPATH` unset, from a directory that is not the
checkout:

```bash
/tmp/pypedal-macos-app/dist/PyPedal.app/Contents/MacOS/PyPedal --version
```
