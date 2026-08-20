# Conformance suite

Language-neutral oracles for the `gprocrustes` kernel and problem grammar.

- Schema: [schema.json](schema.json)
- Fixtures: [fixtures/](fixtures/)
- Binary matrices: [matrices/](matrices/) (little-endian `float64`, column-major)
- Manifest: [manifest.json](manifest.json)
- Spec: [docs/spec/02-conformance-suite.md](../../docs/spec/02-conformance-suite.md)

Regenerate:

```text
python3 tools/generate_conformance.py
```

Do not invent Gower 1975 table digits. `gower-1975-published-history` stays a transcription placeholder until Table 1–2 are copied from Psychometrika 40:33–51.
