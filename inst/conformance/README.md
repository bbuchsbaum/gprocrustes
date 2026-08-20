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

`gower-1975-published-history` transcribes Table 2 (carcass scores) and Table 5 (successive \(S_r\)) from Psychometrika 40:33–51. Table 1 is the ANOVA layout. Do not invent further published digits.
