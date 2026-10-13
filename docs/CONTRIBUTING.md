# Contributing

## Development setup

From the repository root, create the Conda environment from
[`environment.yml`](../environment.yml), then run:

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
python -m compileall -q lib scripts tests
bash -n bin/easypseudogene
```

Additional lint and security audits can be run with separately installed
`ruff` and `bandit`. The GitHub source workflow uses Python compilation,
regression tests, shell syntax and CLI metadata checks; it does not run the
GeneWise installation demo or biological benchmarks.

## Pull requests

1. Keep changes focused and document scientific semantics.
2. Add regression tests before changing coordinates, matching, caching, evidence tiers, or subtype labels.
3. Preserve one-to-one maximum-cardinality benchmark matching.
4. Do not describe sensitive candidates as experimentally validated pseudogenes.
5. Include exact commands and real outputs for new full-scale benchmarks.
6. Update VERSION, README and docs/CHANGELOG.md when schemas or cache semantics change.

## Reporting issues

Include the command, version, input dialect (for example, BRAKER, MAKER,
Ensembl or NCBI), logs, and a minimal non-sensitive fixture. Never attach
credentials or private genomic data without authorization.
