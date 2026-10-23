# Quick start

For a small installation check using bundled synthetic data, activate your
environment and run from the source checkout:

```bash
bash examples/demo/run_demo.sh
```

It runs both comparison modes with two threads and prints `DEMO PASS` after
checking the expected results. See [Demo details](../examples/demo/README.md).

These commands describe the interface in 1.2.0. Use the
matching source checkout and a fresh output directory. See [installation](INSTALL.md)
if the executable or its dependencies are not available.

## Cross

```bash
easypseudogene --mode cross \
  --proteins reference.proteins.fa --genome target.genome.fa \
  --reference-annotation reference.gff3 \
  --species cross_analysis --threads 16 --output cross_run
```

Reference annotation is optional when proteins are supplied. For an annotated
reference genome, replace `--proteins` with `--reference-genome reference.fa`;
reference annotation and gffread are then required.

## Self

```bash
easypseudogene --mode self \
  --proteins species.proteins.fa --genome species.genome.fa \
  --gff species.annotation.gff3 \
  --max-hits-per-protein 20 --miniprot-outs 0.1 \
  --species self_analysis --threads 16 --output self_run
```

The broader search settings retain alternative loci before annotation filtering.
Historical defaults select one hit per query and use `--miniprot-outs 0.99`.
Those defaults can retain only the query's own annotated gene in Self mode.

## Input-only preparation

Append `--check-inputs` to either command. Inspect
`OUTPUT/NAME/work/input_check/inputs/input_qc.json`,
`protein_id_mapping.tsv` and `protein_exclusions.tsv`. Remove the flag to start
analysis in the same output directory.

Plain/gzip FASTA and GTF/GFF3 are accepted. Unique query IDs and assembly-matched
annotation are required. Invalid reference proteins are reported and excluded
by default; `--invalid-proteins error` selects strict rejection.

## View results

Open `OUTPUT/NAME/report.html`, then inspect:

- `results/final_loci.tsv`: final computational candidates.
- `results/final_loci.bed` and `.gff3`: candidate intervals.
- `results/loci.tsv`: all selected loci and decisions.
- `results/queries.tsv`: accepted queries, including queries without selected loci.
- `qc/nonfinal_loci.tsv`: exclusion and unresolved reasons.
- `metadata/run_manifest.json`: workflow and provenance.

TSV/BED are 0-based, half-open; GFF3 is 1-based, inclusive. Evidence grades are
separate from selection and are not biological confidence probabilities.

Use `--help` for common options. Optional search engines and follow-up evidence
annotation are documented in [advanced usage](ADVANCED_USAGE.md).

## Optional cleanup after success

Add `--cleanup-intermediates` to either analysis command to remove generated
work files after successful publication. Default runs retain intermediates.
The final report, result tables, QC, metadata and logs remain available.

For an existing completed current-layout run:

```bash
easypseudogene-cleanup --run-dir OUTPUT/NAME --dry-run
easypseudogene-cleanup --run-dir OUTPUT/NAME
```

Cleaned runs retain `metadata/cleanup.json` and, when applicable,
`results/final_genewise.txt`. Use a new output directory for another analysis.
See [cleanup details](../README.md#intermediate-file-cleanup).
