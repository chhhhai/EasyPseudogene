# Advanced usage

The two comparison modes remain Cross and Self. This document covers optional
execution choices and evidence annotation for the development interface in
1.1.1-dev. Run `easypseudogene --help-advanced` for the complete option list.

## Configuration

Defaults are read from `config/default.conf`; `--config FILE` adds overrides.
Keep changes in a separate configuration file and use a fresh output directory.
The standard workflow records filtering settings in
`metadata/effective_config.txt`; fast workflows record their expanded settings.

For example, a configuration can disable family keyword filtering:

```bash
# experiment.conf
ENABLE_FAMILY_FILTER=0
```

```bash
easypseudogene --mode cross --proteins reference.fa --genome target.fa \
  --config experiment.conf --species cross_analysis --output cross_run
```

A final candidate count reflects the enabled rules. Compare results only with
known input scope and effective settings.

## Self workflow with MMseqs search and local structural validation

```bash
easypseudogene --mode self --workflow fast \
  --proteins species.proteins.fa --genome species.genome.fa \
  --gff species.annotation.gff3 \
  --species self_catalog --threads 16 --output catalog_run
```

The simple interface selects the MMseqs balanced preset and enables miniprot
structural validation. `--sensitivity sensitive` selects the higher-sensitivity
preset. The balanced settings retain longest isoforms, hard-mask CDS before
search, filter source/CDS overlaps and merge candidate loci.

Unified final candidates require validated miniprot disruptions and medium/high
sequence evidence. This path does not apply standard GeneWise family/copy
filters; its final count has a different selection meaning.

## Explicit search engines and presets

| Choice | Search and processing | Structural validation |
| --- | --- | --- |
| `mmseqs-balanced` | MMseqs2 translated search and BLAST6-style clustering | Enabled |
| `mmseqs-sensitive` | Increased MMseqs sensitivity and retained hit budget | Enabled |
| `tblastn-production-candidates-only` | CDS-masked, parallel tblastn and clustering | Not enabled by this preset |
| `tblastn-production-validated` | Same candidate search, with local validation | Enabled |
| `--self-catalog-engine fast --fast-catalog` | Integrated miniprot PAF catalog processing | Enable explicitly |
| `--fast-catalog` with the legacy engine | miniprot hits and heuristic catalog classification | Enable explicitly |

For explicit engines, `--self-catalog-validate-disruptions` enables local
validation. The simple `--workflow fast` interface enables it automatically.
Unvalidated hits remain screening candidates in unified results and do not
enter `results/final_loci.tsv`.

Legacy engine scores and the internal validated tiers retain their own
definitions under `work/fast_catalog/`. Use the common fields in `results/`
for the documented evidence grades; do not equate internal tier names with
those grades.

## Resource controls

- `--threads` is the total CPU budget for the tblastn runner.
- Explicit BLAST job/thread options must fit within that budget.
- Chromosome-parallel search, SQLite processing, alternate chaining, overlap
  bases and merged-span limits are expert controls shown in advanced help.
- More retained miniprot hits and lower `--miniprot-outs` can improve search
  coverage while increasing GeneWise work.

## Additional standard Cross evidence

- `--reference-annotation`: build reference annotations and assess protein
  neighborhoods. `--reference-gtf` remains a compatibility alias.
- `--reference-annotations`: provide an existing query-keyed annotation TSV;
  it is an alternative to `--reference-annotation`.
- `--redundancy-paf`: supply additional genome alignments for copy filtering;
  Cross otherwise automatically searches copies for disrupted queries.
- `--synteny-tsv`: supply independent expected target intervals. Enable
  `ENABLE_SYNTENY_FILTER=1` in the configuration to use them as a filter.

The automatic protein-neighborhood assessment is a proxy annotation after
filtering. It is separate from independent synteny filtering and does not
change the already selected final candidate set.

## Optional unitary evidence annotation

This is a separate command after candidate discovery:

```bash
easypseudogene-unitary \
  --unitary-candidates validated.tsv \
  --proteins reference.fa --target-proteins target.fa \
  --target-genome target.genome.fa --outgroup-proteins outgroup.fa \
  --unitary-synteny-bed synteny.bed \
  --threads 16 --species comparison --output unitary_run
```

The candidate table must satisfy the helper's input schema. Do not assume
that renaming another table creates a compatible input. A unitary-candidate
label requires the relevant intact-copy exclusions, structural disruption,
synteny and outgroup evidence. Missing evidence is recorded with `needs_*`
labels. These labels remain computational evidence, not experimental proof.
The historical `--mode unitary` interface remains a compatibility route.

## Optional GeneWise alignment page

Add `--visualize` to a standard analysis to write
`results/genewise_alignment.html`. The normal `report.html` summary is generated
without this flag. The fast catalog branches exit before the GeneWise-page step.

## Existing intermediate results

A complete run is reused only when its content fingerprints still match.
`--skip-genewise` and `--skip-filter` are expert options for appropriate existing
intermediates; they do not establish a complete result set in a fresh directory.
Do not relocate cached directories: input and intermediate paths are part of
provenance. Use a new output directory after upgrading.

See [methods](METHODS.md) for filter rules and [validation](VALIDATION.md) for
the actual scope of recorded checks.
