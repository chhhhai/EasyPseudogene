# EasyPseudogene

**Pseudogene candidate discovery with two comparison modes: Cross and Self.**

EasyPseudogene compares protein sequences with a target genome, identifies
premature stops and frameshifts, and reports candidates together with their
sequence evidence and filtering decisions. The standard workflow uses
miniprot and GeneWise.

> Version: **1.2.0**. Use this documentation with the matching
> source checkout; earlier versions have different commands and output layouts.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Installation demo](#installation-demo)
- [Quick start](#quick-start)
- [Example: human proteins against a dolphin genome](#example-human-proteins-against-a-dolphin-genome)
- [Input requirements](#input-requirements)
- [Outputs](#outputs)
- [Intermediate-file cleanup](#intermediate-file-cleanup)
- [Workflow and interpretation](#workflow-and-interpretation)
- [Documentation](#documentation)
- [Citation, license and support](#citation-license-and-support)

## Features

- **Cross:** search a target genome using reference proteins from another species.
- **Self:** search a genome using its own proteins and matching annotation.
- **Generic inputs:** plain/gzip FASTA and GTF/GFF3, including BRAKER-style
  coding models and MAKER-style relationships.
- **Traceable results:** retain query IDs, map queries to candidate loci, and
  record input exclusions, filter outcomes and effective parameters.
- **Structured outputs:** candidate TSV/BED/GFF3, an offline summary report,
  input QC and reproducible run metadata.

## Installation

On Linux, clone the repository and create the standard Conda environment:

```bash
git clone https://github.com/chhhhai/EasyPseudogene.git
cd EasyPseudogene
conda env create -f environment.yml
conda activate easypseudogene
export PATH="$PWD/bin:$PATH"
```

The entry point can also be run as `bash bin/easypseudogene` from the checkout.
Bash and standard Unix utilities, including `flock`, must be available.
GeneWise is located through `PATH`, then the active Conda environment's `bin`
directory. Its configuration is detected from the binary installation prefix,
symlink target, active Conda environment and common system locations.
For a custom installation, keep these explicit overrides:

```bash
easypseudogene --mode cross --proteins proteins.fa --genome genome.fa \
  --threads 96 --wiseconfig /path/to/wisecfg \
  --wisepath /path/to/genewise/bin --output cross_run
```

`--wisepath` locates the executable; `--wiseconfig` locates its data files.
Configuration checks happen before analysis and report missing files explicitly.
Further setup and optional dependencies are covered in [INSTALL.md](docs/INSTALL.md).
`gffread` is only needed when translating an annotated reference genome.

## Installation demo

After activating your environment, run this from the source checkout:

```bash
bash examples/demo/run_demo.sh
```

The bundled synthetic dataset has two proteins and approximately 10.5 kb of
genome sequence across the two modes. It runs Cross and Self with two threads,
checks the expected candidates, TSV/BED/GFF3 coordinates and report links, and
prints **`DEMO PASS`** on success. No sequence downloads or GeneWise path flags
are needed for a standard installation.

Open `demo_output/cross/report.html` to inspect the detected `demo_A.t1`
candidate: one verified stop, one frameshift, 100% query coverage and 99.33%
amino-acid identity. A generated [merged example report](examples/demo/expected/report.html) includes
statistics, a searchable candidate table (20 per page), expandable alignments
and result downloads in one page. The GeneWise overview includes five
statistics and seven offline charts, with an all-alignments/final-candidates
scope selector. Chart statistics are independent of table search and paging.

Open `demo_output/self/report.html` for the Self control. Use a
fresh output directory for another run, for example
`bash examples/demo/run_demo.sh demo_output_2`. This checks installation and
execution, not biological accuracy. See [Demo details](examples/demo/README.md).

## Quick start

### Cross: another species' proteins against the target genome

```bash
easypseudogene --mode cross \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --reference-annotation reference.gff3 \
  --species cross_analysis --threads 16 --output cross_run
```

`--reference-annotation` is optional when proteins are supplied. Provide it
when available to support reference annotation filters and neighborhood evidence.

### Self: proteins against their own genome

```bash
easypseudogene --mode self \
  --proteins species.proteins.fa \
  --genome species.genome.fa \
  --gff species.annotation.gff3 \
  --max-hits-per-protein 20 --miniprot-outs 0.1 \
  --species self_analysis --threads 16 --output self_run
```

This example retains additional hits so the query's annotated gene does not
occupy the only selected hit. The historical defaults remain one selected hit
per query and `--miniprot-outs 0.99`; they can miss off-annotation homologs in
Self mode. Broader settings increase downstream work.

### Check inputs before starting an analysis

Append `--check-inputs` to either command. Preparation reports are written to
`OUTPUT/NAME/work/input_check/inputs/`; no search or GeneWise analysis is started.
Run again without the flag to analyze the same inputs.

Threads are optional: the default is at most 16, limited by CPU affinity and
recognized scheduler allocations. Use `--help` for common options and
`--help-advanced` for additional controls.

## Example: human proteins against a dolphin genome

This Cross example uses [Ensembl release 115 human sequences](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/)
and the [NCBI dolphin assembly GCF_011762595.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011762595.2/).
Run the commands from the EasyPseudogene source checkout, with the standard
dependencies installed and `easypseudogene` available on `PATH`.

```bash
epg_source="$PWD"
mkdir -p human_dolphin_example
cd human_dolphin_example

# Download and decompress the reference proteins and CDS.
ensembl_base=https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens
wget -O human.pep.fa.gz \
  "$ensembl_base/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
wget -O human.cds.fa.gz \
  "$ensembl_base/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"
gzip -dc human.pep.fa.gz > human.pep.fa
gzip -dc human.cds.fa.gz > human.cds.fa

# Select one protein per gene, using the longest CDS with a matching protein.
python "$epg_source/examples/ensembl_select_longest_cds_then_protein.py" \
  --cds human.cds.fa --pep human.pep.fa --out-prefix human

# Download and decompress the target genome.
ncbi_base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595
assembly=GCF_011762595.2_mTurTru1.mat.Y
wget -O dolphin.fa.gz "$ncbi_base/$assembly/${assembly}_genomic.fna.gz"
gzip -dc dolphin.fa.gz > dolphin.fa

easypseudogene --mode cross \
  --proteins human.longest_cds.proteins.fa \
  --genome dolphin.fa \
  --species human_dolphin --threads 96 \
  --output result --visualize
```

Adjust `--threads` to your available resources. The selector also writes
`human.longest_cds.tsv`, recording the selected gene, transcript and protein IDs.
Reference annotation is optional for this protein-based Cross run; without it,
annotation-dependent filters and reference-neighborhood evidence are unavailable.

Open `result/human_dolphin/report.html` for statistics and final-candidate
alignments. `--visualize` remains accepted for compatibility; a second viewer
is no longer generated by the standard workflow. Final candidates are in
`result/human_dolphin/results/final_loci.tsv`.
This is a full-data usage example, not a validated biological benchmark; it has
not been rerun for this documentation revision.

## Input requirements

| Input | Cross | Self | Accepted format |
| --- | --- | --- | --- |
| Target genome | Required | Required | DNA FASTA, plain or gzip |
| Protein sequences | Required, or generated from an annotated reference genome | Required | Unaligned protein FASTA, plain or gzip |
| Reference annotation | Optional with proteins; required for reference-genome translation | Not used | GTF/GFF3, plain or gzip |
| Target annotation | Not required for standard Cross | Required | GTF/GFF3, plain or gzip |
| Explicit ID mapping | Optional | Optional | TSV via `--protein-id-map` |

### FASTA and annotation rules

- Use unique first-token FASTA IDs. Ensembl names and `gene_symbol` fields
  are not required.
- Self annotation and genome must come from the same assembly, with matching
  contig IDs and valid coordinates. Reference-genome translation has the same
  requirement for its reference files.
- Coding models should identify transcripts/proteins and their genes through
  GFF3 `ID`/`Parent`, GTF `gene_id`/`transcript_id`, or explicit protein IDs.
  Missing containers can be inferred from supported coding relationships.
- Keep protein IDs associated with their annotation. For renamed IDs, provide
  an explicit mapping; inspect unresolved associations in the mapping report.
- Terminal protein `*` markers are removed. Empty sequences, internal stops
  and invalid amino-acid symbols are excluded and reported by default.
  `--invalid-proteins error` selects strict rejection; duplicate IDs fail.

For example, a BRAKER reference can be supplied directly:

```bash
easypseudogene --mode cross --proteins braker.aa \
  --reference-annotation braker.gtf --genome target.fa \
  --species braker_cross --output cross_run
```

Reference-genome translation is an alternative Cross input route:

```bash
easypseudogene --mode cross --reference-genome reference.fa \
  --reference-annotation reference.gff3 --genome target.fa \
  --species genome_cross --output cross_run
```

Use either `--proteins` or `--reference-genome`. Translation uses gffread and
selects one valid longest representative protein per gene.

## Outputs

Open the English `report.html` first, then inspect `results/final_loci.tsv`.

```text
OUTPUT/NAME/
├── report.html
├── results/
│   ├── final_loci.tsv
│   ├── final_loci.bed
│   ├── final_loci.gff3
│   ├── final_queries.txt
│   ├── loci.tsv
│   ├── queries.tsv
│   └── summary.json
├── qc/
├── metadata/
├── logs/
└── work/
```

| File | Purpose |
| --- | --- |
| `report.html` | Offline statistics, paginated candidate table, alignments, QC and downloads |
| `results/final_loci.tsv` | Final computational candidates and evidence fields |
| `results/final_loci.bed` / `.gff3` | Final candidate coordinates |
| `results/loci.tsv` | All selected loci, decisions and exclusion reasons |
| `results/queries.tsv` | Every accepted input query, including queries without selected loci |
| `qc/input_qc.json` | Input preparation status and statistics |
| `qc/protein_exclusions.tsv` | Invalid reference proteins excluded during preparation |
| `qc/nonfinal_loci.tsv` | Loci outside the final candidate set |
| `metadata/run_manifest.json` | Version, workflow, file paths and counts |
| `metadata/reference_cache.json` | Input, code, tool and artifact fingerprints |
| `logs/run.log` | Analysis log |
| `work/` | Staged inputs and intermediates for the stages actually used |

TSV/BED intervals are **0-based, half-open**; GFF3 intervals are **1-based,
inclusive**. The final GFF3 contains one candidate interval feature per locus.
Optional evidence tables may also appear.

Use a fresh output directory after changing code, inputs or parameters.
A verified complete run can be reused without rewriting its published results;
incomplete, modified or incompatible runs are refused. Zero candidates is a
valid completed result.

## Intermediate-file cleanup

Intermediate files are retained by default. Add `--cleanup-intermediates` to
Cross or Self runs to remove generated work files **after successful analysis
and publication**:

```bash
easypseudogene --mode cross \
  --proteins reference.proteins.fa --genome target.genome.fa \
  --threads 16 --species sample --output result \
  --cleanup-intermediates
```

The report, all published result tables, QC, metadata and logs are retained.
For standard runs with final candidates, their original GeneWise output is
preserved in `results/final_genewise.txt`. Report charts and expandable
alignments are already embedded in `report.html` and remain usable offline.

For an existing completed run, preview cleanup first, then run it:

```bash
easypseudogene-cleanup --run-dir result/sample --dry-run
easypseudogene-cleanup --run-dir result/sample
```

The helper validates the completed run and published-file fingerprints and
acquires both run locks. It refuses active/incomplete runs, modified outputs,
snapshots, unsupported legacy layouts, and work directories containing source
inputs. `--dry-run` writes and deletes nothing. No analysis tools or original
sequence inputs are needed for standalone cleanup.

Cleanup removes known generated stage directories under `work/`, including
staged input copies, search databases, extracted regions, temporary files and
per-locus GeneWise files. Original inputs outside the run are retained, as are
unrecognized entries under `work/`. A minimal `work/` directory remains. Paths,
logical byte counts and cleanup status are recorded in `metadata/cleanup.json`.
Symlinks within selected work directories are unlinked without deleting their
targets; mounted directories are refused.

A cleaned directory is a final-output archive: intermediate-stage continuation
and same-directory analysis reuse are unavailable. Use a new output directory
for another analysis. If cleanup is interrupted or reports `partial`, inspect
`metadata/cleanup.json` and repeat the standalone command after resolving the
reported problem. Completed output integrity is verified again before deletion.

Cleanup reduces storage **after completion**; it does not reduce peak disk
usage while the analysis is running. Logical byte counts are not a guarantee
of physically reclaimed disk space, for example when hard links are present.

Cleanup has been checked with fresh synthetic analyses, 18 runtime checks and
13 dedicated regression tests. See [validation scope](docs/VALIDATION.md#optional-cleanup-revision).

## Workflow and interpretation

The standard workflow follows this order:

1. Normalize inputs and check associations.
2. Run miniprot and select candidate hits.
3. In Self mode, filter hits overlapping the query's own or annotated loci.
4. Prepare per-locus jobs, extract DNA and run GeneWise.
5. Count stops/frameshifts, select representative proteins, and apply the
   enabled family, annotation, copy and alignment-quality filters.
6. In Cross mode, check same-locus intact-alignment conflicts and assess
   reference neighborhoods when the corresponding inputs are available.
7. Publish candidates, QC, metadata and the summary report.

The standard copy filters can remove a disrupted locus when another
full-length homolog exists. This behavior can exclude duplicated pseudogene
candidates; its conditions and alternative workflows are described in
[Methods](docs/METHODS.md) and [Advanced usage](docs/ADVANCED_USAGE.md).

Results use `high`, `medium`, `low` and `unassessed` **sequence-evidence grades**.
Grades are separate from filter decisions and are not calibrated probabilities
of true pseudogene status.



## Documentation

- [Quick start](docs/QUICKSTART.md)
- [Installation and dependencies](docs/INSTALL.md)
- [Methods and filtering rules](docs/METHODS.md)
- [Advanced workflows and configuration](docs/ADVANCED_USAGE.md)
- [Validation scope and limitations](docs/VALIDATION.md)
- [Version history](docs/CHANGELOG.md)

## Citation, license and support

Record the EasyPseudogene version, input assemblies and effective parameters
in your methods. Cite the tools actually used:

- miniprot: Li (2023), *Bioinformatics*.
- GeneWise: Birney et al. (2004), *Genome Research*.
- MMseqs2, when used: Steinegger and Söding (2017), *Nature Biotechnology*.

EasyPseudogene is distributed under the [MIT license](LICENSE). Report issues
through [GitHub Issues](https://github.com/chhhhai/EasyPseudogene/issues),
including the software version, command and relevant logs.
