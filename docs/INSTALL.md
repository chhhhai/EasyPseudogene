# Installation

Use the source version matching the documentation. These instructions describe
the interface in 1.1.0; earlier versions use different commands
and output layouts.

## Platform and core dependencies

Use Linux with Bash and standard Unix tools, including `awk`, `sort`, `comm`,
`xargs`, `tee` and `flock`. Python 3.11 is recommended. Other operating systems
have not been validated in the recorded real-data runs.

| Execution path | Required analysis tools |
| --- | --- |
| Standard Cross/Self | miniprot, GeneWise, samtools, Python |
| Protein-input preparation with `--check-inputs` | Python |
| Annotated reference-genome translation | gffread, in addition to the selected path |
| Optional MMseqs workflow/diagnostics | MMseqs2; miniprot for structural validation |
| Optional tblastn workflow | BLAST+; miniprot when structural validation is enabled |

## Minimal standard environment

```bash
conda create -n easypseudogene -c conda-forge -c bioconda \
  python=3.11 miniprot wise2 samtools gffread
conda activate easypseudogene
```

The `wise2` package provides `genewise`. Optional search tools can be installed
in the same environment:

```bash
conda install -n easypseudogene -c conda-forge -c bioconda mmseqs2 blast
```

The root `environment.yml` contains the minimal standard environment above.
Alternatively, create it from the matching source checkout with
`conda env create -f environment.yml`. Advanced search tools and plotting tools
are optional and are installed separately.

## Source checkout

From the matching checkout:

```bash
cd EasyPseudogene
export PATH="$PWD/bin:$PATH"
```

If executable permissions were lost while transferring the source:

```bash
chmod +x bin/easypseudogene bin/easypseudogene-unitary scripts/*.sh
```

Alternatively, invoke the main entry point as `bash bin/easypseudogene`.
No edits to input FASTA or annotation IDs are needed for installation.

## GeneWise configuration

GeneWise needs both an executable and configuration data. Installing a package
named `wise2` does not mean its executable is named `wise2`: the command is
`genewise`. Activate the environment where the package was installed.

The standard workflow resolves these before staging the input data:

1. Use `--wisepath` or `WISEPATH` when supplied; otherwise search `PATH`, then
   `$CONDA_PREFIX/bin/genewise` in the active environment.
2. Use `--wiseconfig` or `WISECONFIGDIR` when supplied; otherwise inspect the
   executable prefix and symlink target, the active Conda prefix, and
   `/usr/local` and `/usr`. Configuration candidates include
   `share/wise2/wisecfg`, versioned `share/wise*` directories and `wisecfg`.
   Source installations with `src/bin/genewise` are also covered.
3. Require an executable `genewise` and readable, nonempty `codon.table`,
   `BLOSUM62.bla` and `gene.stat` in the selected configuration directory.

Explicit paths take priority. Invalid explicit paths fail with a specific
error instead of silently switching installations. The resolved executable and
configuration are logged; the executable and configuration content are included
in the run cache identity. These file checks do not prove binary compatibility
or validate the contents of the configuration files.

For custom installations, or when automatic discovery fails, provide:

```bash
easypseudogene --mode cross --proteins reference.fa --genome target.fa \
  --wiseconfig /path/to/wisecfg --wisepath /path/to/genewise/bin \
  --species cross_analysis --output cross_run
```

`--wisepath` expects a directory containing `genewise`; `--wiseconfig` expects
the directory containing its configuration files. Either override can be used
alone when the other component is detected automatically. The standalone
`scripts/run_genewise.sh` runner shares the same detection through the environment
variables `WISEPATH` and `WISECONFIGDIR`.

Input-only checks and workflows that skip GeneWise do not require this setup.
Automatic discovery does not search every installed Conda environment. If wise2
was installed in another environment, activate that environment or supply its
paths explicitly.

## Installation checks

```bash
easypseudogene --version
easypseudogene --help
```

These display the interface and version; they do not validate biological
accuracy or exercise every dependency. Use `--check-inputs` with your own files
for preparation, then run a small representative analysis.

## Resources and troubleshooting

Memory, disk space and running time depend on genome size, search sensitivity,
retained hit count and candidate complexity. Specify a CPU budget with
`--threads`; the default is at most 16 and respects recognized allocations.

- Missing executable: activate the intended environment and check its PATH.
- Annotation mismatch: use genome and annotation from the same assembly;
  contig names are not guessed.
- Unresolved protein association: inspect `protein_id_mapping.tsv`; provide
  `--protein-id-map` for renamed records.
- Cache refusal: use a fresh output directory after changes or an incomplete run.
- GeneWise failure: inspect `logs/run.log`, per-locus errors and `genewise_failed.tsv`.
  Missing executable and missing configuration data are reported separately.
  A resolved executable can still fail because of missing shared libraries or
  incompatible configuration data; inspect the per-locus error in that case.

The source genome is staged under `work/inputs/` so indexing can occur locally.
Start with the [quick start](QUICKSTART.md) after installation.
