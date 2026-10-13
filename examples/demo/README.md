# Synthetic installation demo

This small, bundled dataset checks the standard Cross and Self workflows without
downloading external sequences. It exercises input preparation, miniprot,
GeneWise, filtering, result exports, report links and GeneWise path discovery.
It is an installation check, not a biological accuracy benchmark.

## See a detected candidate

The [merged report](expected/report.html) contains statistics, final-candidate
alignments with highlighted disruptions, filtering diagnostics and downloads. The candidate table shows 20 rows per
page, supports search across all candidates, and expands alignments below the
selected row. All alignments start collapsed. Five GeneWise statistics and
seven overview charts are above the table. Use the scope selector to compare
all completed GeneWise alignments with final candidate alignments; changing
the table page or search does not change the charts.
It uses results from an actual execution of this synthetic demo. Open the HTML
in a browser, or run the demo below to reproduce the analysis.

| Candidate | Stops | Frameshifts | Query coverage | Amino-acid identity |
| --- | --- | --- | --- | --- |
| `demo_A.t1` | 1 | 1 | 100% | 99.33% |

The report contains only this final candidate. Click its candidate ID to see the red
stop marker and orange frameshift marker. The intact `demo_B.t1` control is
excluded. The [raw GeneWise candidate output](expected/genewise_candidate.txt)
and TSV/BED/GFF3 files are also included under `expected/`.

## Run both modes

Activate an environment containing Python 3, miniprot, wise2 and samtools.
From the matching EasyPseudogene source checkout:

```bash
conda activate easypseudogene
bash examples/demo/run_demo.sh
```

The default uses two threads and writes a new `demo_output/` directory in your
current directory. No `--wiseconfig` or `--wisepath` flags are required for a
standard installation. Automatic discovery must be able to locate both the
GeneWise executable and configuration. Custom installations can supply
`WISEPATH` and `WISECONFIGDIR` as environment variables.

For another output directory or CPU budget:

```bash
bash examples/demo/run_demo.sh my_demo_output 4
```

Existing output directories are refused. The script uses the entry point in
this checkout, so it will not accidentally call an older version from `PATH`.
It stops on failure and prints **`DEMO PASS`** only after the result checks pass.

## Data and expected results

| Input | Size/content |
| --- | --- |
| `data/proteins.fa` | Two unrelated synthetic proteins, 300 amino acids each |
| `data/cross_genome.fa` | Two contigs, 4,205 bases in total |
| `data/self_genome.fa` | Three contigs, 6,308 bases in total |
| `data/self_annotation.gff3` | Two intact source genes with mRNA/CDS relationships |
| `data/truth.json` | Synthetic origin, seed, engineered events and input SHA256 hashes |

| Mode | Expected result |
| --- | --- |
| Cross | One final candidate for `demo_A.t1`, with one stop and one frameshift; intact `demo_B.t1` is excluded because it has no disruption |
| Self | One off-annotation disrupted locus is analyzed, then excluded by the standard full-length homolog filter; zero final candidates |

The Self result is an intentional negative control. Both reference proteins have
intact annotated source loci in the Self genome. This demonstrates annotation
masking and the standard copy-filter behavior; zero final candidates alone is
not sufficient for the demo to pass.

GeneWise versions may differ in event counts or alignments. The checker requires
the engineered stop and frameshift to be recovered, the expected candidate
queries and decisions, matching export coordinates and complete outputs. A
failure needs inspection; do not relax filters merely to obtain `DEMO PASS`.

The 3.0.13-dev combined demo was checked on 2026-10-13 with an existing Conda
installation, automatic GeneWise discovery and two threads. One recorded run
completed in 7 seconds; 104 regression tests passed. Runtime on other systems
can differ.

## Outputs

```text
demo_output/
├── demo.log
├── demo_check.json
├── cross/
│   ├── report.html
│   ├── results/
│   ├── qc/
│   ├── metadata/
│   ├── logs/
│   └── work/
└── self/
    ├── report.html
    ├── results/
    ├── qc/
    ├── metadata/
    ├── logs/
    └── work/
```

Open `cross/report.html` to inspect the detected candidate and run statistics,
or `self/report.html` for the negative control. A separate final-candidate
viewer is no longer produced. The JSON check record stores `passed` or `failed`.

The report, candidate search and alignment details work offline, with no CDN
or external JavaScript dependency. The standalone `easypseudogene-visualize`
command remains available for advanced inspection of raw GeneWise outputs.

Recheck a completed demo:

```bash
python examples/demo/check_results.py --output demo_output
```

## Run one mode manually

From the source checkout, choose a fresh output directory:

```bash
bash bin/easypseudogene --mode cross \
  --proteins examples/demo/data/proteins.fa \
  --genome examples/demo/data/cross_genome.fa \
  --threads 2 --species cross --output cross_demo --visualize

bash bin/easypseudogene --mode self \
  --proteins examples/demo/data/proteins.fa \
  --genome examples/demo/data/self_genome.fa \
  --gff examples/demo/data/self_annotation.gff3 \
  --max-hits-per-protein 20 --miniprot-outs 0.1 \
  --threads 2 --species self --output self_demo
```

The combined checker expects both modes in a shared output root, as produced by
`run_demo.sh`. These separate manual commands are for inspecting individual runs.

## Provenance and regeneration

All sequences are synthetic, generated with the Python standard library and
seed `20261013`. They are distributed under the repository's MIT license. The
generator inserts an internal stop at amino acid 140 and deletes one base at
coding-sequence position 646 (0-based) in the disrupted A copy. No real species
or experimentally verified pseudogene is represented.

The bundled files are ready to use. Maintainers can regenerate them with:

```bash
python examples/demo/generate_data.py
```

This rewrites the bundled data and truth manifest. Normal users do not need to
run the generator.
