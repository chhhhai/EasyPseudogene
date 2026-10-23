# Validation scope

Execution checks and biological accuracy are separate questions. The records
below apply to their frozen software versions and input scopes. They do not
establish a universal precision, recall or full literature reproduction claim.

## Recorded 3.0.10-dev validation

The frozen 3.0.10-dev source passed **93 regression tests and 20 runtime checks**
on the project compute environment. Runtime analysis used comput4 with 32
threads. Source inputs and older analysis directories were preserved.

| Check | Scope | Recorded result |
| --- | --- | --- |
| T. rubra input preparation | Original BRAKER proteins and GTF/GFF3 | 20,669 records; 20,668 accepted; one internal-stop record excluded |
| T. dohrnii input preparation | Original BRAKER proteins and annotation | 30,536 protein records prepared |
| Cross GTF/GFF3 runs | 200 original reference genes against the whole user-identified T. dohrnii genome | 191 selected loci; 2 final candidates; decisions match frozen prior results |
| Standard Self, broader search | 40 original proteins against the whole target; 20 selected hits/query and `outs=0.1` | 34 off-annotation loci; 0 final candidates after standard filters |
| MMseqs Self catalog | Same 40 proteins against the whole target | 796 screening loci; 20 final structural-evidence candidates representing 11 queries |
| Same-genome translation control | Real contig/GTF; translated representative proteins | 35 alignments; 0 final candidates |
| Empty analyses | Standard and catalog paths | Valid empty results, report and GFF3 produced |
| Directory structure | Eight real-data output inventories | Five visible root directories and report; no unused GeneWise placeholders in catalog output |
| Export consistency | Published candidates | TSV/BED/GFF3 candidate IDs and coordinates agree |
| Offline report | Generated relative file links | All checked links resolve |
| Cache protection | Same invocation, changed parameters, edited report/GFF3 | Complete reuse preserves published content; incompatible or changed artifacts refused |

The standard and catalog final counts differ because their selection rules
differ. The catalog calls are computational candidates requiring review.

## Data and compatibility limits

- Both species' inputs were real BRAKER3 predictions, not synthetic formatting
  examples. The tested target identity and paths were supplied by the user.
- Cross discovery sampled 200 reference genes; Self discovery sampled 40
  proteins. Whole target genomes were searched, but full input preparation
  must not be described as full-proteome discovery.
- GTF/GFF3 relationship parsing and plain/gzip input cases were exercised.
- MAKER-style relationships are supported by the generic parser. A real MAKER
  dataset has not been included in these recorded acceptance runs.
- The same-genome translation run checks execution and sequence handling; it
  is not a cross-species accuracy benchmark.
- The reported candidates lack a curated independent truth set, completed
  manual clearance and experimental validation.

## Recorded documentation/interface revision

**1.1.0** was published from **3.0.11-dev**, which rewrites the documentation and removes the whitelist interface
and tracking helpers. It also expands GeneWise executable/configuration discovery
and validates required configuration files before analysis. Candidate filtering
rules are unchanged. Static syntax and documentation checks are recorded for
this revision.

On 2026-10-13, the 3.0.11-dev checkout passed **98 regression tests**, including
five checks for demo reproducibility and failure handling. The bundled
[synthetic installation demo](../examples/demo/README.md) also completed with
two threads in the existing Conda environment, without explicit `--wisepath` or
`--wiseconfig` options:

- Cross: two selected loci, one final `demo_A.t1` candidate, one recovered stop
  and one frameshift; the intact `demo_B.t1` control was rejected.
- Self: the off-annotation disrupted A copy completed GeneWise, then was
  excluded by the standard full-length homolog filter; zero final candidates.
- TSV/BED/GFF3 IDs and coordinates, expected query IDs, complete cache records,
  report links and the Cross alignment viewer passed the demo checker.
- The combined runner printed `DEMO PASS` in a recorded run taking 25 seconds.
  Runtime depends on tools, CPU and filesystem; this is not a general estimate.

These synthetic checks validate installation and execution. No new full BRAKER
discovery or biological accuracy benchmark was performed for this revision.

## English report revision

**3.0.12-dev** translates the shared HTML report interface into English. Source
syntax and rendered reports from copied 3.0.11-dev demo results were checked.
The previously completed analyses and their cached reports were preserved.
No new regression suite, GeneWise analysis or biological benchmark was run for
this presentation change. Report previews retain the original analysis version.

## Positive demo and viewer revision

On 2026-10-13, **3.0.13-dev** passed **104 regression tests** and a new two-thread
synthetic Cross/Self demo. One recorded combined run took 7 seconds. Cross
produced one final `demo_A.t1` candidate with one verified stop, one frameshift,
100% query coverage and approximately 99.33% amino-acid identity. Self remained
the expected zero-final-candidate control.

The viewer now shares the analysis parser for codon-verified stops and alignment
metrics. In particular, a GeneWise `X` is a stop only when its displayed DNA
codon is TAA, TAG or TGA. The final-only Demo viewer agrees with the TSV counts
and identity and highlights both engineered events. A generated compact
[result snapshot](../examples/demo/expected/README.md) is bundled for viewing.
These are synthetic computational candidates; no biological accuracy value or
new full BRAKER benchmark is claimed.

Version-specific records from the local project are preserved outside the
software checkout under `software_repair_20261012/output_layout_runtime_validation/`.
Those project artifacts are not bundled as portable repository datasets.
Historical release manifests describe their recorded releases and do not
verify this newer development checkout.

## Reporting biological performance

For a reproducible precision/recall result, record:

1. Exact software version, inputs, assembly and annotation releases.
2. Discovery scope and all enabled filters/parameters.
3. An independent truth set and whether it is curated or an annotation-derived
   silver standard.
4. Locus matching rules, strand, reciprocal overlap and treatment of multiple
   loci/queries.
5. Separate candidate counts, precision, recall and runtime/resource measurements.
6. Review of assembly errors, ambiguous sequence and disputed candidates.

Any learned scoring evaluation should separate related gene families and
hold out species. Recovering a literature gene ID is not equivalent to
recovering the same candidate locus or reproducing all manual filtering stages.

## Current public development revision

**1.1.1-dev** incorporates the English report, positive demo and unified-report
changes recorded under **3.0.12-dev** through **3.0.16-dev**. The 104-test and
Cross/Self execution results above apply to **3.0.13-dev**. The later report
revisions were rendered from those existing checked results, with syntax
checks, without a new regression suite or biological analysis. Publication
itself does not establish additional validation.

## Optional cleanup revision

**1.2.0** adds optional post-publication cleanup and a standalone cleanup helper.
On **2026-10-23**, the released runtime code passed **117 regression tests**,
including **13 cleanup-specific tests**, and **18 runtime checks** using new
synthetic analysis outputs and isolated copies of those outputs.

The runtime checks included:

- A fresh Cross/Self installation demo without cleanup: Cross retained one
  final candidate and Self retained zero; intermediate files remained present.
- Successful automatic cleanup for standard Cross, standard Self and Self
  MMseqs-fast workflows, plus standalone cleanup for completed Cross/Self runs.
- A dry-run that changed no file contents or directory entries; repeated cleanup
  was idempotent. Adding the cleanup flag to an otherwise identical complete
  run reused that analysis and then cleaned it.
- Unchanged report, result-table, QC and log fingerprints after standalone
  cleanup; retained report links remained valid. Cross TSV/BED/GFF3 outputs
  matched the uncleaned run, and embedded overview records remained available.
- Preservation of the final Cross GeneWise text with one verified stop and one
  frameshift. Zero-final-candidate and miniprot-only workflows did not require
  a GeneWise archive.
- Refusal of locked/incomplete runs, modified published outputs, missing final
  GeneWise records, source inputs under work/, work-directory symlinks,
  snapshots and unsupported legacy layouts. Symlink targets, original inputs
  and unrecognized work entries were retained.
- A real filesystem permission failure during deletion: cleanup recorded
  `partial`, kept the published outputs, refreshed fingerprints and completed
  successfully after permissions were restored and the command was retried.
- An intentionally failed search command: automatic cleanup did not run, the
  analysis cache remained `running`, and intermediates remained available.
- Refusal to restage or rerun analysis inside a cleaned directory; original
  bundled input hashes remained unchanged. The mount-point guard was exercised
  separately through a mocked unit check, without mounting any filesystem.

The five main cleanup outputs had zero remaining work files. No runtime code
change was needed. Two initial test-harness assumptions were corrected before
all 18 runtime checks passed: fast workflows use miniprot rather than GeneWise,
and GeneWise discovery prepends its executable directory to PATH.

These are small synthetic installation and cleanup checks, using two analysis
threads. They do not establish cleanup performance on a 253 GB intermediate
directory, peak disk usage on large genomes, or biological precision/recall.
Existing scientific analyses and four-node production jobs were not cleaned,
stopped or restarted. Detailed logs and file-fingerprint records are retained
outside the source distribution in the local validation audit.
