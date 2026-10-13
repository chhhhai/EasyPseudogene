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

## Current documentation/interface revision

**1.1.0** (prepared from the **3.0.11-dev** development revision) rewrites the documentation and removes the whitelist interface
and tracking helpers. It also expands GeneWise executable/configuration discovery
and validates required configuration files before analysis. Candidate filtering
rules are unchanged. Static syntax and documentation checks are recorded for
this revision.

On 2026-10-13, the current checkout passed **98 regression tests**, including
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

The GitHub upload copy preserves the runtime code, tests and bundled demo data
from that checked checkout. Packaging checks compare file hashes, parse source
syntax and resolve documentation links. Tests and the demo were not rerun from
the upload copy, and its minimal environment specification was not freshly
installed as part of packaging.

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
