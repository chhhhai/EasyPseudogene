# Methods and result interpretation

This describes the current implementation, not a proposed algorithm. Effective
parameters can override the defaults below. See [advanced usage](ADVANCED_USAGE.md)
for additional evidence and execution choices.

## Inputs and coding-model associations

The entry point normalizes plain/gzip FASTA and GTF/GFF3. Protein first-token
IDs are preserved; resolved gene/transcript associations are added to normalized
headers and stored in `metadata/protein_id_mapping.tsv`. Input quality and
exclusions are published under `qc/`.

The reader supports GFF3 ID/Parent links, explicit protein IDs and GTF coding
relationships. Unresolved or ambiguous joins remain visible. Conflicting
relationships, duplicate FASTA IDs and incompatible assembly coordinates fail.
Reference annotations supplied with protein-only Cross inputs are not expected
to share the target genome's contig names.

The optional ID-map TSV uses `query_id` and `transcript_id` and/or `gene_id`.
Use `--protein-id-map FILE` when sequences have been renamed.

For reference-genome translation, gffread extracts translated CDS using strand
and phase. The helper keeps valid proteins of at least 30 amino acids and one
longest representative per gene. This input route requires matching reference
annotation and genome.

## Standard search and local alignment

1. miniprot searches the whole target genome. Defaults are `outn=1000`,
   `outs=0.99`, and one selected hit per query.
2. Hit selection ranks PAF aligned length. Default minimum PAF aligned length
   is 30 and minimum MAPQ is 0. PAF alignment length is not an amino-acid length
   threshold.
3. Candidate intervals are extended by 1,000 bp on each side within contig
   bounds. In Self mode, annotation filtering occurs before DNA extraction.
4. Stable locus IDs connect each candidate to the original query. samtools
   extracts local DNA; GeneWise runs on the appropriate strand in parallel.
5. Cross automatically searches additional copies for disrupted queries when
   copy filtering is enabled and no additional PAF is supplied.

Optional MMseqs diagnostics in the standard workflow do not restrict miniprot
input queries or select its candidates.

## Self annotation filtering

The default mask features are `gene,CDS`. A candidate is excluded when its
fraction overlapping the inferred own locus or annotated mask reaches 0.8.
`--keep-pseudogene-annotations` can retain candidates overlapping existing
pseudogene annotations; its default overlap threshold is 0.1.

Selection of top hits occurs before this filter. With a perfect annotated
self match and only one selected hit, no off-annotation hit may remain.
The README's broader Self settings are explicit search choices, not new defaults.

## Standard post-GeneWise filters

Protein representative selection is performed after GeneWise and mutation
counting. It uses the normalized original proteins, groups by gene identifier
when available, and retains the longest representative; plain IDs are supported.

| Stage | Current default condition |
| --- | --- |
| Disruption candidate | At least one reported stop or frameshift |
| Representative protein | Keep the longest protein per supported gene/header group |
| Family keywords | Exclude configured olfactory/vomeronasal/zinc-finger and related naming patterns |
| Predicted annotations | Exclude configured keywords and biotype `TEC` |
| Other-copy filter: Cross | Exclude loci with another potential copy covering at least 0.8 of the query |
| Other-copy filter: Self | Exclude loci with another copy covering at least 0.8, match fraction at least 0.5 and zero reported stops/frameshifts |
| Independent synteny | Disabled; when enabled, both interval overlap fractions must reach 0.1 |
| Single-event heuristics | Splicing/ending patterns and low local flank identity; combined flank threshold 0.4, window 10 |
| Late disruption | Exclude when the earliest disruption is at or beyond 0.9 of the parsed query alignment sequence |
| Multi-disruption quality | Coverage at least 0.7, identity at least 0.05, aligned query length at least 200 aa |
| Multi-event local support: Cross | Require a qualifying event-local identity of at least 0.4 when event measurements exist; not enabled for Self |
| Global quality | Coverage at least 0.7, identity at least 0.3, aligned query length at least 200 aa |
| Same-locus conflict: Cross | Exclude GeneWise-only disruptions contradicted by a miniprot hit with complete zero stop/fs tags, coverage at least 0.8 and match fraction at least 0.5 |

The basic copy-hit measurements use minimum aligned length 30, query coverage
0.05, match fraction 0.3 and MAPQ 0. Nearby hits can be grouped within 10,000 bp;
the candidate locus itself is excluded from other-copy counts.

The predicted-annotation stage does not automatically remove every single-exon
gene. The main workflow does not forward that helper's optional GTF-based
intronless exclusion route.

GeneWise pretty outputs verify displayed stop codons against TAA/TAG/TGA.
Compact output has a fallback counting path. Some ending/splicing rules use
legacy text patterns; their behavior depends on output format and should be
interpreted as heuristics.

The standard filters are configurable. The historical
`unitary_pseudogene_ids.txt` is a family/copy-filter stage alias and does not
establish unitary gene loss. Candidates with exactly two frameshifts and no
stops are recorded in a review list without automatic exclusion on that basis.
No human review or RELAX analysis is automatically performed.

## Reference neighborhood evidence

When Cross receives reference GTF/GFF3, the pipeline assesses a protein-neighbor
proxy from reference gene order and target alignment anchors after filtering.
This annotation does not alter the final candidate list. Independent synteny
interval filtering is a separate, disabled-by-default stage.

## Optional catalog selection

The MMseqs balanced Self catalog selects longest isoforms, hard-masks target
CDS, performs translated search, clusters hits and applies source/CDS overlap
filters. Its initial defaults include E-value 0.1, bitscore 40, aligned length
30 aa and cluster score 50. Source overlap at least 0.8 is excluded; allowed CDS
overlap is zero. Locus merging uses same strand, gap 2,000 bp, union coordinates
and maximum span 20,000 bp.

Local validation uses the original genome and normalized proteins. It extracts
candidate windows with a dynamic flank, compares the initial and recorded
alternative parents, and selects the strongest eligible miniprot alignment.
The base flank is 500 bp, maximum 5,000 bp; up to ten alternatives are considered.
This is not a genome-wide orthology inference or an event-consensus algorithm.

Unified final catalog candidates require eligible miniprot validation, complete
stop/fs tags, a valid overlapping validation interval, at least one disruption,
and medium/high common evidence. MMseqs fallback support remains unassessed
for these structural grades. Standard GeneWise family/copy filters are not run
in this catalog path.

## Shared evidence grades

| Grade | Meaning |
| --- | --- |
| high | Stop/frameshift present; query coverage at least 0.8 and identity at least 0.5 |
| medium | Stop/frameshift present; coverage at least 0.2 and identity at least 0.3; below high |
| low | Assessed alignment without a disruption, or below the thresholds |
| unassessed | Structural validation or required measurements are missing/incomplete |

GeneWise identity is aligned-query amino-acid identity. miniprot identity is a
PAF match-fraction proxy. `evidence_query_id` records the protein used for
structural validation and can differ from initial `query_id`.

Grades are independent of standard filter selection and are not probabilities.
Internal catalog tier scores have different definitions. Locus and query counts
must be reported separately because one query can support multiple loci.

## Coordinates and output boundaries

TSV/BED coordinates are 0-based, half-open. GFF3 is 1-based, inclusive. Standard
exports use the selected unflanked miniprot interval; eligible catalog exports
use validated alignment boundaries. GFF3 stores candidate interval features.
Detailed alignments, inputs and provenance remain available for review.
