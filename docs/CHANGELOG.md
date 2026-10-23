# Changelog

## 1.2.0 — 2026-10-23

- Add opt-in `--cleanup-intermediates` to clean known work directories only
  after successful publication and cache completion in every workflow.
- Add `easypseudogene-cleanup --run-dir DIR [--dry-run]` for completed runs.
  Check run locks and published fingerprints before removing intermediates.
- Preserve reports, results, QC, metadata, logs and final-candidate GeneWise
  raw text. Protect source inputs and symlink targets; refuse mounted work.
- Record deletion scope and status, refresh artifact fingerprints and prevent
  accidental restaging or in-place continuation of cleaned output archives.
- This revision has source/shell syntax review only. No cleanup command,
  regression suite or biological analysis has been executed for this change.

## 1.1.1-dev — 2026-10-13 (unreleased)

- Merge English statistics, final-candidate alignments, QC and downloads into
  one report.html; use a searchable candidate table with 20 rows per page.
- Restore five GeneWise statistics and seven offline overview charts with
  selectable all-completed-alignments and final-candidate scopes.
- Correct codon-verified stop display and add a portable positive demo snapshot.
- Incorporate development revisions 3.0.12-dev through 3.0.16-dev. Recorded
  validation remains version-specific; no new regression suite was run for
  this publication step.


## 1.1.0 — 2026-10-13

- Publish the maintained development source as the public v1.1.0 release;
  synchronize CLI, citation and current usage documentation version metadata.
- Support Cross/Self comparison modes, standard/fast workflows, generic
  FASTA/GTF/GFF3 inputs, input checks and reproducible cache metadata.
- Publish organized results, QC, metadata, standalone HTML reports and GFF3
  loci with sequence evidence and filtering decisions.
- Bundle a synthetic installation demo, improve GeneWise discovery, reorganize
  documentation and add Linux source checks.
- Remove the whitelist interface; use a fresh output directory when upgrading.
- Historical 3.x entries below record development provenance, rather than
  public GitHub release numbering. See docs/VALIDATION.md for validation scope.


## 3.0.16-dev — 2026-10-13 (unreleased)

- Restore five GeneWise overview statistics and seven charts in the unified
  report: score, identity and length distributions, disruption event totals,
  mutation counts per alignment, and the two identity scatter plots.
- Let users switch between all completed selected GeneWise alignments and
  final candidate alignments. Label alignment counts independently of final
  candidate counts; table pagination/search do not affect chart statistics.
- Embed compact measurements and offline SVG charts without a CDN dependency.
  Exclude missing measurements, report score availability, and disclose scatter
  downsampling above 2,000 points. Preserve the 20-row candidate table.
- Restore the overview using existing checked demo measurements. No new
  alignment analysis or regression suite was run for this presentation change.

## 3.0.15-dev — 2026-10-13 (unreleased)

- Replace the final-candidate card list with a compact table of IDs, location,
  evidence grade, stops, frameshifts, query coverage and identity.
- Show 20 candidates per page. Search all candidates and reset to page one
  when the search changes; show matching counts and an explicit empty state.
- Expand alignment evidence below the selected row, with all details initially
  collapsed. Keep full TSV/BED/GFF3 exports independent of pagination.
- Regenerate the example report from existing checked demo data. Generate a
  separately labelled 100-row layout illustration, not an analysis result.
  No new alignment analysis or regression suite was run for this UI revision.

## 3.0.14-dev — 2026-10-13 (unreleased)

- Merge statistics, searchable final-candidate details, highlighted GeneWise
  alignments, filtering diagnostics and downloads into a single offline report.
- Match alignments by final locus ID so rejected copies of the same protein
  cannot enter the final-candidate display. Show missing evidence explicitly.
- Stop generating a second viewer in the standard workflow and bundled demo;
  retain --visualize as a compatibility option and the standalone viewer tool.
- Regenerate the portable preview from previously checked 3.0.13-dev demo
  results. No new alignment analysis or regression suite was run for this
  presentation revision.

## 3.0.13-dev — 2026-10-13 (unreleased)

- Add a final-only positive Demo viewer and a generated, portable result
  snapshot with one detected synthetic candidate and its TSV/BED/GFF3 outputs.
- Use the codon-aware analysis parser in the GeneWise viewer so verified X
  stops, frameshifts and identity agree with the result tables. Highlight
  verified stops and retain ambiguous codons as unassessed residues.
- Label generic viewer records as alignments, escape user-supplied text and
  keep tables/alignment details usable when the Chart.js CDN is unavailable.
- Pass 104 regression tests and a fresh synthetic Cross/Self demo using two
  threads. Candidate filtering rules remain unchanged.

## 3.0.12-dev — 2026-10-13 (unreleased)

- Make the shared report.html interface English for both comparison modes and
  standard/catalog workflows: headings, counters, evidence descriptions,
  candidate preview, download labels, QC notes and empty-table messages.
- Keep candidate calculations, filtering thresholds, result schemas and
  original input IDs unchanged. Existing completed reports retain their
  original language unless separately regenerated.
- Record syntax checks and report rendering from copied 3.0.11-dev demo
  results. No new regression suite or alignment analysis was run for this
  presentation change.


## 3.0.11-dev — 2026-10-13 (unreleased)

- Prepare a compact GitHub source layout: detailed documentation under docs/,
  a minimal standard environment, current citation metadata and a source-check
  workflow. Exclude local backups, previous result archives, old release
  records, container recipes and internal benchmark/release helper scripts.
  Runtime code and test contents are preserved from the maintained checkout.
- Bundle a reproducible synthetic installation demo for Cross/Self, with a
  two-thread runner, positive/negative controls and an automated result checker.
- Restore the human/dolphin Cross example with decompression, explicit mode and
  current result paths; retain the manual GeneWise installation overrides.
- Share GeneWise executable/configuration discovery between the CLI and batch
  runner. Add active-Conda binary fallback, symlink and source-prefix handling,
  broader configuration locations, early file checks and specific errors.
  Input-only checks and workflows without GeneWise remain independent of it.
- Rewrite the README, quick start and installation guide in English around
  Cross/Self usage, generic inputs and the current output layout. Move method,
  advanced-workflow and validation details into dedicated documentation.
- Remove the whitelist CLI option, environment forwarding, cache input and
  tracked-ID report. Delete the unused whitelist prefilter and tracking helper.
  Existing candidate filtering rules remain the same.
- Record the frozen 3.0.10-dev validation scope separately from this revision.
  The current checkout passes 98 regression tests and the bundled synthetic
  Cross/Self demo; automatic GeneWise discovery works with the existing tested
  Conda installation. No new full BRAKER discovery run is claimed.

## 3.0.10-dev — 2026-10-12 (unreleased)

- Place CLI intermediates and staged inputs under `work/`; keep published
  `results/`, `qc/`, `metadata/` and `logs/` at the run root. Create only the
  stage directories actually used by the selected workflow.
- Write a standalone `report.html` for cross/self workflows, with input QC,
  stage counts, evidence grades, final candidate previews and file links.
  Reports describe standard and fast selection rules separately.
- Add `results/final_loci.gff3`, using 1-based inclusive coordinates and
  escaped IDs/attributes; TSV/BED coordinates remain 0-based half-open.
- Move CLI cache metadata to `metadata/reference_cache.json`, protect the
  summary report and updated artifact paths, and reject old layouts without
  migrating their existing results. Use fresh output directories when upgrading.
- Keep explicit legacy publishers/stage helpers able to read their old paths.
  Optional GeneWise alignment pages now live under `results/`.
- Pass 93 regression tests and 20 runtime checks on real BRAKER data. Verify
  actual directory inventories, offline report links, GFF3/TSV/BED equivalence,
  complete cache reuse and rejection of modified report/GFF3 artifacts.
  Cross and self result tables match the previous validated baselines; scope
  remains 200 cross reference genes and 40 self proteins, not a full-proteome run.

## 3.0.9-dev — 2026-10-12 (unreleased)

- Keep `cross`/`self` as the comparison modes; expose a short basic help and
  `--workflow standard|fast`, with engine/preset controls in `--help-advanced`.
  The simple self fast workflow selects MMseqs balanced with disruption validation.
- Default to at most 16 CPU threads, bounded by affinity and declared scheduler
  allocation; preserve explicit thread budgets.
- Add `--check-inputs` in a separate `input_check/` directory so preflight does
  not invalidate a later run's result cache.
- Report and skip invalid reference protein sequences by default in the CLI;
  `--invalid-proteins error` retains strict rejection. Duplicate IDs, bad
  annotation relationships and incompatible genomes still fail.
- Stage genomes in the run directory so indexes do not need to be written next
  to source FASTA files. Publish input QC and protein exclusions for both workflows.
- Publish every fast engine through `results/`, `qc/`, and `metadata/`, preserving
  legacy catalogs. Unvalidated hits remain screening candidates; fast final
  candidates require validated disruptions and medium/high sequence evidence.
- Add common evidence fields and explicit high/medium/low/unassessed rules to
  both result publishers. Evidence grades are separate from filtering outcomes
  and are not calibrated probabilities of true pseudogene status.
- Require complete miniprot structural tags for common evidence grades; handle
  empty fast validation inputs without invoking miniprot on empty FASTA files.
- Retain existing GeneWise filters and keep RELAX outside the main workflows.
- Fix self-mode filtering with multiple hits per protein: retain each locus's
  matching flanked interval, and reject mismatched/noncontaining intervals.
- Reuse verified complete runs without rewriting results or metadata; preserve
  refusal of changed inputs/configuration and incomplete or modified results.
- Validate real BRAKER input preflight, cross/self standard and self fast runs,
  GTF/GFF3 and gzip inputs, zero-candidate publication, reference-genome
  translation and cache guards. Cross checks use 200 T. rubra reference genes
  against the whole T. dohrnii target; self checks use 40 original proteins.
  These are execution checks, not a curated biological accuracy benchmark.

## 3.0.8-dev — 2026-10-12 (unreleased)

- Add a shared GTF/GFF3 reader for database annotations and local BRAKER/MAKER
  coding models, including exon/CDS-only GTF, GFF3 ID/Parent relationships,
  shared child features, URL-escaped attributes and embedded FASTA sections.
- Accept `--reference-annotation` (legacy `--reference-gtf` remains an alias),
  `--target-annotation` (alias of `--gff`) and an optional `--protein-id-map` TSV.
- Stage plain/gzip FASTA and annotation inputs, preserve original protein IDs
  and headers in a mapping report, and remove only terminal protein stop markers.
- Join proteins by explicit protein IDs or transcript IDs, covering BRAKER
  `.aa` and MAKER protein exports without requiring Ensembl identifiers.
  Ambiguous or unmatched joins remain visible; custom joins can be supplied.
- Infer missing gene/transcript containers from coding relationships; use
  the same model for reference translation and neighborhood evidence.
- Record input checksums, normalization source and ID associations; reject
  conflicting relationships, duplicate FASTA IDs and incompatible contig bounds.
- Keep scientific filters and the absence of RELAX in the main workflow.
- Fix publication of valid empty self-mode runs: when annotation masking
  removes every hit, create the zero-locus provenance files before filtering.

## 3.0.7-dev — 2026-09-30 (unreleased)

- Match versioned Ensembl GTF gene, transcript, and protein identifiers to
  stable IDs in reference FASTA headers and neighborhood evidence.
- Fix `cross` runs without a reference GTF failing in predicted-annotation
  filtering because a conditional local `re` import shadowed the module import.
- In `cross` mode, automatically run a broad miniprot copy search only for
  GeneWise-disrupted reference queries; accept `--reference-gtf` to generate
  matched annotation and a per-locus flanking-gene-order report.
- Accept an annotated `--reference-genome` with `--reference-gtf` and use
  gffread to translate CDS and select one representative protein per gene.
- Publish the neighborhood proxy status alongside each locus without making
  unresolved neighborhood evidence an automatic exclusion.
- Use a conservative >=80% other-copy policy in `cross` mode, leaving the
  previous intact-copy policy configurable. Exclude strong selected-locus
  miniprot alignments that contradict a GeneWise-only disruption and record
  every decision. Check local evidence around every multi-disruption call.
- Report a copy exclusion as `other_full_length_homolog_filter` because a
  miniprot hit alone does not establish that the other copy is functional.

- Accept a matching reference GTF-derived annotation table so large-family
  filtering can recognize gene symbols even when protein FASTA headers contain
  only Ensembl IDs. Activate the configured predicted-annotation filter.
- Assess additional full-length genome copies relative to each candidate locus;
  a single separate copy now counts even when the disrupted locus is partial.
  The `intact` policy requires a sufficiently similar hit without reported
  stops or frameshifts; the `cross` default also screens potential copies.
  Record other full-length hits separately. Accept a broader search PAF
  through `--redundancy-paf`.
- Approximate the paper's 30-bp coding window on each side of a single
  disruption with 10 aligned amino acids on each side and a 40% combined
  identity threshold. Keep the former `and` and `or` rules as explicit options.
- Wire the optional conserved-synteny gate to an independent interval TSV,
  and fail if the gate is enabled without evidence. Publish an evidence-status
  table recording missing annotation, synteny, and relaxed-selection evidence.
  Final cross results remain computational candidates.

## 3.0.6-dev — 2026-09-29 (unreleased)

- Expose two comparison modes: `cross` (default) and `self`.
- Keep `reference` as a legacy alias for `cross` and retain old `--mode unitary`
  commands for compatibility.
- Present unitary evidence annotation as an optional cross-species follow-up
  through `bin/easypseudogene-unitary`, without changing its algorithm.
- Record comparison mode in organized GeneWise run metadata.
- Name the main-pipeline family/copy-filter stage `family_copy_pass`; keep
  `unitary_pseudogene_ids.txt` as a compatibility alias.

## 3.0.5-dev — 2026-09-28 (unreleased)

- Keep GeneWise and the existing stage paths. Publish reader-facing
  `results/`, `qc/`, and `metadata/` directories with per-locus decisions,
  per-query counts, stage totals, effective thresholds, and a run manifest.
- Clarify that the legacy `unitary_pseudogene_ids.txt` is an operational
  filtering stage, not proof of unitary gene loss.
- Correct the README output tree and whitelist description.

## 3.0.4-dev — 2026-09-28 (unreleased)

- Compute coverage, identity, mutation position and flanking identity from
  query amino acids in GeneWise alignments. Intron spaces and gap columns no
  longer inflate alignment quality or yield coverage above 100%.

## 3.0.3-dev — 2026-09-28 (unreleased)

- Parse GeneWise alignments when a long genomic label touches its coordinate;
  previously these valid alignment blocks were omitted from mutation counts
  and quality filtering. Add a fixed-column regression case.

## 3.0.2-dev — 2026-09-26 (unreleased)

- Correct codon-backed stop parsing, including reverse-strand pretty alignments.
- Remove whitelist injection; export separate tracked-ID status.
- Preserve custom filter configuration and unitary genome-scan arguments.
- Keep independent GeneWise loci and propagate failed child processes.
- Guard reference caches with content fingerprints and completion checks.
- Accept plain FASTA and prefer stable gene IDs for isoform grouping.
- Enforce tblastn CPU budgets and preserve explicit resource options.
- Unify cross-species OOF family folds and correct P-GRe coordinate conversion.
- Skip unused reference MMseqs diagnostics by default.
- Add audit regression tests; old releases and benchmark artifacts remain frozen.

## 3.0.1 — 2026-07-15

### Added
- Frozen TAIR10 and rice RepeatMasker class evidence and provenance.
- MCScanX TAIR10–A. lyrata gene-order synteny from 9,619 cross-species anchors.
- Real-species unitary audit with TAIR10 target, A. lyrata reference, and Capsella outgroup.
- Held-out PR curves, calibration, chromosome-bootstrap forest plot, and evidence attrition.
- Blinded dual-review computational audit workflow and Cohen-kappa summarizer.
- GitHub Actions CI, CITATION CFF, contribution/conduct documents, Conda, Docker, Apptainer, and Zenodo metadata.

### Changed
- Distinguish best homology hits from intact functional copies in unitary output fields.
- Release manifest now includes CI, citation, environment, container, and release metadata.
- Whole-repository Ruff F/E9 and Bandit medium/high checks are CI gates.

### Scientific scope
- Real-species unitary and high-confidence duplicated-like outputs remain computational candidates, not experimental or literature confirmation.

## 3.0.0 — 2026-07-15

### Added
- Explicit sensitive, miniprot-supported, strict-core-disruption, and subtype catalog tiers.
- `tblastn-production-candidates-only` and `tblastn-production-validated` presets.
- Auditable one-to-one prediction–truth matched-pair tables.
- Plotting guards that reject metrics without one-to-one maximum-cardinality matching.
- Stratified manual-review tables with configurable IGV links.
- RepeatMasker class-specific overlap evidence and low-complexity penalties.
- Processed-like evidence scores from intron loss, poly(A), TSD, repeat context, and disruptions.
- Duplicated-like scores from exon conservation, distance, disruption, and optional synteny.
- Conservative `--mode unitary` with target-proteome, whole-target-genome, synteny, disruption, and outgroup gates.
- Experimental DAG HSP chaining; proximity clustering remains the validated default.
- Multi-species manifests, MMseqs2 cross-species orthogroups, orthogroup-safe folds, and leave-one-species-out evaluation with training-only out-of-fold threshold selection, PR-AUC, calibration, and auditable held-out scores.
- Resumable chromosome-parallel tblastn with fingerprints, checksums, retries, disk preflight, and signal cleanup.

### Corrected
- Final benchmark matching is strictly one-to-one; legacy many-to-many values are retained only for provenance.
- Serial BLAST-fast now honors CDS hard masking.
- BLAST-fast now honors disruption validation.
- NCBI `protein_id` mappings are recognized.

### Validation
- Full TAIR10, Oryza sativa IRGSP-1.0, and Arabidopsis lyrata GCF_000004255.2 runs.
- TAIR10 and rice DAG regressions showed no advantage; DAG remains experimental.
- Resume and corrupted-output recovery were exercised on a real tblastn smoke dataset.
