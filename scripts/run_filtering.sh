#!/usr/bin/env bash
set -euo pipefail

# Post-GeneWise filtering pipeline

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EASYPSEUDOGENE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${EASYPSEUDOGENE_DIR}/lib"

# Parse arguments
PROTEINS=""
COMPARISON_MODE="self"
SPECIES=""
OUTPUT_DIR=""
CLI_REFERENCE_ANNOTATIONS=""
CLI_REDUNDANCY_PAF=""
CLI_COPY_SEARCH_SOURCE=""
CLI_SYNTENY_TSV=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG_FILE="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --mode) COMPARISON_MODE="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --reference-annotations) CLI_REFERENCE_ANNOTATIONS="$2"; shift 2 ;;
        --redundancy-paf) CLI_REDUNDANCY_PAF="$2"; shift 2 ;;
        --copy-search-source) CLI_COPY_SEARCH_SOURCE="$2"; shift 2 ;;
        --synteny-tsv) CLI_SYNTENY_TSV="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done
[[ "${COMPARISON_MODE}" == "cross" || "${COMPARISON_MODE}" == "self" ]] || { echo "ERROR: --mode must be cross or self" >&2; exit 1; }

CONFIG_FILE="${CONFIG_FILE:-${EASYPSEUDOGENE_DIR}/config/default.conf}"
[[ -f "${CONFIG_FILE}" ]] || { echo "ERROR: config not found: ${CONFIG_FILE}" >&2; exit 1; }
source "${EASYPSEUDOGENE_DIR}/config/default.conf"
source "${CONFIG_FILE}"
REFERENCE_ANNOTATION_TSV="${CLI_REFERENCE_ANNOTATIONS:-${REFERENCE_ANNOTATION_TSV:-}}"
SYNTENY_TSV="${CLI_SYNTENY_TSV:-${SYNTENY_TSV:-}}"
sp_dir="${EPG_WORK_DIR:-${OUTPUT_DIR}/${SPECIES}}"
gw_dir="${sp_dir}/genewise"
miniprot_dir="${sp_dir}/miniprot"
filter_dir="${sp_dir}/post_genewise_filter"
prefilter_dir="${sp_dir}/prefilter"
ORIGINAL_PROTEINS="${PROTEINS}"
FILTER_PAF="${miniprot_dir}/${SPECIES}.paf"
if [[ -f "${sp_dir}/regions/locus_map.tsv" ]]; then
    PROTEINS="${sp_dir}/regions/locus_proteins.fa"
fi
[[ -z "${CLI_REDUNDANCY_PAF}" ]] || FILTER_PAF="${CLI_REDUNDANCY_PAF}"
[[ -z "${REFERENCE_ANNOTATION_TSV}" || -f "${REFERENCE_ANNOTATION_TSV}" ]] || { echo "ERROR: reference annotation TSV not found: ${REFERENCE_ANNOTATION_TSV}" >&2; exit 1; }
[[ -f "${FILTER_PAF}" ]] || { echo "ERROR: redundancy PAF not found: ${FILTER_PAF}" >&2; exit 1; }
[[ -z "${SYNTENY_TSV}" || -f "${SYNTENY_TSV}" ]] || { echo "ERROR: synteny TSV not found: ${SYNTENY_TSV}" >&2; exit 1; }
if [[ "${ENABLE_SYNTENY_FILTER:-0}" == "1" ]]; then
    [[ -n "${SYNTENY_TSV}" ]] || { echo "ERROR: ENABLE_SYNTENY_FILTER=1 requires --synteny-tsv or SYNTENY_TSV" >&2; exit 1; }
    [[ -f "${sp_dir}/regions/locus_map.tsv" ]] || { echo "ERROR: synteny filter requires locus_map.tsv" >&2; exit 1; }
fi

mkdir -p "${filter_dir}" "${prefilter_dir}"
# Recreate the alignment on every filtering pass; thresholds never reuse stale derived files.
rm -f "${filter_dir}/protein-alignment.txt"
(set -o posix; set) | awk '/^(ENABLE_|GLOBAL_|MINIPROT_|REDUNDANCY_|MULTI_DISRUPTION_|LOW_ID_|END_MUTATION_|FAMILY_PATTERNS|PREDICTED_|REFERENCE_ANNOTATION_|SYNTENY_|CROSS_)/' > "${filter_dir}/effective_config.txt"

# Combine GeneWise outputs if needed
if [[ ! -s "${gw_dir}/pseudogene-genewise.txt" ]]; then
    find "${gw_dir}" -maxdepth 1 -name "*.gw.txt" -type f -size +0c | sort | xargs -r cat \
        > "${gw_dir}/pseudogene-genewise.txt"
fi
cp -f "${gw_dir}/pseudogene-genewise.txt" "${filter_dir}/pseudogene-genewise.txt"

# Count mutations
echo "  Counting mutations..."
python3 "${LIB_DIR}/count_genewise_mutations.py" "${filter_dir}/pseudogene-genewise.txt" \
    > "${filter_dir}/mutation_counts.tsv"

awk '($2+$3)>=1 {print $1}' "${filter_dir}/mutation_counts.tsv" | sort -u \
    > "${filter_dir}/putative_pseudogene_ids.raw.txt"

# Prefilter
VALID_TRANSCRIPTS="${prefilter_dir}/valid_transcript_ids.txt"
if [[ "${ENABLE_PREFILTER:-1}" == "1" ]]; then
    echo "  Prefiltering proteins..."
    prefilter_cmd=(python3 "${LIB_DIR}/prefilter_proteins.py")
    "${prefilter_cmd[@]}" \
        --input "${ORIGINAL_PROTEINS}" \
        --output "${prefilter_dir}/proteins.filtered.fa" \
        --report "${prefilter_dir}/prefilter_report.tsv"
    awk '/^>/ {sub(/^>/, "", $1); print $1}' "${prefilter_dir}/proteins.filtered.fa" | sort -u \
        > "${VALID_TRANSCRIPTS}"
else
    awk '/^>/ {sub(/^>/, "", $1); print $1}' "${ORIGINAL_PROTEINS}" | sort -u \
        > "${VALID_TRANSCRIPTS}"
fi

if [[ -f "${sp_dir}/regions/locus_map.tsv" ]]; then
    awk 'NR==FNR {valid[$1]=1; next} FNR>1 && ($2 in valid) {print $1}' "${VALID_TRANSCRIPTS}" "${sp_dir}/regions/locus_map.tsv" | sort -u > "${VALID_TRANSCRIPTS}.loci"
    VALID_TRANSCRIPTS="${VALID_TRANSCRIPTS}.loci"
fi
comm -12 <(sort "${filter_dir}/putative_pseudogene_ids.raw.txt") <(sort "${VALID_TRANSCRIPTS}") \
    > "${filter_dir}/putative_pseudogene_ids.txt"

unitary_list="${filter_dir}/putative_pseudogene_ids.txt"
printf 'evidence\tstatus\tsource\n' > "${filter_dir}/evidence_status.tsv"
if [[ -n "${REFERENCE_ANNOTATION_TSV}" ]]; then
    printf 'reference_annotation\tprovided\t%s\n' "${REFERENCE_ANNOTATION_TSV}" >> "${filter_dir}/evidence_status.tsv"
else
    printf 'reference_annotation\theader_only\t%s\n' "${ORIGINAL_PROTEINS}" >> "${filter_dir}/evidence_status.tsv"
    echo "  WARNING: no reference annotation TSV; family/predicted filtering uses FASTA headers only" >&2
fi
if [[ "${ENABLE_REDUNDANCY_FILTER:-1}" == "1" ]]; then
    printf 'functional_copy_search\t%s\t%s\n' "${CLI_COPY_SEARCH_SOURCE:-$([[ -n "${CLI_REDUNDANCY_PAF}" ]] && echo supplied || echo screening_paf)}" "${FILTER_PAF}" >> "${filter_dir}/evidence_status.tsv"
else
    printf 'functional_copy_search\tdisabled\t\n' >> "${filter_dir}/evidence_status.tsv"
fi
if [[ "${ENABLE_SYNTENY_FILTER:-0}" == "1" ]]; then
    printf 'conserved_synteny\tprovided\t%s\n' "${SYNTENY_TSV}" >> "${filter_dir}/evidence_status.tsv"
else
    printf 'conserved_synteny\tnot_assessed\t\n' >> "${filter_dir}/evidence_status.tsv"
fi
printf 'relaxed_selection\tnot_assessed\t\n' >> "${filter_dir}/evidence_status.tsv"

# Unitary filters
if [[ "${ENABLE_FAMILY_FILTER:-1}" == "1" ]]; then
    family_cmd=(python3 "${LIB_DIR}/filter_large_family.py" \
        --proteins "${PROTEINS}" \
        --out-list "${filter_dir}/large_family_ids.txt" \
        --out-report "${filter_dir}/large_family_report.tsv" \
        --patterns "${FAMILY_PATTERNS:-olfactory receptor|vomeronasal|zinc finger|^Olfr[0-9]|^Vmn[12]r|^Zfp[0-9]|^OR[0-9]|^V1R|^V2R|^TAAR|^ZNF[0-9]}")
    [[ -z "${REFERENCE_ANNOTATION_TSV}" ]] || family_cmd+=(--annotations "${REFERENCE_ANNOTATION_TSV}")
    [[ ! -f "${sp_dir}/regions/locus_map.tsv" ]] || family_cmd+=(--locus-map "${sp_dir}/regions/locus_map.tsv")
    "${family_cmd[@]}"
    comm -23 <(sort "${unitary_list}") <(sort "${filter_dir}/large_family_ids.txt") \
        > "${filter_dir}/unitary_step1.no_family.txt"
    unitary_list="${filter_dir}/unitary_step1.no_family.txt"
fi

if [[ "${ENABLE_PREDICTED_FILTER:-1}" == "1" ]]; then
    predicted_cmd=(python3 "${LIB_DIR}/filter_predicted_or_intronless.py"
        --proteins "${PROTEINS}"
        --out-list "${filter_dir}/predicted_annotation_ids.txt"
        --out-report "${filter_dir}/predicted_annotation_report.tsv"
        --patterns "${PREDICTED_PATTERNS:-predicted|intronless cDNA|expressed sequence}"
        --exclude-biotypes "${PREDICTED_EXCLUDE_BIOTYPES:-}")
    [[ -z "${REFERENCE_ANNOTATION_TSV}" ]] || predicted_cmd+=(--annotations "${REFERENCE_ANNOTATION_TSV}")
    [[ ! -f "${sp_dir}/regions/locus_map.tsv" ]] || predicted_cmd+=(--locus-map "${sp_dir}/regions/locus_map.tsv")
    "${predicted_cmd[@]}"
    comm -23 <(sort "${unitary_list}") <(sort "${filter_dir}/predicted_annotation_ids.txt") \
        > "${filter_dir}/unitary_step_predicted_annotation.txt"
    unitary_list="${filter_dir}/unitary_step_predicted_annotation.txt"
fi

if [[ "${ENABLE_REDUNDANCY_FILTER:-1}" == "1" ]]; then
    copy_policy="intact"
    [[ "${COMPARISON_MODE}" != "cross" ]] || copy_policy="${REDUNDANCY_COPY_POLICY:-potential}"
    redundancy_cmd=(python3 "${LIB_DIR}/filter_functional_redundancy_miniprot.py" \
        --paf "${FILTER_PAF}" \
        --out-prefix "${filter_dir}/miniprot_redundancy" \
        --min-aln-len "${MINIPROT_MIN_ALN_LEN:-30}" \
        --min-cov "${MINIPROT_MIN_COV:-0.05}" \
        --min-id "${MINIPROT_MIN_ID:-0.3}" \
        --min-mapq "${MINIPROT_MIN_MAPQ:-0}" \
        --copy-cov "${REDUNDANCY_COPY_COV:-0.8}" \
        --copy-min-id "${REDUNDANCY_COPY_MIN_ID:-0.5}" \
        --copy-policy "${copy_policy}" \
        --target-gap "${REDUNDANCY_TARGET_GAP:-10000}")
    [[ ! -f "${sp_dir}/regions/locus_map.tsv" ]] || redundancy_cmd+=(--locus-map "${sp_dir}/regions/locus_map.tsv")
    "${redundancy_cmd[@]}"
    comm -23 <(sort "${unitary_list}") <(sort "${filter_dir}/miniprot_redundancy_redundant.txt") \
        > "${filter_dir}/unitary_step2.no_redundancy.txt"
    unitary_list="${filter_dir}/unitary_step2.no_redundancy.txt"
fi

if [[ "${ENABLE_SYNTENY_FILTER:-0}" == "1" ]]; then
    python3 "${LIB_DIR}/filter_conserved_synteny.py" \
        --synteny-tsv "${SYNTENY_TSV}" \
        --locus-map "${sp_dir}/regions/locus_map.tsv" \
        --in-list "${unitary_list}" \
        --out-pass "${filter_dir}/synteny_pass_ids.txt" \
        --out-fail "${filter_dir}/conserved_synteny_failed.txt" \
        --out-report "${filter_dir}/conserved_synteny_report.tsv" \
        --min-overlap "${SYNTENY_MIN_OVERLAP:-0.1}"
    unitary_list="${filter_dir}/synteny_pass_ids.txt"
fi

cp -f "${unitary_list}" "${filter_dir}/unitary_candidates.txt"

# False-positive filters
if [[ "${ENABLE_FALSE_POSITIVE_FILTER:-1}" == "1" ]]; then
    awk '($2+$3)==1 {print $1}' "${filter_dir}/mutation_counts.tsv" | sort -u \
        > "${filter_dir}/one_mutation_all.txt"
    comm -12 <(sort "${filter_dir}/one_mutation_all.txt") <(sort "${filter_dir}/unitary_candidates.txt") \
        > "${filter_dir}/one_mutation_list.txt"
    awk '($2+$3)>=2 {print $1}' "${filter_dir}/mutation_counts.tsv" | sort -u \
        > "${filter_dir}/multi_mutation_all.txt"
    comm -12 <(sort "${filter_dir}/multi_mutation_all.txt") <(sort "${filter_dir}/unitary_candidates.txt") \
        > "${filter_dir}/multi_mutation_for_manual_review.txt"

    (
        cd "${filter_dir}"
        python3 "${LIB_DIR}/get_ending_splicing_one_stop.py" \
            --one-mutation-list "one_mutation_list.txt" \
            --genewise-file "pseudogene-genewise.txt" \
            --out-splicing "splicing-stop.txt" \
            --out-ending "ending-stop.txt" >/dev/null
        python3 "${LIB_DIR}/get_ending_splicing_one_shift.py" \
            --one-mutation-list "one_mutation_list.txt" \
            --genewise-file "pseudogene-genewise.txt" \
            --out-splicing "splicing-shift.txt" \
            --out-ending "ending-shift.txt" >/dev/null
        python3 "${LIB_DIR}/extract_protein_alignment.py" \
            --input "pseudogene-genewise.txt" \
            --output "protein-alignment.txt"
        python3 "${LIB_DIR}/filter_low_identity_flanks.py" \
            --protein-alignment "protein-alignment.txt" \
            --one-mutation-list "one_mutation_list.txt" \
            --out-stop "low-identity-stop.txt" \
            --out-shift "low-identity-shift.txt" \
            --flank-size "${LOW_ID_FLANK_SIZE:-10}" \
            --flank-th "${LOW_ID_FLANK_TH:-0.4}" \
            --overall-th "${LOW_ID_OVERALL_TH:-0.1}" \
            --mode "${LOW_ID_MODE:-combined}"
        python3 "${LIB_DIR}/filter_end_mutation_90pct.py" \
            --protein-alignment "protein-alignment.txt" \
            --out-list "end_mutation_90pct.txt" \
            --out-report "end_mutation_90pct_report.tsv" \
            --min-frac "${END_MUTATION_FRAC:-0.9}"
    )

    cat "${filter_dir}/splicing-stop.txt" \
        "${filter_dir}/ending-stop.txt" \
        "${filter_dir}/low-identity-stop.txt" \
        "${filter_dir}/splicing-shift.txt" \
        "${filter_dir}/ending-shift.txt" \
        "${filter_dir}/low-identity-shift.txt" \
        "${filter_dir}/end_mutation_90pct.txt" 2>/dev/null | sort -u \
        > "${filter_dir}/false_positive_ids.txt"

    comm -23 <(sort "${filter_dir}/unitary_candidates.txt") <(sort "${filter_dir}/false_positive_ids.txt") \
        > "${filter_dir}/high_confidence_pseudogenes.txt"
else
    cp -f "${filter_dir}/unitary_candidates.txt" "${filter_dir}/high_confidence_pseudogenes.txt"
fi

final_list="${filter_dir}/high_confidence_pseudogenes.txt"

# Always define the multi-mutation list, including when the FP filter is disabled.
awk '($2+$3)>=2 {print $1}' "${filter_dir}/mutation_counts.tsv" | sort -u > "${filter_dir}/multi_mutation_all.txt"
comm -12 <(sort "${filter_dir}/multi_mutation_all.txt") <(sort "${filter_dir}/unitary_candidates.txt") > "${filter_dir}/multi_mutation_for_manual_review.txt"

# Multi-disruption quality filter
if [[ "${ENABLE_MULTI_DISRUPTION_QUALITY_FILTER:-1}" == "1" ]]; then
    if [[ ! -s "${filter_dir}/protein-alignment.txt" ]]; then
        ( cd "${filter_dir}" && python3 "${LIB_DIR}/extract_protein_alignment.py" \
            --input "pseudogene-genewise.txt" \
            --output "protein-alignment.txt" )
    fi
    python3 "${LIB_DIR}/filter_multi_disruption_quality.py" \
        --protein-alignment "${filter_dir}/protein-alignment.txt" \
        --proteins "${PROTEINS}" \
        --mutation-counts "${filter_dir}/mutation_counts.tsv" \
        --in-list "${filter_dir}/multi_mutation_for_manual_review.txt" \
        --min-cov "${MULTI_DISRUPTION_MIN_COV:-0.7}" \
        --min-id "${MULTI_DISRUPTION_MIN_ID:-0.05}" \
        --min-aln-len "${MULTI_DISRUPTION_MIN_ALN_LEN:-200}" \
        --min-event-local-id "$([[ "${COMPARISON_MODE}" == "cross" ]] && echo "${MULTI_DISRUPTION_MIN_EVENT_LOCAL_ID:-0.4}" || echo 0)" \
        --event-flank "${MULTI_DISRUPTION_EVENT_FLANK:-10}" \
        --out-list "${filter_dir}/multi_disruption_low_quality.txt" \
        --out-report "${filter_dir}/multi_disruption_quality_report.tsv"
    comm -23 <(sort "${final_list}") <(sort "${filter_dir}/multi_disruption_low_quality.txt") \
        > "${filter_dir}/high_confidence_pseudogenes_tmp_1.txt"
    final_list="${filter_dir}/high_confidence_pseudogenes_tmp_1.txt"
fi

# Global alignment quality filter
final_list_before_strict="${final_list}"
if [[ "${ENABLE_GLOBAL_ALIGNMENT_FILTER:-1}" == "1" ]]; then
    if [[ ! -s "${filter_dir}/protein-alignment.txt" ]]; then
        ( cd "${filter_dir}" && python3 "${LIB_DIR}/extract_protein_alignment.py" \
            --input "pseudogene-genewise.txt" \
            --output "protein-alignment.txt" )
    fi
    python3 "${LIB_DIR}/filter_alignment_quality.py" \
        --protein-alignment "${filter_dir}/protein-alignment.txt" \
        --proteins "${PROTEINS}" \
        --in-list "${final_list}" \
        --min-cov "${GLOBAL_MIN_COV:-0.7}" \
        --min-id "${GLOBAL_MIN_ID:-0.3}" \
        --min-aln-len "${GLOBAL_MIN_ALN_LEN:-200}" \
        --out-pass "${filter_dir}/high_confidence_pseudogenes_tmp_2.txt" \
        --out-fail "${filter_dir}/high_confidence_pseudogenes_removed.txt" \
        --out-report "${filter_dir}/high_confidence_pseudogenes_report.tsv"
    final_list="${filter_dir}/high_confidence_pseudogenes_tmp_2.txt"
fi

# A high-coverage miniprot alignment to the selected locus that reports no
# disruption contradicts a GeneWise-only mutation call. Other disagreements
# remain in the candidate set for review.
if [[ "${COMPARISON_MODE}" == "cross" && "${ENABLE_CROSS_INTACT_ALIGNMENT_FILTER:-1}" == "1" ]]; then
    python3 "${LIB_DIR}/filter_cross_intact_alignment.py" \
        --screening-paf "${miniprot_dir}/${SPECIES}.paf" \
        --locus-map "${sp_dir}/regions/locus_map.tsv" \
        --mutation-counts "${filter_dir}/mutation_counts.tsv" \
        --in-list "${final_list}" \
        --out-exclude "${filter_dir}/miniprot_intact_conflict.txt" \
        --out-report "${filter_dir}/miniprot_intact_conflict_report.tsv" \
        --min-cov "${CROSS_INTACT_MIN_COV:-0.8}" \
        --min-id "${CROSS_INTACT_MIN_ID:-0.5}"
    comm -23 <(sort "${final_list}") <(sort "${filter_dir}/miniprot_intact_conflict.txt") \
        > "${filter_dir}/high_confidence_pseudogenes_tmp_3.txt"
    final_list="${filter_dir}/high_confidence_pseudogenes_tmp_3.txt"
    printf 'same_locus_aligner_conflict\tassessed\t%s\n' "${filter_dir}/miniprot_intact_conflict_report.tsv" >> "${filter_dir}/evidence_status.tsv"
fi

# Create final output files
cp -f "${filter_dir}/unitary_candidates.txt" "${filter_dir}/family_copy_pass_ids.txt"
# Retain the historical file name for existing downstream scripts.
cp -f "${filter_dir}/unitary_candidates.txt" "${filter_dir}/unitary_pseudogene_ids.txt"
cp -f "${final_list}" "${filter_dir}/final_candidate_ids.txt"
cp -f "${final_list}" "${filter_dir}/high_confidence_pseudogene_ids.txt"

# The paper manually assessed pairs of frameshifts because a compensating
# downstream event can restore the reading frame. Counts alone do not prove
# compensation, so keep these loci in a separate review queue.
awk '$2==0 && $3==2 {print $1}' "${filter_dir}/mutation_counts.tsv" | sort -u \
    > "${filter_dir}/two_frameshift_only_all.txt"
comm -12 <(sort "${final_list}") <(sort "${filter_dir}/two_frameshift_only_all.txt") \
    > "${filter_dir}/two_frameshift_only_review.txt"

echo "  Filtering done for ${SPECIES}."
echo "    putative       : $(wc -l < "${filter_dir}/putative_pseudogene_ids.txt" | tr -d " ")"
echo "    family/copy pass: $(wc -l < "${filter_dir}/family_copy_pass_ids.txt" | tr -d " ")"
echo "    high_confidence: $(wc -l < "${final_list}" | tr -d " ")"

if [[ -f "${sp_dir}/regions/locus_map.tsv" ]]; then
    python3 - "${sp_dir}/regions/locus_map.tsv" "${filter_dir}" <<'PYEXPORT'
import csv,sys
from pathlib import Path
folder=Path(sys.argv[2]);keep=set((folder/'high_confidence_pseudogene_ids.txt').read_text().splitlines())
with open(sys.argv[1]) as handle:rows=[r for r in csv.DictReader(handle,delimiter='\t') if r['locus_id'] in keep]
with (folder/'high_confidence_loci.bed').open('w') as handle:
    csv.writer(handle,delimiter='\t',lineterminator='\n').writerows((r['chrom'],r['start0'],r['end0'],r['locus_id'],0,r['strand']) for r in rows)
(folder/'high_confidence_queries.txt').write_text(''.join(q+'\n' for q in sorted({r['query_id'] for r in rows})))
PYEXPORT
fi
