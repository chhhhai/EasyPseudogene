#!/usr/bin/env bash
set -euo pipefail

# Post-GeneWise filtering pipeline

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EASYPSEUDOGENE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${EASYPSEUDOGENE_DIR}/lib"

# Load default config
CONFIG_FILE="${EASYPSEUDOGENE_DIR}/config/default.conf"
if [[ -f "${CONFIG_FILE}" ]]; then
    source "${CONFIG_FILE}"
fi

# Parse arguments
PROTEINS=""
SPECIES=""
OUTPUT_DIR=""
REF_WHITELIST="${REF_WHITELIST:-}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --proteins) PROTEINS="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --whitelist) REF_WHITELIST="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

sp_dir="${OUTPUT_DIR}/${SPECIES}"
gw_dir="${sp_dir}/genewise"
miniprot_dir="${sp_dir}/miniprot"
filter_dir="${sp_dir}/post_genewise_filter"
prefilter_dir="${OUTPUT_DIR}/prefilter"

mkdir -p "${filter_dir}" "${prefilter_dir}"

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
    prefilter_cmd=(python3 "${LIB_DIR}/prefilter_proteins_with_whitelist.py")
    if [[ -n "${REF_WHITELIST}" && -s "${REF_WHITELIST}" ]]; then
        prefilter_cmd+=(--whitelist "${REF_WHITELIST}")
    else
        prefilter_cmd=(python3 "${LIB_DIR}/prefilter_proteins.py")
    fi
    "${prefilter_cmd[@]}" \
        --input "${PROTEINS}" \
        --output "${prefilter_dir}/proteins.filtered.fa" \
        --report "${prefilter_dir}/prefilter_report.tsv"
    grep "^>" "${prefilter_dir}/proteins.filtered.fa" | sed 's/^>//' | awk '{print $1}' | sort -u \
        > "${VALID_TRANSCRIPTS}"
else
    grep "^>" "${PROTEINS}" | sed 's/^>//' | awk '{print $1}' | sort -u \
        > "${VALID_TRANSCRIPTS}"
fi

comm -12 <(sort "${filter_dir}/putative_pseudogene_ids.raw.txt") <(sort "${VALID_TRANSCRIPTS}") \
    > "${filter_dir}/putative_pseudogene_ids.txt"

if [[ -n "${REF_WHITELIST}" && -s "${REF_WHITELIST}" ]]; then
    sort -u "${filter_dir}/putative_pseudogene_ids.txt" "${REF_WHITELIST}" \
        > "${filter_dir}/putative_pseudogene_ids.txt.tmp"
    mv "${filter_dir}/putative_pseudogene_ids.txt.tmp" "${filter_dir}/putative_pseudogene_ids.txt"
fi

unitary_list="${filter_dir}/putative_pseudogene_ids.txt"

# Helper function to merge whitelist
merge_whitelist() {
    local list_file="$1"
    if [[ -n "${REF_WHITELIST}" && -s "${REF_WHITELIST}" && -s "${list_file}" ]]; then
        sort -u "${list_file}" "${REF_WHITELIST}" > "${list_file}.tmp"
        mv "${list_file}.tmp" "${list_file}"
    fi
}

# Unitary filters
if [[ "${ENABLE_FAMILY_FILTER:-1}" == "1" ]]; then
    python3 "${LIB_DIR}/filter_large_family.py" \
        --proteins "${PROTEINS}" \
        --out-list "${filter_dir}/large_family_ids.txt" \
        --out-report "${filter_dir}/large_family_report.tsv" \
        --patterns "${FAMILY_PATTERNS:-olfactory receptor|vomeronasal|zinc finger|^OR|^V1R|^V2R|^TAAR|^ZNF|^ZFP}"
    comm -23 <(sort "${unitary_list}") <(sort "${filter_dir}/large_family_ids.txt") \
        > "${filter_dir}/unitary_step1.no_family.txt"
    unitary_list="${filter_dir}/unitary_step1.no_family.txt"
    merge_whitelist "${unitary_list}"
fi

if [[ "${ENABLE_REDUNDANCY_FILTER:-1}" == "1" ]]; then
    python3 "${LIB_DIR}/filter_functional_redundancy_miniprot.py" \
        --paf "${miniprot_dir}/${SPECIES}.paf" \
        --out-prefix "${filter_dir}/miniprot_redundancy" \
        --min-aln-len "${MINIPROT_MIN_ALN_LEN:-30}" \
        --min-cov "${MINIPROT_MIN_COV:-0.3}" \
        --min-id "${MINIPROT_MIN_ID:-0.3}" \
        --min-mapq "${MINIPROT_MIN_MAPQ:-0}" \
        --copy-cov "${REDUNDANCY_COPY_COV:-0.8}" \
        --target-gap "${REDUNDANCY_TARGET_GAP:-10000}"
    comm -23 <(sort "${unitary_list}") <(sort "${filter_dir}/miniprot_redundancy_redundant.txt") \
        > "${filter_dir}/unitary_step2.no_redundancy.txt"
    unitary_list="${filter_dir}/unitary_step2.no_redundancy.txt"
    merge_whitelist "${unitary_list}"
fi

cp -f "${unitary_list}" "${filter_dir}/unitary_candidates.txt"
merge_whitelist "${filter_dir}/unitary_candidates.txt"

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
            --flank-size "${LOW_ID_FLANK_SIZE:-30}" \
            --flank-th "${LOW_ID_FLANK_TH:-0.2}" \
            --overall-th "${LOW_ID_OVERALL_TH:-0.1}" \
            --mode "${LOW_ID_MODE:-and}"
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
    merge_whitelist "${filter_dir}/high_confidence_pseudogenes.txt"
else
    cp -f "${filter_dir}/unitary_candidates.txt" "${filter_dir}/high_confidence_pseudogenes.txt"
    merge_whitelist "${filter_dir}/high_confidence_pseudogenes.txt"
fi

final_list="${filter_dir}/high_confidence_pseudogenes.txt"

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
        --out-list "${filter_dir}/multi_disruption_low_quality.txt" \
        --out-report "${filter_dir}/multi_disruption_quality_report.tsv"
    comm -23 <(sort "${final_list}") <(sort "${filter_dir}/multi_disruption_low_quality.txt") \
        > "${filter_dir}/high_confidence_pseudogenes_tmp_1.txt"
    final_list="${filter_dir}/high_confidence_pseudogenes_tmp_1.txt"
    merge_whitelist "${final_list}"
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
    merge_whitelist "${final_list}"
fi

# Create final output files
cp -f "${filter_dir}/unitary_candidates.txt" "${filter_dir}/unitary_pseudogene_ids.txt"
cp -f "${final_list}" "${filter_dir}/high_confidence_pseudogene_ids.txt"

echo "  Filtering done for ${SPECIES}."
echo "    putative       : $(wc -l < "${filter_dir}/putative_pseudogene_ids.txt" | tr -d " ")"
echo "    unitary        : $(wc -l < "${filter_dir}/unitary_candidates.txt" | tr -d " ")"
echo "    high_confidence: $(wc -l < "${final_list}" | tr -d " ")"
