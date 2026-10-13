#!/usr/bin/env bash
set -euo pipefail

# Batch GeneWise runner

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EASYPSEUDOGENE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Parse arguments
PROTEINS=""
GENOME=""
SPECIES=""
OUTPUT_DIR=""
THREADS=""
JOBS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --proteins) PROTEINS="$2"; shift 2 ;;
        --genome) GENOME="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

JOBS="${JOBS:-${THREADS}}"

sp_dir="${EPG_WORK_DIR:-${OUTPUT_DIR}/${SPECIES}}"
miniprot_dir="${sp_dir}/miniprot"
regions_dir="${sp_dir}/regions"
genewise_dir="${sp_dir}/genewise"
logs_dir="${EPG_LOG_DIR:-${sp_dir}/logs}"
tmp_dir="${sp_dir}/tmp"

source "${EASYPSEUDOGENE_DIR}/lib/genewise_env.sh"
epg_setup_genewise

mkdir -p "${regions_dir}" "${genewise_dir}" "${logs_dir}"

samtools faidx "${GENOME}"
samtools faidx "${PROTEINS}"

best_hits="${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv"
if [[ ! -f "${best_hits}" ]]; then
    echo "ERROR: missing ${best_hits}. Run screening first." >&2
    exit 1
fi

jobs_tsv="${genewise_dir}/${SPECIES}.jobs.tsv"
id_map="${regions_dir}/id_map.tsv"

python3 "${EASYPSEUDOGENE_DIR}/lib/prepare_genewise_jobs.py" \
    --proteins "${PROTEINS}" --hits "${miniprot_dir}/${SPECIES}.best_hits.tsv" \
    --flanked "${best_hits}" --regions "${regions_dir}" --jobs "${jobs_tsv}" \
    --paf "${miniprot_dir}/${SPECIES}.paf"

: > "${logs_dir}/genewise_failed.tsv"

export PROTEINS GENOME regions_dir genewise_dir logs_dir

echo "  Running GeneWise on $(wc -l < "${jobs_tsv}" | tr -d ' ') regions..."

# Export all needed variables (xargs should inherit them, but be explicit)
export PROTEINS GENOME regions_dir genewise_dir logs_dir
if [[ -n "${WISECONFIGDIR:-}" ]]; then
    export WISECONFIGDIR
fi

xargs_status=0
tr '\t\n' '\0' < "${jobs_tsv}" | xargs -0 -r -P "${JOBS}" -n 6 bash -c '
  set -euo pipefail
  # Re-export in sub-shell to ensure they are available
  export PROTEINS="${PROTEINS}" GENOME="${GENOME}" regions_dir="${regions_dir}" genewise_dir="${genewise_dir}" logs_dir="${logs_dir}"
  if [[ -n "${WISECONFIGDIR:-}" ]]; then
    export WISECONFIGDIR="${WISECONFIGDIR}"
  fi
  
  gene="$1"; strand="$2"; tname="$3"; start1="$4"; end1="$5"; safe_id="$6";
  out="${genewise_dir}/${safe_id}.gw.txt"
  err="${logs_dir}/${safe_id}.gw.err"
  pep="${regions_dir}/${safe_id}.pep.fa"
  dna="${regions_dir}/${safe_id}.dna.fa"

  # Standalone invocations recompute jobs. CLI reuse is guarded by content hashes.
  samtools faidx "${GENOME}" "${tname}:${start1}-${end1}" > "${dna}" 2>> "${logs_dir}/genewise_failed.tsv" || { printf "%s\tdna_missing\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; exit 1; }

  if [[ "${strand}" == "-" ]]; then
    "${GENEWISE_BIN}" -trev "${pep}" "${dna}" > "${out}.partial" 2> "${err}" || { printf "%s\tgenewise_failed\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; rm -f "${out}.partial"; exit 1; }
  else
    "${GENEWISE_BIN}" -tfor "${pep}" "${dna}" > "${out}.partial" 2> "${err}" || { printf "%s\tgenewise_failed\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; rm -f "${out}.partial"; exit 1; }
  fi
  if [[ ! -s "${out}.partial" ]] || ! grep -q "^Query protein:" "${out}.partial" || ! grep -q "^//" "${out}.partial"; then
    printf "%s\tinvalid_genewise_output\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"
    rm -f "${out}.partial"
    exit 1
  fi
  mv "${out}.partial" "${out}"
  touch "${out}.done"
' _ || xargs_status=$?

if (( xargs_status != 0 )) || [[ -s "${logs_dir}/genewise_failed.tsv" ]]; then
    echo "ERROR: GeneWise jobs failed; no combined result published. See ${logs_dir}/genewise_failed.tsv" >&2
    exit 1
fi

# Concatenate only this job set, never unrelated files left in a directory.
: > "${genewise_dir}/pseudogene-genewise.txt.partial"
while IFS=$'\t' read -r gene strand chrom start end locus; do
    cat "${genewise_dir}/${locus}.gw.txt" >> "${genewise_dir}/pseudogene-genewise.txt.partial"
done < "${jobs_tsv}"
mv "${genewise_dir}/pseudogene-genewise.txt.partial" "${genewise_dir}/pseudogene-genewise.txt"

success_count=$(wc -l < "${jobs_tsv}" | tr -d ' ')
failed_count=$(wc -l < "${logs_dir}/genewise_failed.tsv" 2>/dev/null || echo "0")
echo "  GeneWise done: ${success_count} successful alignments"
if [[ "${failed_count}" -gt 0 ]]; then
    echo "  Failed: ${failed_count} regions (see ${logs_dir}/genewise_failed.tsv for details)"
    if [[ "${failed_count}" -lt 20 ]]; then
        echo "  Failure reasons:"
        cut -f2 "${logs_dir}/genewise_failed.tsv" 2>/dev/null | sort | uniq -c | sort -rn || true
    fi
fi
