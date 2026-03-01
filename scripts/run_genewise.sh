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

sp_dir="${OUTPUT_DIR}/${SPECIES}"
miniprot_dir="${sp_dir}/miniprot"
regions_dir="${sp_dir}/regions"
genewise_dir="${sp_dir}/genewise"
logs_dir="${sp_dir}/logs"
tmp_dir="${sp_dir}/tmp"

mkdir -p "${regions_dir}" "${genewise_dir}" "${logs_dir}" "${tmp_dir}"

samtools faidx "${GENOME}"
samtools faidx "${PROTEINS}"

best_hits="${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv"
if [[ ! -s "${best_hits}" ]]; then
    echo "ERROR: missing ${best_hits}. Run screening first." >&2
    exit 1
fi

jobs_tsv="${genewise_dir}/${SPECIES}.jobs.tsv"
id_map="${regions_dir}/id_map.tsv"

awk -v OFS="\t" '{
  gene=$1; safe=gene; gsub(/[^A-Za-z0-9_.-]/,"_",safe);
  print safe, gene;
}' "${best_hits}" | sort -u > "${id_map}"

awk -v OFS="\t" '{
  gene=$1; strand=$2; tname=$3; start0=$4; end0=$5;
  safe=gene; gsub(/[^A-Za-z0-9_.-]/,"_",safe);
  start1=start0+1; end1=end0;
  if (start1 < 1 || end1 < start1) next;
  print gene, strand, tname, start1, end1, safe;
}' "${best_hits}" > "${jobs_tsv}"

: > "${logs_dir}/genewise_failed.tsv"

export PROTEINS GENOME regions_dir genewise_dir logs_dir

# Set GeneWise environment if provided
if [[ -n "${WISECONFIGDIR:-}" ]]; then
    export WISECONFIGDIR
    echo "  Using WISECONFIGDIR=${WISECONFIGDIR}"
fi
if [[ -n "${WISEPATH:-}" ]]; then
    export PATH="${WISEPATH}:${PATH}"
    echo "  Added ${WISEPATH} to PATH"
fi

# Check if genewise is available
if ! command -v genewise >/dev/null 2>&1; then
    echo "ERROR: genewise command not found in PATH" >&2
    echo "  Please set --wiseconfig and --wisepath when running easypseudogene" >&2
    exit 1
fi

echo "  Running GeneWise on $(wc -l < "${jobs_tsv}" | tr -d ' ') regions..."

# Export all needed variables (xargs should inherit them, but be explicit)
export PROTEINS GENOME regions_dir genewise_dir logs_dir
if [[ -n "${WISECONFIGDIR:-}" ]]; then
    export WISECONFIGDIR
fi

cat "${jobs_tsv}" | xargs -P "${JOBS}" -n 6 bash -c '
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

  if [[ -s "${out}" ]]; then
    exit 0
  fi

  if [[ ! -s "${pep}" ]]; then
    samtools faidx "${PROTEINS}" "${gene}" > "${pep}" 2>> "${logs_dir}/genewise_failed.tsv" || { printf "%s\tpep_missing\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; exit 0; }
  fi
  if [[ ! -s "${dna}" ]]; then
    samtools faidx "${GENOME}" "${tname}:${start1}-${end1}" > "${dna}" 2>> "${logs_dir}/genewise_failed.tsv" || { printf "%s\tdna_missing\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; exit 0; }
  fi

  if [[ "${strand}" == "-" ]]; then
    genewise -trev "${pep}" "${dna}" > "${out}" 2> "${err}" || { printf "%s\tgenewise_failed\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; rm -f "${out}"; exit 0; }
  else
    genewise -tfor "${pep}" "${dna}" > "${out}" 2> "${err}" || { printf "%s\tgenewise_failed\n" "${gene}" >> "${logs_dir}/genewise_failed.tsv"; rm -f "${out}"; exit 0; }
  fi
' _

find "${genewise_dir}" -maxdepth 1 -name "*.gw.txt" -type f -size +0c | sort | xargs -r cat \
  > "${genewise_dir}/pseudogene-genewise.txt"

success_count=$(find "${genewise_dir}" -maxdepth 1 -name "*.gw.txt" -type f -size +0c | wc -l | tr -d ' ')
failed_count=$(wc -l < "${logs_dir}/genewise_failed.tsv" 2>/dev/null || echo "0")
echo "  GeneWise done: ${success_count} successful alignments"
if [[ "${failed_count}" -gt 0 ]]; then
    echo "  Failed: ${failed_count} regions (see ${logs_dir}/genewise_failed.tsv for details)"
    if [[ "${failed_count}" -lt 20 ]]; then
        echo "  Failure reasons:"
        cut -f2 "${logs_dir}/genewise_failed.tsv" 2>/dev/null | sort | uniq -c | sort -rn || true
    fi
fi
