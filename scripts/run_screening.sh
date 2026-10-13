#!/usr/bin/env bash
set -euo pipefail

# Pre-GeneWise screening using mmseqs2 and miniprot

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EASYPSEUDOGENE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${EASYPSEUDOGENE_DIR}/lib"

# Parse arguments
PROTEINS=""
GENOME=""
SPECIES=""
OUTPUT_DIR=""
THREADS=""
FLANK="${FLANK:-1000}"
MIN_ALN_LEN="${MIN_ALN_LEN:-30}"
MIN_MAPQ="${MIN_MAPQ:-0}"
MMSEQS_SEARCH_TYPE="${MMSEQS_SEARCH_TYPE:-3}"
MMSEQS_SENS="${MMSEQS_SENS:-7.5}"
MMSEQS_EVALUE="${MMSEQS_EVALUE:-1e-3}"
MMSEQS_MAX_SEQS="${MMSEQS_MAX_SEQS:-5}"
SKIP_MMSEQS="${SKIP_MMSEQS:-0}"
MINIPROT_MAX_HITS="${MINIPROT_MAX_HITS:-1}"
MINIPROT_OUTN="${MINIPROT_OUTN:-1000}"
MINIPROT_OUTS="${MINIPROT_OUTS:-0.99}"
SKIP_REGION_EXTRACT="${SKIP_REGION_EXTRACT:-0}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --proteins) PROTEINS="$2"; shift 2 ;;
        --genome) GENOME="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --flank) FLANK="$2"; shift 2 ;;
        --skip-mmseqs) SKIP_MMSEQS=1; shift ;;
        --max-hits-per-protein) MINIPROT_MAX_HITS="$2"; shift 2 ;;
        --miniprot-outn) MINIPROT_OUTN="$2"; shift 2 ;;
        --miniprot-outs) MINIPROT_OUTS="$2"; shift 2 ;;
        --skip-region-extract) SKIP_REGION_EXTRACT=1; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

sp_dir="${EPG_WORK_DIR:-${OUTPUT_DIR}/${SPECIES}}"
mmseqs_dir="${sp_dir}/mmseqs"
miniprot_dir="${sp_dir}/miniprot"
regions_dir="${sp_dir}/regions"
logs_dir="${EPG_LOG_DIR:-${sp_dir}/logs}"
tmp_dir="${sp_dir}/tmp"

mkdir -p "${miniprot_dir}" "${logs_dir}"

samtools faidx "${PROTEINS}"
samtools faidx "${GENOME}"

# 1) mmseqs2: optional fast protein-to-genome screening. Downstream selection
# currently uses miniprot, so users can skip mmseqs to improve runtime.
if [[ "${SKIP_MMSEQS}" == "1" ]]; then
  echo "  Skipping mmseqs2 (--skip-mmseqs)..."
else
  echo "  Running mmseqs2..."
  mkdir -p "${mmseqs_dir}" "${tmp_dir}"
  mmseqs createdb "${PROTEINS}" "${mmseqs_dir}/query_db"
  mmseqs createdb "${GENOME}" "${mmseqs_dir}/target_db"
  mmseqs search "${mmseqs_dir}/query_db" "${mmseqs_dir}/target_db" "${mmseqs_dir}/result_db" "${tmp_dir}" \
    --search-type "${MMSEQS_SEARCH_TYPE}" -s "${MMSEQS_SENS}" -e "${MMSEQS_EVALUE}" \
    --max-seqs "${MMSEQS_MAX_SEQS}" --threads "${THREADS}"
  mmseqs convertalis "${mmseqs_dir}/query_db" "${mmseqs_dir}/target_db" "${mmseqs_dir}/result_db" \
    "${mmseqs_dir}/${SPECIES}.m8" \
    --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits"
fi

# 2) miniprot: spliced protein-to-genome alignment
echo "  Running miniprot..."
miniprot -t "${THREADS}" --outn="${MINIPROT_OUTN}" --outs="${MINIPROT_OUTS}" "${GENOME}" "${PROTEINS}" > "${miniprot_dir}/${SPECIES}.paf"

# 3) Pick top miniprot hit(s) per protein. Default keeps historical behavior
# (1 hit/query); --max-hits-per-protein increases recall for cataloging mode.
echo "  Selecting top ${MINIPROT_MAX_HITS} miniprot hit(s) per protein..."
awk -v min_len="${MIN_ALN_LEN}" -v min_mapq="${MIN_MAPQ}" 'BEGIN{OFS="\t"}{
    q=$1; qlen=$2; qstart=$3; qend=$4; strand=$5; t=$6; tlen=$7; tstart=$8; tend=$9; nmatch=$10; alen=$11; mapq=$12;
    AS="NA"; ms="NA"; np="NA"; fs="0"; st="0"; cg="NA"; cs="NA";
    for (i=13; i<=NF; i++) {
      if ($i ~ /^AS:i:/) {split($i,a,":"); AS=a[3]}
      else if ($i ~ /^ms:i:/) {split($i,a,":"); ms=a[3]}
      else if ($i ~ /^np:i:/) {split($i,a,":"); np=a[3]}
      else if ($i ~ /^fs:i:/) {split($i,a,":"); fs=a[3]}
      else if ($i ~ /^st:i:/) {split($i,a,":"); st=a[3]}
      else if ($i ~ /^cg:Z:/) {cg=substr($i,6)}
      else if ($i ~ /^cs:Z:/) {cs=substr($i,6)}
    }
    if (alen < min_len || mapq < min_mapq) next;
    n++;
    payload = q OFS strand OFS t OFS tstart OFS tend OFS tlen OFS alen OFS mapq OFS qlen OFS qstart OFS qend OFS nmatch OFS AS OFS ms OFS np OFS fs OFS st OFS cg OFS cs;
    print q, alen, n, payload;
  }' "${miniprot_dir}/${SPECIES}.paf" \
  | sort -k1,1 -k2,2nr -k3,3n \
  | awk -v max_hits="${MINIPROT_MAX_HITS}" 'BEGIN{FS="\t"; OFS="\t"}{
      q=$1; count[q]++; if (count[q] <= max_hits) {for (i=4;i<=NF;i++) printf "%s%s", $i, (i<NF?OFS:ORS)}
    }' \
  | sort -k1,1 -k3,3 -k4,4n > "${miniprot_dir}/${SPECIES}.best_hits.tsv"

# 4) Add flanks
awk -v flank="${FLANK}" 'BEGIN{OFS="\t"}{
    start=$4 - flank; if (start < 0) start=0;
    end=$5 + flank; if (end > $6) end=$6;
    print $1,$2,$3,start,end,$6,$7,$8;
    # Preserve optional evidence columns by appending them to the just-written line via a second pass below is not needed for GeneWise; fast classifier reads unflanked best_hits.tsv.
  }' \
  "${miniprot_dir}/${SPECIES}.best_hits.tsv" > "${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv"
# Region extraction is performed once, by run_genewise.sh, after self filtering.

echo "  Pre-screening done: $(wc -l < "${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv" | tr -d ' ') candidate regions"
