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

while [[ $# -gt 0 ]]; do
    case $1 in
        --proteins) PROTEINS="$2"; shift 2 ;;
        --genome) GENOME="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --flank) FLANK="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

sp_dir="${OUTPUT_DIR}/${SPECIES}"
mmseqs_dir="${sp_dir}/mmseqs"
miniprot_dir="${sp_dir}/miniprot"
regions_dir="${sp_dir}/regions"
logs_dir="${sp_dir}/logs"
tmp_dir="${sp_dir}/tmp"

mkdir -p "${mmseqs_dir}" "${miniprot_dir}" "${regions_dir}" "${logs_dir}" "${tmp_dir}"

samtools faidx "${PROTEINS}"
samtools faidx "${GENOME}"

# 1) mmseqs2: fast protein-to-genome screening
echo "  Running mmseqs2..."
mmseqs createdb "${PROTEINS}" "${mmseqs_dir}/query_db"
mmseqs createdb "${GENOME}" "${mmseqs_dir}/target_db"
mmseqs search "${mmseqs_dir}/query_db" "${mmseqs_dir}/target_db" "${mmseqs_dir}/result_db" "${tmp_dir}" \
  --search-type "${MMSEQS_SEARCH_TYPE}" -s "${MMSEQS_SENS}" -e "${MMSEQS_EVALUE}" \
  --max-seqs "${MMSEQS_MAX_SEQS}" --threads "${THREADS}"
mmseqs convertalis "${mmseqs_dir}/query_db" "${mmseqs_dir}/target_db" "${mmseqs_dir}/result_db" \
  "${mmseqs_dir}/${SPECIES}.m8" \
  --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits"

# 2) miniprot: spliced protein-to-genome alignment
echo "  Running miniprot..."
miniprot -t "${THREADS}" "${GENOME}" "${PROTEINS}" > "${miniprot_dir}/${SPECIES}.paf"

# 3) Pick best miniprot hit per protein
echo "  Selecting best hits..."
awk -v min_len="${MIN_ALN_LEN}" -v min_mapq="${MIN_MAPQ}" 'BEGIN{OFS="\t"}{
    q=$1; strand=$5; t=$6; tlen=$7; tstart=$8; tend=$9; alen=$11; mapq=$12;
    if (alen < min_len || mapq < min_mapq) next;
    score=alen;
    if (!(q in best) || score > best[q]) {
      best[q]=score;
      line[q]=q OFS strand OFS t OFS tstart OFS tend OFS tlen OFS alen OFS mapq;
    }
  }
  END{for (q in line) print line[q];}' \
  "${miniprot_dir}/${SPECIES}.paf" | sort -k1,1 > "${miniprot_dir}/${SPECIES}.best_hits.tsv"

# 4) Add flanks
awk -v flank="${FLANK}" 'BEGIN{OFS="\t"}{
    start=$4 - flank; if (start < 0) start=0;
    end=$5 + flank; if (end > $6) end=$6;
    print $1,$2,$3,start,end,$6,$7,$8;
  }' \
  "${miniprot_dir}/${SPECIES}.best_hits.tsv" > "${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv"

# 5) Extract region/protein FASTA
echo "  Extracting regions..."
: > "${regions_dir}/id_map.tsv"
: > "${logs_dir}/skipped.tsv"

while IFS=$'\t' read -r gene strand tname start0 end0 tlen alen mapq; do
  safe_id="$(echo "${gene}" | sed 's/[^A-Za-z0-9_.-]/_/g')"
  start1=$((start0 + 1))
  end1=$((end0))

  if (( start1 < 1 || end1 < start1 )); then
    printf "%s\t%s\t%s\t%s\t%s\n" "${gene}" "${tname}" "${start1}" "${end1}" "invalid_region" >> "${logs_dir}/skipped.tsv"
    continue
  fi

  echo -e "${safe_id}\t${gene}" >> "${regions_dir}/id_map.tsv"
  samtools faidx "${GENOME}" "${tname}:${start1}-${end1}" > "${regions_dir}/${safe_id}.dna.fa"
  samtools faidx "${PROTEINS}" "${gene}" > "${regions_dir}/${safe_id}.pep.fa"
done < "${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv"

echo "  Pre-screening done: $(wc -l < "${miniprot_dir}/${SPECIES}.best_hits.flanked.tsv" | tr -d ' ') candidate regions"
