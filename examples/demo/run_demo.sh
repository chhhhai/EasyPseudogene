#!/usr/bin/env bash
set -euo pipefail

demo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${demo_dir}/../.." && pwd)"
if (( $# > 2 )); then
    echo "Usage: bash examples/demo/run_demo.sh [NEW_OUTPUT_DIRECTORY] [THREADS]" >&2
    exit 2
fi
output="${1:-demo_output}"
threads="${2:-2}"
[[ "${threads}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be positive" >&2; exit 2; }
if [[ -e "${output}" || -L "${output}" ]]; then
    echo "ERROR: demo output already exists: ${output}. Choose a fresh directory." >&2
    exit 1
fi
command -v python3 >/dev/null 2>&1 || { echo "ERROR: Python 3 is required" >&2; exit 1; }
command -v tee >/dev/null 2>&1 || { echo "ERROR: tee is required" >&2; exit 1; }
output="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "${output}")"
mkdir -p "$(dirname "${output}")"
mkdir "${output}"
exec > >(tee "${output}/demo.log") 2>&1
echo "Synthetic installation demo. Output: ${output}"
start="${SECONDS}"

bash "${repo_dir}/bin/easypseudogene" --mode cross \
    --proteins "${demo_dir}/data/proteins.fa" \
    --genome "${demo_dir}/data/cross_genome.fa" \
    --threads "${threads}" --species cross --output "${output}" --visualize

bash "${repo_dir}/bin/easypseudogene" --mode self \
    --proteins "${demo_dir}/data/proteins.fa" \
    --genome "${demo_dir}/data/self_genome.fa" \
    --gff "${demo_dir}/data/self_annotation.gff3" \
    --max-hits-per-protein 20 --miniprot-outs 0.1 \
    --threads "${threads}" --species self --output "${output}"

python3 "${demo_dir}/check_results.py" --output "${output}"
echo "Demo completed in $(( SECONDS - start )) seconds."
echo "Cross report: ${output}/cross/report.html"
echo "Self report:  ${output}/self/report.html"
