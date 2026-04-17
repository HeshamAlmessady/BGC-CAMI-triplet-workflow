#!/usr/bin/env bash
set -u
set -o pipefail

#########################################
# Usage
#########################################

if [[ $# -lt 1 ]]; then
    echo "Usage:"
    echo "bash run_quast_bgcquast_top5.sh <triplet_workflow_outdir> [TOP_N]"
    exit 1
fi

BASE="$(realpath "$1")"
TOP_N="${2:-5}"

#########################################
# Input files inside workflow output dir
#########################################

TOP_FILE="${BASE}/summary_ranked_hq_bgcrich.tsv"

HYB_DIR="${BASE}/good_hybrid_assembly_based_MAGs_nanosim"
SHORT_DIR="${BASE}/poor_short_assembly_based_MAGs"

#########################################
# External tools
#########################################

BGCQUAST_REPO="/net/projects/halmessa/BGC_project/bgc-quast"

#########################################
# Prepare
#########################################

mkdir -p "${BASE}/logs"

source ~/miniforge3/etc/profile.d/conda.sh

#########################################
# Helper function
#########################################

has_real_clusters() {
    local f="$1"
    [[ -f "$f" ]] && [[ $(wc -l < "$f") -gt 1 ]]
}

#########################################
# Check required files
#########################################

if [[ ! -f "${TOP_FILE}" ]]; then
    echo "ERROR: Cannot find ${TOP_FILE}"
    exit 1
fi

#########################################
# Main loop
#########################################

tail -n +1 "${TOP_FILE}" | head -n "${TOP_N}" | while IFS=$'\t' read -r \
OTU TAXID REF HYB SHORT REF_CONTIGS REF_TOTAL REF_N50 HYB_CONTIGS HYB_TOTAL HYB_N50 SHORT_CONTIGS SHORT_TOTAL SHORT_N50 REFHYB HYBSHORT
do
(
set -euo pipefail

SAFE_OTU="${OTU//./_}"

echo "=================================================="
echo "Processing ${OTU}"
echo "=================================================="

REF_OUTDIR="${BASE}/gecco_run_on_${OTU}_ref"
HYB_OUTDIR="${BASE}/gecco_run_on_${OTU}_hybrid_MAG"
SHORT_OUTDIR="${BASE}/gecco_run_on_${OTU}_short_MAG"

QUAST_OUT="${BASE}/assembly_to_ref_align_quast_${SAFE_OTU}"
BGCQUAST_OUT="${BASE}/bgcquast_${SAFE_OTU}_compare_to_reference"

HYB_SAFE_FASTA="${HYB_DIR}/hybrid_${SAFE_OTU}.fasta"
SHORT_SAFE_FASTA="${SHORT_DIR}/short_${SAFE_OTU}.fasta"

HYB_SAFE_CLUSTER="${HYB_OUTDIR}/hybrid_${SAFE_OTU}.clusters.tsv"
SHORT_SAFE_CLUSTER="${SHORT_OUTDIR}/short_${SAFE_OTU}.clusters.tsv"

REF_BASENAME="$(basename "${REF}" .fasta)"
REF_CLUSTER="${REF_OUTDIR}/${REF_BASENAME}.clusters.tsv"

mkdir -p "${REF_OUTDIR}" "${HYB_OUTDIR}" "${SHORT_OUTDIR}"

#########################################
# GECCO
#########################################

conda activate gecco_env

if [[ ! -f "${REF_CLUSTER}" ]]; then
    gecco run --genome "${REF}" -o "${REF_OUTDIR}"
fi

if [[ ! -f "${HYB_SAFE_CLUSTER}" ]]; then
    gecco run --genome "${HYB_SAFE_FASTA}" -o "${HYB_OUTDIR}" || true
fi

if [[ ! -f "${SHORT_SAFE_CLUSTER}" ]]; then
    gecco run --genome "${SHORT_SAFE_FASTA}" -o "${SHORT_OUTDIR}" || true
fi

conda deactivate

#########################################
# Skip empty BGC outputs
#########################################

if ! has_real_clusters "${REF_CLUSTER}"; then
    echo "[SKIP] ${OTU}: no ref BGCs"
    exit 0
fi

if ! has_real_clusters "${HYB_SAFE_CLUSTER}"; then
    echo "[SKIP] ${OTU}: no hybrid BGCs"
    exit 0
fi

if ! has_real_clusters "${SHORT_SAFE_CLUSTER}"; then
    echo "[SKIP] ${OTU}: no short BGCs"
    exit 0
fi

#########################################
# QUAST
#########################################

conda activate quast

if [[ ! -f "${QUAST_OUT}/report.tsv" ]]; then
    quast.py \
      "${HYB_SAFE_FASTA}" \
      "${SHORT_SAFE_FASTA}" \
      -r "${REF}" \
      -o "${QUAST_OUT}" \
      -t 16
fi

conda deactivate

#########################################
# BGC-QUAST
#########################################

conda activate bgc-quast

if [[ ! -f "${BGCQUAST_OUT}/report.tsv" ]]; then
python "${BGCQUAST_REPO}/bgc-quast.py" \
"${HYB_SAFE_CLUSTER}" \
"${SHORT_SAFE_CLUSTER}" \
--mode compare-to-reference \
--reference-mining-result "${REF_CLUSTER}" \
--quast-output-dir "${QUAST_OUT}" \
--output-dir "${BGCQUAST_OUT}"
fi

conda deactivate

echo "Done ${OTU}"

) || {
echo "[FAILED] ${OTU}"
continue
}

done
