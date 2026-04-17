#!/usr/bin/env bash
set -u
set -o pipefail

BASE="/net/projects/halmessa/BGC_project/plant_BGC_paper/run_from_script"
TOP_FILE="${BASE}/summary_ranked_hq_bgcrich.tsv"
TOP_N=5

HYB_DIR="${BASE}/good_hybrid_assembly_based_MAGs_nanosim"
SHORT_DIR="${BASE}/poor_short_assembly_based_MAGs"

LOGDIR="/net/projects/halmessa/BGC_project/logs_gecco_hybrid"
BGCQUAST_REPO="/net/projects/halmessa/BGC_project/bgc-quast"

mkdir -p "${LOGDIR}" "${BASE}/logs"

source ~/miniforge3/etc/profile.d/conda.sh

has_real_clusters() {
    local f="$1"
    [[ -f "$f" ]] && [[ $(wc -l < "$f") -gt 1 ]]
}

tail -n +2 "${TOP_FILE}" | head -n "${TOP_N}" | while IFS=$'\t' read -r \
    OTU TAXID REF HYB SHORT REF_CONTIGS REF_TOTAL REF_N50 HYB_CONTIGS HYB_TOTAL HYB_N50 SHORT_CONTIGS SHORT_TOTAL SHORT_N50 REFHYB HYBSHORT
do
    (
        set -euo pipefail

        SAFE_OTU="${OTU//./_}"

        echo "=================================================="
        echo "Processing ${OTU}"
        echo "Reference: ${REF}"
        echo "Hybrid MAG: ${HYB}"
        echo "Short MAG: ${SHORT}"
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

        # Prepare sanitized MAG filenames
        if [[ ! -f "${HYB_SAFE_FASTA}" ]]; then
            cp "${HYB}" "${HYB_SAFE_FASTA}"
        fi

        if [[ ! -f "${SHORT_SAFE_FASTA}" ]]; then
            cp "${SHORT}" "${SHORT_SAFE_FASTA}"
        fi

        # -----------------------------
        # 1. GECCO
        # -----------------------------
        conda activate gecco_env

        echo "[GECCO] reference ${OTU}"
        if [[ ! -f "${REF_CLUSTER}" ]]; then
            gecco run \
              --genome "${REF}" \
              -o "${REF_OUTDIR}" \
              > "${BASE}/logs/${SAFE_OTU}.gecco_ref.stdout.log" \
              2> "${BASE}/logs/${SAFE_OTU}.gecco_ref.stderr.log"
        else
            echo "  Reference GECCO output exists, skipping."
        fi

        echo "[GECCO] hybrid ${OTU}"
        if [[ ! -f "${HYB_SAFE_CLUSTER}" ]]; then
            gecco run \
              --genome "${HYB_SAFE_FASTA}" \
              -o "${HYB_OUTDIR}" \
              > "${BASE}/logs/${SAFE_OTU}.gecco_hybrid.stdout.log" \
              2> "${BASE}/logs/${SAFE_OTU}.gecco_hybrid.stderr.log"

            if [[ -f "${HYB_OUTDIR}/${OTU}.clusters.tsv" && ! -f "${HYB_SAFE_CLUSTER}" ]]; then
                cp "${HYB_OUTDIR}/${OTU}.clusters.tsv" "${HYB_SAFE_CLUSTER}"
            fi
        else
            echo "  Hybrid GECCO output exists, skipping."
        fi

        echo "[GECCO] short ${OTU}"
        if [[ ! -f "${SHORT_SAFE_CLUSTER}" ]]; then
            gecco run \
              --genome "${SHORT_SAFE_FASTA}" \
              -o "${SHORT_OUTDIR}" \
              > "${BASE}/logs/${SAFE_OTU}.gecco_short.stdout.log" \
              2> "${BASE}/logs/${SAFE_OTU}.gecco_short.stderr.log" || true

            if [[ -f "${SHORT_OUTDIR}/${OTU}.clusters.tsv" && ! -f "${SHORT_SAFE_CLUSTER}" ]]; then
                cp "${SHORT_OUTDIR}/${OTU}.clusters.tsv" "${SHORT_SAFE_CLUSTER}"
            fi
        else
            echo "  Short GECCO output exists, skipping."
        fi

        conda deactivate

        # -----------------------------
        # Check whether GECCO found BGCs
        # -----------------------------
        if ! has_real_clusters "${REF_CLUSTER}"; then
            echo "[SKIP] ${OTU}: reference has no detectable BGCs" | tee -a "${BASE}/logs/skipped_candidates.log"
            exit 0
        fi

        if ! has_real_clusters "${HYB_SAFE_CLUSTER}"; then
            echo "[SKIP] ${OTU}: hybrid MAG has no detectable BGCs" | tee -a "${BASE}/logs/skipped_candidates.log"
            exit 0
        fi

        if ! has_real_clusters "${SHORT_SAFE_CLUSTER}"; then
            echo "[SKIP] ${OTU}: short MAG has no detectable BGCs" | tee -a "${BASE}/logs/skipped_candidates.log"
            exit 0
        fi

        # -----------------------------
        # 2. QUAST
        # -----------------------------
        conda activate quast

        echo "[QUAST] ${OTU}"
        if [[ ! -f "${QUAST_OUT}/report.tsv" ]]; then
            quast.py \
              "${HYB_SAFE_FASTA}" \
              "${SHORT_SAFE_FASTA}" \
              -r "${REF}" \
              -o "${QUAST_OUT}" \
              -t 16 \
              > "${BASE}/logs/${SAFE_OTU}.quast.stdout.log" \
              2> "${BASE}/logs/${SAFE_OTU}.quast.stderr.log"
        else
            echo "  QUAST output exists, skipping."
        fi

        conda deactivate

        # -----------------------------
        # 3. BGC-QUAST
        # -----------------------------
        conda activate bgc-quast

        echo "[BGC-QUAST] ${OTU}"
        if [[ ! -f "${BGCQUAST_OUT}/report.tsv" ]]; then
            python "${BGCQUAST_REPO}/bgc-quast.py" \
              "${HYB_SAFE_CLUSTER}" \
              "${SHORT_SAFE_CLUSTER}" \
              --mode compare-to-reference \
              --reference-mining-result "${REF_CLUSTER}" \
              --quast-output-dir "${QUAST_OUT}" \
              --output-dir "${BGCQUAST_OUT}" \
              > "${BASE}/logs/${SAFE_OTU}.bgcquast.stdout.log" \
              2> "${BASE}/logs/${SAFE_OTU}.bgcquast.stderr.log"
        else
            echo "  BGC-QUAST output exists, skipping."
        fi

        conda deactivate

        echo "Done: ${OTU}"
        echo
    ) || {
        echo "[ERROR] ${OTU} failed; continuing to next candidate" | tee -a "${BASE}/logs/failed_candidates.log"
        continue
    }
done
