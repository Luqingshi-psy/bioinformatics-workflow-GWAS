#!/usr/bin/env bash
# Retry all failed bivariate fit2 tasks one at a time
set -euo pipefail

cd ${PROJECT_ROOT}

THREADS=1
MIXER_REF_DIR="/data"
SUMSTATS_DIR="./output/MiXeR/sumstats"
OUT_DIR="./output/MiXeR/results"

function mixer_py {
    MSYS_NO_PATHCONV=1 docker run --rm \
        -v "${PROJECT_ROOT}/TRAIT_PAIR:/home" \
        -v "C:\\1000G_EUR_Phase3_plink:/data" \
        -w /home \
        ghcr.io/precimed/gsa-mixer:2.2.1 \
        python /tools/mixer/precimed/mixer.py "$@"
}

COMMON_FLAGS="--bim-file ${MIXER_REF_DIR}/1000G.EUR.QC.@.bim"
COMMON_FLAGS="${COMMON_FLAGS} --ld-file  ${MIXER_REF_DIR}/1000G.EUR.QC.@.run4.ld"
COMMON_FLAGS="${COMMON_FLAGS} --threads ${THREADS}"
COMMON_FLAGS="${COMMON_FLAGS} --exclude-ranges MHC"
COMMON_FLAGS="${COMMON_FLAGS} --fit-sequence diffevo-fast"

TRAITS=("TRAIT_A" "de Lange_IBD" "de Lange_CD" "de Lange_UC" "Liu_IBD" "Liu_CD" "Liu_UC")

echo "=== Retrying all failed fit2 tasks (serial mode) ==="
for TRAIT2 in "${TRAITS[@]}"; do
    LABEL="ASD_vs_${TRAIT2}"
    for REP in $(seq 1 20); do
        EXTRACT="--extract ${MIXER_REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"
        OUT="${OUT_DIR}/${LABEL}.fit.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "[$(date +%H:%M:%S)] [skip] ${LABEL} rep${REP} (already exists)"
        else
            echo "[$(date +%H:%M:%S)] [run]  ${LABEL} rep${REP}"
            if mixer_py fit2 ${COMMON_FLAGS} ${EXTRACT} \
                --trait1-file "${SUMSTATS_DIR}/TRAIT_B.sumstats.gz" \
                --trait2-file "${SUMSTATS_DIR}/${TRAIT2}.sumstats.gz" \
                --trait1-params "${OUT_DIR}/ASD.fit.rep${REP}.json" \
                --trait2-params "${OUT_DIR}/${TRAIT2}.fit.rep${REP}.json" \
                --out "${OUT}"; then
                echo "[$(date +%H:%M:%S)] [done] ${LABEL} rep${REP}"
            else
                echo "[$(date +%H:%M:%S)] [FAIL] ${LABEL} rep${REP}"
            fi
        fi
    done
done
echo "Retry finished. Check for any remaining FAILs."