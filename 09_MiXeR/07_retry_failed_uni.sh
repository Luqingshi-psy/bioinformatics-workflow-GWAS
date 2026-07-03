#!/usr/bin/env bash
# Retry all failed univariate fit1 and test1 tasks serially
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

# fit1 flags: include diffevo-fast
FIT_FLAGS="--bim-file ${MIXER_REF_DIR}/1000G.EUR.QC.@.bim"
FIT_FLAGS="${FIT_FLAGS} --ld-file  ${MIXER_REF_DIR}/1000G.EUR.QC.@.run4.ld"
FIT_FLAGS="${FIT_FLAGS} --threads ${THREADS}"
FIT_FLAGS="${FIT_FLAGS} --exclude-ranges MHC"
FIT_FLAGS="${FIT_FLAGS} --fit-sequence diffevo-fast"

# test1 flags: no fit-sequence
TEST_FLAGS="--bim-file ${MIXER_REF_DIR}/1000G.EUR.QC.@.bim"
TEST_FLAGS="${TEST_FLAGS} --ld-file  ${MIXER_REF_DIR}/1000G.EUR.QC.@.run4.ld"
TEST_FLAGS="${TEST_FLAGS} --threads ${THREADS}"
TEST_FLAGS="${TEST_FLAGS} --exclude-ranges MHC"

TRAITS=("TRAIT_A" "de Lange_IBD" "de Lange_CD" "de Lange_UC" "Liu_IBD" "Liu_CD" "Liu_UC")
ALL_TRAITS=("TRAIT_B" "${TRAITS[@]}")

echo "=== Retrying missing fit1 tasks ==="
for TRAIT in "${ALL_TRAITS[@]}"; do
    for REP in $(seq 1 20); do
        EXTRACT="--extract ${MIXER_REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"
        OUT="${OUT_DIR}/${TRAIT}.fit.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "[$(date +%H:%M:%S)] [skip] fit1 ${TRAIT} rep${REP}"
        else
            echo "[$(date +%H:%M:%S)] [run]  fit1 ${TRAIT} rep${REP}"
            if mixer_py fit1 ${FIT_FLAGS} ${EXTRACT} --trait1-file "${SUMSTATS_DIR}/${TRAIT}.sumstats.gz" --out "${OUT}"; then
                echo "[$(date +%H:%M:%S)] [done] fit1 ${TRAIT} rep${REP}"
            else
                echo "[$(date +%H:%M:%S)] [FAIL] fit1 ${TRAIT} rep${REP}"
            fi
        fi
    done
done

echo "=== Retrying missing test1 tasks ==="
for TRAIT in "${ALL_TRAITS[@]}"; do
    for REP in $(seq 1 20); do
        OUT="${OUT_DIR}/${TRAIT}.test.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "[$(date +%H:%M:%S)] [skip] test1 ${TRAIT} rep${REP}"
        else
            echo "[$(date +%H:%M:%S)] [run]  test1 ${TRAIT} rep${REP}"
            if mixer_py test1 ${TEST_FLAGS} --trait1-file "${SUMSTATS_DIR}/${TRAIT}.sumstats.gz" --load-params "${OUT_DIR}/${TRAIT}.fit.rep${REP}.json" --out "${OUT}"; then
                echo "[$(date +%H:%M:%S)] [done] test1 ${TRAIT} rep${REP}"
            else
                echo "[$(date +%H:%M:%S)] [FAIL] test1 ${TRAIT} rep${REP}"
            fi
        fi
    done
done
echo "Univariate retry completed."
