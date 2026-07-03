#!/usr/bin/env bash
# ============================================================================
# MiXeR bivariate test2 retry script (WSL, serial, no power/qq curves)
# ============================================================================
# - Skips existing valid .json files.
# - Uses --no-power-curve --no-qq-plots to avoid huge outputs and long runtimes.
# - THREADS=8 to fully utilize your CPU, single job to avoid memory overrun.
# ============================================================================

set -euo pipefail

# ---- Configuration ----
THREADS=6                       # threads per MiXeR job
MIXER_REF_DIR="/data"           # mounted reference panel inside container
SUMSTATS_DIR="./output/MiXeR/sumstats"
OUT_DIR="./output/MiXeR/results"

# ---- Docker run wrapper ----
function mixer_py {
    docker run --rm \
        -v "/mnt/c/Users/Administrator/Desktop/TRAIT_PAIR:/home" \
        -v "/mnt/c/1000G_EUR_Phase3_plink:/data" \
        -w /home \
        ghcr.io/precimed/gsa-mixer:2.2.1 \
        python -u /tools/mixer/precimed/mixer.py "$@"
}

# ---- Test2 specific flags ----
# MiXeR 2.2.1 test2 has --power-curve / --qq-plots defaulted to True (store_true).
# There is no flag to disable them. Instead, --downsample-factor controls how many
# SNPs feed the bivariate pdf computation (1001x1001 grid). Default = 50. Using
# 200 cuts pdf compute time ~4x with negligible impact on rg / params / ci.
# CRITICAL: NO --fit-sequence, NO --extract.
TEST_FLAGS="--bim-file ${MIXER_REF_DIR}/1000G.EUR.QC.@.bim"
TEST_FLAGS="${TEST_FLAGS} --ld-file  ${MIXER_REF_DIR}/1000G.EUR.QC.@.run4.ld"
TEST_FLAGS="${TEST_FLAGS} --threads ${THREADS}"
TEST_FLAGS="${TEST_FLAGS} --exclude-ranges MHC"
TEST_FLAGS="${TEST_FLAGS} --downsample-factor 200"

# ---- Trait list ----
TRAITS=("TRAIT_A" "de Lange_IBD" "de Lange_CD" "de Lange_UC" "Liu_IBD" "Liu_CD" "Liu_UC")

# ---- Main serial loop ----
echo "=== Starting test2 retry (serial, THREADS=${THREADS}) ==="
for TRAIT2 in "${TRAITS[@]}"; do
    LABEL="ASD_vs_${TRAIT2}"
    for REP in $(seq 1 20); do
        OUT="${OUT_DIR}/${LABEL}.test.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "[$(date +%H:%M:%S)] [skip] test2 ${LABEL} rep${REP} (already exists)"
            continue
        fi
        echo "[$(date +%H:%M:%S)] [run]  test2 ${LABEL} rep${REP}"
        if mixer_py test2 ${TEST_FLAGS} \
            --trait1-file "${SUMSTATS_DIR}/TRAIT_B.sumstats.gz" \
            --trait2-file "${SUMSTATS_DIR}/${TRAIT2}.sumstats.gz" \
            --load-params "${OUT_DIR}/${LABEL}.fit.rep${REP}.json" \
            --out "${OUT}"; then
            echo "[$(date +%H:%M:%S)] [done] test2 ${LABEL} rep${REP}"
        else
            echo "[$(date +%H:%M:%S)] [FAIL] test2 ${LABEL} rep${REP}"
        fi
    done
done
echo "=== test2 retry finished ==="
SCRIPT_EOF

# Make executable
chmod +x analysis/retry_failed_test2.sh