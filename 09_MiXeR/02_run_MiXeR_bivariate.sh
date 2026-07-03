#!/usr/bin/env bash
# MiXeR bivariate analysis: TRAIT_B × gut diseases (Bash parallel)
# Run from Git Bash in ${PROJECT_ROOT} pair (A x B)
# Prerequisites: Docker Desktop, MSYS_NO_PATHCONV=1

set -euo pipefail

# ---- Configuration ----
N_REP=20
THREADS=2          # threads per MiXeR job
PARALLEL_JOBS=2    # max simultaneous jobs

# Paths inside the container
MIXER_REF_DIR="/data"
SUMSTATS_DIR="./output/MiXeR/sumstats"
OUT_DIR="./output/MiXeR/results"
mkdir -p "${OUT_DIR}"

# Docker wrapper functions with MSYS_NO_PATHCONV=1
function mixer_py {
    MSYS_NO_PATHCONV=1 docker run --rm \
        -v "${PROJECT_ROOT}/TRAIT_PAIR:/home" \
        -v "C:\\1000G_EUR_Phase3_plink:/data" \
        -w /home \
        ghcr.io/precimed/gsa-mixer:2.2.1 \
        python /tools/mixer/precimed/mixer.py "$@"
}

function mixer_figs_py {
    MSYS_NO_PATHCONV=1 docker run --rm \
        -v "${PROJECT_ROOT}/TRAIT_PAIR:/home" \
        -v "C:\\1000G_EUR_Phase3_plink:/data" \
        -w /home \
        ghcr.io/precimed/gsa-mixer:2.2.1 \
        python /tools/mixer/precimed/mixer_figures.py "$@"
}

COMMON_FLAGS="--bim-file ${MIXER_REF_DIR}/1000G.EUR.QC.@.bim"
COMMON_FLAGS="${COMMON_FLAGS} --ld-file  ${MIXER_REF_DIR}/1000G.EUR.QC.@.run4.ld"
COMMON_FLAGS="${COMMON_FLAGS} --threads ${THREADS}"
COMMON_FLAGS="${COMMON_FLAGS} --exclude-ranges MHC"
COMMON_FLAGS="${COMMON_FLAGS} --fit-sequence diffevo-fast"

TRAITS=("TRAIT_A" "de Lange_IBD" "de Lange_CD" "de Lange_UC" "Liu_IBD" "Liu_CD" "Liu_UC")
ALL_TRAITS=("TRAIT_B" "${TRAITS[@]}")

# ---- Parallel job runner ----
function run_parallel {
    local MAX_JOBS=$1
    local TASK_LIST=$2
    local RUNNING=0

    while IFS= read -r CMD; do
        while [ $RUNNING -ge $MAX_JOBS ]; do
            wait -n 2>/dev/null || true
            RUNNING=$((RUNNING - 1))
        done
        eval "$CMD" &
        RUNNING=$((RUNNING + 1))
    done < "$TASK_LIST"
    wait
}

# ---- Step 1: Univariate fits ----
echo "=== Step 1: Univariate fits ==="
TASK_FILE="${OUT_DIR}/tasks_fit1.txt"
> "$TASK_FILE"

for TRAIT in "${ALL_TRAITS[@]}"; do
    for REP in $(seq 1 ${N_REP}); do
        EXTRACT="--extract ${MIXER_REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"
        OUT="${OUT_DIR}/${TRAIT}.fit.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "echo \"[$(date +%H:%M:%S)] [skip] fit1 ${TRAIT} rep${REP}\"" >> "$TASK_FILE"
        else
            echo "echo \"[$(date +%H:%M:%S)] [run]  fit1 ${TRAIT} rep${REP}\" && mixer_py fit1 ${COMMON_FLAGS} ${EXTRACT} --trait1-file \"${SUMSTATS_DIR}/${TRAIT}.sumstats.gz\" --out \"${OUT}\" && echo \"[$(date +%H:%M:%S)] [done] fit1 ${TRAIT} rep${REP}\" || echo \"[$(date +%H:%M:%S)] [FAIL] fit1 ${TRAIT} rep${REP}\"" >> "$TASK_FILE"
        fi
    done
done
run_parallel ${PARALLEL_JOBS} "$TASK_FILE"

# ---- Step 1b: Univariate tests ----
echo "=== Step 1b: Univariate tests ==="
TASK_FILE="${OUT_DIR}/tasks_test1.txt"
> "$TASK_FILE"

for TRAIT in "${ALL_TRAITS[@]}"; do
    for REP in $(seq 1 ${N_REP}); do
        OUT="${OUT_DIR}/${TRAIT}.test.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "echo \"[$(date +%H:%M:%S)] [skip] test1 ${TRAIT} rep${REP}\"" >> "$TASK_FILE"
        else
            echo "echo \"[$(date +%H:%M:%S)] [run]  test1 ${TRAIT} rep${REP}\" && mixer_py test1 ${COMMON_FLAGS} --trait1-file \"${SUMSTATS_DIR}/${TRAIT}.sumstats.gz\" --load-params \"${OUT_DIR}/${TRAIT}.fit.rep${REP}.json\" --out \"${OUT}\" && echo \"[$(date +%H:%M:%S)] [done] test1 ${TRAIT} rep${REP}\" || echo \"[$(date +%H:%M:%S)] [FAIL] test1 ${TRAIT} rep${REP}\"" >> "$TASK_FILE"
        fi
    done
done
run_parallel ${PARALLEL_JOBS} "$TASK_FILE"

# ---- Step 2: Bivariate fits ----
echo "=== Step 2: Bivariate fits ==="
TASK_FILE="${OUT_DIR}/tasks_fit2.txt"
> "$TASK_FILE"

for TRAIT2 in "${TRAITS[@]}"; do
    LABEL="ASD_vs_${TRAIT2}"
    for REP in $(seq 1 ${N_REP}); do
        EXTRACT="--extract ${MIXER_REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"
        OUT="${OUT_DIR}/${LABEL}.fit.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "echo \"[$(date +%H:%M:%S)] [skip] fit2 ${LABEL} rep${REP}\"" >> "$TASK_FILE"
        else
            echo "echo \"[$(date +%H:%M:%S)] [run]  fit2 ${LABEL} rep${REP}\" && mixer_py fit2 ${COMMON_FLAGS} ${EXTRACT} --trait1-file \"${SUMSTATS_DIR}/TRAIT_B.sumstats.gz\" --trait2-file \"${SUMSTATS_DIR}/${TRAIT2}.sumstats.gz\" --trait1-params \"${OUT_DIR}/ASD.fit.rep${REP}.json\" --trait2-params \"${OUT_DIR}/${TRAIT2}.fit.rep${REP}.json\" --out \"${OUT}\" && echo \"[$(date +%H:%M:%S)] [done] fit2 ${LABEL} rep${REP}\" || echo \"[$(date +%H:%M:%S)] [FAIL] fit2 ${LABEL} rep${REP}\"" >> "$TASK_FILE"
        fi
    done
done
run_parallel ${PARALLEL_JOBS} "$TASK_FILE"

# ---- Step 2b: Bivariate tests ----
echo "=== Step 2b: Bivariate tests ==="
TASK_FILE="${OUT_DIR}/tasks_test2.txt"
> "$TASK_FILE"

for TRAIT2 in "${TRAITS[@]}"; do
    LABEL="ASD_vs_${TRAIT2}"
    for REP in $(seq 1 ${N_REP}); do
        OUT="${OUT_DIR}/${LABEL}.test.rep${REP}"
        if [ -f "${OUT}.json" ]; then
            echo "echo \"[$(date +%H:%M:%S)] [skip] test2 ${LABEL} rep${REP}\"" >> "$TASK_FILE"
        else
            echo "echo \"[$(date +%H:%M:%S)] [run]  test2 ${LABEL} rep${REP}\" && mixer_py test2 ${COMMON_FLAGS} --trait1-file \"${SUMSTATS_DIR}/TRAIT_B.sumstats.gz\" --trait2-file \"${SUMSTATS_DIR}/${TRAIT2}.sumstats.gz\" --load-params \"${OUT_DIR}/${LABEL}.fit.rep${REP}.json\" --out \"${OUT}\" && echo \"[$(date +%H:%M:%S)] [done] test2 ${LABEL} rep${REP}\" || echo \"[$(date +%H:%M:%S)] [FAIL] test2 ${LABEL} rep${REP}\"" >> "$TASK_FILE"
        fi
    done
done
run_parallel ${PARALLEL_JOBS} "$TASK_FILE"

# ---- Step 3: Combine replicates ----
echo "=== Step 3: Combining replicates ==="

for TRAIT in TRAIT_B "${TRAITS[@]}"; do
    echo "  combining ${TRAIT}..."
    mixer_figs_py combine --json "${OUT_DIR}/${TRAIT}.fit.rep@.json" --out "${OUT_DIR}/${TRAIT}.fit"
    mixer_figs_py combine --json "${OUT_DIR}/${TRAIT}.test.rep@.json" --out "${OUT_DIR}/${TRAIT}.test"
done

for TRAIT2 in "${TRAITS[@]}"; do
    LABEL="ASD_vs_${TRAIT2}"
    echo "  combining ${LABEL}..."
    mixer_figs_py combine --json "${OUT_DIR}/${LABEL}.fit.rep@.json" --out "${OUT_DIR}/${LABEL}.fit"
    mixer_figs_py combine --json "${OUT_DIR}/${LABEL}.test.rep@.json" --out "${OUT_DIR}/${LABEL}.test"
done

echo "[$(date +%H:%M:%S)] All done. Results in ${OUT_DIR}"
echo "Next: run 12c_analyze_mixer.R to extract and visualize estimates."