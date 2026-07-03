#!/usr/bin/env python3
"""
Convert GWAS sumstats to pleiotropyFDR .mat format.

pleiotropyFDR requires MATLAB .mat files containing:
  logpvec  — vector of -log10(p), length = n_ref_snps (9,545,380)
  zvec     — vector of Z scores, same length

SNPs not in the reference are set to NaN.
Reference file: 9545380.ref  (tab-sep: SNP CHR BP A1 A2)

Prerequisites:
  pip install scipy numpy pandas
  Download reference:
    wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/9545380.ref

Usage:
  python 13a_format_sumstats_pleiofdr.py
"""

import os, sys, glob
import numpy as np
import pandas as pd
from scipy.io import savemat

# ─── Paths ────────────────────────────────────────────────────────────────────
IBD_DIR  = "${PROJECT_ROOT}"
REF_FILE = "${PROJECT_ROOT}"
OUT_DIR  = "${PROJECT_ROOT}"
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Load reference ───────────────────────────────────────────────────────────
if not os.path.exists(REF_FILE):
    sys.exit(f"Reference file not found: {REF_FILE}")

print(f"Loading reference ({REF_FILE})...")
ref = pd.read_csv(REF_FILE, sep="\t",
                  names=["SNP", "CHR", "BP", "A1", "A2"],
                  dtype={"SNP": str, "CHR": str, "BP": str, "A1": str, "A2": str})
N_REF = len(ref)
print(f"  {N_REF:,} reference SNPs")

snp2idx = {snp: i for i, snp in enumerate(ref["SNP"])}

AMBIGUOUS = {frozenset(["A", "T"]), frozenset(["C", "G"])}

def is_ambiguous(a1, a2):
    return frozenset([str(a1).upper(), str(a2).upper()]) in AMBIGUOUS

# ─── Core conversion function ─────────────────────────────────────────────────
def sumstats_to_mat(label, df, snp_col="SNP", beta_col="BETA", se_col="SE",
                    p_col="P", a1_col="A1", a2_col="A2"):
    """Convert a sumstats DataFrame to pleiotropyFDR .mat file."""
    logpvec = np.full(N_REF, np.nan)
    zvec    = np.full(N_REF, np.nan)

    matched = missed = ambig = 0
    for _, row in df.iterrows():
        snp = str(row[snp_col])
        if snp not in snp2idx:
            missed += 1
            continue
        a1 = str(row[a1_col]).upper()
        a2 = str(row[a2_col]).upper()
        if is_ambiguous(a1, a2):
            ambig += 1
            continue
        try:
            p_val = float(row[p_col])
            beta  = float(row[beta_col])
            se    = float(row[se_col])
        except (ValueError, TypeError):
            continue
        if p_val <= 0 or np.isnan(p_val):
            continue
        idx = snp2idx[snp]
        logpvec[idx] = -np.log10(p_val)
        zvec[idx]    = beta / se
        matched += 1

    print(f"  {label}: matched={matched:,}  missed={missed:,}  ambig={ambig:,}")

    # Check allele alignment with reference
    ref_a1 = ref["A1"].str.upper().values
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    # Flip Z for SNPs where alleles are reversed
    for _, row in df.iterrows():
        snp = str(row[snp_col])
        if snp not in snp2idx:
            continue
        idx = snp2idx[snp]
        a1 = str(row[a1_col]).upper()
        if a1 == complement.get(ref_a1[idx], ""):
            if not np.isnan(zvec[idx]):
                zvec[idx] = -zvec[idx]   # flip sign

    out = os.path.join(OUT_DIR, f"{label}.mat")
    savemat(out, {"logpvec": logpvec, "zvec": zvec})
    print(f"  → saved: {out}")

# ─── TRAIT_B ──────────────────────────────────────────────────────────────────────
print("\n[ASD]")
asd = pd.read_csv(f"{IBD_DIR}/ieu-a-1185_extracted_ASD.txt", sep="\t")
sumstats_to_mat("TRAIT_B", asd, p_col="P")

# ─── de Lange ─────────────────────────────────────────────────────────────────
for dset in ["TRAIT_C", "CD", "UC"]:
    files = glob.glob(f"{IBD_DIR}/delange/*{dset}*.txt")
    if not files:
        print(f"[de Lange {dset}] WARNING: file not found, skipping.")
        continue
    print(f"\n[de Lange {dset}]")
    df = pd.read_csv(files[0], sep="\t")
    sumstats_to_mat(f"de Lange_{dset}", df)

# ─── Liu ──────────────────────────────────────────────────────────────────────
for dset in ["TRAIT_C", "CD", "UC"]:
    files = glob.glob(f"{IBD_DIR}/liu/*{dset}*.txt")
    if not files:
        print(f"[Liu {dset}] WARNING: file not found, skipping.")
        continue
    print(f"\n[Liu {dset}]")
    df = pd.read_csv(files[0], sep="\t")
    sumstats_to_mat(f"Liu_{dset}", df)

# ─── TRAIT_A ──────────────────────────────────────────────────────────────────────
IBS_CANDIDATES = [
    "${PROJECT_ROOT}",
    "${PROJECT_ROOT}",
]
for ibs_f in IBS_CANDIDATES:
    if os.path.exists(ibs_f):
        print(f"\n[IBS] {ibs_f}")
        ibs = pd.read_csv(ibs_f, sep="\t")
        ibs.columns = [c.upper() for c in ibs.columns]
        # Compute P from BETA and SE if not present
        if "P" not in ibs.columns:
            from scipy import stats
            z = ibs["BETA"].astype(float) / ibs["SE"].astype(float)
            ibs["P"] = 2 * stats.norm.sf(np.abs(z))
        sumstats_to_mat("TRAIT_A", ibs)
        break
else:
    print("\n[IBS] WARNING: TRAIT_A GWAS not found; skipping.")

print(f"\nAll .mat files in: {OUT_DIR}")
print("Next: run 13b_run_pleiofdr.m in MATLAB")
