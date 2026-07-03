#!/usr/bin/env python3
"""
Format GWAS sumstats for MiXeR input.
Required output columns: SNP  CHR  BP  A1  A2  N  Z
- Ambiguous (A/T, C/G) SNPs are removed.
- MHC region (chr6:26-34 Mb) is removed.
- Z is computed as BETA/SE (sign kept).
"""

import pandas as pd
import numpy as np
import os, sys

IBD_DIR = "${PROJECT_ROOT}"
OUT_DIR = "${PROJECT_ROOT}"
os.makedirs(OUT_DIR, exist_ok=True)

MHC = (6, 26_000_000, 34_000_000)

AMBIGUOUS = {frozenset(["A", "T"]), frozenset(["C", "G"])}

def is_ambiguous(a1, a2):
    return frozenset([a1.upper(), a2.upper()]) in AMBIGUOUS

def in_mhc(chrom, pos):
    try:
        return int(chrom) == MHC[0] and MHC[1] <= int(pos) <= MHC[2]
    except Exception:
        return False

def compute_z(df):
    """Compute Z from BETA and SE."""
    return df["BETA"].astype(float) / df["SE"].astype(float)

def save_mixer(df, label):
    out = os.path.join(OUT_DIR, f"{label}.sumstats.gz")
    df[["SNP", "CHR", "BP", "A1", "A2", "N", "Z"]].to_csv(out, sep="\t", index=False)
    print(f"  → {out}  ({len(df):,} SNPs)")

def filter_common(df, chr_col, pos_col, a1_col, a2_col):
    # Drop ambiguous
    amb = df.apply(lambda r: is_ambiguous(r[a1_col], r[a2_col]), axis=1)
    df = df[~amb].copy()
    # Drop MHC
    mhc_mask = df.apply(lambda r: in_mhc(r[chr_col], r[pos_col]), axis=1)
    df = df[~mhc_mask].copy()
    # Upper-case alleles
    df[a1_col] = df[a1_col].str.upper()
    df[a2_col] = df[a2_col].str.upper()
    return df

# ─── TRAIT_B ─────────────────────────────────────────────────────────────────────
print("[ASD]")
asd = pd.read_csv(f"{IBD_DIR}/ieu-a-1185_extracted_ASD.txt", sep="\t", dtype=str)
# cols: SNP CHR position A1 A2 BETA SE EAF P N
asd = asd.rename(columns={"position": "BP"})
asd["Z"] = compute_z(asd)
asd["N"] = asd["N"].astype(float)
asd = filter_common(asd, "CHR", "BP", "A1", "A2")
save_mixer(asd, "TRAIT_B")

# ─── de Lange datasets ───────────────────────────────────────────────────────
# cols: SNP A1 A2 BETA SE EAF P N  (no CHR/BP → need to add from TRAIT_B or ref)
# de Lange files lack CHR and BP — we merge on SNP from TRAIT_B to get coordinates.
asd_coords = asd[["SNP", "CHR", "BP"]].drop_duplicates("SNP")

for dset in ["TRAIT_C", "CD", "UC"]:
    print(f"[de Lange {dset}]")
    # find file
    import glob
    files = glob.glob(f"{IBD_DIR}/delange/*{dset}*.txt")
    if not files:
        print(f"  WARNING: no file found for de Lange {dset}, skipping.")
        continue
    dl = pd.read_csv(files[0], sep="\t", dtype=str)
    # cols: SNP A1 A2 BETA SE EAF P N
    dl = dl.merge(asd_coords, on="SNP", how="inner")
    dl["Z"] = compute_z(dl)
    dl["N"] = dl["N"].astype(float)
    dl = filter_common(dl, "CHR", "BP", "A1", "A2")
    save_mixer(dl, f"de Lange_{dset}")

# ─── Liu datasets ────────────────────────────────────────────────────────────
# cols: SNP CHR position A1 A2 BETA SE EAF P N
for dset in ["TRAIT_C", "CD", "UC"]:
    print(f"[Liu {dset}]")
    files = glob.glob(f"{IBD_DIR}/liu/*{dset}*.txt")
    if not files:
        print(f"  WARNING: no file found for Liu {dset}, skipping.")
        continue
    liu = pd.read_csv(files[0], sep="\t", dtype=str)
    liu = liu.rename(columns={"position": "BP"})
    liu["Z"] = compute_z(liu)
    liu["N"] = liu["N"].astype(float)
    liu = filter_common(liu, "CHR", "BP", "A1", "A2")
    save_mixer(liu, f"Liu_{dset}")

# ─── TRAIT_A ─────────────────────────────────────────────────────────────────────
# TRAIT_A GWAS: CHR SNP A1 A2 BETA SE N (no EAF, no P, no BP)
# BP needs to come from TRAIT_B coords (merge on SNP).
# Locate TRAIT_A file
IBS_FILES = [
    "${PROJECT_ROOT}",
    "${PROJECT_ROOT}",
    "${PROJECT_ROOT}",
]
ibs_file = None
for f in IBS_FILES:
    if os.path.exists(f):
        ibs_file = f
        break

if ibs_file:
    print(f"[IBS] using {ibs_file}")
    ibs = pd.read_csv(ibs_file, sep="\t", dtype=str)
    # Standardise column names
    col_map = {c.upper(): c for c in ibs.columns}
    ibs.columns = [c.upper() for c in ibs.columns]
    if "POSITION" in ibs.columns:
        ibs = ibs.rename(columns={"POSITION": "BP"})
    if "BP" not in ibs.columns:
        # TRAIT_A already has CHR; only merge BP from TRAIT_B coords to avoid CHR_x/CHR_y conflict
        ibs = ibs.merge(asd_coords[["SNP", "BP"]], on="SNP", how="inner")
    ibs["Z"] = compute_z(ibs)
    ibs["N"] = ibs["N"].astype(float)
    ibs = filter_common(ibs, "CHR", "BP", "A1", "A2")
    save_mixer(ibs, "TRAIT_A")
else:
    print("[IBS] WARNING: TRAIT_A GWAS file not found. Please place gwas_TRAIT_A.txt in the data directory.")

print("\nAll done. Formatted files in:", OUT_DIR)
print("Next step: run 12b_run_mixer.sh (requires MiXeR reference panel).")
