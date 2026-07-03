#!/usr/bin/env python3
"""
06a_prepare_smr_Anxiety.py

SMR GWAS format：SNP  A1  A2  freq  b  se  p  N
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

DATA_DIR = "${PROJECT_ROOT}"
OUT_DIR  = "${PROJECT_ROOT}"
os.makedirs(OUT_DIR, exist_ok=True)

IBS_N = 486601
ANX_N = 418399

# ── TRAIT_A GWAS ──────────────────────────────────────────────────
print("Process TRAIT_A GWAS...")
ibs = pd.read_csv(os.path.join(DATA_DIR, "gwas_TRAIT_A.txt"), sep="\t", low_memory=False)
print(f  msg：{len(ibs):,})

ibs['beta'] = pd.to_numeric(ibs['BETA'], errors='coerce')
ibs['se']   = pd.to_numeric(ibs['SE'],   errors='coerce')
ibs = ibs[(ibs['se'] > 0) & ibs['se'].notna() & ibs['beta'].notna()].copy()

ibs['zscore'] = ibs['beta'] / ibs['se']
ibs = ibs[np.isfinite(ibs['zscore'])].copy()
ibs['pvalue'] = 2 * stats.norm.sf(np.abs(ibs['zscore']))
ibs['freq']   = 0.5

ibs_out = pd.DataFrame({
    'SNP':  ibs['SNP'],
    'A1':   ibs['A1'].str.upper(),
    'A2':   ibs['A2'].str.upper(),
    'freq': ibs['freq'],
    'b':    ibs['beta'],
    'se':   ibs['se'],
    'p':    ibs['pvalue'],
    'N':    IBS_N
})
ibs_out = ibs_out.dropna(subset=['SNP','A1','A2','b','se','p'])
ibs_out = ibs_out[ibs_out['p'] > 0]
out_ibs = os.path.join(OUT_DIR, "IBS_smr.txt")
ibs_out.to_csv(out_ibs, sep="\t", index=False)
print(f  msg：{out_ibs}  ({len(ibs_out):,} SNPs))

# ── disease_anxiety GWAS ─────────────────────────────────────────────────
print("\nProcess disease_anxiety GWAS...")
anx = pd.read_csv(os.path.join(DATA_DIR, "gwas_anxiety_fmt.txt"), sep="\t", low_memory=False)
print(f  msg：{len(anx):,})

anx['beta'] = pd.to_numeric(anx['BETA'], errors='coerce')
anx['se']   = pd.to_numeric(anx['SE'],   errors='coerce')
anx['pval'] = pd.to_numeric(anx['P'],    errors='coerce')
anx = anx[(anx['se'] > 0) & anx['se'].notna() & anx['beta'].notna() & anx['pval'].notna()].copy()
anx = anx[np.isfinite(anx['beta']) & np.isfinite(anx['se'])].copy()

if 'EAF' in anx.columns:
    anx['freq'] = pd.to_numeric(anx['EAF'], errors='coerce').fillna(0.5)
else:
    anx['freq'] = 0.5

anx_out = pd.DataFrame({
    'SNP':  anx['SNP'],
    'A1':   anx['A1'].str.upper(),
    'A2':   anx['A2'].str.upper(),
    'freq': anx['freq'],
    'b':    anx['beta'],
    'se':   anx['se'],
    'p':    anx['pval'],
    'N':    ANX_N
})
anx_out = anx_out.dropna(subset=['SNP','A1','A2','b','se','p'])
anx_out = anx_out[anx_out['p'] > 0]
out_anx = os.path.join(OUT_DIR, "Anxiety_smr.txt")
anx_out.to_csv(out_anx, sep="\t", index=False)
print(f  msg：{out_anx}  ({len(anx_out):,} SNPs))

print(\nmsg。msg：bash 06b_run_smr_Anxiety.sh)
