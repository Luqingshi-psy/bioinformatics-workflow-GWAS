#!/usr/bin/env python3
"""
06_step1_prepare_smr_gwas.py

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
ASD_N = 46351

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

# ── TRAIT_B GWAS ──────────────────────────────────────────────────
print("\nProcess TRAIT_B GWAS...")
asd = pd.read_csv(os.path.join(DATA_DIR, "gwas_TRAIT_B.txt"), sep="\t", low_memory=False)
print(f  msg：{len(asd):,})

asd['beta'] = pd.to_numeric(asd['BETA'], errors='coerce')
asd['se']   = pd.to_numeric(asd['SE'],   errors='coerce')
asd['pval'] = pd.to_numeric(asd['P'],    errors='coerce')
asd = asd[(asd['se'] > 0) & asd['se'].notna() & asd['beta'].notna() & asd['pval'].notna()].copy()
asd = asd[np.isfinite(asd['beta']) & np.isfinite(asd['se'])].copy()

if 'EAF' in asd.columns:
    asd['freq'] = pd.to_numeric(asd['EAF'], errors='coerce').fillna(0.5)
else:
    asd['freq'] = 0.5

asd_out = pd.DataFrame({
    'SNP':  asd['SNP'],
    'A1':   asd['A1'].str.upper(),
    'A2':   asd['A2'].str.upper(),
    'freq': asd['freq'],
    'b':    asd['beta'],
    'se':   asd['se'],
    'p':    asd['pval'],
    'N':    ASD_N
})
asd_out = asd_out.dropna(subset=['SNP','A1','A2','b','se','p'])
asd_out = asd_out[asd_out['p'] > 0]
out_asd = os.path.join(OUT_DIR, "ASD_smr.txt")
asd_out.to_csv(out_asd, sep="\t", index=False)
print(f  msg：{out_asd}  ({len(asd_out):,} SNPs))

print(\nmsg。msg：bash 06_step2_run_smr.sh)
