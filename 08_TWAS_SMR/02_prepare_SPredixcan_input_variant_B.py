#!/usr/bin/env python3
"""
05a_prepare_spredixcan_Anxiety.py

input：
  disease_anxiety：    gwas_anxiety_fmt.txt columns：SNP, CHR, position, A1, A2, BETA, SE, EAF, P, N
"""

import pandas as pd
import numpy as np
import os, sys, time
from pyliftover import LiftOver

DATA_DIR = "${PROJECT_ROOT}"
BIM_FILE = "${PROJECT_ROOT}"
OUT_DIR  = "${PROJECT_ROOT}"
os.makedirs(OUT_DIR, exist_ok=True)

IBS_N = 486601
ANX_N = 418399

print(msg liftover（hg19 → hg38）...)
_lo = LiftOver('hg19', 'hg38')

def build_pos_map(df_with_chr_bp):
    unique_pos = df_with_chr_bp[['CHR', 'BP']].drop_duplicates()
    n = len(unique_pos)
    print(f  msg liftover msg：{n:,})
    pos_map = {}
    t0 = time.time()
    for i, row in enumerate(unique_pos.itertuples(index=False), 1):
        if i % 200000 == 0:
            elapsed = time.time() - t0
            eta = elapsed / i * (n - i)
            print(f    {i:,}/{n:,}  ({100*i/n:.1f}%)  msg {elapsed:.0f}s  msg {eta:.0f}s)
        try:
            res = _lo.convert_coordinate(f'chr{int(row.CHR)}', int(row.BP))
            if res:
                pos_map[(int(row.CHR), int(row.BP))] = int(res[0][1])
        except Exception:
            pass
    n_ok = len(pos_map)
    print(f  liftover msg：{n_ok:,}/{n:,} ({100*n_ok/n:.1f}%))
    return pos_map

def apply_pos_map(df, pos_map):
    keys  = list(zip(df['CHR'].astype(int), df['BP'].astype(int)))
    bp38  = [pos_map.get(k) for k in keys]
    df    = df.copy()
    df['BP_b38'] = bp38
    df = df[pd.notnull(df['BP_b38'])].copy()
    df['BP_b38'] = df['BP_b38'].astype(int)
    df['variant_id'] = (
        'chr' + df['CHR'].astype(str) + '_' +
        df['BP_b38'].astype(str) + '_' +
        df['A2'].str.upper() + '_' + df['A1'].str.upper() + '_b38'
    )
    return df

def save_spredixcan(df, outname):
    out_cols = ["variant_id", "chromosome", "position",
                "effect_allele", "non_effect_allele",
                "zscore", "beta", "standard_error", "pvalue", "sample_size"]
    if "frequency" in df.columns:
        out_cols.insert(5, "frequency")
    df = df[out_cols].drop_duplicates(subset=["variant_id"])
    outpath = os.path.join(OUT_DIR, outname)
    df.to_csv(outpath, sep="\t", index=False, compression="gzip")
    print(f  msg → {outpath}  ({len(df):,} SNPs))

print(\nmsg 1000G BIM（hg19 msg）...)
bim = pd.read_csv(BIM_FILE, sep="\t", header=None,
                  names=["CHR_bim","SNP","CM","BP","A1_bim","A2_bim"])
bim = bim[["SNP","CHR_bim","BP"]].drop_duplicates(subset="SNP")
bim['CHR_bim'] = pd.to_numeric(bim['CHR_bim'], errors='coerce')
print(f"  BIM：{len(bim):,} SNPs")

# ============================================================
# ============================================================
print("\n" + "="*60)
print(msg IBS_GWAS...)

ibs = pd.read_csv(os.path.join(DATA_DIR, "gwas_TRAIT_A.txt"), sep="\t", low_memory=False)
ibs['CHR'] = pd.to_numeric(ibs['CHR'], errors='coerce')
print(f  msg：{len(ibs):,})

ibs = ibs.merge(bim, on="SNP", how="inner")
ibs = ibs[ibs['CHR'] == ibs['CHR_bim']].copy()
print(f  BIM msg：{len(ibs):,} SNPs)
if len(ibs) == 0:
    sys.exit([msg] TRAIT_A × BIM msg 0)

ibs['beta']           = pd.to_numeric(ibs['BETA'], errors='coerce')
ibs['standard_error'] = pd.to_numeric(ibs['SE'],   errors='coerce')
ibs = ibs[(ibs['standard_error'] > 0) & ibs['standard_error'].notna()]
ibs['zscore'] = ibs['beta'] / ibs['standard_error']
ibs = ibs[np.isfinite(ibs['zscore'])]
try:
    from scipy import stats as _stats
    ibs['pvalue'] = 2 * _stats.norm.sf(ibs['zscore'].abs())
except ImportError:
    ibs['pvalue'] = 2 * np.exp(-0.717 * ibs['zscore'].abs() - 0.416 * ibs['zscore']**2)
ibs['sample_size'] = IBS_N
print(f  msg SNPs：{len(ibs):,})

print(  msg hg19→hg38 msg...)
ibs['CHR'] = ibs['CHR'].astype(int)
ibs['BP']  = ibs['BP'].astype(int)
pos_map_ibs = build_pos_map(ibs[['CHR','BP']])
ibs = apply_pos_map(ibs, pos_map_ibs)
print(f  liftover msg：{len(ibs):,} SNPs)

ibs = ibs.rename(columns={'CHR':'chromosome', 'BP':'position',
                           'A1':'effect_allele', 'A2':'non_effect_allele'})
save_spredixcan(ibs, "IBS_GWAS_spredixcan.txt.gz")

# ============================================================
# ============================================================
print("\n" + "="*60)
print("Process disease_anxiety GWAS...")

anx = pd.read_csv(os.path.join(DATA_DIR, "gwas_anxiety_fmt.txt"), sep="\t", low_memory=False)
anx = anx.rename(columns={'position':'BP'})
anx['CHR'] = pd.to_numeric(anx['CHR'], errors='coerce')
anx['BP']  = pd.to_numeric(anx['BP'],  errors='coerce')
anx = anx.dropna(subset=['CHR', 'BP']).copy()
anx['CHR'] = anx['CHR'].astype(int)
anx['BP']  = anx['BP'].astype(int)
print(f  msg：{len(anx):,})

anx['beta']           = pd.to_numeric(anx['BETA'], errors='coerce')
anx['standard_error'] = pd.to_numeric(anx['SE'],   errors='coerce')
anx = anx[(anx['standard_error'] > 0) & anx['standard_error'].notna()]
anx['zscore'] = anx['beta'] / anx['standard_error']
anx = anx[np.isfinite(anx['zscore'])]
anx['pvalue'] = pd.to_numeric(anx['P'], errors='coerce')
anx = anx[anx['pvalue'].notna() & (anx['pvalue'] > 0)]
anx['sample_size'] = ANX_N

if 'EAF' in anx.columns:
    anx = anx.rename(columns={'EAF':'frequency'})

print(f  msg SNPs：{len(anx):,})

print(  msg hg19→hg38 msg...)
pos_map_anx = build_pos_map(anx[['CHR','BP']])
anx = apply_pos_map(anx, pos_map_anx)
print(f  liftover msg：{len(anx):,} SNPs)

anx = anx.rename(columns={'CHR':'chromosome', 'BP':'position',
                           'A1':'effect_allele', 'A2':'non_effect_allele'})
save_spredixcan(anx, "Anxiety_spredixcan.txt.gz")

print(\nmsg。msg：bash 05b_run_spredixcan_Anxiety.sh)
