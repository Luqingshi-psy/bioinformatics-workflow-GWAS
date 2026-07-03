#!/usr/bin/env python3
"""
Equivalent of 12c_analyze_mixer.R
Parses combined MiXeR bivariate results, prints summary, saves CSV + 2 PDF plots.

Usage (from project root):
    docker run --rm -v /mnt/c/Users/Administrator/Desktop/trait pair (A x B):/home -w /home \
        ghcr.io/precimed/gsa-mixer:2.2.1 \
        python /home/analysis/12c_analyze_mixer.py
"""
import json
import os
import sys
import csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Resolve project root from script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.dirname(SCRIPT_DIR)
MIXER_DIR = os.path.join(PROJ_DIR, "output", "MiXeR", "results")
OUT_DIR = os.path.join(PROJ_DIR, "output", "MiXeR")

PAIRS = [
    ("TRAIT_B", "TRAIT_A",         "TRAIT_B vs TRAIT_A"),
    ("TRAIT_B", "de Lange_IBD", "TRAIT_B vs TRAIT_C"),
    ("TRAIT_B", "de Lange_CD",  "TRAIT_B vs TRAIT_E"),
    ("TRAIT_B", "de Lange_UC",  "TRAIT_B vs TRAIT_D"),
    ("TRAIT_B", "Liu_IBD",     "TRAIT_B vs TRAIT_C_Liu"),
    ("TRAIT_B", "Liu_CD",      "TRAIT_B vs TRAIT_E_Liu"),
    ("TRAIT_B", "Liu_UC",      "TRAIT_B vs TRAIT_D_Liu"),
]


def parse_bivar_json(json_file: str, label: str):
    if not os.path.exists(json_file):
        print(f"  [missing] {json_file}")
        return None
    with open(json_file) as f:
        j = json.load(f)
    if "ci" not in j:
        print(f"  [no 'ci'] {json_file}")
        return None

    def g(field):
        c = j.get("ci", {}).get(field, {})
        return {
            "mean":   c.get("mean", np.nan),
            "std":    c.get("std", np.nan),
            "median": c.get("median", np.nan),
        }

    return {
        "label":      label,
        "nc1":        g("nc1")["mean"],
        "nc1_se":     g("nc1")["std"],
        "nc2":        g("nc2")["mean"],
        "nc2_se":     g("nc2")["std"],
        "nc12":       g("nc12")["mean"],
        "nc12_se":    g("nc12")["std"],
        "dice":       g("dice")["mean"],
        "dice_se":    g("dice")["std"],
        "rho_ge":     g("rho_beta")["mean"],   # genetic-effect correlation
        "rho_ge_se":  g("rho_beta")["std"],
        "rg":         g("rg")["mean"],
        "rg_se":      g("rg")["std"],
    }


def main():
    print(f"Project root: {PROJ_DIR}")
    print(f"MiXeR dir   : {MIXER_DIR}")
    print(f"Output dir  : {OUT_DIR}\n")

    rows = []
    for t1, t2, lab in PAIRS:
        jf = os.path.join(MIXER_DIR, f"{t1}_vs_{t2}.fit.json")
        r = parse_bivar_json(jf, lab)
        if r is not None:
            rows.append(r)
        else:
            print(f"  [skip] {lab}: file not found")

    df = pd.DataFrame(rows)

    print("\n=== MiXeR Bivariate Results (combined across 20 reps) ===")
    show_cols = ["label", "nc12", "nc12_se", "dice", "dice_se",
                 "rg", "rg_se", "rho_ge", "rho_ge_se"]
    print(df[show_cols].to_string(index=False))

    # ── Save summary CSV ────────────────────────────────────────────────────
    summary_csv = os.path.join(OUT_DIR, "mixer_summary.csv")
    df.to_csv(summary_csv, index=False)
    print(f"\nSaved: {summary_csv}")

    if df.empty:
        print("No data; exiting.")
        return

    # ── Dice coefficient bar chart ─────────────────────────────────────────
    df_sorted = df.sort_values("dice", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    y = np.arange(len(df_sorted))
    ax.barh(y, df_sorted["dice"], xerr=df_sorted["dice_se"],
            color="#4E79A7", edgecolor="black", alpha=0.85, capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df_sorted["label"])
    ax.set_xlabel("Dice coefficient (polygenic overlap)")
    ax.set_title("ASD × Gut Diseases — Shared Genetic Architecture (MiXeR)")
    ax.axvline(0, color="black", linewidth=0.6)
    plt.tight_layout()
    out_pdf1 = os.path.join(OUT_DIR, "mixer_dice_barplot.pdf")
    fig.savefig(out_pdf1)
    plt.close(fig)
    print(f"Saved: {out_pdf1}")

    # ── Venn-like summary CSV ──────────────────────────────────────────────
    # MiXeR semantics: nc1 = trait1 unique, nc2 = trait2 unique, nc12 = shared
    # (verified: nc1 + nc12 = univariate total for trait1, etc.)
    venn = pd.DataFrame({
        "label":      df["label"],
        "unique_t1":  df["nc1"],                # TRAIT_B unique causal variants (×1000)
        "shared":     df["nc12"],               # shared causal variants (×1000)
        "unique_t2":  df["nc2"],                # gut-disease unique (×1000)
        "nc12":       df["nc12"],
        "nc1u":       df["nc1"] + df["nc12"],   # TRAIT_B total (univariate)
        "nc2u":       df["nc2"] + df["nc12"],   # gut total (univariate)
        "dice":       df["dice"],
    })
    print("\n=== Shared vs unique causal variants (in thousands of variants) ===")
    print(venn.to_string(index=False))
    venn_csv = os.path.join(OUT_DIR, "mixer_venn_summary.csv")
    venn.to_csv(venn_csv, index=False)
    print(f"Saved: {venn_csv}")

    # ── Stacked bar: proportion shared ─────────────────────────────────────
    # Use the verified-unique values directly
    components = pd.DataFrame({
        "label":         df["label"],
        "ASD only":      df["nc1"],
        "Shared":        df["nc12"],
        "Gut only":      df["nc2"],
    }).set_index("label")

    fig, ax = plt.subplots(figsize=(8, 5))
    components.plot.barh(
        stacked=True, ax=ax,
        color=["#4E79A7", "#F28E2B", "#59A14F"],
        edgecolor="black", width=0.7
    )
    ax.set_xlabel("Proportion of causal variants")
    ax.set_ylabel("")
    ax.set_title("Polygenic overlap composition")
    ax.legend(loc="lower right", framealpha=0.9)
    plt.tight_layout()
    out_pdf2 = os.path.join(OUT_DIR, "mixer_overlap_proportion.pdf")
    fig.savefig(out_pdf2)
    plt.close(fig)
    print(f"Saved: {out_pdf2}")


if __name__ == "__main__":
    main()
