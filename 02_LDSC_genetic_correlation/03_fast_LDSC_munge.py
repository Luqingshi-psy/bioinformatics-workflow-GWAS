#!/usr/bin/env python3
import csv
import gzip
import math
from collections import Counter
from pathlib import Path


STANDARD_DIR = Path("${PROJECT_ROOT}")
HM3_PATH = Path("${PROJECT_ROOT}")
OUT_DIR = Path("${PROJECT_ROOT}")

TRAITS = [
    "TRAIT_B",
    "TRAIT_A",
    "ADHD",
    "SCZ",
    "disease_anxiety",
    "IBD_deLange",
    "CD_deLange",
    "UC_deLange",
]

N_FALLBACK = {
    "TRAIT_A": 486601.0,
    "IBD_deLange": 58331.0,
    "CD_deLange": 58331.0,
    "UC_deLange": 58331.0,
}

AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def load_hm3():
    snps = {}
    with HM3_PATH.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            snp = row["SNP"]
            a1 = row["A1"].upper()
            a2 = row["A2"].upper()
            if is_valid_allele(a1, a2):
                snps[snp] = (a1, a2)
    return snps


def is_valid_allele(a1, a2):
    return (
        len(a1) == 1
        and len(a2) == 1
        and a1 in "ACGT"
        and a2 in "ACGT"
        and (a1, a2) not in AMBIGUOUS
    )


def munge_trait(trait, hm3):
    input_path = STANDARD_DIR / f"{trait}.standardized.tsv.gz"
    output_path = OUT_DIR / f"{trait}.sumstats.gz"
    log_path = OUT_DIR / f"{trait}.fast_munge_qc.csv"
    counts = Counter()
    with gzip.open(input_path, "rt") as src, gzip.open(output_path, "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=["SNP", "A1", "A2", "Z", "N"], delimiter="\t")
        writer.writeheader()
        seen = set()
        for row in reader:
            counts["input"] += 1
            snp = row["SNP"]
            hm3_alleles = hm3.get(snp)
            if hm3_alleles is None:
                counts["drop_not_hm3"] += 1
                continue
            if snp in seen:
                counts["drop_duplicate"] += 1
                continue
            a1 = row["EA"].upper()
            a2 = row["OA"].upper()
            if not is_valid_allele(a1, a2):
                counts["drop_bad_or_ambiguous_allele"] += 1
                continue
            z = safe_float(row["Z"])
            hm3_a1, hm3_a2 = hm3_alleles
            if a1 == hm3_a1 and a2 == hm3_a2:
                out_z = z
                counts["hm3_same_orientation"] += 1
            elif a1 == hm3_a2 and a2 == hm3_a1:
                out_z = -z
                counts["hm3_swap_orientation"] += 1
            else:
                counts["drop_hm3_allele_mismatch"] += 1
                continue
            n = safe_float(row["N"])
            if (not math.isfinite(n) or n <= 0) and trait in N_FALLBACK:
                n = N_FALLBACK[trait]
                counts["used_n_fallback"] += 1
            if not math.isfinite(out_z) or not math.isfinite(n) or n <= 0:
                counts["drop_bad_z_or_n"] += 1
                continue
            writer.writerow({"SNP": snp, "A1": hm3_a1, "A2": hm3_a2, "Z": out_z, "N": n})
            seen.add(snp)
            counts["kept"] += 1
    with log_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait"] + sorted(counts.keys()))
        writer.writeheader()
        writer.writerow({"trait": trait, **counts})
    print(f"[done] {trait}: kept {counts['kept']:,} -> {output_path}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    hm3 = load_hm3()
    print(f"Loaded HapMap3 SNPs: {len(hm3):,}")
    for trait in TRAITS:
        munge_trait(trait, hm3)
    summary_path = OUT_DIR / "fast_munge_qc_summary.csv"
    rows = []
    for trait in TRAITS:
        path = OUT_DIR / f"{trait}.fast_munge_qc.csv"
        with path.open() as handle:
            rows.append(next(csv.DictReader(handle)))
    fieldnames = sorted({key for row in rows for key in row})
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
