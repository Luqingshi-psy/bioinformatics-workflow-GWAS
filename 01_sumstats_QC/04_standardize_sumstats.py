#!/usr/bin/env python3
import argparse
import csv
import math
from collections import Counter

from analysis_config import MAX_ABS_Z, RESULTS_DIR, SNP_ONLY, TRAITS
from sumstats_utils import (
    dict_writer,
    finite,
    is_snp,
    open_text,
    p_from_lp,
    parse_vcf_sample,
    safe_float,
    z_from_beta_se,
)


OUT_DIR = RESULTS_DIR / "standardized"
FIELDS = ["SNP", "CHR", "BP", "EA", "OA", "BETA", "SE", "Z", "P", "N", "SOURCE"]


def valid_record(record, counters):
    counters["parsed"] += 1
    if not record["SNP"] or record["SNP"] == ".":
        counters["drop_missing_snp"] += 1
        return False
    if SNP_ONLY and not is_snp(record["EA"], record["OA"]):
        counters["drop_non_snp"] += 1
        return False
    if not finite(record["BETA"]) or not finite(record["SE"]) or record["SE"] <= 0:
        counters["drop_bad_beta_se"] += 1
        return False
    if not finite(record["Z"]):
        counters["drop_bad_z"] += 1
        return False
    if abs(record["Z"]) > MAX_ABS_Z:
        counters["drop_extreme_z"] += 1
        return False
    if not finite(record["P"]) or record["P"] <= 0 or record["P"] > 1:
        counters["drop_bad_p"] += 1
        return False
    counters["kept"] += 1
    return True


def iter_vcf_records(trait, path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, rsid, ref, alt = fields[:5]
            sample = parse_vcf_sample(fields[8], fields[9])
            beta = safe_float(sample.get("ES"))
            se = safe_float(sample.get("SE"))
            z = safe_float(sample.get("EZ"))
            if not finite(z):
                z = z_from_beta_se(beta, se)
            p = p_from_lp(sample.get("LP"))
            n = safe_float(sample.get("SS"))
            if not finite(n):
                n = safe_float(sample.get("NC"))
            yield {
                "SNP": rsid,
                "CHR": chrom,
                "BP": safe_float(pos),
                "EA": alt.upper(),
                "OA": ref.upper(),
                "BETA": beta,
                "SE": se,
                "Z": z,
                "P": p,
                "N": n,
                "SOURCE": trait,
            }


def iter_anxiety_records(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            odds_ratio = safe_float(row.get("OR"))
            beta = math.log(odds_ratio) if finite(odds_ratio) and odds_ratio > 0 else math.nan
            se = safe_float(row.get("SE"))
            yield {
                "SNP": row.get("SNP"),
                "CHR": row.get("CHR"),
                "BP": safe_float(row.get("BP")),
                "EA": str(row.get("A1", "")).upper(),
                "OA": str(row.get("A2", "")).upper(),
                "BETA": beta,
                "SE": se,
                "Z": z_from_beta_se(beta, se),
                "P": safe_float(row.get("P")),
                "N": safe_float(row.get("Neff_half")),
                "SOURCE": "disease_anxiety",
            }


def standardize_trait(trait, path, limit=None):
    counters = Counter()
    out_path = OUT_DIR / f"{trait}.standardized.tsv.gz"
    handle, writer = dict_writer(out_path, FIELDS, gz=True)
    iterator = iter_vcf_records(trait, path) if str(path).endswith(".vcf.gz") else iter_anxiety_records(path)
    with handle:
        for record in iterator:
            if limit is not None and counters["parsed"] >= limit:
                counters["stopped_at_limit"] = limit
                break
            if valid_record(record, counters):
                writer.writerow(record)
    counters["trait"] = trait
    counters["output"] = str(out_path)
    return dict(counters)


def main():
    parser = argparse.ArgumentParser(description="Standardize GWAS summary statistics.")
    parser.add_argument("--limit", type=int, default=None, help="Optional per-trait parsed-row cap for smoke tests.")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = [standardize_trait(trait, path, limit=args.limit) for trait, path in TRAITS.items()]
    fieldnames = sorted({key for row in summaries for key in row})
    summary_path = OUT_DIR / "standardization_qc_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
