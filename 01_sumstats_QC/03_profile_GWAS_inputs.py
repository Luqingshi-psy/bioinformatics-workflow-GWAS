#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path

from analysis_config import RESULTS_DIR, TRAITS
from sumstats_utils import open_text, parse_vcf_sample, safe_float, write_json


OUT_DIR = RESULTS_DIR / "input_profile"


def profile_vcf(path, max_rows=None):
    counters = Counter()
    format_keys_seen = Counter()
    sample_name = None
    header = None
    first_variants = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                counters["meta_lines"] += 1
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                sample_name = header[-1]
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            counters["rows"] += 1
            if max_rows is not None and counters["rows"] > max_rows:
                counters["stopped_at_max_rows"] = max_rows
                break
            chrom, pos, rsid, ref, alt = fields[:5]
            counters["snp_rows" if len(ref) == 1 and len(alt) == 1 else "non_snp_rows"] += 1
            if rsid == ".":
                counters["missing_rsid"] += 1
            format_keys_seen[fields[8]] += 1
            sample = parse_vcf_sample(fields[8], fields[9])
            for key in ("ES", "SE", "LP", "SS", "NC", "AF", "EZ"):
                if key in sample and safe_float(sample.get(key)) == safe_float(sample.get(key)):
                    counters[f"has_{key}"] += 1
            if len(first_variants) < 5:
                first_variants.append(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "rsid": rsid,
                        "ref": ref,
                        "alt": alt,
                        "format": fields[8],
                        "sample": fields[9],
                    }
                )

    return {
        "path": str(path),
        "file_type": "vcf.gz",
        "header": header,
        "sample_name": sample_name,
        "counts": dict(counters),
        "format_keys_top": format_keys_seen.most_common(10),
        "first_variants": first_variants,
        "effect_allele_assumption": "ES is relative to ALT; standardized EA=ALT and OA=REF.",
    }


def profile_table(path, max_rows=None):
    counters = Counter()
    first_rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        for row in reader:
            counters["rows"] += 1
            if max_rows is not None and counters["rows"] > max_rows:
                counters["stopped_at_max_rows"] = max_rows
                break
            if row.get("SNP") in {"", ".", None}:
                counters["missing_rsid"] += 1
            if row.get("A1") and row.get("A2") and len(row["A1"]) == 1 and len(row["A2"]) == 1:
                counters["snp_rows"] += 1
            else:
                counters["non_snp_rows"] += 1
            for key in ("OR", "SE", "P", "Neff_half", "Nca", "Nco", "INFO"):
                if safe_float(row.get(key)) == safe_float(row.get(key)):
                    counters[f"has_{key}"] += 1
            if len(first_rows) < 5:
                first_rows.append(row)
    return {
        "path": str(path),
        "file_type": "tabular",
        "header": header,
        "counts": dict(counters),
        "first_rows": first_rows,
        "effect_allele_assumption": "OR is treated as relative to A1; standardized EA=A1 and OA=A2.",
    }


def main():
    parser = argparse.ArgumentParser(description="Profile GWAS input files.")
    parser.add_argument("--max-rows", type=int, default=None, help="Optional per-file row cap for fast QC.")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []
    for trait, path in TRAITS.items():
        if not Path(path).exists():
            summary = {"path": str(path), "error": "missing file"}
        elif str(path).endswith(".vcf.gz"):
            summary = profile_vcf(path, max_rows=args.max_rows)
        else:
            summary = profile_table(path, max_rows=args.max_rows)
        summary["trait"] = trait
        write_json(OUT_DIR / f"{trait}.profile.json", summary)
        counts = summary.get("counts", {})
        summaries.append(
            {
                "trait": trait,
                "path": str(path),
                "file_type": summary.get("file_type", "missing"),
                "rows": counts.get("rows", 0),
                "snp_rows": counts.get("snp_rows", 0),
                "non_snp_rows": counts.get("non_snp_rows", 0),
                "missing_rsid": counts.get("missing_rsid", 0),
                "effect_allele_assumption": summary.get("effect_allele_assumption", ""),
            }
        )

    out_csv = OUT_DIR / "input_profile_summary.csv"
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
