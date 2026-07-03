#!/usr/bin/env python3
import csv
import gzip
import math
from collections import Counter
from pathlib import Path


BASE_DIR = Path("${PROJECT_ROOT}")
DATA_DIR = Path("${PROJECT_ROOT}")
STANDARD_DIR = BASE_DIR / "results" / "standardized"
PAIR_DIR = BASE_DIR / "results" / "harmonized_asd_pairs"
OUT_DIR = BASE_DIR / "results" / "diagnostics"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def parse_sample(format_keys, sample):
    return dict(zip(format_keys.split(":"), sample.split(":")))


def summarize_products(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return {"n": 0}
    positive = sum(1 for v in finite if v > 0)
    finite_sorted = sorted(finite)
    return {
        "n": len(finite),
        "positive_fraction": positive / len(finite),
        "negative_fraction": 1 - positive / len(finite),
        "mean": sum(finite) / len(finite),
        "median": finite_sorted[len(finite_sorted) // 2],
    }


def pair_file_summary(trait, max_examples=30):
    path = PAIR_DIR / f"ASD_vs_{trait}.harmonized.tsv.gz"
    classes = Counter()
    z_products_by_class = {}
    examples_by_class = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cls = row["HARMONIZATION_CLASS"]
            classes[cls] += 1
            z_products_by_class.setdefault(cls, []).append(safe_float(row["ASD_Z"]) * safe_float(row["TRAIT_Z"]))
            examples_by_class.setdefault(cls, [])
            if len(examples_by_class[cls]) < max_examples:
                examples_by_class[cls].append(
                    {
                        "SNP": row["SNP"],
                        "CHR": row["CHR"],
                        "BP": row["BP"],
                        "EA": row["EA"],
                        "OA": row["OA"],
                        "ASD_Z": row["ASD_Z"],
                        "TRAIT_Z_after_harmonization": row["TRAIT_Z"],
                    }
                )
    z_summary = {cls: summarize_products(vals) for cls, vals in z_products_by_class.items()}
    return classes, z_summary, examples_by_class


def scz_raw_vcf_id_review(max_examples=50):
    counts = Counter()
    examples = []
    with gzip.open(DATA_DIR / "gwas_SCZ.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            sample = parse_sample(fields[8], fields[9])
            vcf_id = fields[2]
            fmt_id = sample.get("ID", "")
            counts["rows"] += 1
            if vcf_id == ".":
                counts["vcf_id_missing"] += 1
            if fmt_id in {"", "."}:
                counts["format_id_missing"] += 1
            if fmt_id and fmt_id != "." and fmt_id != vcf_id:
                counts["format_id_differs_from_vcf_id"] += 1
                if len(examples) < max_examples:
                    examples.append(
                        {
                            "CHR": fields[0],
                            "BP": fields[1],
                            "VCF_ID": vcf_id,
                            "FORMAT_ID": fmt_id,
                            "REF": fields[3],
                            "ALT": fields[4],
                        }
                    )
            else:
                counts["format_id_same_or_missing"] += 1
    return counts, examples


def coord_fallback_examples_from_standardized(trait, max_examples=50):
    trait_by_coord = {}
    trait_snps = set()
    with open_text(STANDARD_DIR / f"{trait}.standardized.tsv.gz") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            trait_snps.add(row["SNP"])
            key = (row["CHR"], str(int(float(row["BP"]))))
            trait_by_coord.setdefault(key, row)

    counts = Counter()
    examples = []
    with open_text(STANDARD_DIR / "TRAIT_B.standardized.tsv.gz") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["SNP"] in trait_snps:
                counts["direct_snp_match"] += 1
                continue
            key = (row["CHR"], str(int(float(row["BP"]))))
            other = trait_by_coord.get(key)
            if other is None:
                counts["missing_by_snp_and_coord"] += 1
                continue
            counts["coord_fallback"] += 1
            if row["EA"] == other["EA"] and row["OA"] == other["OA"]:
                counts["fallback_same_orientation"] += 1
            elif row["EA"] == other["OA"] and row["OA"] == other["EA"]:
                counts["fallback_swap_orientation"] += 1
            else:
                counts["fallback_other_orientation"] += 1
            if len(examples) < max_examples:
                examples.append(
                    {
                        "CHR": row["CHR"],
                        "BP": row["BP"],
                        "ASD_SNP": row["SNP"],
                        f"{trait}_SNP": other["SNP"],
                        "ASD_EA": row["EA"],
                        "ASD_OA": row["OA"],
                        f"{trait}_EA": other["EA"],
                        f"{trait}_OA": other["OA"],
                    }
                )
    return counts, examples


def write_csv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    scz_pair_classes, scz_z_summary, scz_pair_examples = pair_file_summary("SCZ")
    anxiety_pair_classes, anxiety_z_summary, anxiety_pair_examples = pair_file_summary("disease_anxiety")
    scz_raw_counts, scz_raw_examples = scz_raw_vcf_id_review()
    scz_fallback_counts, scz_fallback_examples = coord_fallback_examples_from_standardized("SCZ")

    write_csv(OUT_DIR / "scz_raw_vcf_id_mismatch_examples.csv", scz_raw_examples)
    write_csv(OUT_DIR / "scz_coord_fallback_examples.csv", scz_fallback_examples)
    for cls, rows in anxiety_pair_examples.items():
        write_csv(OUT_DIR / f"anxiety_harmonized_{cls}_examples.csv", rows)
    for cls, rows in scz_pair_examples.items():
        write_csv(OUT_DIR / f"scz_harmonized_{cls}_examples.csv", rows)

    report = OUT_DIR / "scz_anxiety_p0_risk_review.md"
    report.write_text(
        "\n".join(
            [
                "# SCZ and anxiety P0 risk review",
                "",
                "## SCZ coord fallback",
                "",
                f"- Raw SCZ VCF ID review: {dict(scz_raw_counts)}",
                f"- Standardized ASD-SCZ fallback review: {dict(scz_fallback_counts)}",
                f"- Harmonized ASD-SCZ class counts: {dict(scz_pair_classes)}",
                f"- Harmonized ASD-SCZ Z-product summary by class: {scz_z_summary}",
                "",
                "Interpretation: SCZ has many coordinate fallback matches because many variants share CHR:BP with TRAIT_B but not rsID. This is consistent with rsID/version or identifier-source differences rather than immediate genome-build failure, because fallback sites are coordinate matched and overwhelmingly have same allele orientation in the standardized comparison. Keep the coord fallback flag in QC and inspect examples before using SCZ as a benchmark.",
                "",
                "## disease_anxiety allele and effect-direction review",
                "",
                f"- Harmonized ASD-anxiety class counts: {dict(anxiety_pair_classes)}",
                f"- Harmonized ASD-anxiety Z-product summary by class after sign flipping: {anxiety_z_summary}",
                "",
                "Interpretation: anxiety has lower SNP coverage than the VCF traits, which explains the larger missing_in_trait. The large swap count means many anxiety rows present the same biallelic SNP in the opposite allele order relative to ASD; the harmonization step flips those effects. This does not by itself prove OR is relative to A1. That assumption still needs source-document confirmation, but the post-harmonization Z-product summaries do not show an obvious catastrophic sign inversion.",
                "",
                "## Generated example files",
                "",
                f"- {OUT_DIR / 'scz_raw_vcf_id_mismatch_examples.csv'}",
                f"- {OUT_DIR / 'scz_coord_fallback_examples.csv'}",
                f"- {OUT_DIR / 'anxiety_harmonized_swap_examples.csv'}",
                f"- {OUT_DIR / 'anxiety_harmonized_match_examples.csv'}",
            ]
        )
    )
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
