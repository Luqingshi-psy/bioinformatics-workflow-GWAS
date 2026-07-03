#!/usr/bin/env python3
import csv
from collections import Counter

from analysis_config import ASD_CENTERED_TRAITS, DROP_AMBIGUOUS_IN_HARMONIZATION, RESULTS_DIR
from sumstats_utils import complement, dict_writer, is_ambiguous, open_text, safe_float


STANDARD_DIR = RESULTS_DIR / "standardized"
OUT_DIR = RESULTS_DIR / "harmonized_asd_pairs"
FIELDS = [
    "SNP",
    "CHR",
    "BP",
    "EA",
    "OA",
    "ASD_BETA",
    "ASD_SE",
    "ASD_Z",
    "ASD_P",
    "ASD_N",
    "TRAIT",
    "TRAIT_BETA",
    "TRAIT_SE",
    "TRAIT_Z",
    "TRAIT_P",
    "TRAIT_N",
    "HARMONIZATION_CLASS",
]


def coord_key(row):
    try:
        bp = str(int(float(row["BP"])))
    except (TypeError, ValueError):
        bp = str(row.get("BP", ""))
    return (str(row.get("CHR", "")), bp)


def load_standardized(path):
    by_snp = {}
    by_coord = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            by_snp[row["SNP"]] = row
            key = coord_key(row)
            if key not in by_coord:
                by_coord[key] = row
    return by_snp, by_coord


def classify(ref_ea, ref_oa, trait_ea, trait_oa):
    ref_ea, ref_oa = ref_ea.upper(), ref_oa.upper()
    trait_ea, trait_oa = trait_ea.upper(), trait_oa.upper()
    if trait_ea == ref_ea and trait_oa == ref_oa:
        return "match", 1
    if trait_ea == ref_oa and trait_oa == ref_ea:
        return "swap", -1
    comp_ea = complement(trait_ea)
    comp_oa = complement(trait_oa)
    if comp_ea == ref_ea and comp_oa == ref_oa:
        return "strand_match", 1
    if comp_ea == ref_oa and comp_oa == ref_ea:
        return "strand_swap", -1
    return "incompatible", 0


def harmonize_pair(asd_by_snp, asd_by_coord, trait_name):
    trait_path = STANDARD_DIR / f"{trait_name}.standardized.tsv.gz"
    trait_by_snp, trait_by_coord = load_standardized(trait_path)
    counters = Counter({"trait": trait_name})
    out_path = OUT_DIR / f"ASD_vs_{trait_name}.harmonized.tsv.gz"
    handle, writer = dict_writer(out_path, FIELDS, gz=True)
    with handle:
        for snp, ref in asd_by_snp.items():
            other = trait_by_snp.get(snp)
            if other is None:
                other = trait_by_coord.get(coord_key(ref))
                if other is not None:
                    counters["coord_fallback"] += 1
            if other is None:
                counters["missing_in_trait"] += 1
                continue
            counters["overlap"] += 1
            if DROP_AMBIGUOUS_IN_HARMONIZATION and (
                is_ambiguous(ref["EA"], ref["OA"]) or is_ambiguous(other["EA"], other["OA"])
            ):
                counters["drop_ambiguous"] += 1
                continue
            label, sign = classify(ref["EA"], ref["OA"], other["EA"], other["OA"])
            counters[label] += 1
            if sign == 0:
                continue
            trait_beta = safe_float(other["BETA"]) * sign
            trait_z = safe_float(other["Z"]) * sign
            writer.writerow(
                {
                    "SNP": snp,
                    "CHR": ref["CHR"],
                    "BP": ref["BP"],
                    "EA": ref["EA"],
                    "OA": ref["OA"],
                    "ASD_BETA": ref["BETA"],
                    "ASD_SE": ref["SE"],
                    "ASD_Z": ref["Z"],
                    "ASD_P": ref["P"],
                    "ASD_N": ref["N"],
                    "TRAIT": trait_name,
                    "TRAIT_BETA": trait_beta,
                    "TRAIT_SE": other["SE"],
                    "TRAIT_Z": trait_z,
                    "TRAIT_P": other["P"],
                    "TRAIT_N": other["N"],
                    "HARMONIZATION_CLASS": label,
                }
            )
            counters["kept"] += 1
    counters["output"] = str(out_path)
    return dict(counters)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    asd_by_snp, asd_by_coord = load_standardized(STANDARD_DIR / "TRAIT_B.standardized.tsv.gz")
    summaries = [harmonize_pair(asd_by_snp, asd_by_coord, trait) for trait in ASD_CENTERED_TRAITS]
    fieldnames = sorted({key for row in summaries for key in row})
    qc_path = OUT_DIR / "asd_centered_harmonization_qc.csv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {qc_path}")


if __name__ == "__main__":
    main()
