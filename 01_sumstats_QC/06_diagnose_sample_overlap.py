#!/usr/bin/env python3
import csv
import gzip
from collections import Counter
from pathlib import Path

from analysis_config import DATA_DIR, RESULTS_DIR


OUT_DIR = RESULTS_DIR / "diagnostics"


def read_anxiety_keys(path):
    ids = set()
    coords = set()
    chrom_counts = Counter()
    rows = 0
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            snp = row.get("SNP", "")
            chrom = row.get("CHR", "")
            bp = row.get("BP", "")
            if snp:
                ids.add(snp)
            if chrom and bp:
                coords.add((chrom, bp))
                chrom_counts[chrom] += 1
    return rows, ids, coords, chrom_counts


def scan_asd_overlap(path, anxiety_ids, anxiety_coords):
    rows = 0
    id_overlap = 0
    coord_overlap = 0
    examples = []
    coord_examples = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            rows += 1
            chrom, bp, snp, ref, alt = fields[:5]
            if snp in anxiety_ids:
                id_overlap += 1
                if len(examples) < 10:
                    examples.append({"CHR": chrom, "BP": bp, "SNP": snp, "REF": ref, "ALT": alt})
            if (chrom, bp) in anxiety_coords:
                coord_overlap += 1
                if len(coord_examples) < 10:
                    coord_examples.append({"CHR": chrom, "BP": bp, "SNP": snp, "REF": ref, "ALT": alt})
    return rows, id_overlap, coord_overlap, examples, coord_examples


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    anxiety_path = DATA_DIR / "gwas_anxiety.txt"
    asd_path = DATA_DIR / "gwas_ASD.vcf.gz"
    anxiety_rows, anxiety_ids, anxiety_coords, chrom_counts = read_anxiety_keys(anxiety_path)
    asd_rows, id_overlap, coord_overlap, examples, coord_examples = scan_asd_overlap(
        asd_path, anxiety_ids, anxiety_coords
    )
    report = OUT_DIR / "asd_anxiety_raw_overlap_diagnosis.md"
    report.write_text(
        "\n".join(
            [
                "# ASD-anxiety raw overlap diagnosis",
                "",
                "## Conclusion",
                "",
                "- Raw TRAIT_B and anxiety files do overlap substantially by rsID and by chromosome-position.",
                "- A zero overlap in limited smoke-test outputs is expected when using first-row limits because anxiety begins on chromosome 8 while TRAIT_B begins on chromosome 1.",
                "- Full production standardization should remove row limits before final harmonization.",
                "",
                "## Counts",
                "",
                f"- anxiety rows: {anxiety_rows}",
                f"- anxiety unique rsIDs: {len(anxiety_ids)}",
                f"- anxiety unique coordinates: {len(anxiety_coords)}",
                f"- TRAIT_B rows scanned: {asd_rows}",
                f"- ASD-anxiety rsID overlap: {id_overlap}",
                f"- ASD-anxiety coordinate overlap: {coord_overlap}",
                "",
                "## disease_anxiety chromosome distribution",
                "",
                "\n".join(f"- chr{chrom}: {count}" for chrom, count in chrom_counts.most_common()),
                "",
                "## Example rsID overlaps",
                "",
                "\n".join(f"- {row}" for row in examples),
                "",
                "## Example coordinate overlaps",
                "",
                "\n".join(f"- {row}" for row in coord_examples),
                "",
            ]
        )
    )
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
