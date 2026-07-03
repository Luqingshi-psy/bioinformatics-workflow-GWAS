#!/usr/bin/env python3
from pathlib import Path

from analysis_config import RESULTS_DIR


def csv_to_markdown(path):
    lines = path.read_text().splitlines()
    if not lines:
        return ""
    header = lines[0].split(",")
    rows = [line.split(",") for line in lines[1:]]
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |"]
    for row in rows:
        row = row + [""] * (len(header) - len(row))
        out.append("| " + " | ".join(row[: len(header)]) + " |")
    return "\n".join(out)


def main():
    report = ["# IBS-ASD P0 QC report", ""]
    for rel in [
        "input_profile/input_profile_summary.csv",
        "standardized/standardization_qc_summary.csv",
        "harmonized_asd_pairs/asd_centered_harmonization_qc.csv",
    ]:
        path = RESULTS_DIR / rel
        report.append(f"## {rel}")
        report.append("")
        if path.exists():
            report.append(csv_to_markdown(path))
        else:
            report.append("Not generated yet.")
        report.append("")
    out = RESULTS_DIR / "p0_qc_report.md"
    out.write_text("\n".join(report))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

