#!/usr/bin/env python3
import csv
import re
from pathlib import Path


RG_DIR = Path("${PROJECT_ROOT}")
OUT_PATH = RG_DIR / "ASD_centered_ldsc_rg_summary.csv"


def parse_log(path):
    rows = []
    capture = False
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("p1 ") and "rg" in stripped and "p" in stripped:
                capture = True
                continue
            if capture:
                if not stripped or stripped.startswith("Analysis finished"):
                    break
                parts = re.split(r"\s+", stripped)
                if len(parts) >= 6 and parts[0].endswith(".sumstats.gz"):
                    rows.append(
                        {
                            "p1": parts[0],
                            "p2": parts[1],
                            "rg": parts[2],
                            "se": parts[3],
                            "z": parts[4],
                            "p": parts[5],
                            "log": str(path),
                        }
                    )
    return rows


def main():
    rows = []
    for log in sorted(RG_DIR.glob("ASD_vs_*.log")):
        if log.name.endswith(".stdout.log") or log.name.endswith(".stderr.log"):
            continue
        rows.extend(parse_log(log))
    with OUT_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["p1", "p2", "rg", "se", "z", "p", "log"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
