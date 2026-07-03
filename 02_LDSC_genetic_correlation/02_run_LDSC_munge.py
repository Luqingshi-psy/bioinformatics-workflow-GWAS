#!/usr/bin/env python3
import subprocess
from pathlib import Path


LDSC_PYTHON = Path("${PROJECT_ROOT}")
LDSC_DIR = Path("${PROJECT_ROOT}")
STANDARD_DIR = Path("${PROJECT_ROOT}")
REF_DIR = Path("${PROJECT_ROOT}")
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


def run_trait(trait):
    input_path = STANDARD_DIR / f"{trait}.standardized.tsv.gz"
    out_prefix = OUT_DIR / trait
    log_path = OUT_DIR / f"{trait}.munge.stdout.log"
    err_path = OUT_DIR / f"{trait}.munge.stderr.log"
    expected = OUT_DIR / f"{trait}.sumstats.gz"
    if expected.exists() and expected.stat().st_size > 0:
        print(f"[skip] {trait}: {expected} exists")
        return
    cmd = [
        str(LDSC_PYTHON),
        str(LDSC_DIR / "munge_sumstats.py"),
        "--sumstats",
        str(input_path),
        "--snp",
        "SNP",
        "--a1",
        "EA",
        "--a2",
        "OA",
        "--p",
        "P",
        "--N-col",
        "N",
        "--signed-sumstats",
        "Z,0",
        "--merge-alleles",
        str(REF_DIR / "w_hm3.snplist"),
        "--out",
        str(out_prefix),
    ]
    print(f"[run] {trait}")
    with log_path.open("w") as stdout, err_path.open("w") as stderr:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for trait in TRAITS:
        run_trait(trait)
    print(f"Done. Outputs: {OUT_DIR}")


if __name__ == "__main__":
    main()
