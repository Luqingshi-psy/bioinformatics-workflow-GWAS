#!/usr/bin/env python3
import subprocess
from pathlib import Path


LDSC_PYTHON = Path("${PROJECT_ROOT}")
LDSC_DIR = Path("${PROJECT_ROOT}")
REF_DIR = Path("${PROJECT_ROOT}")
MUNGE_DIR = Path("${PROJECT_ROOT}")
OUT_DIR = Path("${PROJECT_ROOT}")

SECONDARY_TRAITS = [
    "TRAIT_A",
    "ADHD",
    "SCZ",
    "disease_anxiety",
    "IBD_deLange",
    "CD_deLange",
    "UC_deLange",
]


def run_pair(trait):
    asd = MUNGE_DIR / "TRAIT_B.sumstats.gz"
    other = MUNGE_DIR / f"{trait}.sumstats.gz"
    out_prefix = OUT_DIR / f"ASD_vs_{trait}"
    log_path = OUT_DIR / f"ASD_vs_{trait}.stdout.log"
    err_path = OUT_DIR / f"ASD_vs_{trait}.stderr.log"
    cmd = [
        str(LDSC_PYTHON),
        str(LDSC_DIR / "ldsc.py"),
        "--rg",
        f"{asd},{other}",
        "--ref-ld-chr",
        str(REF_DIR) + "/",
        "--w-ld-chr",
        str(REF_DIR) + "/",
        "--out",
        str(out_prefix),
    ]
    print(f"[run] ASD_vs_{trait}")
    with log_path.open("w") as stdout, err_path.open("w") as stderr:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for trait in SECONDARY_TRAITS:
        run_pair(trait)
    print(f"Done. Outputs: {OUT_DIR}")


if __name__ == "__main__":
    main()
