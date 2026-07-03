#!/usr/bin/env python3
from analysis_config import RESULTS_DIR, TRAITS


STANDARD_DIR = RESULTS_DIR / "standardized"
OUT_DIR = RESULTS_DIR / "ldsc_munge"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Edit LDSC paths before running.",
        "# These commands assume standardized columns: SNP, EA, OA, BETA, SE, P, N.",
        "",
        "LDSC_DIR=/path/to/ldsc",
        "REF_DIR=/path/to/eur_w_ld_chr",
        "",
    ]
    for trait in TRAITS:
        input_path = STANDARD_DIR / f"{trait}.standardized.tsv.gz"
        out_prefix = OUT_DIR / trait
        lines.append(
            "python ${LDSC_DIR}/munge_sumstats.py "
            f"--sumstats {input_path} "
            "--snp SNP --a1 EA --a2 OA --p P --N-col N "
            f"--out {out_prefix} --merge-alleles ${{REF_DIR}}/w_hm3.snplist"
        )
    script = OUT_DIR / "munge_ldsc_inputs.zsh"
    script.write_text("\n".join(lines) + "\n")
    print(f"Wrote {script}")


if __name__ == "__main__":
    main()

