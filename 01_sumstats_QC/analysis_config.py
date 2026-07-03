from pathlib import Path

DATA_DIR = Path("${PROJECT_ROOT}")
BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"

TRAITS = {
    "TRAIT_B": DATA_DIR / "gwas_ASD.vcf.gz",
    "TRAIT_A": DATA_DIR / "gwas_IBS.vcf.gz",
    "IBD_deLange": DATA_DIR / "gwas_IBD_deLange.vcf.gz",
    "CD_deLange": DATA_DIR / "gwas_CD_deLange.vcf.gz",
    "UC_deLange": DATA_DIR / "gwas_UC_deLange.vcf.gz",
    "ADHD": DATA_DIR / "gwas_ADHD.vcf.gz",
    "SCZ": DATA_DIR / "gwas_SCZ.vcf.gz",
    "disease_anxiety": DATA_DIR / "gwas_anxiety.txt",
}

ASD_CENTERED_TRAITS = [
    "TRAIT_A",
    "ADHD",
    "SCZ",
    "disease_anxiety",
    "IBD_deLange",
    "CD_deLange",
    "UC_deLange",
]

MAX_ABS_Z = 80.0
SNP_ONLY = True
DROP_AMBIGUOUS_IN_HARMONIZATION = True

