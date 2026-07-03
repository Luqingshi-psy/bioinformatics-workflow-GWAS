# GWAS Integrative Workflow

For GWAS summary-stats-based integrative analysis. Built around a
two-trait pair but parameterised so it generalises to any GWAS pair.

> Placeholders used throughout this repo:
> - `TRAIT_A`, `TRAIT_B`, `TRAIT_C`, ... : traits in a GWAS pair (e.g.
>   trait A = a neuropsychiatric phenotype, trait B = a GI phenotype)
> - `compound_X` : a single TCM / herbal compound (relevant to the
>   `07_MR_analysis/` mediator-MR scripts)
> - `gene_A`, `gene_B`, ... : any candidate gene referenced in a script
> - `metabolite_A`, `metabolite_B`, ... : any metabolite / mediator
> - `${PROJECT_ROOT}` : the user-local project root directory; replace
>   per environment, or set the corresponding environment variables
>   (`GWAS_DATA_DIR`, `LD_REF_DIR`, `LD_REF_CLEAN`, `HDL_OUT_DIR`, ...).
>
> Concrete trait / gene / metabolite names appear only where they are
> load-bearing for the script to run (e.g. `gwas_TRAIT_A.vcf.gz` inside
> a string literal). Rename those locally when you reuse the script
> for a different study.

> All scripts here are the author's own code. Third-party tools are
> listed in [`third_party/REFERENCES.md`](third_party/REFERENCES.md) —
> install them from upstream; their full source is **not** redistributed.

---

## Repository layout

```
GWAS_workflow/
├── README.md
├── LICENSE
├── third_party/
│   └── REFERENCES.md
├── 01_sumstats_QC/                       # IEU fetch, standardise, harmonise, sample-overlap diagnosis
├── 02_LDSC_genetic_correlation/          # Global genetic correlation (LD Score regression)
├── 03_HDL_local_rg/                      # Local genetic correlation (HDL, HDL.L)
├── 04_LAVA/                              # Local genetic correlation in pre-defined LD blocks
├── 05_PLACO_pleiotropy/                  # PLACO / PLACO+ shared-loci test
├── 06_colocalization/                    # coloc GWAS-GWAS + GWAS-eQTL SNP lookup
├── 07_MR_analysis/                       # Two-sample / bidirectional / mediator MR
├── 08_TWAS_SMR/                          # S-PrediXcan + SMR tissue-wide integration
├── 09_MiXeR/                             # MiXeR bivariate Gaussian copula mixture
├── 10_pleioFDR_conjFDR/                  # pleioFDR / conjFDR / condFDR (MATLAB + R)
└── 11_mediation_MR/                      # Two-step mediation MR (mediator = metabolite / microbiota)
```

Each numbered folder = one analytical stage. Read top-to-bottom = run
order.

---

## Pipeline

| Step | Folder | What it does |
|------|--------|--------------|
| 01   | `01_sumstats_QC/`                       | Fetch from IEU OpenGWAS, standardise, harmonise, diagnose sample overlap |
| 02   | `02_LDSC_genetic_correlation/`          | Global genetic correlation (LD Score regression) |
| 03   | `03_HDL_local_rg/`                      | Local genetic correlation (HDL, HDL.L), gene annotation of significant loci |
| 04   | `04_LAVA/`                              | Local genetic correlation in pre-defined LD blocks |
| 05   | `05_PLACO_pleiotropy/`                  | PLACO / PLACO+ shared-loci test across two traits |
| 06   | `06_colocalization/`                    | `coloc` GWAS-GWAS + GWAS-eQTL SNP lookup |
| 07   | `07_MR_analysis/`                       | Two-sample / bidirectional / mediator Mendelian Randomisation |
| 08   | `08_TWAS_SMR/`                          | S-PrediXcan + SMR tissue-wide integration |
| 09   | `09_MiXeR/`                             | MiXeR bivariate Gaussian copula mixture (shared / unique SNP counts) |
| 10   | `10_pleioFDR_conjFDR/`                  | pleioFDR / conjFDR / condFDR (MATLAB + R post-processing) |
| 11   | `11_mediation_MR/`                      | Two-step mediation MR (mediator = metabolite / microbiota) |

---

## Quick start

```bash
# 1) Fetch & QC
Rscript GWAS_workflow/01_sumstats_QC/01_fetch_GWAS_from_IEU.R
python  GWAS_workflow/01_sumstats_QC/03_profile_GWAS_inputs.py
python  GWAS_workflow/01_sumstats_QC/04_standardize_sumstats.py
python  GWAS_workflow/01_sumstats_QC/05_harmonize_TRAIT_B_centered_pairs.py

# 2) Global / local genetic correlation
python  GWAS_workflow/02_LDSC_genetic_correlation/02_run_LDSC_munge.py
Rscript GWAS_workflow/03_HDL_local_rg/01_HDL_global_rg.R
Rscript GWAS_workflow/03_HDL_local_rg/03_HDL_L_local_rg.R

# 3) Shared loci
Rscript GWAS_workflow/05_PLACO_pleiotropy/PLACO_shared_loci.R
Rscript GWAS_workflow/06_colocalization/coloc_GWAS_GWAS.R

# 4) MR / TWAS / SMR
Rscript GWAS_workflow/07_MR_analysis/01_mediator_microbiota_to_TRAIT_pair_MR.R
Rscript GWAS_workflow/08_TWAS_SMR/03_analyze_TWAS_results.R

# 5) Advanced shared-heritability models
bash    GWAS_workflow/09_MiXeR/02_run_MiXeR_bivariate.sh
matlab -batch "run('GWAS_workflow/10_pleioFDR_conjFDR/02_run_conjFDR_condFDR.m')"
```

---

## File-naming convention

Each script is prefixed with a **2-digit ordinal** reflecting its order
within its analytical stage. Examples:

- `01_fetch_GWAS_from_IEU.R`        ← stage 1, script 1
- `02_HDL_global_rg.R`              ← stage 3 (HDL), second in the global-rg pair
- `mediation_MR_mediator_metabolite.R` ← stage 11, mediator = metabolite MR

Naming is content-driven, not project-driven: names describe the
analysis (`coloc_GWAS_GWAS.R`, `LAVA_local_rg.R`), so the same files
can be reused for different phenotype pairs by renaming the placeholder
strings inside.

---

## Configuration & data layout

The pipeline expects a layout like:

```
<project_root>/
├── data/
│   ├── gwas_TRAIT_A.vcf.gz
│   ├── gwas_TRAIT_B.vcf.gz
│   └── ...
├── ref/
│   ├── UKB_imputed_SVD_eigen99_extraction/
│   └── clean_no_duplicates/      # symlink dir without "(1).bim" duplicates
├── results/
│   ├── standardized/
│   ├── ldsc_munge/
│   ├── ldsc_rg/
│   ├── HDL_global/
│   └── ...
└── scripts/                # <- this repo lives here
```

Paths default to `${PROJECT_ROOT}` placeholders. Override via env vars:

| Variable          | Used in                       | Default            |
|-------------------|-------------------------------|--------------------|
| `GWAS_DATA_DIR`   | most GWAS scripts             | `${PROJECT_ROOT}/data/` |
| `LD_REF_DIR`      | LDSC / HDL                    | `${PROJECT_ROOT}/ref/UKB_imputed_SVD_eigen99_extraction/...` |
| `LD_REF_CLEAN`    | HDL clean-symlink step        | `${PROJECT_ROOT}/ref/clean_no_duplicates` |
| `HDL_OUT_DIR`     | `03_HDL_local_rg/*`           | `${PROJECT_ROOT}/output/HDL_global/` |
| ... (per-script)  | see top-of-file constants     | ...                |

---

## Requirements

- **R >= 4.2** with `BiocManager`, plus packages listed in
  `third_party/REFERENCES.md` (CRAN + Bioconductor).
- **Python >= 3.10** with `pandas`, plus `ldsc` (clone from
  github.com/bulik/ldsc into a virtualenv).
- **MATLAB** (only for `10_pleioFDR_conjFDR/02_run_conjFDR_condFDR.m`).
- **R packages not on CRAN** (`HDL`, `LAVA`, `placo`, `TwoSampleMR`,
  `ieugwasr`, `MR-PRESSO`, `RadialMR`, `MRMix`, `prop-coloc`) — install
  via `remotes::install_github(...)` per `third_party/REFERENCES.md`.
- **MiXeR** — Docker image (`precimed/mixer`) is the easiest install.

---

## Reproducibility caveats

- Many scripts default to `${PROJECT_ROOT}` placeholders. Set the
  matching environment variables (see table above) before running.
- Sample-size-dependent thresholds (FDR, R-squared, magic-number
  p-cutoffs) are hard-coded for the trait pair used during
  development; revisit before applying to a new pair.
- External reference panels (LD scores, GTEx eQTL, microbiota
  sumstats) must be downloaded separately; URLs and installation
  notes are in `third_party/REFERENCES.md`.

---

## Contributing

This is a personal workflow. Issues and PRs are welcome, especially:
- Path-config refactors (replace `${PROJECT_ROOT}` with config files).
- New analytical stages (e.g. fine-mapping, eCAVIAR, transcriptome-wide
  MR).
- Documentation of expected input / output for any script.

---

## License

MIT — see [`LICENSE`](LICENSE).

## Citation

If you use any of these workflows in published work, please cite the
underlying methodological papers (LDSC, HDL, LAVA, PLACO, coloc, TWAS,
SMR, MiXeR, pleioFDR, MR-PRESSO, etc.) — links in
`third_party/REFERENCES.md`.