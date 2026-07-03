# Third-Party Tools & Repositories

This repository contains only the **author's own workflow scripts**. The
third-party open-source tools listed below were used *as dependencies* by
the scripts in `GEO_data_processing/` and `GWAS_data_analysis/`. Their full
source code is **NOT** redistributed in this repository to respect each
project's license — install them directly from the upstream repositories
linked below.

If a script in this repo fails because one of these tools is missing,
follow the install command(s) listed for that tool.

Throughout this repository, placeholders are used in place of concrete
trait / gene / compound names:

| Placeholder        | Meaning                                              |
|--------------------|------------------------------------------------------|
| `TRAIT_A`, `TRAIT_B`, ... | Any GWAS trait in a pair (e.g. neuropsychiatric, GI) |
| `compound_X`       | A TCM / herbal compound under study                  |
| `gene_A`, `gene_B`, ... | Any candidate gene referenced in a script        |
| `metabolite_A`, ... | Any metabolite / mediator in an MR pipeline         |
| `${PROJECT_ROOT}`  | User-local project root directory                    |

When you reuse the workflow for a different study, substitute the
placeholders with your concrete names inside the script string literals.

---

## GWAS / Statistical Genetics

### R packages (CRAN / Bioconductor)

| Tool          | Used in scripts                 | Install                                     | License      |
|---------------|--------------------------------|---------------------------------------------|--------------|
| `ieugwasr`    | `01_fetch_GWAS_from_IEU.R`     | `install.packages("ieugwasr")`             | MIT          |
| `TwoSampleMR` | `01_fetch_GWAS_from_IEU.R`     | `install.packages("TwoSampleMR")`          | MIT          |
| `coloc`       | `06_colocalization/*`          | `install.packages("coloc")`                 | GPL-3        |
| `LAVA`        | `04_LAVA/*`                    | https://github.com/josefin-werme/LAVA       | GPL-3        |
| `HDL`         | `03_HDL_local_rg/*`            | https://github.com/zhenin/HDL               | GPL-3        |
| `placo`       | `05_PLACO_pleiotropy/*`        | https://github.com/RayDebashree/PLACO       | GPL-3        |
| `MRPRESSO`    | MR sensitivity (referenced)    | https://github.com/rondolab/MR-PRESSO       | GPL-3        |
| `RadialMR`    | MR sensitivity (referenced)    | https://github.com/WSpiller/RadialMR        | GPL-3        |
| `MRMix`       | MR mixture model (referenced)  | https://github.com/gqi/MRMix                | Artistic-2.0 |
| `MendelianRandomization` | MR core    | `install.packages("MendelianRandomization")`| GPL-2        |
| `MRInstruments` / `mrcieu` | helpers | `install.packages(...)` | see CRAN |

### Python packages / tools

| Tool          | Used in scripts                     | Install                                     | License |
|---------------|-------------------------------------|---------------------------------------------|---------|
| `ldsc`        | `02_LDSC_genetic_correlation/*`     | https://github.com/bulik/ldsc               | GPL-3   |
| `MiXeR`       | `09_MiXeR/*`                        | https://github.com/precimed/mixer           | MIT     |
| `pleioFDR`    | `10_pleioFDR_conjFDR/*`             | https://github.com/precimed/pleioFDR        | MIT     |
| `prop-coloc`  | colocalization references           | https://github.com/jsu27/prop-coloc        | MIT     |

### MATLAB

| Tool           | Used in scripts                       | Install / source                                          | License |
|----------------|---------------------------------------|-----------------------------------------------------------|---------|
| `pleiotropyFDR`| `10_pleioFDR_conjFDR/02_run_*.m`     | https://github.com/precimed/pleiotropyFDR (MATLAB scripts) | MIT     |

---

## GEO / Transcriptomics

### R packages (CRAN / Bioconductor)

| Tool              | Used in scripts                | Install                                                          | License |
|-------------------|--------------------------------|------------------------------------------------------------------|---------|
| `limma`           | `01_DEG_analysis/*`            | `BiocManager::install("limma")`                                  | GPL     |
| `edgeR`           | `01_DEG_analysis/*`            | `BiocManager::install("edgeR")`                                  | GPL     |
| `biomaRt`         | `04_gene_id_conversion/*`      | `BiocManager::install("biomaRt")`                                | LGPL    |
| `clusterProfiler` | `02_enrichment/*`              | `BiocManager::install("clusterProfiler")`                        | MIT     |
| `org.Hs.eg.db`    | enrichment / id-conversion     | `BiocManager::install("org.Hs.eg.db")`                           | Artistic-2.0 |
| `ggplot2`         | plotting across GEO scripts    | `install.packages("ggplot2")`                                    | MIT     |
| `ggvenn` / `venn` | `03_Venn_intersection/*`       | `install.packages("ggvenn")` / `install.packages("venn")`       | MIT     |
| `readxl` / `writexl` | I/O in 02/03                  | `install.packages(...)`                                          | MIT     |
| `openxlsx`        | id-conversion                  | `install.packages("openxlsx")`                                   | MIT     |

### Python packages

| Tool           | Used in scripts                       | Install                          | License |
|----------------|---------------------------------------|----------------------------------|---------|
| `pandas`       | `04_gene_id_conversion/*`, np scripts | `pip install pandas`             | BSD-3   |
| `openpyxl`     | GEO Excel I/O                          | `pip install openpyxl`           | MIT     |
| `requests`     | GEO / PubChem downloads                | `pip install requests`           | Apache-2|

### Perl libraries (Network pharmacology)

| Tool             | Used in scripts                                              | Install                                                     | License |
|------------------|--------------------------------------------------------------|-------------------------------------------------------------|---------|
| `Text::CSV`      | `05_network_pharmacology/*score*.pl`                        | `cpan Text::CSV`                                            | Artistic-1.0 |
| `List::Util`     | `05_network_pharmacology/*score*.pl`                        | core                                                        | GPL/Artistic |
| `POSIX`          | `05_network_pharmacology/*score*.pl`                        | core                                                        | GPL/Artistic |

---

## Reference data sources

These are referenced by `01_sumstats_QC/` and other scripts:

- **IEU OpenGWAS** (https://gwas.mrcieu.ac.uk/) — used by
  `01_fetch_GWAS_from_IEU.R`. Requires an OpenGWAS JWT (see
  `OPENGWAS_JWT` env var in the script).
- **NCBI Gene / HomoloGene** (`homo_sapiens.gene_info`) — used by
  `04_gene_id_conversion/entrez_to_symbol.py`.
- **Ensembl BioMart** — used by `04_gene_id_conversion/biomart_id_translation.Rmd`.
- **TCM compound databases** (TCMSP / HERB / SwissTargetPrediction) — used
  by `05_network_pharmacology/*`. Each has its own download procedure;
  see the local install/notes file in your project.
- **PubChem** — used by `05_network_pharmacology/03_pubchem_*.R`.
- **GTEx** eQTL reference — used by `08_TWAS_SMR/*` and
  `06_colocalization/coloc_eQTL_SNP_lookup.R`. Pre-downloaded models
  required.
- **Microbiota GWAS summary statistics** (e.g. MiBioGen) — used by
  `07_MR_analysis/01_mediator_microbiota_to_TRAIT_pair_MR.R`. Download
  from https://mibiogen.gcc.rug.nl/.

---

## License of *this* repository

The author's own scripts in `GEO_data_processing/` and `GWAS_data_analysis/`
are released under the **MIT License** (see `LICENSE`). They are research
workflows — provided as-is, no warranty. Pull requests welcome.

Third-party software listed above retains its own license. By using any
script here you also agree to comply with the upstream licenses of the
tools it depends on.