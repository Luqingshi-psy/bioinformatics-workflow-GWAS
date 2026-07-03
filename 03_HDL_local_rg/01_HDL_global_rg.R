# ============================================================
# 01_HDL_global_rg.R
# Genome-wide genetic correlation (HDL.rg) for a trait pair.
#
# Reference panel: official UKB imputed SVD eigen99 extraction
#   (see third_party/REFERENCES.md for download)
#   The original directory may contain a macOS duplicate
#   "(1).bim" file that makes HDL.rg grep match two BIM files
#   and fail — hence the clean symlink step below.
#
# Trait pairs (configure for your project):
#   1. TRAIT_A_GWAS    (large N)  × TRAIT_B (moderate N)
#   2. TRAIT_A_FinnGen-style cohort (FinnGen-style cohort-style cohort)  × TRAIT_B (moderate N)
# ============================================================

library(HDL)
library(data.table)

# ── Paths — override via env vars or edit here ─────────────────
DATA_DIR     <- Sys.getenv("GWAS_DATA_DIR", "${PROJECT_ROOT}/data/")
OFFICIAL_REF <- Sys.getenv("LD_REF_DIR",    "${PROJECT_ROOT}/ref/UKB_imputed_SVD_eigen99_extraction/UKB_imputed_SVD_eigen99_extraction")
# Clean symlink dir without "(1).bim" duplicates (no trailing slash)
CLEAN_REF    <- Sys.getenv("LD_REF_CLEAN",  "${PROJECT_ROOT}/ref/clean_no_duplicates")
OUT_DIR      <- Sys.getenv("HDL_OUT_DIR",   "${PROJECT_ROOT}/output/HDL_global/")

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CLEAN_REF, showWarnings = FALSE, recursive = TRUE)

# ── Step 1：Build clean symlink dir excluding duplicate BIM ─────────────
setup_clean_ref <- function(official.path, clean.dir, overwrite = TRUE) {
  all_files  <- list.files(official.path)
  keep_files <- all_files[!grepl("\\(1\\)", all_files, fixed = TRUE)]

  n_rda <- sum(grepl("\\.rda$", keep_files))
  n_bim <- sum(grepl("\\.bim$", keep_files))

  if (!overwrite) {
    existing_rda <- length(list.files(clean.dir, pattern = "\\.rda$"))
    if (existing_rda >= n_rda) {
      cat(sprintf(Clean ref msg（%d .rda），msg。\n, existing_rda))
      return(clean.dir)
    }
  }

  old <- list.files(clean.dir, full.names = TRUE)
  if (length(old) > 0) unlink(old)

  # rebuild symlink
  for (f in keep_files) {
    file.symlink(file.path(official.path, f),
                 file.path(clean.dir,     f))
  }

  cat(sprintf(Clean ref msg：%d .rda，%d .bim（msg %d msg）\n,
              n_rda, n_bim, length(all_files) - length(keep_files)))
  return(clean.dir)
}

cat(msg...\n)
ref_path <- setup_clean_ref(OFFICIAL_REF, CLEAN_REF, overwrite = TRUE)

# verify
n_rda <- length(list.files(ref_path, pattern = "\\.rda$"))
n_bim <- length(list.files(ref_path, pattern = "\\.bim$"))
cat(sprintf(msg：%d .rda，%d .bim, n_rda, n_bim))
if (n_rda < 10 || n_bim < 10) stop(msg，msg T7MAC msg。)
if (n_rda != n_bim) cat(sprintf(  [msg] .rda msg .bim msg\n)) else cat("  OK\n")

prep_hdl <- function(file) {
  dt <- fread(file, showProgress = FALSE)
  if (!("Z" %in% names(dt))) dt[, Z := BETA / SE]
  dt[, A1 := toupper(A1)]
  dt[, A2 := toupper(A2)]
  dt <- dt[!is.na(Z) & is.finite(Z) & !is.na(N) & N > 0]
  dt[, .(SNP, A1, A2, N, Z)]
}

cat(\nmsg GWAS msg...\n)
gwas_ibs  <- prep_hdl(file.path(DATA_DIR, "gwas_TRAIT_A.txt"))
gwas_asd  <- prep_hdl(file.path(DATA_DIR, "gwas_TRAIT_B.txt"))
gwas_finn <- prep_hdl(file.path(DATA_DIR, "finn_ibs.txt"))
cat(sprintf(  IBS_GWAS:    %d SNPs (N msg = %d)\n,
            nrow(gwas_ibs),  as.integer(median(gwas_ibs$N))))
cat(sprintf(  TRAIT_B:         %d SNPs (N msg = %d)\n,
            nrow(gwas_asd),  as.integer(median(gwas_asd$N))))
cat(sprintf(  IBS_FinnGen: %d SNPs (N msg = %d)\n,
            nrow(gwas_finn), as.integer(median(gwas_finn$N))))

# ── Step 3：run HDL.rg ───────────────────────────────────────
run_HDL_rg <- function(gwas1, gwas2, name1, name2, LD.path.rg, out.dir) {
  tag      <- paste0(name1, "_x_", name2)
  log_file <- file.path(out.dir, paste0("HDL_", tag, ".log"))
  rds_file <- file.path(out.dir, paste0("HDL_", tag, "_result.rds"))

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("HDL.rg:", name1, "×", name2, "\n")
  cat(strrep("=", 60), "\n")

  res <- tryCatch({
    HDL.rg(gwas1.df    = gwas1,
            gwas2.df    = gwas2,
            LD.path     = LD.path.rg,
            Nref        = 335272,
            N0          = 0,
            output.file = log_file,
            eigen.cut   = "automatic",
            lim         = exp(-18))
  }, error = function(e) {
    cat(msg：, conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res)) {
    saveRDS(res, rds_file)
    est    <- res$estimates.df
    rg_est <- est["Genetic_Correlation", "Estimate"]
    rg_se  <- est["Genetic_Correlation", "se"]
    h1_est <- est["Heritability_1",      "Estimate"]
    h2_est <- est["Heritability_2",      "Estimate"]
    P_val  <- res$P
    cat(sprintf(\nmsg：\n))
    cat(sprintf("  h²(%s) = %.4f\n",                 name1, h1_est))
    cat(sprintf("  h²(%s) = %.4f\n",                 name2, h2_est))
    cat(sprintf("  rg             = %.4f (SE = %.4f)\n", rg_est, rg_se))
    cat(sprintf(  Pmsg            = %.3e\n,          P_val))
    cat(sprintf(  msg：%s\n,               rds_file))
  }
  return(res)
}

# Analysis pairs 1：TRAIT_A_GWAS × TRAIT_B
res1 <- run_HDL_rg(
  gwas1      = gwas_ibs,
  gwas2      = gwas_asd,
  name1      = "TRAIT_A_GWAS",
  name2      = "TRAIT_B",
  LD.path.rg = ref_path,
  out.dir    = OUT_DIR
)

# Analysis pairs 2：IBS_FinnGen-style cohort × TRAIT_B
res2 <- run_HDL_rg(
  gwas1      = gwas_finn,
  gwas2      = gwas_asd,
  name1      = "TRAIT_A_FinnGen",
  name2      = "TRAIT_B",
  LD.path.rg = ref_path,
  out.dir    = OUT_DIR
)

cat("\n\n", strrep("=", 60), "\n", sep = "")
cat(HDL msg\n)
cat(strrep("=", 60), "\n")

summary_rows <- list()
for (pair_info in list(
  list(res = res1, n1 = "TRAIT_A_GWAS",    n2 = "TRAIT_B"),
  list(res = res2, n1 = "TRAIT_A_FinnGen", n2 = "TRAIT_B")
)) {
  if (!is.null(pair_info$res)) {
    est <- pair_info$res$estimates.df
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Trait1 = pair_info$n1,
      Trait2 = pair_info$n2,
      h2_T1  = round(est["Heritability_1",      "Estimate"], 4),
      h2_T2  = round(est["Heritability_2",      "Estimate"], 4),
      rg     = round(est["Genetic_Correlation", "Estimate"], 4),
      rg_SE  = round(est["Genetic_Correlation", "se"],       4),
      P      = formatC(pair_info$res$P, format = "e", digits = 2),
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  print(summary_df, row.names = FALSE)
  write.csv(summary_df,
            file.path(OUT_DIR, "HDL_global_summary.csv"),
            row.names = FALSE)
  cat(\nmsg:, file.path(OUT_DIR, "HDL_global_summary.csv"), "\n")
}

cat(\nmsg HDL msg。\n)
