#!/usr/bin/env Rscript
# ============================================================
# 01_HDL_global_Anxiety.R
# trait pair (A x B) Genome-wide genetic correlation analysis（HDL.rg）
# Analysis pairs：TRAIT_A_GWAS × Anxiety，IBS_FinnGen-style cohort × disease_anxiety # ============================================================

library(HDL)
library(data.table)

DATA_DIR     <- "${PROJECT_ROOT}"
OFFICIAL_REF <- "${PROJECT_ROOT}"
CLEAN_REF    <- "${PROJECT_ROOT}"
OUT_DIR      <- "${PROJECT_ROOT}"

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CLEAN_REF, showWarnings = FALSE, recursive = TRUE)

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
  for (f in keep_files)
    file.symlink(file.path(official.path, f), file.path(clean.dir, f))
  cat(sprintf(Clean ref msg：%d .rda，%d .bim\n, n_rda, n_bim))
  return(clean.dir)
}

cat(msg...\n)
ref_path <- setup_clean_ref(OFFICIAL_REF, CLEAN_REF, overwrite = FALSE)

n_rda <- length(list.files(ref_path, pattern = "\\.rda$"))
n_bim <- length(list.files(ref_path, pattern = "\\.bim$"))
cat(sprintf(msg：%d .rda，%d .bim, n_rda, n_bim))
if (n_rda < 10 || n_bim < 10) stop(msg，msg T7MAC msg。)
cat("  OK\n")

prep_hdl <- function(file) {
  dt <- fread(file, showProgress = FALSE)
  if (!("Z" %in% names(dt))) dt[, Z := BETA / SE]
  dt[, A1 := toupper(A1)]
  dt[, A2 := toupper(A2)]
  dt <- dt[!is.na(Z) & is.finite(Z) & !is.na(N) & N > 0]
  dt[, .(SNP, A1, A2, N, Z)]
}

cat(\nmsg GWAS msg...\n)
gwas_ibs     <- prep_hdl(file.path(DATA_DIR, "gwas_TRAIT_A.txt"))
gwas_anxiety <- prep_hdl(file.path(DATA_DIR, "gwas_anxiety_fmt.txt"))
gwas_finn    <- prep_hdl(file.path(DATA_DIR, "finn_ibs.txt"))
cat(sprintf(  IBS_GWAS:    %d SNPs (N msg = %d)\n,
            nrow(gwas_ibs),     as.integer(median(gwas_ibs$N))))
cat(sprintf(  disease_anxiety:     %d SNPs (N msg = %d)\n,
            nrow(gwas_anxiety), as.integer(median(gwas_anxiety$N))))
cat(sprintf(  IBS_FinnGen: %d SNPs (N msg = %d)\n,
            nrow(gwas_finn),    as.integer(median(gwas_finn$N))))

# ── run HDL.rg ───────────────────────────────────────────────
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
  }, error = function(e) { cat(msg：, conditionMessage(e), "\n"); NULL })

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

res1 <- run_HDL_rg(gwas_ibs,  gwas_anxiety, "TRAIT_A_GWAS",    "disease_anxiety", ref_path, OUT_DIR)
res2 <- run_HDL_rg(gwas_finn, gwas_anxiety, "TRAIT_A_FinnGen", "disease_anxiety", ref_path, OUT_DIR)

# ── summary ──────────────────────────────────────────────────────
summary_rows <- list()
for (pair_info in list(
  list(res = res1, n1 = "TRAIT_A_GWAS",    n2 = "disease_anxiety"),
  list(res = res2, n1 = "TRAIT_A_FinnGen", n2 = "disease_anxiety")
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
            file.path(OUT_DIR, "HDL_global_Anxiety_summary.csv"),
            row.names = FALSE)
  cat(\nmsg:, file.path(OUT_DIR, "HDL_global_Anxiety_summary.csv"), "\n")
}
cat(\nglobal HDL analysis (trait_pair_A_anxiety)msg。\n)
