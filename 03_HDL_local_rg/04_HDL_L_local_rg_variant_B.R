#!/usr/bin/env Rscript
# ============================================================
# 04_HDL_L_local_rg_Anxiety.R
# Local genetic correlation (HDL.L) — chromosome-by-chromosome
# with checkpoint restart support.
# ============================================================

if (!requireNamespace("pbmcapply", quietly = TRUE))
  install.packages("pbmcapply")

library(HDL)
library(data.table)
library(parallel)
library(pbmcapply)

DATA_DIR <- "${PROJECT_ROOT}"
LD_PATH  <- "${PROJECT_ROOT}"
BIM_PATH <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"
CKPT_DIR <- file.path(OUT_DIR, "checkpoints")

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

prep_hdll <- function(file) {
  dt <- fread(file, showProgress = FALSE)
  if (!("Z" %in% names(dt))) dt[, Z := BETA / SE]
  dt[, A1 := toupper(A1)]
  dt[, A2 := toupper(A2)]
  dt <- dt[!is.na(Z) & is.finite(Z) & !is.na(N) & N > 0]
  dt[, .(SNP, A1, A2, N, Z)]
}

cat(msg GWAS msg...\n)
gwas_ibs     <- prep_hdll(file.path(DATA_DIR, "gwas_TRAIT_A.txt"))
gwas_anxiety <- prep_hdll(file.path(DATA_DIR, "gwas_anxiety_fmt.txt"))
gwas_finn    <- prep_hdll(file.path(DATA_DIR, "finn_ibs.txt"))
cat(sprintf("  IBS_GWAS:    %d SNPs\n", nrow(gwas_ibs)))
cat(sprintf("  disease_anxiety:     %d SNPs\n", nrow(gwas_anxiety)))
cat(sprintf("  IBS_FinnGen: %d SNPs\n", nrow(gwas_finn)))

load(file.path(LD_PATH, "HDLL_LOC_snps.RData"))
cat(sprintf(msg %d msg（pieces）\n, nrow(NEWLOC)))

log_msg <- function(log_file, msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

run_hdll_pair <- function(gwas1, gwas2, name1, name2, cores = 2) {
  tag      <- paste0(name1, "_x_", name2)
  rds_file <- file.path(OUT_DIR,  paste0("res_HDLL_", tag, ".rds"))
  csv_file <- file.path(OUT_DIR,  paste0("res_HDLL_", tag, ".csv"))
  log_file <- file.path(OUT_DIR,  paste0("progress_", tag, ".log"))

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("HDL.L:", name1, "×", name2, "\n")
  cat(strrep("=", 60), "\n")

  if (file.exists(rds_file)) {
    cat(msg，msg：, rds_file, "\n")
    return(readRDS(rds_file))
  }

  cat(sprintf(HDL.L %s × %s  msg %s\n,
              name1, name2, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = log_file)

  process_piece <- function(chr, piece) {
    tryCatch(
      HDL.L(gwas1.df   = gwas1,   gwas2.df  = gwas2,
            Trait1name = name1,   Trait2name = name2,
            LD.path    = LD_PATH, bim.path   = BIM_PATH,
            chr        = as.character(chr),
            piece      = as.integer(piece),
            N0 = 0, Nref = 335272, eigen.cut = 0.99, lim = exp(-18)),
      error = function(e) {
        message(sprintf(  [msg] chr%s piece%s: %s, chr, piece, conditionMessage(e)))
        NULL
      }
    )
  }

  chrs      <- sort(unique(NEWLOC$CHR))
  n_total   <- nrow(NEWLOC)
  n_done    <- 0L
  all_results <- list()

  log_msg(log_file, sprintf(msg，msg %d msg，%d msg，%d msg,
                             length(chrs), n_total, cores))

  for (chr in chrs) {
    ckpt_file <- file.path(CKPT_DIR, sprintf("%s_chr%s.rds", tag, chr))

    if (file.exists(ckpt_file)) {
      chr_res <- readRDS(ckpt_file)
      n_done  <- n_done + nrow(chr_res)
      all_results <- c(all_results, list(chr_res))
      log_msg(log_file, sprintf(chr%s  msg（%d pieces）  msg %d/%d (%.1f%%),
                                chr, nrow(chr_res), n_done, n_total,
                                100 * n_done / n_total))
      next
    }

    pieces <- NEWLOC$piece[NEWLOC$CHR == chr]
    log_msg(log_file, sprintf(chr%s  msg  %d pieces ..., chr, length(pieces)))

    args_list <- lapply(pieces, function(p) list(chr = chr, piece = p))
    chr_results_raw <- pbmclapply(
      args_list,
      function(args) process_piece(args$chr, args$piece),
      mc.cores           = cores,
      ignore.interactive = TRUE
    )

    chr_res <- rbindlist(chr_results_raw, fill = TRUE, use.names = TRUE)
    n_done  <- n_done + nrow(chr_res)
    all_results <- c(all_results, list(chr_res))
    saveRDS(chr_res, ckpt_file)

    n_sig <- if ("p" %in% names(chr_res)) sum(!is.na(chr_res$p) & chr_res$p < 0.05) else 0L
    log_msg(log_file, sprintf(chr%s  msg  %d/%d pieces  msg(p<0.05):%d  msg %d/%d (%.1f%%),
                              chr, nrow(chr_res), length(pieces),
                              n_sig, n_done, n_total, 100 * n_done / n_total))
  }

  res_df <- rbindlist(all_results, fill = TRUE, use.names = TRUE)

  if (nrow(res_df) > 0) {
    saveRDS(res_df, rds_file)
    fwrite(res_df, csv_file)
    n_sig_total <- if ("p" %in% names(res_df))
      sum(!is.na(res_df$p) & res_df$p < 0.05) else 0L
    log_msg(log_file, sprintf(msg！msg %d msg，msg(p<0.05)：%d msg, nrow(res_df), n_sig_total))
    cat(sprintf(msg：%s\n, csv_file))
    if ("p" %in% names(res_df) && n_sig_total > 0) {
      cat(\nTop msg：\n)
      cols <- intersect(c("chr", "piece", "rg", "rg_lower", "rg_upper", "p"), names(res_df))
      sig  <- res_df[!is.na(res_df$p) & res_df$p < 0.05, ]
      print(as.data.frame(sig[order(sig$p), .SD, .SDcols = cols][seq_len(min(10, nrow(sig)))]),
            row.names = FALSE)
    }
  }
  return(res_df)
}

NCORES <- 2

# Pair 1: TRAIT_A_GWAS × Anxiety（primary analysis）
res1 <- run_hdll_pair(gwas_ibs,  gwas_anxiety, "TRAIT_A_GWAS",    "disease_anxiety", cores = NCORES)
rm(res1); gc()

# Pair 2: IBS_FinnGen-style cohort × Anxiety（replicate check）
res2 <- run_hdll_pair(gwas_finn, gwas_anxiety, "TRAIT_A_FinnGen", "disease_anxiety", cores = NCORES)
rm(res2); gc()

cat("\n\n", strrep("=", 60), "\n", sep = "")
cat(HDL.L msg（p < 0.05）\n)
cat(strrep("=", 60), "\n")

for (tag in c("IBS_GWAS_x_Anxiety", "IBS_FinnGen_x_Anxiety")) {
  csv_file <- file.path(OUT_DIR, paste0("res_HDLL_", tag, ".csv"))
  if (!file.exists(csv_file)) next
  d <- fread(csv_file)
  if (!"p" %in% names(d)) next
  sig <- d[!is.na(p) & p < 0.05]
  cat(sprintf(\n%s：%d msg / %d msg\n, tag, nrow(sig), nrow(d)))
  if (nrow(sig) > 0) {
    cols <- intersect(c("chr", "piece", "rg", "rg_lower", "rg_upper", "p"), names(sig))
    print(sig[order(p), ..cols][seq_len(min(10, nrow(sig)))], row.names = FALSE)
  }
}
cat(\nHDL.L analysis (trait_pair_A_anxiety)msg。\n)
