#!/usr/bin/env Rscript
# ============================================================
# 08_p0b_bidirectional_mr.R
# Phase 0-B:Bidirectional MR（TRAIT_B→TRAIT_A，trait pair (A x B)）+ Steiger filtering
#
# Purpose：
#
#   TRAIT_B GWAS：SNP CHR position A1 A2 BETA SE EAF P N
# ============================================================

library(TwoSampleMR)
library(data.table)
library(ggplot2)

# ── pathconfig ─────────────────────────────────────────────────
GWAS_DIR <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IBS_N <- 486601
ASD_N <- 46351

P_GWS   <- 5e-8
P_RELAX <- 1e-5

cat(msg GWAS msg...\n)

# TRAIT_A：CHR SNP A1 A2 BETA SE N
ibs_raw <- fread(file.path(GWAS_DIR, "gwas_TRAIT_A.txt"), showProgress = FALSE)
setnames(ibs_raw, names(ibs_raw), toupper(names(ibs_raw)))
ibs_raw[, BETA := as.numeric(BETA)]
ibs_raw[, SE   := as.numeric(SE)]
ibs_raw[, P    := 2 * pnorm(-abs(BETA / SE))]
ibs_raw[, N    := IBS_N]
cat(sprintf("  IBS：%d SNPs\n", nrow(ibs_raw)))

# TRAIT_B：SNP CHR position A1 A2 BETA SE EAF P N
asd_raw <- fread(file.path(GWAS_DIR, "gwas_TRAIT_B.txt"), showProgress = FALSE)
setnames(asd_raw, names(asd_raw), toupper(names(asd_raw)))
asd_raw[, BETA := as.numeric(BETA)]
asd_raw[, SE   := as.numeric(SE)]
asd_raw[, P    := as.numeric(P)]
asd_raw[, EAF  := as.numeric(EAF)]
asd_raw[, N    := ASD_N]
cat(sprintf("  ASD：%d SNPs\n", nrow(asd_raw)))

cat(\nmsg...\n)

asd_iv <- asd_raw[P < P_GWS]
p_thresh_asd <- P_GWS
if (nrow(asd_iv) < 3) {
  asd_iv <- asd_raw[P < P_RELAX]
  p_thresh_asd <- P_RELAX
  cat(sprintf(  [WARN] TRAIT_B p<5e-8 msg %d SNPs，msg p<1e-5\n, nrow(asd_iv)))
}
cat(sprintf(  TRAIT_B msg（p<%s）：%d SNPs\n,
            formatC(p_thresh_asd, format="e", digits=0), nrow(asd_iv)))

ibs_iv <- ibs_raw[P < P_GWS]
p_thresh_ibs <- P_GWS
if (nrow(ibs_iv) < 3) {
  ibs_iv <- ibs_raw[P < P_RELAX]
  p_thresh_ibs <- P_RELAX
  cat(sprintf(  [WARN] TRAIT_A p<5e-8 msg %d SNPs，msg p<1e-5\n, nrow(ibs_iv)))
}
cat(sprintf(  TRAIT_A msg（p<%s）：%d SNPs\n,
            formatC(p_thresh_ibs, format="e", digits=0), nrow(ibs_iv)))

cat(\nmsg...\n)

exp_asd <- format_data(
  as.data.frame(asd_iv),
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "BETA",
  se_col            = "SE",
  pval_col          = "P",
  effect_allele_col = "A1",
  other_allele_col  = "A2",
  eaf_col           = "EAF",
  samplesize_col    = "N"
)
exp_asd$exposure <- "TRAIT_B"

ibs_as_outcome_raw <- ibs_raw[SNP %in% exp_asd$SNP]
out_ibs <- format_data(
  as.data.frame(ibs_as_outcome_raw),
  type              = "outcome",
  snp_col           = "SNP",
  beta_col          = "BETA",
  se_col            = "SE",
  pval_col          = "P",
  effect_allele_col = "A1",
  other_allele_col  = "A2",
  samplesize_col    = "N"
)
out_ibs$outcome            <- "TRAIT_A"
out_ibs$samplesize.outcome <- IBS_N
out_ibs$eaf.outcome        <- 0.5

exp_ibs <- format_data(
  as.data.frame(ibs_iv),
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "BETA",
  se_col            = "SE",
  pval_col          = "P",
  effect_allele_col = "A1",
  other_allele_col  = "A2",
  samplesize_col    = "N"
)
exp_ibs$exposure    <- "TRAIT_A"
exp_ibs$eaf.exposure <- 0.5

asd_as_outcome_raw <- asd_raw[SNP %in% exp_ibs$SNP]
out_asd <- format_data(
  as.data.frame(asd_as_outcome_raw),
  type              = "outcome",
  snp_col           = "SNP",
  beta_col          = "BETA",
  se_col            = "SE",
  pval_col          = "P",
  effect_allele_col = "A1",
  other_allele_col  = "A2",
  eaf_col           = "EAF",
  samplesize_col    = "N"
)
out_asd$outcome            <- "TRAIT_B"
out_asd$samplesize.outcome <- ASD_N

cat(\nmsg...\n)

har_asd_ibs <- harmonise_data(exp_asd, out_ibs, action = 2)
har_asd_ibs <- har_asd_ibs[har_asd_ibs$mr_keep == TRUE, ]
cat(sprintf(  ASD→IBS msg：%d SNPs\n, nrow(har_asd_ibs)))

har_ibs_asd <- harmonise_data(exp_ibs, out_asd, action = 2)
har_ibs_asd <- har_ibs_asd[har_ibs_asd$mr_keep == TRUE, ]
cat(sprintf(  IBS→ASD msg：%d SNPs\n, nrow(har_ibs_asd)))

fwrite(as.data.table(har_asd_ibs),
       file.path(OUT_DIR, "harmonised_ASD_to_IBS.txt"), sep = "\t")
fwrite(as.data.table(har_ibs_asd),
       file.path(OUT_DIR, "harmonised_IBS_to_ASD.txt"), sep = "\t")

cat(\nMR msg...\n)

run_mr_analysis <- function(har, label) {
  n <- nrow(har)
  if (n < 2) {
    cat(sprintf(  [SKIP] %s：SNP msg（%d msg）\n, label, n))
    return(NULL)
  }
  methods <- c("mr_ivw", "mr_egger_regression",
               "mr_weighted_median", "mr_weighted_mode")
  if (n < 3) {
    methods <- c("mr_ivw", "mr_weighted_median")
    cat(sprintf(  [INFO] %s：msg %d SNPs，msg Egger msg\n, label, n))
  }
  res <- mr(har, method_list = methods)
  res <- as.data.table(res)
  res[, `:=`(direction = label,
              OR    = exp(b),
              CI_lo = exp(b - 1.96 * se),
              CI_hi = exp(b + 1.96 * se))]
  res
}

mr_asd_ibs <- run_mr_analysis(har_asd_ibs, "ASD → IBS")
mr_ibs_asd <- run_mr_analysis(har_ibs_asd, "IBS → ASD")

mr_all <- rbindlist(Filter(Negate(is.null), list(mr_asd_ibs, mr_ibs_asd)))

if (nrow(mr_all) > 0) {
  setorder(mr_all, direction, method)
  fwrite(mr_all, file.path(OUT_DIR, "bidirectional_MR_results.txt"), sep = "\t")

  cat(\n=== msg MR msg（IVW）===\n)
  ivw <- mr_all[method == "Inverse variance weighted"]
  print(ivw[, .(direction, nsnp,
                b     = round(b, 4),
                OR    = round(OR, 3),
                CI_lo = round(CI_lo, 3),
                CI_hi = round(CI_hi, 3),
                pval  = formatC(pval, format = "e", digits = 2))],
        row.names = FALSE)
}

cat(\nmsg...\n)

for (item in list(
  list(har = har_asd_ibs, label = "ASD_to_IBS"),
  list(har = har_ibs_asd, label = "IBS_to_ASD")
)) {
  h   <- item$har
  lbl <- item$label
  if (nrow(h) < 3) next

  het <- tryCatch(mr_heterogeneity(h), error = function(e) NULL)
  if (!is.null(het) && nrow(het) > 0) {
    fwrite(as.data.table(het),
           file.path(OUT_DIR, sprintf("heterogeneity_%s.txt", lbl)), sep = "\t")
    cat(sprintf(  %s msg（IVW Q-pval）：%.3f\n,
                lbl, het$Q_pval[het$method == "Inverse variance weighted"]))
  }

  ple <- tryCatch(mr_pleiotropy_test(h), error = function(e) NULL)
  if (!is.null(ple) && nrow(ple) > 0) {
    cat(sprintf(  %s MR-Egger msg：%.5f（p = %.4f）\n,
                lbl, ple$egger_intercept, ple$pval))
    fwrite(as.data.table(ple),
           file.path(OUT_DIR, sprintf("pleiotropy_%s.txt", lbl)), sep = "\t")
  }
}

# ── 7. Steiger filtering ──────────────────────────────────────────
cat(\nSteiger msg...\n)


run_steiger <- function(har, label, n_exp, n_out, eaf_placeholder = FALSE) {
  if (nrow(har) < 2) return(NULL)

  har2 <- as.data.frame(har)
  har2$samplesize.exposure <- n_exp
  har2$samplesize.outcome  <- n_out

  if (!"eaf.exposure" %in% names(har2) || all(is.na(har2$eaf.exposure)))
    har2$eaf.exposure <- 0.5
  if (!"eaf.outcome" %in% names(har2) || all(is.na(har2$eaf.outcome)))
    har2$eaf.outcome <- 0.5

  st <- tryCatch(steiger_filtering(har2), error = function(e) {
    cat(sprintf(  [WARN] %s Steiger msg：%s\n, label, conditionMessage(e))
    )
    NULL
  })
  if (is.null(st)) return(NULL)
  st <- as.data.table(st)
  st[, direction := label]

  n_correct <- sum(st$steiger_dir == TRUE,  na.rm = TRUE)
  n_wrong   <- sum(st$steiger_dir == FALSE, na.rm = TRUE)
  n_total   <- n_correct + n_wrong

  eaf_note <- if (eaf_placeholder) （EAF msg，msg） else ""
  cat(sprintf(  %s：%d/%d SNPs msg（%.0f%%）%s\n,
              label, n_correct, n_total,
              100 * n_correct / max(n_total, 1), eaf_note))
  st
}

st_asd_ibs <- run_steiger(har_asd_ibs, "ASD → IBS",
                           n_exp = ASD_N, n_out = IBS_N,
                           eaf_placeholder = TRUE)
st_ibs_asd <- run_steiger(har_ibs_asd, "IBS → ASD",
                           n_exp = IBS_N, n_out = ASD_N,
                           eaf_placeholder = TRUE)

steiger_all <- rbindlist(Filter(Negate(is.null), list(st_asd_ibs, st_ibs_asd)))
if (nrow(steiger_all) > 0)
  fwrite(steiger_all,
         file.path(OUT_DIR, "steiger_filtering_results.txt"), sep = "\t")

cat(\nSteiger msg...\n)

mr_steiger_list <- list()

if (!is.null(st_asd_ibs)) {
  keep <- st_asd_ibs[steiger_dir == TRUE]$SNP
  har_f <- har_asd_ibs[har_asd_ibs$SNP %in% keep, ]
  cat(sprintf(  ASD→IBS Steiger msg：%d / %d SNPs msg\n,
              nrow(har_f), nrow(har_asd_ibs)))
  res_f <- run_mr_analysis(har_f, "ASD → TRAIT_A (Steiger filtered)")
  if (!is.null(res_f)) mr_steiger_list[["ASD_IBS"]] <- res_f
}

if (!is.null(st_ibs_asd)) {
  keep <- st_ibs_asd[steiger_dir == TRUE]$SNP
  har_f <- har_ibs_asd[har_ibs_asd$SNP %in% keep, ]
  cat(sprintf(  IBS→ASD Steiger msg：%d / %d SNPs msg\n,
              nrow(har_f), nrow(har_ibs_asd)))
  res_f <- run_mr_analysis(har_f, "IBS → TRAIT_B (Steiger filtered)")
  if (!is.null(res_f)) mr_steiger_list[["TRAIT_PAIR"]] <- res_f
}

if (length(mr_steiger_list) > 0) {
  mr_steiger_dt <- rbindlist(mr_steiger_list)
  fwrite(mr_steiger_dt,
         file.path(OUT_DIR, "bidirectional_MR_steiger_filtered.txt"), sep = "\t")

  cat(\n=== Steiger msg MR（IVW）===\n)
  ivw_st <- mr_steiger_dt[method == "Inverse variance weighted"]
  print(ivw_st[, .(direction, nsnp,
                   b     = round(b, 4),
                   OR    = round(OR, 3),
                   CI_lo = round(CI_lo, 3),
                   CI_hi = round(CI_hi, 3),
                   pval  = formatC(pval, format = "e", digits = 2))],
        row.names = FALSE)
}

cat(\nmsg...\n)

all_mr_for_plot <- rbindlist(Filter(Negate(is.null), list(
  if (nrow(mr_all) > 0) mr_all[method == "Inverse variance weighted"] else NULL,
  if (length(mr_steiger_list) > 0)
    rbindlist(mr_steiger_list)[method == "Inverse variance weighted"] else NULL
)))

if (nrow(all_mr_for_plot) > 0) {
  color_map <- c(
    "ASD → IBS"                    = "#D6604D",
    "IBS → ASD"                    = "#2166AC",
    "ASD → TRAIT_A (Steiger filtered)" = "#F4A582",
    "IBS → TRAIT_B (Steiger filtered)" = "#92C5DE"
  )

  p_forest <- ggplot(all_mr_for_plot,
                     aes(x = OR, y = direction, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                   height = 0.25, linewidth = 0.7) +
    geom_point(size = 4) +
    scale_x_log10() +
    scale_color_manual(values = color_map, guide = "none") +
    theme_bw(base_size = 11) +
    labs(title = msg MR msg（IVW，msg Steiger msg）,
         subtitle = OR>1 msg p<0.05 msg,
         x = "OR (95% CI，log scale)", y = NULL) +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(OUT_DIR, "bidirectional_MR_forest.pdf"),
         p_forest, width = 9, height = max(3, nrow(all_mr_for_plot) * 0.6))
  cat(  msg\n)
}

cat(\n=== msg ===\n)
cat("┌─────────────────────────────────────────────────────────┐\n")
cat(│  msg        OR > 1, p < 0.05  │  msg                  │\n)
cat("├──────────────┤────────────────────────────────────────── │\n")
cat(│  TRAIT_B → TRAIT_A   YES               │  TRAIT_B msg→IBS msg  │\n)
cat(│  TRAIT_B → TRAIT_A   NO (p ≥ 0.05)     │  msg（msg）│\n)
cat(│  TRAIT_A → TRAIT_B   YES               │  TRAIT_A msg→ASD msg  │\n)
cat(│  TRAIT_A → TRAIT_B   NO (p ≥ 0.05)     │  msg（msg）│\n)
cat("├──────────────────────────────────────────────────────────┤\n")
cat(│  Steiger msg = msg                 │\n)
cat(│  TRAIT_A EAF msg（msg 0.5 msg），Steiger R² msg       │\n)
cat("└──────────────────────────────────────────────────────────┘\n")

cat(sprintf(\n=== P0-B msg ===\n))
cat(sprintf(msg：%s\n, OUT_DIR))
cat(  bidirectional_MR_results.txt           msg MR msg\n)
cat(  steiger_filtering_results.txt          msg SNP msg Steiger msg\n)
cat(  bidirectional_MR_steiger_filtered.txt  Steiger msg\n)
cat(  bidirectional_MR_forest.pdf            msg\n)
