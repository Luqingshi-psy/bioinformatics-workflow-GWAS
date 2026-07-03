#!/usr/bin/env Rscript
# ============================================================
# 07_bidirectional_mr_Anxiety.R
# Bidirectional MR（disease_anxiety→TRAIT_A，trait pair (A x B)）+ Steiger filtering
#
# Purpose：
#
#   disease_anxiety GWAS：SNP CHR position A1 A2 BETA SE EAF P N
# ============================================================

library(TwoSampleMR)
library(data.table)
library(ggplot2)

# ── pathconfig ─────────────────────────────────────────────────
GWAS_DIR <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IBS_N <- 486601
ANX_N <- 418399

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

# disease_anxiety：SNP CHR position A1 A2 BETA SE EAF P N
anx_raw <- fread(file.path(GWAS_DIR, "gwas_anxiety_fmt.txt"), showProgress = FALSE)
setnames(anx_raw, names(anx_raw), toupper(names(anx_raw)))
anx_raw[, BETA := as.numeric(BETA)]
anx_raw[, SE   := as.numeric(SE)]
anx_raw[, P    := as.numeric(P)]
anx_raw[, EAF  := as.numeric(EAF)]
anx_raw[, N    := ANX_N]
cat(sprintf(  msg：%d SNPs\n, nrow(anx_raw)))

cat(\nmsg...\n)

anx_iv <- anx_raw[P < P_GWS]
p_thresh_anx <- P_GWS
if (nrow(anx_iv) < 3) {
  anx_iv <- anx_raw[P < P_RELAX]
  p_thresh_anx <- P_RELAX
  cat(sprintf(  [WARN] disease_anxiety p<5e-8 msg %d SNPs，msg p<1e-5\n, nrow(anx_iv)))
}
cat(sprintf(  msg（p<%s）：%d SNPs\n,
            formatC(p_thresh_anx, format="e", digits=0), nrow(anx_iv)))

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

exp_anx <- format_data(
  as.data.frame(anx_iv),
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
exp_anx$exposure <- "disease_anxiety"

ibs_as_outcome_raw <- ibs_raw[SNP %in% exp_anx$SNP]
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
exp_ibs$exposure     <- "TRAIT_A"
exp_ibs$eaf.exposure <- 0.5

anx_as_outcome_raw <- anx_raw[SNP %in% exp_ibs$SNP]
out_anx <- format_data(
  as.data.frame(anx_as_outcome_raw),
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
out_anx$outcome            <- "disease_anxiety"
out_anx$samplesize.outcome <- ANX_N

cat(\nmsg...\n)

har_anx_ibs <- harmonise_data(exp_anx, out_ibs, action = 2)
har_anx_ibs <- har_anx_ibs[har_anx_ibs$mr_keep == TRUE, ]
cat(sprintf(  msg→IBS msg：%d SNPs\n, nrow(har_anx_ibs)))

har_ibs_anx <- harmonise_data(exp_ibs, out_anx, action = 2)
har_ibs_anx <- har_ibs_anx[har_ibs_anx$mr_keep == TRUE, ]
cat(sprintf(  IBS→msg msg：%d SNPs\n, nrow(har_ibs_anx)))

fwrite(as.data.table(har_anx_ibs),
       file.path(OUT_DIR, "harmonised_Anxiety_to_IBS.txt"), sep = "\t")
fwrite(as.data.table(har_ibs_anx),
       file.path(OUT_DIR, "harmonised_IBS_to_Anxiety.txt"), sep = "\t")

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

mr_anx_ibs <- run_mr_analysis(har_anx_ibs, "Anxiety → IBS")
mr_ibs_anx <- run_mr_analysis(har_ibs_anx, "IBS → Anxiety")

mr_all <- rbindlist(Filter(Negate(is.null), list(mr_anx_ibs, mr_ibs_anx)))

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
  list(har = har_anx_ibs, label = "Anxiety_to_IBS"),
  list(har = har_ibs_anx, label = "IBS_to_Anxiety")
)) {
  h   <- item$har
  lbl <- item$label
  if (nrow(h) < 3) next

  het <- tryCatch(mr_heterogeneity(h), error = function(e) NULL)
  if (!is.null(het) && nrow(het) > 0) {
    fwrite(as.data.table(het),
           file.path(OUT_DIR, sprintf("heterogeneity_%s.txt", lbl)), sep = "\t")
    ivw_het <- het[het$method == "Inverse variance weighted", ]
    if (nrow(ivw_het) > 0)
      cat(sprintf(  %s msg（IVW Q-pval）：%.3f\n, lbl, ivw_het$Q_pval))
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

  if (!("eaf.exposure" %in% names(har2)) || all(is.na(har2$eaf.exposure)))
    har2$eaf.exposure <- 0.5
  if (!("eaf.outcome" %in% names(har2)) || all(is.na(har2$eaf.outcome)))
    har2$eaf.outcome <- 0.5

  st <- tryCatch(steiger_filtering(har2), error = function(e) {
    cat(sprintf(  [WARN] %s Steiger msg：%s\n, label, conditionMessage(e)))
    NULL
  })
  if (is.null(st)) return(NULL)
  st <- as.data.table(st)
  st[, direction := label]

  n_correct <- sum(st$steiger_dir == TRUE,  na.rm = TRUE)
  n_wrong   <- sum(st$steiger_dir == FALSE, na.rm = TRUE)
  n_total   <- n_correct + n_wrong

  eaf_note <- if (eaf_placeholder) （IBS EAF msg，msg） else ""
  cat(sprintf(  %s：%d/%d SNPs msg（%.0f%%）%s\n,
              label, n_correct, n_total,
              100 * n_correct / max(n_total, 1), eaf_note))
  st
}

st_anx_ibs <- run_steiger(har_anx_ibs, "Anxiety → IBS",
                           n_exp = ANX_N, n_out = IBS_N,
                           eaf_placeholder = TRUE)
st_ibs_anx <- run_steiger(har_ibs_anx, "IBS → Anxiety",
                           n_exp = IBS_N, n_out = ANX_N,
                           eaf_placeholder = TRUE)

steiger_all <- rbindlist(Filter(Negate(is.null), list(st_anx_ibs, st_ibs_anx)))
if (nrow(steiger_all) > 0)
  fwrite(steiger_all,
         file.path(OUT_DIR, "steiger_filtering_results.txt"), sep = "\t")

cat(\nSteiger msg...\n)

mr_steiger_list <- list()

if (!is.null(st_anx_ibs)) {
  keep  <- st_anx_ibs[steiger_dir == TRUE]$SNP
  har_f <- har_anx_ibs[har_anx_ibs$SNP %in% keep, ]
  cat(sprintf(  msg→IBS Steiger msg：%d / %d SNPs msg\n,
              nrow(har_f), nrow(har_anx_ibs)))
  res_f <- run_mr_analysis(har_f, "Anxiety → TRAIT_A (Steiger filtered)")
  if (!is.null(res_f)) mr_steiger_list[["ANX_IBS"]] <- res_f
}

if (!is.null(st_ibs_anx)) {
  keep  <- st_ibs_anx[steiger_dir == TRUE]$SNP
  har_f <- har_ibs_anx[har_ibs_anx$SNP %in% keep, ]
  cat(sprintf(  IBS→msg Steiger msg：%d / %d SNPs msg\n,
              nrow(har_f), nrow(har_ibs_anx)))
  res_f <- run_mr_analysis(har_f, "IBS → disease_anxiety (Steiger filtered)")
  if (!is.null(res_f)) mr_steiger_list[["IBS_ANX"]] <- res_f
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
    "Anxiety → IBS"                    = "#D6604D",
    "IBS → Anxiety"                    = "#2166AC",
    "Anxiety → TRAIT_A (Steiger filtered)" = "#F4A582",
    "IBS → disease_anxiety (Steiger filtered)" = "#92C5DE"
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
         subtitle = Anxiety ↔ IBS；OR>1 msg p<0.05 msg,
         x = "OR (95% CI，log scale)", y = NULL) +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(OUT_DIR, "bidirectional_MR_forest.pdf"),
         p_forest, width = 9, height = max(3, nrow(all_mr_for_plot) * 0.6))
  cat(  msg\n)
}

cat(\n=== msg ===\n)
cat("┌──────────────────────────────────────────────────────────┐\n")
cat(│  msg              OR > 1, p < 0.05  │  msg              │\n)
cat("├────────────────────┤─────────────────────────────────────┤\n")
cat(│  disease_anxiety → TRAIT_A   YES                 │  msg→IBS msg │\n)
cat(│  disease_anxiety → TRAIT_A   NO (p >= 0.05)      │  msg        │\n)
cat(│  TRAIT_A → disease_anxiety   YES                 │  TRAIT_A msg→msg │\n)
cat(│  TRAIT_A → disease_anxiety   NO (p >= 0.05)      │  msg        │\n)
cat("├──────────────────────────────────────────────────────────┤\n")
cat(│  Steiger msg = msg               │\n)
cat(│  TRAIT_A EAF msg（msg 0.5 msg），Steiger R² msg     │\n)
cat("└──────────────────────────────────────────────────────────┘\n")

cat(sprintf(\n=== msg ===\n))
cat(sprintf(msg：%s\n, OUT_DIR))
cat(  bidirectional_MR_results.txt           msg MR msg\n)
cat(  steiger_filtering_results.txt          msg SNP msg Steiger msg\n)
cat(  bidirectional_MR_steiger_filtered.txt  Steiger msg\n)
cat(  bidirectional_MR_forest.pdf            msg\n)
