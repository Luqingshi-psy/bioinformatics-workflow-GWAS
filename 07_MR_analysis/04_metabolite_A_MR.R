#!/usr/bin/env Rscript
# ============================================================
# 08_metabolite_B_metabolite_mr.R  （metabolite_Bmetabolite met-a-476 data version）
# metabolite → trait pair (A x B) Two-sample MR
#
# Exposure：metabolite_Bmetabolite（Butyrate, met-a-476, Shin 2014, N=7824）
# Outcome：local trait pair (A x B) GWAS
# LD clumping：plink1.9（local，${PROJECT_ROOT}
#
# dependency：TwoSampleMR, data.table, ggplot2
#   remotes::install_github("MRCIEU/TwoSampleMR")
# ============================================================
install.packages("${PROJECT_ROOT}", repos = NULL, type = "source")
install.packages("${PROJECT_ROOT}", repos = NULL, type = "source")
install.packages("${PROJECT_ROOT}", repos = NULL, type = "source")
install.packages("${PROJECT_ROOT}", repos = NULL, type = "source")
install.packages("${PROJECT_ROOT}", repos = NULL, type = "source")

library(TwoSampleMR)
library(data.table)
library(ggplot2)

# ── pathconfig ─────────────────────────────────────────────────
METABOLITE_VCF <- "${PROJECT_ROOT}"
GWAS_DIR       <- "${PROJECT_ROOT}"
OUT_DIR        <- "${PROJECT_ROOT}"
PLINK          <- "${PROJECT_ROOT}"
REF_PLINK      <- "${PROJECT_ROOT}"
TMP_DIR        <- file.path(OUT_DIR, "tmp")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

IBS_N <- 486601
ASD_N <- 46351

# ════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════
parsed_file <- file.path(TMP_DIR, "met_a476_parsed.tsv")
if(!file.exists(parsed_file)) {
  stop(msg met_a476_parsed.tsv)
}

cat(msg...\n)
met <- fread(parsed_file, sep = "\t", header = TRUE)
cat(sprintf(  msg SNP msg：%d\n, nrow(met)))

# compute P value = 10^(-LP)
met[, pval := 10^(-LP)]

P_THRESHOLD <- 5e-8
vcf_sig <- met[pval < P_THRESHOLD]
cat(sprintf(\nmsg p < 5e-8，msg %d msg SNP\n, nrow(vcf_sig)))

if (nrow(vcf_sig) == 0) {
  stop(msg p<5e-8 msg SNP，msg)
}

METABOLITE_N <- 7824
instruments_raw <- data.table(
  SNP      = vcf_sig$ID,
  chr      = as.character(vcf_sig$CHROM),
  bp       = as.integer(vcf_sig$POS),
  A1       = toupper(vcf_sig$ALT),
  A2       = toupper(vcf_sig$REF),
  beta     = vcf_sig$ES,
  se       = vcf_sig$SE,
  pval     = vcf_sig$pval,
  eaf      = vcf_sig$AF,
  N        = METABOLITE_N,
  exposure = "Butyrate"
)

cat(sprintf(msg SNP msg (clump msg)：%d\n, nrow(instruments_raw)))

# ════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════
cat(\nplink LD msg...\n)

clump_with_plink <- function(dt, exposure_name, plink_bin, ref, tmp_dir) {
  tmp_assoc <- file.path(tmp_dir, paste0(exposure_name, "_clump_input.txt"))
  tmp_out   <- file.path(tmp_dir, exposure_name)
  fwrite(dt[, .(SNP, P = pval)], tmp_assoc, sep = "\t")
  
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000 --out %s --silent',
    plink_bin, ref, tmp_assoc, tmp_out
  )
  ret <- system(cmd, intern = FALSE)
  
  clump_file <- paste0(tmp_out, ".clumped")
  if (!file.exists(clump_file)) {
    cat(sprintf(  [WARN] %s：clump msg\n, exposure_name))
    return(NULL)
  }
  clumped <- fread(clump_file, showProgress = FALSE)
  keep_snps <- clumped$SNP
  dt_keep <- dt[SNP %in% keep_snps]
  cat(sprintf(  %s：%d → %d SNP（clump msg）\n, exposure_name, nrow(dt), nrow(dt_keep)))
  return(dt_keep)
}

instruments_clumped <- clump_with_plink(instruments_raw, "Butyrate", PLINK, REF_PLINK, TMP_DIR)

if (is.null(instruments_clumped) || nrow(instruments_clumped) == 0) {
  stop(clump msg，msg plink msg LD msg)
}
cat(sprintf(\nclump msg SNP：%d\n, nrow(instruments_clumped)))
fwrite(instruments_clumped, file.path(OUT_DIR, "instruments_clumped.txt"), sep = "\t")

cat(\nmsg...\n)
exposure_dat <- format_data(
  as.data.frame(instruments_clumped),
  type                 = "exposure",
  snp_col              = "SNP",
  beta_col             = "beta",
  se_col               = "se",
  pval_col             = "pval",
  effect_allele_col    = "A1",
  other_allele_col     = "A2",
  eaf_col              = "eaf",
  samplesize_col       = "N",
  phenotype_col        = "exposure"
)
cat(sprintf(  msg：%d msg\n, nrow(exposure_dat)))

cat(\nmsg GWAS msg...\n)

read_local_outcome <- function(file, outcome_name, n_total,
                               snp_col, beta_col, se_col, pval_col,
                               ea_col, oa_col, eaf_col = NULL, compute_p = FALSE) {
  dt <- fread(file, showProgress = FALSE)
  if (compute_p) {
    dt[, pval_tmp := 2 * pnorm(-abs(as.numeric(get(beta_col)) / as.numeric(get(se_col))))]
    pval_col <- "pval_tmp"
  }
  if (!is.null(eaf_col) && eaf_col %in% names(dt)) {
    dt[, eaf_tmp := as.numeric(get(eaf_col))]
    eaf_col_use <- "eaf_tmp"
  } else {
    dt[, eaf_tmp := 0.5]
    eaf_col_use <- "eaf_tmp"
  }
  
  iv_snps <- unique(instruments_clumped$SNP)
  dt_sub <- dt[get(snp_col) %in% iv_snps]
  cat(sprintf(  %s：msg %d / %d msg SNP\n, outcome_name, nrow(dt_sub), length(iv_snps)))
  if (nrow(dt_sub) == 0) return(NULL)
  
  out <- format_data(
    as.data.frame(dt_sub),
    type                 = "outcome",
    snp_col              = snp_col,
    beta_col             = beta_col,
    se_col               = se_col,
    pval_col             = pval_col,
    effect_allele_col    = ea_col,
    other_allele_col     = oa_col,
    eaf_col              = eaf_col_use
  )
  out$outcome           <- outcome_name
  out$samplesize.outcome <- n_total
  as.data.table(out)
}

outcome_ibs <- read_local_outcome(
  file.path(GWAS_DIR, "gwas_TRAIT_A.txt"),
  outcome_name = "TRAIT_A", n_total = IBS_N,
  snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  pval_col = NULL, ea_col = "A1", oa_col = "A2",
  compute_p = TRUE
)
outcome_asd <- read_local_outcome(
  file.path(GWAS_DIR, "gwas_TRAIT_B.txt"),
  outcome_name = "TRAIT_B", n_total = ASD_N,
  snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  pval_col = "P", ea_col = "A1", oa_col = "A2", eaf_col = "EAF"
)

outcomes <- rbindlist(Filter(Negate(is.null), list(outcome_ibs, outcome_asd)))
if (nrow(outcomes) == 0) stop(msg)

cat(\nmsg...\n)
harmonised <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat  = as.data.frame(outcomes),
  action       = 2
)
cat(sprintf(msg：%d msg（%d msg × msg）\n,
            nrow(harmonised), length(unique(paste(harmonised$exposure, harmonised$outcome)))))
fwrite(as.data.table(harmonised), file.path(OUT_DIR, "harmonised_data.txt"), sep = "\t")

harmonised <- harmonised[harmonised$mr_keep == TRUE, ]
cat(sprintf(mr_keep==TRUE：%d msg\n, nrow(harmonised)))

# ── 6. MR analysis ────────────────────────────────────────────────
cat(\nMR msg...\n)

mr_results <- mr(harmonised, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode"
))
mr_dt <- as.data.table(mr_results)
mr_dt[, `:=`(OR = exp(b), CI_lo = exp(b - 1.96*se), CI_hi = exp(b + 1.96*se))]
setorder(mr_dt, outcome, exposure, method)
fwrite(mr_dt, file.path(OUT_DIR, "MR_results_all.txt"), sep = "\t")

cat(\n=== IVW msg ===\n)
ivw <- mr_dt[method == "Inverse variance weighted"][order(pval)]
print(ivw[, .(exposure, outcome, nsnp,
              b = round(b, 4), OR = round(OR, 3),
              CI_lo = round(CI_lo, 3), CI_hi = round(CI_hi, 3),
              pval = formatC(pval, format = "e", digits = 2))], row.names = FALSE)

cat(\nmsg（MR-Egger msg）...\n)
pleio <- mr_pleiotropy_test(harmonised)
if (!is.null(pleio) && nrow(pleio) > 0) {
  print(as.data.frame(pleio[, c("exposure","outcome","egger_intercept","se","pval")]),
        row.names = FALSE)
  fwrite(as.data.table(pleio), file.path(OUT_DIR, "MR_pleiotropy_test.txt"), sep = "\t")
}

het <- mr_heterogeneity(harmonised)
if (!is.null(het) && nrow(het) > 0) {
  fwrite(as.data.table(het), file.path(OUT_DIR, "MR_heterogeneity.txt"), sep = "\t")
}

overlap_ibs <- instruments_clumped$SNP %in% outcome_ibs$SNP
overlap_asd <- instruments_clumped$SNP %in% outcome_asd$SNP
cat(IBS msg SNP msg：, sum(overlap_ibs), "/", nrow(instruments_clumped), "\n")
cat(ASD msg SNP msg：, sum(overlap_asd), "/", nrow(instruments_clumped), "\n")
cat(msg SNP：\n)
print(instruments_clumped$SNP[!overlap_ibs | !overlap_asd])

# ── 8. visualisation ─────────────────────────────────────────────────
cat(\nmsg...\n)

if (nrow(ivw) > 0) {
  ivw[, sig := pval < 0.05]
  p_forest <- ggplot(ivw, aes(x = OR, 
                              y = reorder(paste0(exposure, "\n→", outcome), OR),
                              color = outcome, shape = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.5) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "p<0.05", "FALSE" = "p≥0.05")) +
    scale_color_manual(values = c("TRAIT_A" = "#2166AC", "TRAIT_B" = "#D6604D")) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7)) +
    labs(title = msg → msg MR（IVW）,
         x = "OR (95% CI，log scale)", y = NULL,
         color = msg, shape = msg)
  ggsave(file.path(OUT_DIR, "MR_forest_IVW.pdf"), p_forest,
         width = 10, height = max(4, nrow(ivw) * 0.5))
}

sig_pairs <- unique(mr_dt[method == "Inverse variance weighted" & pval < 0.05,
                          paste(exposure, outcome)])
if (length(sig_pairs) > 0) {
  har_sig <- harmonised[paste(harmonised$exposure, harmonised$outcome) %in% sig_pairs, ]
  if (nrow(har_sig) > 0) {
    res_sig <- mr(har_sig, method_list = c("mr_ivw", "mr_egger_regression",
                                           "mr_weighted_median"))
    p_scatter <- mr_scatter_plot(res_sig, har_sig)
    pdf(file.path(OUT_DIR, "MR_scatter_significant.pdf"), width = 7, height = 6)
    for (plt in p_scatter) print(plt)
    dev.off()
    cat(  msg（msg）\n)
  }
}

cat(\n=== msg MR msg ===\n)
cat(msg：, OUT_DIR, "\n")
cat(\nmsg：\n)
cat(  instruments_clumped.txt   clump msg\n)
cat(  MR_results_all.txt        msg × msg\n)
cat(  MR_pleiotropy_test.txt    msg（msg≠0=msg）\n)
cat(  MR_forest_IVW.pdf         msg\n)
cat(  MR_scatter_significant.pdf msg（p<0.05msg）\n)