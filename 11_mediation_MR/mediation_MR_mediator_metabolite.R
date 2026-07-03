#!/usr/bin/env Rscript
# ============================================================
# mediation_MR_PCS.R
# 
# Data sources：
#   TRAIT_A GWAS: gwas_TRAIT_A.txt
#   TRAIT_B GWAS: gwas_TRAIT_B.txt
# ============================================================

suppressMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

IBS_FILE  <- "${PROJECT_ROOT}"
ASD_FILE  <- "${PROJECT_ROOT}"
MET_FILE  <- "${PROJECT_ROOT}"
OUT_DIR   <- "${PROJECT_ROOT}"
PLINK     <- "${PROJECT_ROOT}"
REF_PLINK <- "${PROJECT_ROOT}"
TMP_DIR   <- file.path(OUT_DIR, "tmp")

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP_DIR, showWarnings=FALSE, recursive=TRUE)

IBS_N    <- 486601
ASD_N    <- 46351
P_IBS    <- 1e-5
P_MED    <- 1e-5

cat(msg MR：IBS → msg (met-a-xxx) → ASD\n\n)

cat("Read TRAIT_A GWAS...\n")
ibs_raw <- fread(IBS_FILE, showProgress=FALSE)
setnames(ibs_raw, names(ibs_raw), tolower(names(ibs_raw)))
ibs_raw[, p := 2 * pnorm(-abs(as.numeric(beta) / as.numeric(se)))]
cat(sprintf(  msg：%d\n\n, nrow(ibs_raw)))

cat("Read TRAIT_B GWAS...\n")
asd_raw <- fread(ASD_FILE, showProgress=FALSE)
setnames(asd_raw, names(asd_raw), tolower(names(asd_raw)))
if ("position" %in% names(asd_raw) && !"bp" %in% names(asd_raw))
  setnames(asd_raw, "position", "bp")
cat(sprintf(  msg：%d\n\n, nrow(asd_raw)))

cat(msg GWAS (met-a-718)...\n)
met_raw <- fread(MET_FILE, showProgress=FALSE)
setnames(met_raw, names(met_raw), tolower(names(met_raw)))
if (!"pval" %in% names(met_raw) && "p" %in% names(met_raw))
  setnames(met_raw, "p", "pval")
cat(sprintf(  msg：%d\n\n, nrow(met_raw)))

clump_plink <- function(dt_snp_p, label, p_thresh) {
  tmp_assoc <- file.path(TMP_DIR, paste0(label, "_clump.txt"))
  tmp_out   <- file.path(TMP_DIR, label)
  fwrite(dt_snp_p[, .(SNP, P)], tmp_assoc, sep="\t")
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 %g --clump-r2 0.001 --clump-kb 10000 --out %s --silent',
    PLINK, REF_PLINK, tmp_assoc, p_thresh, tmp_out)
  system(cmd, intern=FALSE)
  clump_file <- paste0(tmp_out, ".clumped")
  if (!file.exists(clump_file)) return(character(0))
  fread(clump_file, showProgress=FALSE)$SNP
}

cat(sprintf(msg TRAIT_A IV（p < %g）msg clumping...\n, P_IBS))
ibs_sig <- ibs_raw[p < P_IBS]
if (nrow(ibs_sig) == 0) stop(IBS GWAS msg SNP msg p < , P_IBS, ，msg)
ibs_for_clump <- ibs_sig[, .(SNP = snp, P = p)]
ibs_clumped_snps <- clump_plink(ibs_for_clump, "TRAIT_A", P_IBS)
if (length(ibs_clumped_snps) == 0) stop(IBS clumping msg：msg SNP。msg plink msg，msg p msg)
cat(sprintf(IBS msg（clumped）：%d msg\n, length(ibs_clumped_snps)))

ibs_iv <- ibs_raw[snp %in% ibs_clumped_snps,
                  .(rsID=snp, P=p, A1=toupper(a1), A2=toupper(a2),
                    beta, se, N=n, eaf=0.5)]
fwrite(ibs_iv, file.path(OUT_DIR, "ibs_instruments_clumped.txt"), sep="\t")

ibs_exp_total <- format_data(
  as.data.frame(ibs_iv),
  type="exposure", snp_col="rsID", beta_col="beta", se_col="se",
  pval_col="P", effect_allele_col="A1", other_allele_col="A2",
  eaf_col="eaf", samplesize_col="N")
ibs_exp_total$exposure <- "TRAIT_A"

cat(\n========== msg MR：IBS → TRAIT_B ==========\n)

asd_for_total <- asd_raw[snp %in% ibs_iv$rsID]
cat(sprintf(  TRAIT_B GWAS msg TRAIT_A IV：%d / %d\n, nrow(asd_for_total), nrow(ibs_iv)))
if (nrow(asd_for_total) < 3) stop(ASD GWAS msg TRAIT_A IV msg 3 msg)

asd_out_total <- format_data(
  as.data.frame(asd_for_total),
  type="outcome", snp_col="snp", beta_col="beta", se_col="se",
  pval_col="p", effect_allele_col="a1", other_allele_col="a2")
asd_out_total$outcome <- "TRAIT_B"

harm_total <- harmonise_data(ibs_exp_total, asd_out_total, action=2)
harm_total  <- harm_total[harm_total$mr_keep, ]
cat(sprintf(  msg SNP：%d\n, nrow(harm_total)))

mr_total <- mr(harm_total,
               method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
cat(\n--- γ msg（IVW）---\n)
print(mr_total[, c("method","nsnp","b","se","pval")])

gamma    <- mr_total$b[grepl("Inverse variance", mr_total$method)]
gamma_se <- mr_total$se[grepl("Inverse variance", mr_total$method)]

pleio_total <- mr_pleiotropy_test(harm_total)
cat(sprintf(  Egger msg = %.4f (p=%.4f)\n\n,
            pleio_total$egger_intercept, pleio_total$pval))
write.csv(mr_total, file.path(OUT_DIR, "00_total_effect_TRAIT_PAIR.csv"), row.names=FALSE)

cat(========== msg MR ==========\n)

# ─── Step α：TRAIT_A IV → metabolite ──────────────────────────────────
cat("Step α: TRAIT_A IV → metabolite_A...\n")
met_for_alpha <- met_raw[snp %in% ibs_iv$rsID]
cat(sprintf(  msg GWAS msg TRAIT_A IV：%d / %d\n, nrow(met_for_alpha), nrow(ibs_iv)))
if (nrow(met_for_alpha) < 3) {
  cat(⚠️  Step α msg（<3），msg，msg。\n)
  quit(status=0)
}

met_out_alpha <- format_data(
  as.data.frame(met_for_alpha),
  type="outcome", snp_col="snp", beta_col="beta", se_col="se",
  pval_col="pval", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf",
  samplesize_col="samplesize")
met_out_alpha$outcome <- "metabolite_A"

harm_alpha <- harmonise_data(ibs_exp_total, met_out_alpha, action=2)
harm_alpha <- harm_alpha[harm_alpha$mr_keep, ]
n_alpha    <- nrow(harm_alpha)
cat(sprintf(  Step α msg SNP：%d\n, n_alpha))
if (n_alpha < 3) {
  cat(⚠️  msg SNP<3，msg。\n)
  quit(status=0)
}

mr_alpha <- mr(harm_alpha,
               method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
alpha    <- mr_alpha$b[grepl("Inverse variance", mr_alpha$method)]
alpha_se <- mr_alpha$se[grepl("Inverse variance", mr_alpha$method)]
alpha_p  <- mr_alpha$pval[grepl("Inverse variance", mr_alpha$method)]
cat(sprintf("  α = %.4f (SE=%.4f, p=%.3g, nSNP=%d)\n",
            alpha, alpha_se, alpha_p, n_alpha))
write.csv(mr_alpha, file.path(OUT_DIR, "alpha_IBS_to_PCS.csv"), row.names=FALSE)

# ─── Step β：metabolite IV → TRAIT_B ─────────────────────────────────
cat("\nStep β: metabolite_A → ASD...\n")
met_sig <- met_raw[pval < P_MED]
cat(sprintf(  msg p < %.0e msg SNP msg：%d\n, P_MED, nrow(met_sig)))
if (nrow(met_sig) < 3) {
  P_MED2 <- 5e-5
  met_sig <- met_raw[pval < P_MED2]
  cat(sprintf(  msg p < %.0e：%d msg\n, P_MED2, nrow(met_sig)))
}
if (nrow(met_sig) < 3) {
  cat(⚠️  msg IV msg，msg Step β。\n)
  quit(status=0)
}

# Clumping metabolite IV
met_clump_snps <- clump_plink(
  data.table(SNP = met_sig$snp, P = met_sig$pval),
  "met_PCS", P_MED)
if (length(met_clump_snps) > 0) {
  met_sig <- met_sig[snp %in% met_clump_snps]
  cat(sprintf(  Clumping msg：%d msg SNP\n, nrow(met_sig)))
}
if (nrow(met_sig) < 3) {
  cat(⚠️  Clumping msg IV msg，msg Step β。\n)
  quit(status=0)
}

met_exp_beta <- format_data(
  as.data.frame(met_sig),
  type="exposure", snp_col="snp", beta_col="beta", se_col="se",
  pval_col="pval", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf",
  samplesize_col="samplesize")
met_exp_beta$exposure <- "metabolite_A"

asd_for_beta <- asd_raw[snp %in% met_sig$snp]
cat(sprintf(  TRAIT_B GWAS msg IV：%d / %d\n, nrow(asd_for_beta), nrow(met_sig)))
if (nrow(asd_for_beta) < 3) {
  cat(⚠️  TRAIT_B msg IV msg，msg。\n)
  quit(status=0)
}

asd_out_beta <- format_data(
  as.data.frame(asd_for_beta),
  type="outcome", snp_col="snp", beta_col="beta", se_col="se",
  pval_col="p", effect_allele_col="a1", other_allele_col="a2")
asd_out_beta$outcome <- "TRAIT_B"

harm_beta <- harmonise_data(met_exp_beta, asd_out_beta, action=2)
harm_beta <- harm_beta[harm_beta$mr_keep, ]
n_beta    <- nrow(harm_beta)
cat(sprintf(  Step β msg SNP：%d\n, n_beta))
if (n_beta < 3) {
  cat(⚠️  SNP<3，msg。\n)
  quit(status=0)
}

mr_beta <- mr(harm_beta,
              method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
beta    <- mr_beta$b[grepl("Inverse variance", mr_beta$method)]
beta_se <- mr_beta$se[grepl("Inverse variance", mr_beta$method)]
beta_p  <- mr_beta$pval[grepl("Inverse variance", mr_beta$method)]
cat(sprintf("  β = %.4f (SE=%.4f, p=%.3g, nSNP=%d)\n",
            beta, beta_se, beta_p, n_beta))
write.csv(mr_beta, file.path(OUT_DIR, "beta_PCS_to_ASD.csv"), row.names=FALSE)

indirect    <- alpha * beta
se_indirect <- sqrt(beta^2 * alpha_se^2 + alpha^2 * beta_se^2)
ci_lo       <- indirect - 1.96 * se_indirect
ci_hi       <- indirect + 1.96 * se_indirect
p_indirect  <- 2 * pnorm(-abs(indirect / se_indirect))
prop_med    <- if (abs(gamma) > 1e-10) indirect / gamma else NA

cat(\n========== msg ==========\n)
cat(sprintf("  α（IBS→metabolite_A）: %.4f (SE=%.4f, p=%.3g)\n", alpha, alpha_se, alpha_p))
cat(sprintf("  β（metabolite_A→ASD）: %.4f (SE=%.4f, p=%.3g)\n", beta, beta_se, beta_p))
cat(sprintf(  msg: %.6f (95%%CI: %.6f ~ %.6f, p=%.3g)\n,
            indirect, ci_lo, ci_hi, p_indirect))
cat(sprintf(  msg: %.1f%%\n, prop_med * 100))

result <- data.frame(
  mediator = "metabolite_A (met-a-718)",
  alpha = round(alpha, 5), alpha_se = round(alpha_se, 5),
  alpha_p = alpha_p, alpha_nsnp = n_alpha,
  beta = round(beta, 5), beta_se = round(beta_se, 5),
  beta_p = beta_p, beta_nsnp = n_beta,
  indirect = round(indirect, 6), indirect_se = round(se_indirect, 6),
  indirect_ci_lo = round(ci_lo, 6), indirect_ci_hi = round(ci_hi, 6),
  indirect_p = p_indirect,
  gamma = round(gamma, 5), gamma_se = round(gamma_se, 5),
  prop_mediated = round(prop_med, 4),
  stringsAsFactors = FALSE
)
write.csv(result, file.path(OUT_DIR, "01_mediation_summary.csv"), row.names = FALSE)

p_forest <- ggplot(result, aes(x = indirect, y = mediator)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = indirect_ci_lo, xmax = indirect_ci_hi), height = 0.2) +
  geom_point(size = 3, color = "#D6604D") +
  theme_bw(base_size = 12) +
  labs(title = msg MR：IBS → metabolite_A (met-a-718) → ASD,
       subtitle = sprintf(msg γ = %.4f (p=%.2g)；msg %.1f%%,
                          gamma, mr_total$pval[grepl("Inverse variance", mr_total$method)],
                          prop_med * 100),
       x = msg (β), y = NULL)
ggsave(file.path(OUT_DIR, "02_mediation_forest.pdf"), p_forest,
       width = 7, height = 3)
cat(\nmsg。\n)

cat(\n=== msg ===\nmsg：, OUT_DIR, "\n")