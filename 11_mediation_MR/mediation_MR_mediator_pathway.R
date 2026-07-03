#!/usr/bin/env Rscript
# ============================================================
# 09_mediation_MR_metabolite_C.R  v3
#
# Revision notes（v3）：
# ============================================================

suppressMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

MIBIOGEN_ALLHITS <- "${PROJECT_ROOT}"
MIBIOGEN_DIR     <- "${PROJECT_ROOT}"
IBS_FILE         <- "${PROJECT_ROOT}"
ASD_FILE         <- "${PROJECT_ROOT}"
OUT_DIR          <- "${PROJECT_ROOT}"
PLINK            <- "${PROJECT_ROOT}"
REF_PLINK        <- "${PROJECT_ROOT}"
TMP_DIR          <- file.path(OUT_DIR, "tmp")
IBS_CLUMPED_FILE <- file.path(TMP_DIR, "IBS.clumped")

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP_DIR, showWarnings=FALSE, recursive=TRUE)

IBS_N  <- 486601
ASD_N  <- 46351
P_IBS  <- 1e-5
P_MBIO <- 1e-5

EPS_TAXA <- list(
  Peptostreptococcaceae  = list(
    id  = "family.Peptostreptococcaceae.id.2042",
    zip = "MiBioGen_QmbQTL_summary_family.zip",
    gz  = "family.Peptostreptococcaceae.id.2042.summary.txt.gz"),
  Clostridiales_order    = list(
    id  = "order.Clostridiales.id.1863",
    zip = "MiBioGen_QmbQTL_summary_order.zip",
    gz  = "order.Clostridiales.id.1863.summary.txt.gz"),
  Lachnospiraceae        = list(
    id  = "family.Lachnospiraceae.id.1987",
    zip = "MiBioGen_QmbQTL_summary_family.zip",
    gz  = "family.Lachnospiraceae.id.1987.summary.txt.gz"),
  Ruminococcaceae        = list(
    id  = "family.Ruminococcaceae.id.2050",
    zip = "MiBioGen_QmbQTL_summary_family.zip",
    gz  = "family.Ruminococcaceae.id.2050.summary.txt.gz"),
  Erysipelotrichaceae    = list(
    id  = "family.Erysipelotrichaceae.id.2149",
    zip = "MiBioGen_QmbQTL_summary_family.zip",
    gz  = "family.Erysipelotrichaceae.id.2149.summary.txt.gz"),
  Blautia                = list(
    id  = "genus.Blautia.id.1992",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Blautia.id.1992.summary.txt.gz"),
  Clostridium_s1         = list(
    id  = "genus.Clostridiumsensustricto1.id.1873",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Clostridiumsensustricto1.id.1873.summary.txt.gz"),
  Ruminococcus1          = list(
    id  = "genus.Ruminococcus1.id.11373",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Ruminococcus1.id.11373.summary.txt.gz"),
  Ruminococcus2          = list(
    id  = "genus.Ruminococcus2.id.11374",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Ruminococcus2.id.11374.summary.txt.gz"),
  RuminococcusGnavus     = list(
    id  = "genus..Ruminococcusgnavusgroup.id.14376",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus..Ruminococcusgnavusgroup.id.14376.summary.txt.gz"),
  RuminococcusTorques    = list(
    id  = "genus..Ruminococcustorquesgroup.id.14377",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus..Ruminococcustorquesgroup.id.14377.summary.txt.gz"),
  Lachnoclostridium      = list(
    id  = "genus.Lachnoclostridium.id.11308",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Lachnoclostridium.id.11308.summary.txt.gz"),
  Erysipelatoclostridium = list(
    id  = "genus.Erysipelatoclostridium.id.11381",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Erysipelatoclostridium.id.11381.summary.txt.gz"),
  RuminococcaceaeUCG002  = list(
    id  = "genus.RuminococcaceaeUCG002.id.11360",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.RuminococcaceaeUCG002.id.11360.summary.txt.gz"),
  RuminococcaceaeUCG013  = list(
    id  = "genus.RuminococcaceaeUCG013.id.11370",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.RuminococcaceaeUCG013.id.11370.summary.txt.gz"),
  Ruminiclostridium5     = list(
    id  = "genus.Ruminiclostridium5.id.11355",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Ruminiclostridium5.id.11355.summary.txt.gz"),
  Ruminiclostridium9     = list(
    id  = "genus.Ruminiclostridium9.id.11357",
    zip = "MiBioGen_QmbQTL_summary_genus.zip",
    gz  = "genus.Ruminiclostridium9.id.11357.summary.txt.gz")
)
cat(sprintf(msg：%d msg\n\n, length(EPS_TAXA)))

extract_snps_from_zip_gz <- function(zip_path, gz_name, snp_vec, tmp_out) {
  snp_json <- paste0('["', paste(snp_vec, collapse='","'), '"]')
  py_script <- sprintf('
import zipfile, gzip, sys
zip_path = "%s"
gz_name  = "%s"
snps = set(%s)
out_path = "%s"
with zipfile.ZipFile(zip_path) as zf:
    with zf.open(gz_name) as gz_bytes:
        with gzip.open(gz_bytes) as f:
            with open(out_path, "w") as out:
                header = f.readline().decode("utf-8", errors="replace")
                out.write(header)
                for line in f:
                    row = line.decode("utf-8", errors="replace")
                    fields = row.split("\\t")
                    if len(fields) > 3 and fields[3] in snps:
                        out.write(row)
', zip_path, gz_name, snp_json, tmp_out)

  py_file <- tempfile(fileext=".py")
  writeLines(py_script, py_file)
  ret <- system(sprintf("python3 %s", py_file), intern=FALSE)
  unlink(py_file)
  if (ret != 0 || !file.exists(tmp_out) || file.size(tmp_out) < 10) return(NULL)
  dt <- tryCatch(fread(tmp_out, showProgress=FALSE), error=function(e) NULL)
  dt
}

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

cat(msg MiBioGen allHits（Step β msg）...\n)
mbg <- fread(MIBIOGEN_ALLHITS, quote='"', showProgress=FALSE)
setnames(mbg, names(mbg),
         c("bac","chr","bp","rsID","ref_allele","eff_allele",
           "beta","SE","Z","P","N","Ncohorts"))
mbg[, bac        := gsub('"', '', bac)]
mbg[, rsID       := gsub('"', '', rsID)]
mbg[, eff_allele := toupper(gsub('"', '', eff_allele))]
mbg[, ref_allele := toupper(gsub('"', '', ref_allele))]
mbg <- mbg[!is.na(rsID) & rsID != "NA"]
cat(sprintf(  msg：%d\n\n, nrow(mbg)))

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

if (file.exists(IBS_CLUMPED_FILE)) {
  cat(msg IBS.clumped...\n)
  ibs_clumped_snps <- fread(IBS_CLUMPED_FILE, showProgress=FALSE)$SNP
} else {
  cat(sprintf(msg TRAIT_A IV（p<%g）msg clumping...\n, P_IBS))
  ibs_sig  <- ibs_raw[p < P_IBS]
  ibs_for_clump <- ibs_sig[, .(SNP=snp, P=p)]
  ibs_clumped_snps <- clump_plink(ibs_for_clump, "TRAIT_A", P_IBS)
  if (length(ibs_clumped_snps) == 0) stop(IBS clumping msg)
}
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

if (nrow(asd_for_total) < 3) stop(ASD GWAS msg TRAIT_A IV msg 3 msg，msg SNP msg)

asd_out_total <- format_data(
  as.data.frame(asd_for_total),
  type="outcome", snp_col="snp", beta_col="beta", se_col="se",
  pval_col="p", effect_allele_col="a1", other_allele_col="a2",
  eaf_col=if("eaf" %in% names(asd_for_total)) "eaf" else NULL)
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
mediation_rows <- list()

for (i in seq_along(EPS_TAXA)) {
  tn  <- names(EPS_TAXA)[i]
  info <- EPS_TAXA[[i]]
  cat(sprintf("\n[%d/%d] %s\n", i, length(EPS_TAXA), tn))

  zip_path <- file.path(MIBIOGEN_DIR, info$zip)

  cat(  Step α: msg summary msg TRAIT_A IV msg...\n)
  tmp_alpha <- file.path(TMP_DIR, sprintf("alpha_lookup_%s.txt", tn))

  alpha_dt <- extract_snps_from_zip_gz(
    zip_path  = zip_path,
    gz_name   = info$gz,
    snp_vec   = ibs_iv$rsID,
    tmp_out   = tmp_alpha)

  if (is.null(alpha_dt) || nrow(alpha_dt) == 0) {
    cat(  ⚠️  msg summary msg TRAIT_A IV，msg\n); next
  }
  setnames(alpha_dt, names(alpha_dt),
           c("bac","chr","bp","rsID","ref_allele","eff_allele",
             "beta","SE","Z","P","N","Ncohorts"))
  alpha_dt[, eff_allele := toupper(eff_allele)]
  alpha_dt[, ref_allele := toupper(ref_allele)]
  cat(sprintf(  msg TRAIT_A IV msg summary msg：%d / %d msg\n,
              nrow(alpha_dt), nrow(ibs_iv)))

  mbio_out_alpha <- format_data(
    as.data.frame(alpha_dt),
    type="outcome", snp_col="rsID", beta_col="beta", se_col="SE",
    pval_col="P", effect_allele_col="eff_allele", other_allele_col="ref_allele")
  mbio_out_alpha$outcome <- tn

  harm_alpha <- tryCatch(
    harmonise_data(ibs_exp_total, mbio_out_alpha, action=2),
    error=function(e) { cat(  ❌ Step α msg\n); NULL })
  if (is.null(harm_alpha)) next
  harm_alpha <- harm_alpha[harm_alpha$mr_keep, ]
  n_alpha    <- nrow(harm_alpha)
  cat(sprintf(  Step α msg SNP：%d\n, n_alpha))
  if (n_alpha < 3) { cat(  ⚠️  msg SNP < 3，msg\n); next }

  mr_alpha <- tryCatch(
    mr(harm_alpha, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression")),
    error=function(e) { cat(  ❌ Step α MR msg\n); NULL })
  if (is.null(mr_alpha)) next

  alpha    <- mr_alpha$b[grepl("Inverse variance", mr_alpha$method)]
  alpha_se <- mr_alpha$se[grepl("Inverse variance", mr_alpha$method)]
  alpha_p  <- mr_alpha$pval[grepl("Inverse variance", mr_alpha$method)]
  cat(sprintf("  α = %.4f (SE=%.4f, p=%.3g, nSNP=%d)\n",
              alpha, alpha_se, alpha_p, n_alpha))

  cat(  Step β: msg → ASD...\n)
  mbio_taxon <- mbg[bac == info$id]
  cat(sprintf(  allHits msg SNP msg：%d\n, nrow(mbio_taxon)))

  mbio_sig <- mbio_taxon[P < P_MBIO]
  if (nrow(mbio_sig) < 3) {
    mbio_sig <- mbio_taxon[P < 5e-5]
    cat(sprintf(  msg p<5e-5：%d msg\n, nrow(mbio_sig)))
  }
  if (nrow(mbio_sig) < 3) { cat(  ⚠️  msg IV msg，msg\n); next }

  # plink clumping（microbiota IV）
  keep_mbio <- clump_plink(
    data.table(SNP=mbio_sig$rsID, P=mbio_sig$P),
    paste0("mbio_", tn), P_MBIO)
  if (length(keep_mbio) > 0) {
    mbio_sig <- mbio_sig[rsID %in% keep_mbio]
    cat(sprintf(  Clumping msg：%d msg\n, nrow(mbio_sig)))
  }
  if (nrow(mbio_sig) < 3) { cat(  ⚠️  Clumping msg IV msg，msg\n); next }

  mbio_exp_beta <- format_data(
    as.data.frame(mbio_sig),
    type="exposure", snp_col="rsID", beta_col="beta", se_col="SE",
    pval_col="P", effect_allele_col="eff_allele", other_allele_col="ref_allele")
  mbio_exp_beta$exposure <- tn

  asd_for_beta <- asd_raw[snp %in% mbio_sig$rsID]
  cat(sprintf(  TRAIT_B GWAS msg IV：%d / %d\n, nrow(asd_for_beta), nrow(mbio_sig)))
  if (nrow(asd_for_beta) < 3) { cat(  ⚠️  msg，msg\n); next }

  asd_out_beta <- format_data(
    as.data.frame(asd_for_beta),
    type="outcome", snp_col="snp", beta_col="beta", se_col="se",
    pval_col="p", effect_allele_col="a1", other_allele_col="a2",
    eaf_col=if("eaf" %in% names(asd_for_beta)) "eaf" else NULL)
  asd_out_beta$outcome <- "TRAIT_B"

  harm_beta <- tryCatch(
    harmonise_data(mbio_exp_beta, asd_out_beta, action=2),
    error=function(e) { cat(  ❌ Step β msg\n); NULL })
  if (is.null(harm_beta)) next
  harm_beta <- harm_beta[harm_beta$mr_keep, ]
  n_beta    <- nrow(harm_beta)
  cat(sprintf(  Step β msg SNP：%d\n, n_beta))
  if (n_beta < 3) { cat(  ⚠️  SNP < 3，msg\n); next }

  mr_beta <- tryCatch(
    mr(harm_beta, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression")),
    error=function(e) { cat(  ❌ Step β MR msg\n); NULL })
  if (is.null(mr_beta)) next

  beta    <- mr_beta$b[grepl("Inverse variance", mr_beta$method)]
  beta_se <- mr_beta$se[grepl("Inverse variance", mr_beta$method)]
  beta_p  <- mr_beta$pval[grepl("Inverse variance", mr_beta$method)]
  cat(sprintf("  β = %.4f (SE=%.4f, p=%.3g, nSNP=%d)\n",
              beta, beta_se, beta_p, n_beta))

  indirect    <- alpha * beta
  se_indirect <- sqrt(beta^2 * alpha_se^2 + alpha^2 * beta_se^2)
  ci_lo       <- indirect - 1.96 * se_indirect
  ci_hi       <- indirect + 1.96 * se_indirect
  p_indirect  <- 2 * pnorm(-abs(indirect / se_indirect))
  prop_med    <- if (abs(gamma) > 1e-10) indirect / gamma else NA

  cat(sprintf(  ★ msg = %.5f (95%%CI: %.5f~%.5f, p=%.3g)\n,
              indirect, ci_lo, ci_hi, p_indirect))
  cat(sprintf(    msg = %.1f%%\n, prop_med * 100))

  mediation_rows[[tn]] <- data.frame(
    taxon=tn, taxon_id=info$id,
    alpha=round(alpha,5), alpha_se=round(alpha_se,5),
    alpha_p=alpha_p, alpha_nsnp=n_alpha,
    beta=round(beta,5), beta_se=round(beta_se,5),
    beta_p=beta_p, beta_nsnp=n_beta,
    indirect=round(indirect,6), indirect_se=round(se_indirect,6),
    indirect_ci_lo=round(ci_lo,6), indirect_ci_hi=round(ci_hi,6),
    indirect_p=p_indirect,
    gamma=round(gamma,5), gamma_se=round(gamma_se,5),
    prop_mediated=round(prop_med,4),
    stringsAsFactors=FALSE)

  write.csv(mr_alpha, file.path(OUT_DIR, sprintf("alpha_IBS_to_%s.csv", tn)), row.names=FALSE)
  write.csv(mr_beta,  file.path(OUT_DIR, sprintf("beta_%s_to_ASD.csv",  tn)), row.names=FALSE)
}

# ── 9. summaryoutput ───────────────────────────────────────────────
if (length(mediation_rows) == 0) {
  cat(\n⚠️  msg。msg：\n)
  cat(   1. TRAIT_A IV msg summary gz msg（SNP msg GWAS panel）\n)
  cat(   2. msg allHits msg SNP（Step β IV msg）\n)
  cat(   msg tmp/ msg alpha_lookup_*.txt msg\n)
} else {
  res <- bind_rows(mediation_rows)
  res$p_bonferroni <- p.adjust(res$indirect_p, method="bonferroni")
  res$p_fdr        <- p.adjust(res$indirect_p, method="fdr")
  res <- res[order(res$indirect_p), ]

  cat(\n========== msg MR msg ==========\n)
  print(res[, c("taxon","alpha","alpha_p","beta","beta_p",
                "indirect","indirect_ci_lo","indirect_ci_hi",
                "indirect_p","p_fdr","prop_mediated")],
        row.names=FALSE, digits=3)

  write.csv(res, file.path(OUT_DIR, "01_mediation_summary.csv"), row.names=FALSE)

  if (any(res$indirect_p < 0.05, na.rm=TRUE)) {
    sig <- res[!is.na(res$indirect_p) & res$indirect_p < 0.05, ]
    sig$sig_label <- ifelse(sig$p_fdr < 0.05, "FDR<0.05", "p<0.05")
    p_forest <- ggplot(sig, aes(x=indirect, y=reorder(taxon,-indirect_p), color=sig_label)) +
      geom_vline(xintercept=0, linetype="dashed", color="gray50") +
      geom_errorbarh(aes(xmin=indirect_ci_lo, xmax=indirect_ci_hi), height=0.25) +
      geom_point(size=3) +
      scale_color_manual(values=c("FDR<0.05"="#D6604D","p<0.05"="#F4A582")) +
      theme_bw(base_size=10) +
      labs(title=msg MR：IBS → metabolite_Cmsg → ASD,
           subtitle=sprintf(msg γ=%.4f (p=%.2g), gamma,
                            mr_total$pval[mr_total$method=="mr_ivw"]),
           x=msg (β), y=NULL, color=msg)
    ggsave(file.path(OUT_DIR, "02_mediation_forest.pdf"), p_forest,
           width=9, height=max(4, nrow(sig)*0.55))
    cat(\nmsg\n)
  }
}

cat(\n=== msg ===\nmsg：, OUT_DIR, "\n")
