#!/usr/bin/env Rscript
# ============================================================
# 07_microbiota_mr.R  （local MiBioGen data version）
# microbiota → trait pair (A x B) Two-sample MR
#
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
MIBIOGEN_FILE <- "${PROJECT_ROOT}"
GWAS_DIR      <- "${PROJECT_ROOT}"
OUT_DIR       <- "${PROJECT_ROOT}"
PLINK         <- "${PROJECT_ROOT}"
REF_PLINK     <- "${PROJECT_ROOT}"
TMP_DIR       <- file.path(OUT_DIR, "tmp")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP_DIR, showWarnings=FALSE, recursive=TRUE)

IBS_N <- 486601
ASD_N <- 46351

TARGET_TAXA <- list(
  Faecalibacterium   = "genus.Faecalibacterium.id.2057",
  Roseburia          = "genus.Roseburia.id.2012",
  Bifidobacterium    = "genus.Bifidobacterium.id.436",
  Coprococcus1       = "genus.Coprococcus1.id.11301",
  Coprococcus2       = "genus.Coprococcus2.id.11302",
  Coprococcus3       = "genus.Coprococcus3.id.11303",
  Lachnospira        = "genus.Lachnospira.id.2004",
  RuminococcaceaeUCG002 = "genus.RuminococcaceaeUCG002.id.11360"
)

cat(msg MiBioGen allHits msg...\n)
mbg <- fread(MIBIOGEN_FILE, quote='"', showProgress=FALSE)
setnames(mbg, names(mbg),
         c("bac","chr","bp","rsID","ref_allele","eff_allele","beta","SE","Z","P","N","Ncohorts"))
mbg[, bac := gsub('"', '', bac)]
cat(sprintf(  msg：%d\n, nrow(mbg)))

cat(\nmsg...\n)

p_thresh_strict <- 5e-8
p_thresh_relax  <- 1e-5

extract_instruments_local <- function(taxon_id, taxon_name, p_thresh) {
  dt <- mbg[bac == taxon_id & P < p_thresh]
  if (nrow(dt) == 0) {
    cat(sprintf(  [SKIP] %s：p<%s msg SNP\n, taxon_name, formatC(p_thresh,format="e",digits=0)))
    return(NULL)
  }
  dt2 <- data.table(
    SNP    = gsub('"', '', dt$rsID),
    chr    = dt$chr,
    bp     = dt$bp,
    A1     = toupper(gsub('"', '', dt$eff_allele)),
    A2     = toupper(gsub('"', '', dt$ref_allele)),
    beta   = dt$beta,
    se     = dt$SE,
    pval   = dt$P,
    N      = dt$N,
    exposure = taxon_name
  )
  cat(sprintf(  %s：%d msg SNP（p<%s）\n, taxon_name, nrow(dt2),
              formatC(p_thresh,format="e",digits=0)))
  dt2
}

all_instruments <- list()
for (nm in names(TARGET_TAXA)) {
  tid <- TARGET_TAXA[[nm]]
  p_use <- if (nm == "Bifidobacterium") p_thresh_strict else p_thresh_relax
  res <- extract_instruments_local(tid, nm, p_use)
  if (!is.null(res)) all_instruments[[nm]] <- res
}

if (length(all_instruments) == 0) stop(msg SNP，msg)
instruments_raw <- rbindlist(all_instruments)
cat(sprintf(\nmsg SNP：%d（msg %d msg）\n,
            nrow(instruments_raw), length(all_instruments)))

cat(\nplink LD msg...\n)

clump_with_plink <- function(dt, taxon_name, plink_bin, ref, tmp_dir) {
  tmp_assoc <- file.path(tmp_dir, paste0(taxon_name, "_clump_input.txt"))
  tmp_out   <- file.path(tmp_dir, taxon_name)
  fwrite(dt[, .(SNP, P=pval)], tmp_assoc, sep="\t")

  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 %s --clump-r2 0.001 --clump-kb 10000 --out %s --silent',
    plink_bin, ref, tmp_assoc,
    ifelse(taxon_name=="Bifidobacterium", "5e-8", "1e-5"),
    tmp_out
  )
  ret <- system(cmd, intern=FALSE)
  clump_file <- paste0(tmp_out, ".clumped")
  if (!file.exists(clump_file)) {
    cat(sprintf(  [WARN] %s：clump msg\n, taxon_name))
    return(NULL)
  }
  clumped <- fread(clump_file, showProgress=FALSE)
  keep_snps <- clumped$SNP
  dt_keep <- dt[SNP %in% keep_snps]
  cat(sprintf(  %s：%d → %d SNP（clump msg）\n, taxon_name, nrow(dt), nrow(dt_keep)))
  dt_keep
}

all_clumped <- list()
for (nm in names(all_instruments)) {
  res <- clump_with_plink(all_instruments[[nm]], nm, PLINK, REF_PLINK, TMP_DIR)
  if (!is.null(res) && nrow(res) > 0) all_clumped[[nm]] <- res
}

if (length(all_clumped) == 0) stop(clump msg，msg plink msg LD msg)
instruments_clumped <- rbindlist(all_clumped)
cat(sprintf(\nclump msg SNP：%d\n, nrow(instruments_clumped)))
fwrite(instruments_clumped, file.path(OUT_DIR, "instruments_clumped.txt"), sep="\t")

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
  samplesize_col       = "N",
  phenotype_col        = "exposure"
)
cat(sprintf(  msg：%d msg\n, nrow(exposure_dat)))

cat(\nmsg GWAS msg...\n)

read_local_outcome <- function(file, outcome_name, n_total,
                               snp_col, beta_col, se_col, pval_col,
                               ea_col, oa_col, eaf_col=NULL, compute_p=FALSE) {
  dt <- fread(file, showProgress=FALSE)

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
  out$outcome        <- outcome_name
  out$samplesize.outcome <- n_total
  as.data.table(out)
}

outcome_ibs <- read_local_outcome(
  file.path(GWAS_DIR, "gwas_TRAIT_A.txt"),
  outcome_name = "TRAIT_A", n_total = IBS_N,
  snp_col="SNP", beta_col="BETA", se_col="SE",
  pval_col=NULL, ea_col="A1", oa_col="A2",
  compute_p=TRUE
)
outcome_asd <- read_local_outcome(
  file.path(GWAS_DIR, "gwas_TRAIT_B.txt"),
  outcome_name = "TRAIT_B", n_total = ASD_N,
  snp_col="SNP", beta_col="BETA", se_col="SE",
  pval_col="P", ea_col="A1", oa_col="A2", eaf_col="EAF"
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
fwrite(as.data.table(harmonised), file.path(OUT_DIR, "harmonised_data.txt"), sep="\t")

harmonised <- harmonised[harmonised$mr_keep == TRUE, ]
cat(sprintf(mr_keep==TRUE：%d msg\n, nrow(harmonised)))

# ── 7. MR analysis ────────────────────────────────────────────────
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
fwrite(mr_dt, file.path(OUT_DIR, "MR_results_all.txt"), sep="\t")

cat(\n=== IVW msg ===\n)
ivw <- mr_dt[method == "Inverse variance weighted"][order(pval)]
print(ivw[, .(exposure, outcome, nsnp,
               b=round(b,4), OR=round(OR,3),
               CI_lo=round(CI_lo,3), CI_hi=round(CI_hi,3),
               pval=formatC(pval,format="e",digits=2))], row.names=FALSE)

cat(\nmsg（MR-Egger msg）...\n)
pleio <- mr_pleiotropy_test(harmonised)
if (!is.null(pleio) && nrow(pleio) > 0) {
  print(as.data.frame(pleio[, c("exposure","outcome","egger_intercept","se","pval")]),
        row.names=FALSE)
  fwrite(as.data.table(pleio), file.path(OUT_DIR, "MR_pleiotropy_test.txt"), sep="\t")
}

het <- mr_heterogeneity(harmonised)
if (!is.null(het) && nrow(het) > 0)
  fwrite(as.data.table(het), file.path(OUT_DIR, "MR_heterogeneity.txt"), sep="\t")

# ── 9. visualisation ─────────────────────────────────────────────────
cat(\nmsg...\n)

if (nrow(ivw) > 0) {
  ivw[, sig := pval < 0.05]
  p_forest <- ggplot(ivw, aes(x = OR, y = reorder(paste0(exposure,"\n→",outcome), OR),
                               color = outcome, shape = sig)) +
    geom_vline(xintercept=1, linetype="dashed", color="gray50") +
    geom_errorbarh(aes(xmin=CI_lo, xmax=CI_hi), height=0.25, linewidth=0.5) +
    geom_point(size=3) +
    scale_x_log10() +
    scale_shape_manual(values=c("TRUE"=16, "FALSE"=1),
                       labels=c("TRUE"="p<0.05","FALSE"="p≥0.05")) +
    scale_color_manual(values=c("TRAIT_A"="#2166AC","TRAIT_B"="#D6604D")) +
    theme_bw(base_size=9) +
    theme(axis.text.y=element_text(size=7)) +
    labs(title=msg → msg MR（IVW）,
         x="OR (95% CI，log scale)", y=NULL,
         color=msg, shape=msg)
  ggsave(file.path(OUT_DIR, "MR_forest_IVW.pdf"), p_forest,
         width=10, height=max(4, nrow(ivw)*0.45))
}

sig_pairs <- unique(mr_dt[method=="Inverse variance weighted" & pval<0.05,
                            paste(exposure, outcome)])
if (length(sig_pairs) > 0) {
  har_sig <- harmonised[paste(harmonised$exposure, harmonised$outcome) %in% sig_pairs, ]
  res_sig  <- mr(har_sig, method_list=c("mr_ivw","mr_egger_regression",
                                         "mr_weighted_median"))
  p_scatter <- mr_scatter_plot(res_sig, har_sig)
  pdf(file.path(OUT_DIR, "MR_scatter_significant.pdf"), width=7, height=6)
  for (plt in p_scatter) print(plt)
  dev.off()
  cat(  msg（msg）\n)
}

cat(\n=== msg MR msg ===\n)
cat(msg：, OUT_DIR, "\n")
cat(\nmsg：\n)
cat(  instruments_clumped.txt   clump msg\n)
cat(  MR_results_all.txt        msg × msg × msg\n)
cat(  MR_pleiotropy_test.txt    msg（msg≠0=msg）\n)
cat(  MR_forest_IVW.pdf         msg\n)
cat(  MR_scatter_significant.pdf msg（p<0.05msg）\n)
