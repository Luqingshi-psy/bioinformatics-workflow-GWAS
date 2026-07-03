#!/usr/bin/env Rscript
# ============================================================
# 08_p0a_eqtl_snp_lookup.R
# Phase 0-A:extract top eQTL SNP for gene，
#        validate effect direction in GWAS，
#        build table。
#
# ============================================================

library(data.table)
library(ggplot2)

# ── pathconfig ─────────────────────────────────────────────────
SMR_DIR  <- "${PROJECT_ROOT}"
GWAS_DIR <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_GENES   <- c("gene_A", "gene_B")
P_SMR_THRESH   <- 0.05
P_HEIDI_THRESH <- 0.05

cat(msg SMR msg...\n)

smr_files <- list.files(SMR_DIR, pattern = "\\.smr$", full.names = TRUE)
cat(sprintf(  msg %d msg .smr msg\n, length(smr_files)))

all_rows <- lapply(smr_files, function(f) {
  dt <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  fname   <- tools::file_path_sans_ext(basename(f))
  parts   <- strsplit(fname, "_")[[1]]
  dataset <- tail(parts, 1)

  dt_target <- dt[Gene %in% TARGET_GENES]
  if (nrow(dt_target) == 0) return(NULL)
  dt_target[, `:=`(tissue = tissue, dataset = dataset)]
  dt_target
})

smr_combined <- rbindlist(Filter(Negate(is.null), all_rows))
cat(sprintf(  gene_A/gene_B msg %d msg\n, nrow(smr_combined)))

for (col in c("b_GWAS", "se_GWAS", "p_GWAS",
              "b_eQTL", "se_eQTL", "p_eQTL",
              "b_SMR",  "se_SMR",  "p_SMR", "p_HEIDI")) {
  if (col %in% names(smr_combined))
    smr_combined[, (col) := as.numeric(get(col))]
}

sig <- smr_combined[p_SMR < P_SMR_THRESH &
                    (is.na(p_HEIDI) | p_HEIDI > P_HEIDI_THRESH)]
cat(sprintf(  pSMR<0.05 & pHEIDI>0.05：%d msg\n, nrow(sig)))

fwrite(sig, file.path(OUT_DIR, "smr_mtor_ctsb_significant.txt"), sep = "\t")

cat(\nmsg SNP（msg×msg pSMR）...\n)

sig_ordered <- sig[order(p_SMR)]
anchor <- sig_ordered[, .SD[which.min(p_SMR)], by = .(Gene, dataset)]
anchor_snps <- unique(anchor$topSNP)

cat(sprintf(  msg SNP：%d msg\n, length(anchor_snps)))
cat(\n  msg：\n)
print(anchor[, .(Gene, dataset, tissue, topSNP,
                 b_eQTL = round(b_eQTL, 4),
                 b_SMR  = round(b_SMR, 4),
                 p_SMR  = formatC(p_SMR, format = "e", digits = 2),
                 p_HEIDI = round(p_HEIDI, 3))],
      row.names = FALSE)

cat(\nmsg SNP msg GWAS msg...\n)

ibs_raw <- fread(file.path(GWAS_DIR, "gwas_TRAIT_A.txt"), showProgress = FALSE)
setnames(ibs_raw, names(ibs_raw), toupper(names(ibs_raw)))
ibs_raw[, P_COMPUTED := 2 * pnorm(-abs(as.numeric(BETA) / as.numeric(SE)))]

ibs_sub <- ibs_raw[SNP %in% anchor_snps,
                   .(SNP,
                     A1_ibs   = A1, A2_ibs = A2,
                     beta_ibs = as.numeric(BETA),
                     se_ibs   = as.numeric(SE),
                     p_ibs    = P_COMPUTED)]

asd_raw <- fread(file.path(GWAS_DIR, "gwas_TRAIT_B.txt"), showProgress = FALSE)
setnames(asd_raw, names(asd_raw), toupper(names(asd_raw)))

asd_sub <- asd_raw[SNP %in% anchor_snps,
                   .(SNP,
                     A1_asd   = A1, A2_asd = A2,
                     beta_asd = as.numeric(BETA),
                     se_asd   = as.numeric(SE),
                     p_asd    = as.numeric(P))]

cat(sprintf(  TRAIT_A msg %d / %d msg SNP\n, nrow(ibs_sub), length(anchor_snps)))
cat(sprintf(  TRAIT_B msg %d / %d msg SNP\n, nrow(asd_sub), length(anchor_snps)))

cat(\nmsg...\n)

base <- anchor[, .(Gene, dataset, tissue, topSNP,
                   A1_eqtl = A1, A2_eqtl = A2,
                   b_eQTL, se_eQTL, p_eQTL,
                   b_SMR, se_SMR, p_SMR, p_HEIDI,
                   b_GWAS_smr = b_GWAS,
                   p_GWAS_smr = p_GWAS)]

result <- merge(base, ibs_sub, by.x = "topSNP", by.y = "SNP", all.x = TRUE)
result <- merge(result, asd_sub, by.x = "topSNP", by.y = "SNP", all.x = TRUE)

result[, beta_ibs_aligned := fifelse(
  !is.na(A1_ibs) & toupper(A1_ibs) != toupper(A1_eqtl),
  -beta_ibs, beta_ibs
)]
result[, beta_asd_aligned := fifelse(
  !is.na(A1_asd) & toupper(A1_asd) != toupper(A1_eqtl),
  -beta_asd, beta_asd
)]

result[, dir_consistent_IBS := sign(b_SMR * b_eQTL) == sign(beta_ibs_aligned)]
result[, dir_consistent_ASD := sign(b_SMR * b_eQTL) == sign(beta_asd_aligned)]

result[, gwas_consistent_IBS := sign(b_GWAS_smr) == sign(beta_ibs_aligned)]
result[, gwas_consistent_ASD := sign(b_GWAS_smr) == sign(beta_asd_aligned)]

out_tab <- result[, .(
  Gene, topSNP, A1_eqtl, A2_eqtl, tissue, dataset,
  b_eQTL        = round(b_eQTL, 4),
  p_eQTL        = formatC(p_eQTL, format = "e", digits = 2),
  b_SMR         = round(b_SMR, 4),
  p_SMR         = formatC(p_SMR, format = "e", digits = 2),
  p_HEIDI       = round(p_HEIDI, 3),
  beta_IBS      = round(beta_ibs_aligned, 5),
  p_IBS         = formatC(p_ibs, format = "e", digits = 2),
  beta_ASD      = round(beta_asd_aligned, 5),
  p_ASD         = formatC(p_asd, format = "e", digits = 2),
  dir_IBS_consistent = dir_consistent_IBS,
  dir_ASD_consistent = dir_consistent_ASD
)]

setorder(out_tab, Gene, dataset, p_SMR)
fwrite(out_tab, file.path(OUT_DIR, "anchor_snp_direction_table.txt"), sep = "\t")

cat(\n=== msg SNP msg ===\n)
print(out_tab, row.names = FALSE)

cat(\n=== msg ===\n)

for (gene in TARGET_GENES) {
  rows <- result[Gene == gene]
  if (nrow(rows) == 0) next

  snps <- unique(rows$topSNP)
  for (snp in snps) {
    r <- rows[topSNP == snp][1]
    cat(sprintf(\n[%s] msg SNP：%s（A1 = %s）\n, gene, snp, r$A1_eqtl))
    cat(sprintf(  eQTL msg：b = %+.4f（p = %s）→ A1 msg%smsg\n,
                r$b_eQTL,
                formatC(r$p_eQTL, format = "e", digits = 2),
                ifelse(r$b_eQTL > 0, ↑ msg, ↓ msg)))
    cat(sprintf(  SMR msg：b = %+.4f（p = %s，pHEIDI = %.3f）→ msg↑msg%s\n,
                r$b_SMR,
                formatC(r$p_SMR, format = "e", digits = 2),
                r$p_HEIDI,
                ifelse(r$b_SMR > 0, ↑ msg, ↓ msg)))

    if (!is.na(r$beta_ibs_aligned)) {
      cat(sprintf(  TRAIT_A GWAS msg：β = %+.5f（p = %s）→ msg%s\n,
                  r$beta_ibs_aligned,
                  formatC(r$p_ibs, format = "e", digits = 2),
                  ifelse(isTRUE(r$dir_consistent_IBS), ✓ msg SMR msg, ✗ msg SMR msg)))
    } else {
      cat(  TRAIT_A GWAS：msg SNP\n)
    }

    if (!is.na(r$beta_asd_aligned)) {
      cat(sprintf(  TRAIT_B GWAS msg：β = %+.5f（p = %s）→ msg%s\n,
                  r$beta_asd_aligned,
                  formatC(r$p_asd, format = "e", digits = 2),
                  ifelse(isTRUE(r$dir_consistent_ASD), ✓ msg SMR msg, ✗ msg SMR msg)))
    } else {
      cat(  TRAIT_B GWAS：msg SNP\n)
    }
  }
}

cat(\nmsg SNP msg...\n)

plot_data <- list()
for (i in seq_len(nrow(result))) {
  r <- result[i]
  snp_label <- sprintf("%s\n(%s)", r$Gene, r$topSNP)
  tissue_short <- gsub("GTEx_|BrainMeta_|_IBS|_ASD", "", r$tissue)
  tissue_short <- substr(tissue_short, 1, 25)

  if (!is.na(r$beta_ibs_aligned)) {
    plot_data[[length(plot_data) + 1]] <- data.table(
      label   = snp_label,
      tissue  = tissue_short,
      phenotype = "TRAIT_A",
      beta    = r$beta_ibs_aligned,
      se      = r$se_ibs,
      p_val   = r$p_ibs,
      consistent = r$dir_consistent_IBS
    )
  }
  if (!is.na(r$beta_asd_aligned)) {
    plot_data[[length(plot_data) + 1]] <- data.table(
      label   = snp_label,
      tissue  = tissue_short,
      phenotype = "TRAIT_B",
      beta    = r$beta_asd_aligned,
      se      = r$se_asd,
      p_val   = r$p_asd,
      consistent = r$dir_consistent_ASD
    )
  }
}

if (length(plot_data) > 0) {
  pd <- rbindlist(plot_data)
  pd[, ci_lo := beta - 1.96 * se]
  pd[, ci_hi := beta + 1.96 * se]

  p_anchor <- ggplot(pd, aes(x = beta, y = label,
                              color = phenotype, shape = consistent)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                   linewidth = 0.6, position = position_dodge(0.4)) +
    geom_point(size = 3, position = position_dodge(0.4)) +
    scale_color_manual(values = c("TRAIT_A" = "#2166AC", "TRAIT_B" = "#D6604D")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = msg, "FALSE" = msg)) +
    theme_bw(base_size = 10) +
    labs(title = msg SNP msg TRAIT_A / TRAIT_B GWAS msg,
         subtitle = msg SMR msg，msg→msg→msg,
         x = GWAS β（msg A1 msg）,
         y = NULL,
         color = msg, shape = SMR msg)
  ggsave(file.path(OUT_DIR, "anchor_snp_effect_plot.pdf"),
         p_anchor, width = 8, height = max(3, nrow(pd) * 0.5))
  cat(  msg\n)
}

cat(sprintf(\n=== P0-A msg ===\n))
cat(sprintf(msg：%s\n, OUT_DIR))
cat(  anchor_snp_direction_table.txt   msg SNP msg（msg）\n)
cat(  smr_mtor_ctsb_significant.txt    gene_A/gene_B msg SMR msg\n)
cat(  anchor_snp_effect_plot.pdf       msg\n)
