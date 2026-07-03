#!/usr/bin/env Rscript
# ============================================================
# 06_step3_analyze_smr.R
#
# ============================================================

library(data.table)
library(ggplot2)

SMR_DIR <- "${PROJECT_ROOT}"
OUT_DIR <- "${PROJECT_ROOT}"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_GENES <- c(
  "gene_A", "gene_B", "gene_I", "EP300", "gene_L",
  "gene_K", "gene_J", "EHMT2", "BAG6",
  "gene_M", "gene_N", "APOM", "FKBP5"
)

cat(msg SMR msg...\n)

smr_files <- list.files(SMR_DIR, pattern = "\\.smr$", full.names = TRUE)
cat(sprintf(  msg %d msg .smr msg\n, length(smr_files)))

if (length(smr_files) == 0) {
  stop(msg .smr msg，msg 06_step2_run_smr.sh)
}

all_smr <- lapply(smr_files, function(f) {
  dt <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  #   GTEx_{tissue}_{DS}.smr          → source=GTEx_V8
  #   BrainMeta_chr{n}_{DS}.smr       → source=BrainMeta
  bname <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(bname, "_")[[1]]

  if (grepl("^GTEx", bname)) {
    source_db <- "GTEx_V8"
    dataset   <- tail(parts, 1)
  } else {
    # BrainMeta_chr{n}_{DS}
    source_db <- "BrainMeta"
    dataset   <- tail(parts, 1)
    tissue    <- paste(parts[1:(length(parts)-1)], collapse = "_")
  }

  dt[, `:=`(source = source_db, tissue = tissue, dataset = dataset)]
  dt
})

combined_smr <- rbindlist(Filter(Negate(is.null), all_smr), fill = TRUE)
cat(sprintf(  msg：%d msg-msg\n, nrow(combined_smr)))

#            A1 A2 Freq b_GWAS se_GWAS p_GWAS b_eQTL se_eQTL p_eQTL
#            b_SMR se_SMR p_SMR p_HEIDI nsnp_HEIDI
required <- c("Gene", "p_SMR", "p_HEIDI", "b_SMR", "se_SMR")
missing  <- setdiff(required, names(combined_smr))
if (length(missing) > 0)
  stop(sprintf(msg：%s, paste(missing, collapse=", ")))

for (col in c("p_SMR","p_HEIDI","b_SMR","se_SMR","b_GWAS","p_GWAS","b_eQTL","p_eQTL","nsnp_HEIDI"))
  if (col %in% names(combined_smr))
    combined_smr[, (col) := as.numeric(get(col))]

# ── 3. Filter target genes ──────────────────────────────────────────
target_smr <- combined_smr[Gene %in% TARGET_GENES]
cat(sprintf(  msg（%d msg）msg：%d msg\n,
            length(TARGET_GENES), nrow(target_smr)))

if (nrow(target_smr) == 0) {
  cat(\n[WARNING] msg SMR msg。\n)
  cat(  msg：msg cis-eQTL（p<5e-8）\n)
  cat(  msg：msg cis-eQTL msg（--peqtl-smr 1e-5）msg step2\n)
  cat(  msg top msg。\n\n)

  top_all <- combined_smr[!is.na(p_SMR)][order(p_SMR)][1:min(20,.N)]
  print(top_all[, .(Gene, dataset, tissue, b_SMR=round(b_SMR,3),
                     p_SMR=formatC(p_SMR,format="e",digits=2),
                     p_HEIDI=round(p_HEIDI,3))], row.names=FALSE)
  fwrite(combined_smr, file.path(OUT_DIR, "SMR_all_genes_raw.txt"), sep="\t")
  stop(msg，msg。)
}

fwrite(target_smr, file.path(OUT_DIR, "SMR_target_genes_raw.txt"), sep="\t")

cat(\n=== SMR msg ===\n)
cat(sprintf(msg：%d msg\n, nrow(target_smr)))
cat(sprintf("pSMR < 0.05：%d\n",   sum(target_smr$p_SMR  < 0.05, na.rm=TRUE)))
cat(sprintf("pHEIDI > 0.05：%d\n", sum(target_smr$p_HEIDI > 0.05, na.rm=TRUE)))

smr_pos <- target_smr[p_SMR < 0.05 & (is.na(p_HEIDI) | p_HEIDI > 0.05)]
cat(sprintf(SMR msg（pSMR<0.05 msg pHEIDI>0.05）：%d\n, nrow(smr_pos)))
smr_pos[, heidi_note := ifelse(is.na(p_HEIDI), HEIDImsg(SNP<10), HEIDImsg)]

smr_relax <- target_smr[p_SMR < 0.1 & (is.na(p_HEIDI) | p_HEIDI > 0.05)]

# ── 5. summary ─────────────────────────────────────────────────
plot_data <- if (nrow(smr_pos) > 0) smr_pos else smr_relax

if (nrow(smr_pos) > 0) {
  cat(\n=== SMR msg ===\n)
  summary_tbl <- smr_pos[, .(
    n_hits    = .N,
    tissues   = paste(unique(tissue), collapse="; "),
    best_pSMR = min(p_SMR, na.rm=TRUE),
    b_range   = paste0("[",round(min(b_SMR,na.rm=TRUE),3),",",round(max(b_SMR,na.rm=TRUE),3),"]")
  ), by = .(Gene, dataset)]
  setorder(summary_tbl, dataset, best_pSMR)
  print(summary_tbl, row.names=FALSE)
  fwrite(smr_pos,     file.path(OUT_DIR, "SMR_positive_all.txt"),    sep="\t")
  fwrite(summary_tbl, file.path(OUT_DIR, "SMR_positive_summary.txt"), sep="\t")
} else {
  cat(sprintf(\nmsg pSMR<0.05 msg，msg pSMR<0.1：%d msg\n, nrow(smr_relax)))
  if (nrow(smr_relax) > 0) {
    print(smr_relax[order(p_SMR),
                    .(Gene, dataset, tissue,
                      b_SMR=round(b_SMR,3),
                      p_SMR=formatC(p_SMR,format="e",digits=1),
                      p_HEIDI=round(p_HEIDI,3), heidi_note)], row.names=FALSE)
    fwrite(smr_relax, file.path(OUT_DIR, "SMR_relaxed_p01.txt"), sep="\t")
  }
}

cat(\n=== TRAIT_A × TRAIT_B msg ===\n)
if (nrow(plot_data) >= 2) {
  ibs_d <- plot_data[dataset == "TRAIT_A", .(Gene, tissue, source, b_IBS=b_SMR, p_IBS=p_SMR)]
  asd_d <- plot_data[dataset == "TRAIT_B", .(Gene, tissue, source, b_ASD=b_SMR, p_ASD=p_SMR)]
  shared <- merge(ibs_d, asd_d, by=c("Gene","tissue","source"))
  if (nrow(shared) > 0) {
    shared[, direction := ifelse(sign(b_IBS)==sign(b_ASD), msg, msg（msg）)]
    setorder(shared, p_IBS)
    print(shared[, .(Gene, tissue, source,
                     b_IBS=round(b_IBS,3), b_ASD=round(b_ASD,3),
                     p_IBS=formatC(p_IBS,format="e",digits=1),
                     p_ASD=formatC(p_ASD,format="e",digits=1),
                     direction)], row.names=FALSE)
    fwrite(shared, file.path(OUT_DIR, "SMR_shared_TRAIT_PAIR.txt"), sep="\t")
  } else {
    cat(  msg SMR msg×msg\n)
  }
}

if (nrow(plot_data) > 0) {
  plot_data[, ci_lo := b_SMR - 1.96 * se_SMR]
  plot_data[, ci_hi := b_SMR + 1.96 * se_SMR]
  plot_data[, tissue_short := gsub("Brain_|_basal_ganglia|_Terminal_Ileum|BrainMeta_","",tissue)]
  plot_data[, label := paste0(Gene, " | ", tissue_short)]
  plot_data[, sig := ifelse(p_SMR < 0.05 & (is.na(p_HEIDI)|p_HEIDI>0.05),
                             SMRmsg, msg(p<0.1))]

  p <- ggplot(plot_data, aes(x = b_SMR,
                              y = reorder(label, b_SMR),
                              color = dataset, shape = sig)) +
    geom_vline(xintercept = 0, linetype="dashed", color="gray50") +
    geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.25, linewidth=0.5) +
    geom_point(size=2.5) +
    scale_shape_manual(values = c(SMRmsg=16, msg(p<0.1)=1)) +
    scale_color_manual(values = c("TRAIT_A"="#2166AC", "TRAIT_B"="#D6604D")) +
    facet_wrap(~source, scales="free_y", ncol=1) +
    theme_bw(base_size=9) +
    theme(axis.text.y = element_text(size=7),
          panel.grid.minor = element_blank()) +
    labs(title   = IBS × TRAIT_B SMR msg（msg）,
         subtitle = msg：pSMR<0.05 msg pHEIDI>0.05,
         x = SMR msg b (95% CI), y = NULL,
         color=msg, shape=msg)

  h <- max(5, nrow(plot_data) * 0.35)
  ggsave(file.path(OUT_DIR, "SMR_forest_plot.pdf"), p,
         width=11, height=h, limitsize=FALSE)
  cat(\nmsg：, file.path(OUT_DIR, "SMR_forest_plot.pdf"), "\n")
}

cat(\nSMR msg。msg：, OUT_DIR, "\n")
