#!/usr/bin/env Rscript
# ============================================================
# 05c_analyze_twas_Anxiety.R
# trait pair (A x B) Genome-wide TWAS result integration
#
#   5. visualisation（heatmap）
# ============================================================

library(data.table)
library(ggplot2)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── path ─────────────────────────────────────────────────
TWAS_DIR <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("TRAIT_A_GWAS", "disease_anxiety")
TISSUES  <- c("Brain_Frontal_Cortex_BA9",
              "Brain_Caudate_basal_ganglia",
              "Brain_Cerebellum",
              "Colon_Transverse",
              "Colon_Sigmoid",
              "Small_Intestine_Terminal_Ileum",
              "Whole_Blood")

GENES_OF_INTEREST <- c(
  "gene_M",
  "FKBP5",
  "SLC6A4",
  "BDNF",
  "COMT",
  "gene_N",
  "gene_B",
  "TPH1",
  "HTR3A",
  "HTR4",
  "NRXN1",
  "CNTNAP2"
)

cat(msg TWAS msg...\n)
all_res <- list()

for (ds in DATASETS) {
  for (tis in TISSUES) {
    fname <- file.path(TWAS_DIR, paste0(ds, "__", tis, ".csv"))
    if (!file.exists(fname)) {
      cat(sprintf("  [SKIP] %s × %s\n", ds, tis))
      next
    }
    dt <- tryCatch(
      fread(fname, showProgress = FALSE),
      error = function(e) { cat(sprintf("  [ERROR] %s × %s: %s\n", ds, tis, e$message)); NULL }
    )
    if (is.null(dt) || nrow(dt) == 0) next

    dt[, dataset := ds]
    dt[, tissue  := tis]
    all_res[[paste(ds, tis, sep = "__")]] <- dt
    cat(sprintf(  [OK] %-12s × %-40s → %d msg\n, ds, tis, nrow(dt)))
  }
}

if (length(all_res) == 0) stop(msg TWAS msg，msg step2（05b_run_spredixcan_Anxiety.sh）。)

combined <- rbindlist(all_res, fill = TRUE, use.names = TRUE)
cat(sprintf(\nmsg：%d msg-msg\n, nrow(combined)))

if ("gene_name"   %in% names(combined)) combined[, GENE := gene_name]
if ("zscore"      %in% names(combined)) combined[, Z    := as.numeric(zscore)]
if ("pvalue"      %in% names(combined)) combined[, P    := as.numeric(pvalue)]
if ("effect_size" %in% names(combined)) combined[, ES   := as.numeric(effect_size)]

combined[, FDR := p.adjust(P, method = "BH"), by = .(dataset, tissue)]

cat(sprintf(\nmsg (P < 0.05):  %d\n, sum(combined$P < 0.05, na.rm=TRUE)))
cat(sprintf(FDR msg (q < 0.05):  %d\n, sum(combined$FDR < 0.05, na.rm=TRUE)))

cat(\nmsg FDR msg：\n)
stats_dt <- combined[, .(n_fdr = sum(FDR < 0.05, na.rm=TRUE),
                          n_nom = sum(P  < 0.05, na.rm=TRUE),
                          n_tot = .N),
                      by = .(dataset, tissue)]
print(stats_dt[order(-n_fdr)], row.names=FALSE)

# ── 3. FDR significantsummary ──────────────────────────────────────
sig_fdr <- combined[FDR < 0.05]
if (nrow(sig_fdr) > 0) {
  cat(sprintf(\n=== FDR msg（msg %d msg gene×tissue msg）===\n, nrow(sig_fdr)))
  fdr_summary <- sig_fdr[, .(
    n_tissues = .N,
    tissues   = paste(unique(tissue), collapse="; "),
    min_P     = min(P, na.rm=TRUE),
    min_FDR   = min(FDR, na.rm=TRUE)
  ), by = .(dataset, GENE)]
  print(fdr_summary[order(dataset, min_P)], row.names=FALSE)
  fwrite(sig_fdr, file.path(OUT_DIR, "TWAS_FDR_sig_all.txt"), sep="\t")
}

cat("\n\n", strrep("=",60), "\n", sep="")
cat(msg：IBS msg TWAS msg\n)
cat(strrep("=",60), "\n")

find_shared_genes <- function(dt, thresh_P, thresh_label) {
  ibs_sig <- dt[dataset == "TRAIT_A_GWAS" & P < thresh_P]
  anx_sig <- dt[dataset == "disease_anxiety"  & P < thresh_P]

  if (nrow(ibs_sig) == 0 || nrow(anx_sig) == 0) {
    cat(sprintf(  [%s] TRAIT_A msg，msg\n, thresh_label))
    return(NULL)
  }

  shared <- merge(
    ibs_sig[, .(GENE, tissue, Z_ibs=Z, P_ibs=P, FDR_ibs=FDR)],
    anx_sig[, .(GENE, tissue, Z_anx=Z, P_anx=P, FDR_anx=FDR)],
    by = c("GENE", "tissue")
  )

  if (nrow(shared) == 0) return(NULL)

  shared[, direction := ifelse(sign(Z_ibs) == sign(Z_anx),
                                msg（concordant）, msg（antagonistic）)]
  shared[, joint_P := P_ibs * P_anx]
  setorder(shared, joint_P)
  cat(sprintf(\n[%s] msg：%d msg gene×tissue msg\n, thresh_label, nrow(shared)))
  print(shared[, .(GENE, tissue, Z_ibs=round(Z_ibs,3), Z_anx=round(Z_anx,3),
                    P_ibs=formatC(P_ibs,format="e",digits=1),
                    P_anx=formatC(P_anx,format="e",digits=1),
                    direction)], row.names=FALSE)
  return(shared)
}

n_ibs <- nrow(combined[dataset == "TRAIT_A_GWAS"])
shared_fdr  <- find_shared_genes(combined, 0.05 / max(n_ibs, 1), "FDR~Bonf")
shared_nom  <- find_shared_genes(combined, 0.05, "P<0.05")
shared_01   <- if (is.null(shared_nom)) find_shared_genes(combined, 0.1, "P<0.1") else NULL

best_shared <- shared_fdr %||% shared_nom %||% shared_01
if (!is.null(best_shared) && nrow(best_shared) > 0) {
  fwrite(best_shared, file.path(OUT_DIR, "TWAS_shared_IBS_Anxiety.txt"), sep="\t")
  cat(\nmsg。\n)

  antagonistic <- best_shared[direction == msg（antagonistic）]
  if (nrow(antagonistic) > 0) {
    cat(sprintf(\nmsg（msg，msg %d msg）：\n, nrow(antagonistic)))
    print(antagonistic[, .(GENE, tissue, Z_ibs=round(Z_ibs,2), Z_anx=round(Z_anx,2),
                            P_ibs=formatC(P_ibs,format="e",digits=1),
                            P_anx=formatC(P_anx,format="e",digits=1))], row.names=FALSE)
    fwrite(antagonistic, file.path(OUT_DIR, "TWAS_antagonistic_IBS_Anxiety.txt"), sep="\t")
  }
}

cat("\n\n", strrep("=",60), "\n", sep="")
cat(msg\n)
cat(strrep("=",60), "\n")

for (gene in GENES_OF_INTEREST) {
  rows <- combined[GENE == gene]
  if (nrow(rows) == 0) {
    cat(sprintf(  %-12s：msg TWAS msg（eQTL msg）\n, gene))
    next
  }
  cat(sprintf(\n  %s（%d msg）：\n, gene, nrow(rows)))
  show <- rows[order(P), .(dataset, tissue, Z=round(Z,3),
                             P=formatC(P,format="e",digits=2), FDR=round(FDR,3))]
  print(show, row.names=FALSE)
}

cat(\nmsg...\n)
wide <- dcast(combined[P < 0.1],
              GENE ~ dataset + tissue,
              value.var = "Z", fun.aggregate = mean)
fwrite(wide, file.path(OUT_DIR, "TWAS_wide_table_P01.txt"), sep="\t")
cat(msg（P<0.1）：, file.path(OUT_DIR, "TWAS_wide_table_P01.txt"), "\n")

if (!is.null(best_shared) && nrow(best_shared) > 0) {
  heat_genes <- unique(best_shared$GENE)[1:min(40, length(unique(best_shared$GENE)))]
  heat_data  <- combined[GENE %in% heat_genes & P < 0.5]
  heat_data[, label := paste0(dataset, "\n(", gsub("Brain_|_basal_ganglia|_Terminal_Ileum","",tissue), ")")]
  heat_data[, Z_cap := pmax(pmin(Z, 5), -5)]
  heat_data[, sig_mark := ifelse(FDR < 0.05, "**", ifelse(P < 0.05, "*", ""))]

  p <- ggplot(heat_data, aes(x = label, y = GENE, fill = Z_cap)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig_mark), size = 3, vjust = 0.5) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D",
                         midpoint = 0, name = "Z-score\n(±5 cap)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank()) +
    labs(title    = TRAIT_A x disease_anxiety msg TWAS msg,
         subtitle = "* P<0.05；** FDR<0.05",
         x = msg × msg, y = msg)

  ggsave(file.path(OUT_DIR, "TWAS_shared_heatmap.pdf"),
         p, width = 14, height = max(6, length(heat_genes) * 0.45),
         limitsize = FALSE)
  cat(msg：, file.path(OUT_DIR, "TWAS_shared_heatmap.pdf"), "\n")
}

cat("\n\n", strrep("=",60), "\n", sep="")
cat(TRAIT_A x disease_anxiety TWAS msg\n)
cat(strrep("=",60), "\n")

pleiotropic_genes <- if (exists("antagonistic") && !is.null(antagonistic) && nrow(antagonistic) > 0) {
  unique(antagonistic$GENE)
} else {
  character(0)
}
cat(sprintf(Antagonistic pleiotropy genes (TRAIT_A x disease_anxietymsg P<0.05 msg，msg %d msg）：\n,
            length(pleiotropic_genes)))
if (length(pleiotropic_genes) > 0)
  cat(" ", paste(sort(pleiotropic_genes), collapse=", "), "\n")

sig_specific  <- sig_fdr[!GENE %in% pleiotropic_genes]
ibs_spec  <- sig_specific[dataset == "TRAIT_A_GWAS"]
anx_spec  <- sig_specific[dataset == "disease_anxiety"]

cat(sprintf(\nTRAIT_A specific FDR msg（msg %d msg gene×tissue，%d msg）：\n,
            nrow(ibs_spec), uniqueN(ibs_spec$GENE)))
if (nrow(ibs_spec) > 0) {
  ibs_summary <- ibs_spec[, .(
    n_tissues = .N,
    tissues   = paste(gsub("Brain_|_basal_ganglia|_Terminal_Ileum","",unique(tissue)), collapse="; "),
    best_P    = min(P, na.rm=TRUE),
    best_Z    = Z[which.min(P)]
  ), by = GENE][order(best_P)]
  print(ibs_summary[, .(GENE, n_tissues, tissues,
                         best_P = formatC(best_P, format="e", digits=1),
                         best_Z = round(best_Z, 2))], row.names=FALSE)
  fwrite(ibs_summary, file.path(OUT_DIR, "TWAS_IBS_specific_FDR.txt"), sep="\t")
}

cat(sprintf(\ndisease_anxiety specific FDR msg（msg %d msg gene×tissue，%d msg）：\n,
            nrow(anx_spec), uniqueN(anx_spec$GENE)))
if (nrow(anx_spec) > 0) {
  anx_summary <- anx_spec[, .(
    n_tissues = .N,
    tissues   = paste(gsub("Brain_|_basal_ganglia|_Terminal_Ileum","",unique(tissue)), collapse="; "),
    best_P    = min(P, na.rm=TRUE),
    best_Z    = Z[which.min(P)]
  ), by = GENE][order(best_P)]
  print(anx_summary[, .(GENE, n_tissues, tissues,
                         best_P = formatC(best_P, format="e", digits=1),
                         best_Z = round(best_Z, 2))], row.names=FALSE)
  fwrite(anx_summary, file.path(OUT_DIR, "TWAS_Anxiety_specific_FDR.txt"), sep="\t")
}

cat(sprintf(\nmsg：%d msg TRAIT_A specificmsg，%d msgdisease_anxiety specificmsg（FDR<0.05，msg）\n,
            uniqueN(ibs_spec$GENE), uniqueN(anx_spec$GENE)))
cat(msg：, OUT_DIR, "\n")
cat(TWAS msg。\n)
