#!/usr/bin/env Rscript
# Parse pleiotropyFDR conjFDR results for TRAIT_B × gut diseases
# Input:  output/pleiofdr/results/*/  (CSV/mat files from MATLAB run)
# Output: summary table + Manhattan plots

library(dplyr)
library(ggplot2)
library(data.table)

PLEIOFDR_DIR <- "${PROJECT_ROOT}"
OUT_DIR      <- "${PROJECT_ROOT}"

PAIRS <- list(
  list(label = "TRAIT_B_vs_TRAIT_A",           display = "TRAIT_B vs TRAIT_A",          color = "#E63946"),
  list(label = "TRAIT_B_vs_TRAIT_C",   display = "TRAIT_B vs TRAIT_C",color = "#457B9D"),
  list(label = "TRAIT_B_vs_TRAIT_E",    display = "TRAIT_B vs TRAIT_E", color = "#2A9D8F"),
  list(label = "TRAIT_B_vs_TRAIT_D",    display = "TRAIT_B vs TRAIT_D", color = "#E9C46A"),
  list(label = "TRAIT_B_vs_TRAIT_C_Liu",       display = "TRAIT_B vs TRAIT_C_Liu",     color = "#264653"),
  list(label = "TRAIT_B_vs_TRAIT_E_Liu",        display = "TRAIT_B vs TRAIT_E_Liu",      color = "#F4A261"),
  list(label = "TRAIT_B_vs_TRAIT_D_Liu",        display = "TRAIT_B vs TRAIT_D_Liu",      color = "#A8DADC")
)

# ── Try to load pleiofdr CSV output ──────────────────────────────────────────
# pleiofdr saves CSV files in the output directory; typical name pattern:
#   *_conjfdr_loci.csv  or  *_results.csv
load_conj_results <- function(pair_dir, label) {
  csvs <- list.files(pair_dir, pattern = "(?i)(loci|result|conj).*\\.csv$",
                     recursive = TRUE, full.names = TRUE)
  if (length(csvs) == 0) {
    message("No CSV found in: ", pair_dir)
    return(NULL)
  }
  # Use the one most likely to be the loci table
  f <- csvs[1]
  message("  Reading: ", f)
  df <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$pair_label <- label
  df
}

# ── Collect across pairs ──────────────────────────────────────────────────────
all_loci <- list()
summary_rows <- list()

for (p in PAIRS) {
  conj_dir <- file.path(PLEIOFDR_DIR, paste0(p$label, "_conjFDR"))
  if (!dir.exists(conj_dir)) {
    message("Missing: ", conj_dir)
    next
  }
  df <- load_conj_results(conj_dir, p$display)
  if (!is.null(df)) {
    all_loci[[p$label]] <- df
    n_loci <- nrow(df)
    message(sprintf("  %s: %d loci", p$display, n_loci))
    summary_rows[[p$label]] <- data.frame(
      pair    = p$display,
      n_loci  = n_loci,
      stringsAsFactors = FALSE
    )
  }
}

# ── Summary table ─────────────────────────────────────────────────────────────
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  cat("\n=== conjFDR<0.05 Loci Summary ===\n")
  print(summary_df, row.names = FALSE)
  write.csv(summary_df, file.path(OUT_DIR, "conjfdr_loci_counts.csv"), row.names = FALSE)
}

# ── Locus count barplot ───────────────────────────────────────────────────────
if (length(summary_rows) > 0 && nrow(summary_df) > 0) {
  colors <- setNames(
    sapply(PAIRS, `[[`, "color"),
    sapply(PAIRS, `[[`, "display")
  )
  p_bar <- ggplot(summary_df, aes(x = reorder(pair, n_loci), y = n_loci, fill = pair)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n_loci), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_fill_manual(values = colors, na.value = "grey60") +
    labs(
      x = NULL, y = "Number of conjFDR < 0.05 loci",
      title = "Shared genetic loci: TRAIT_B × gut diseases (conjFDR)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none") +
    expand_limits(y = max(summary_df$n_loci, na.rm = TRUE) * 1.15)

  ggsave(file.path(OUT_DIR, "conjfdr_loci_barplot.pdf"), p_bar, width = 8, height = 5)
  cat("Saved: conjfdr_loci_barplot.pdf\n")
}

# ── Manhattan-like plot for one pair (TRAIT_B vs TRAIT_A as example) ────────────────
if ("TRAIT_B_vs_TRAIT_A" %in% names(all_loci)) {
  df_ibs <- all_loci[["TRAIT_B_vs_TRAIT_A"]]

  # Detect CHR and BP columns (pleiofdr uses 'CHR' and 'BP' or 'pos')
  chr_col <- grep("^chr", names(df_ibs), ignore.case = TRUE, value = TRUE)[1]
  bp_col  <- grep("^(bp|pos|position)", names(df_ibs), ignore.case = TRUE, value = TRUE)[1]
  fdr_col <- grep("^conj", names(df_ibs), ignore.case = TRUE, value = TRUE)[1]

  if (!is.na(chr_col) && !is.na(bp_col) && !is.na(fdr_col)) {
    df_ibs[[chr_col]] <- as.integer(df_ibs[[chr_col]])
    df_ibs[[bp_col]]  <- as.numeric(df_ibs[[bp_col]])
    df_ibs$neglog_fdr <- -log10(df_ibs[[fdr_col]])
    df_ibs$sig        <- df_ibs[[fdr_col]] < 0.05

    p_manh <- ggplot(df_ibs, aes(x = .data[[bp_col]], y = neglog_fdr,
                                  color = factor(.data[[chr_col]] %% 2))) +
      geom_point(size = 0.6, alpha = 0.6) +
      geom_point(data = subset(df_ibs, sig), aes(x = .data[[bp_col]], y = neglog_fdr),
                 color = "red", size = 1.2) +
      facet_grid(. ~ .data[[chr_col]], scales = "free_x", space = "free_x") +
      scale_color_manual(values = c("0" = "#457B9D", "1" = "#A8DADC")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      labs(x = "Chromosome", y = "-log10(conjFDR)",
           title = "Shared loci: TRAIT_B × TRAIT_A (conjFDR < 0.05 in red)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none",
            panel.spacing   = unit(0.1, "lines"),
            axis.text.x     = element_blank())

    ggsave(file.path(OUT_DIR, "conjfdr_manhattan_ASD_IBS.pdf"), p_manh,
           width = 14, height = 4)
    cat("Saved: conjfdr_manhattan_ASD_IBS.pdf\n")
  }
}

cat("\nDone. Results in:", OUT_DIR, "\n")
