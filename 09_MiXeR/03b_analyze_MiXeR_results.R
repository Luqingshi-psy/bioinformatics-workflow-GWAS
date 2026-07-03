#!/usr/bin/env Rscript
# Parse and visualize MiXeR bivariate results for TRAIT_B × gut diseases
# Input:  output/MiXeR/results/*.json (combined replicates from 12b_run_mixer.sh)
# Output: output/MiXeR/mixer_summary.csv  +  venn diagrams

library(jsonlite)
library(ggplot2)
library(dplyr)

MIXER_DIR <- "${PROJECT_ROOT}"
OUT_DIR   <- "${PROJECT_ROOT}"

PAIRS <- list(
  c("TRAIT_B", "TRAIT_A",         "TRAIT_B vs TRAIT_A"),
  c("TRAIT_B", "de Lange_IBD", "TRAIT_B vs TRAIT_C"),
  c("TRAIT_B", "de Lange_CD",  "TRAIT_B vs TRAIT_E"),
  c("TRAIT_B", "de Lange_UC",  "TRAIT_B vs TRAIT_D"),
  c("TRAIT_B", "Liu_IBD",     "TRAIT_B vs TRAIT_C_Liu"),
  c("TRAIT_B", "Liu_CD",      "TRAIT_B vs TRAIT_E_Liu"),
  c("TRAIT_B", "Liu_UC",      "TRAIT_B vs TRAIT_D_Liu")
)

# ── Parse a combined bivariate JSON ──────────────────────────────────────────
parse_bivar_json <- function(json_file, label) {
  if (!file.exists(json_file)) {
    message("Missing: ", json_file)
    return(NULL)
  }
  j <- fromJSON(json_file)

  # MiXeR key estimates
  n1   <- j$ci$n1$mean;       n1_se  <- j$ci$n1$std
  n2   <- j$ci$n2$mean;       n2_se  <- j$ci$n2$std
  n12  <- j$ci$n12$mean;      n12_se <- j$ci$n12$std
  dice <- j$ci$dice$mean;     dice_se <- j$ci$dice$std
  rho  <- j$ci$rho$mean;      rho_se  <- j$ci$rho$std

  data.frame(
    label   = label,
    n1      = n1,   n1_se  = n1_se,
    n2      = n2,   n2_se  = n2_se,
    n12     = n12,  n12_se = n12_se,
    dice    = dice, dice_se = dice_se,
    rho_ge  = rho,  rho_ge_se = rho_se,  # genetic effect correlation
    stringsAsFactors = FALSE
  )
}

# ── Collect results ───────────────────────────────────────────────────────────
results <- list()
for (pair in PAIRS) {
  t1 <- pair[1]; t2 <- pair[2]; lab <- pair[3]
  json_f <- file.path(MIXER_DIR, paste0(t1, "_vs_", t2, ".fit.json"))
  res <- parse_bivar_json(json_f, lab)
  if (!is.null(res)) results[[lab]] <- res
}
df <- bind_rows(results)

# Print to console
cat("\n=== MiXeR Bivariate Results ===\n")
print(df[, c("label", "n12", "n12_se", "dice", "dice_se", "rho_ge", "rho_ge_se")],
      row.names = FALSE)

# Save
write.csv(df, file.path(OUT_DIR, "mixer_summary.csv"), row.names = FALSE)
cat("\nSaved: ", file.path(OUT_DIR, "mixer_summary.csv"), "\n")

# ── Dice coefficient bar chart ────────────────────────────────────────────────
if (nrow(df) > 0) {
  p <- ggplot(df, aes(x = reorder(label, dice), y = dice, fill = label)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = dice - dice_se, ymax = dice + dice_se), width = 0.3) +
    coord_flip() +
    labs(
      x = NULL, y = "Dice coefficient (polygenic overlap)",
      title = "ASD × Gut Diseases — Shared Genetic Architecture (MiXeR)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "mixer_dice_barplot.pdf"), p, width = 7, height = 5)
  cat("Saved: mixer_dice_barplot.pdf\n")
}

# ── Venn-like summary: shared vs unique causal variants ──────────────────────
venn_data <- df %>%
  mutate(
    unique_t1  = n1 - n12,
    unique_t2  = n2 - n12,
    shared     = n12,
    total_t1   = n1,
    total_t2   = n2
  ) %>%
  select(label, unique_t1, shared, unique_t2, n12, dice)

cat("\n=== Shared vs unique causal variants ===\n")
print(venn_data, row.names = FALSE)
write.csv(venn_data, file.path(OUT_DIR, "mixer_venn_summary.csv"), row.names = FALSE)

# ── Stacked bar: proportion shared ───────────────────────────────────────────
if (nrow(df) > 0) {
  long_df <- df %>%
    mutate(trait2_only = n2 - n12, asd_only = n1 - n12) %>%
    select(label, asd_only, n12, trait2_only) %>%
    tidyr::pivot_longer(-label, names_to = "component", values_to = "count") %>%
    mutate(component = factor(component,
      levels = c("asd_only", "n12", "trait2_only"),
      labels = c("ASD only", "Shared", "Gut disease only")))

  p2 <- ggplot(long_df, aes(x = label, y = count, fill = component)) +
    geom_bar(stat = "identity", position = "fill") +
    coord_flip() +
    scale_fill_manual(values = c("ASD only" = "#4E79A7",
                                 "Shared"   = "#F28E2B",
                                 "Gut disease only" = "#59A14F")) +
    labs(x = NULL, y = "Proportion of causal variants", fill = NULL,
         title = "Polygenic overlap composition") +
    theme_bw(base_size = 13)
  ggsave(file.path(OUT_DIR, "mixer_overlap_proportion.pdf"), p2, width = 8, height = 5)
  cat("Saved: mixer_overlap_proportion.pdf\n")
}
