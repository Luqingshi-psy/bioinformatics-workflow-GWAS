#!/usr/bin/env Rscript
# Parse and visualize MiXeR bivariate results for TRAIT_B × gut diseases
# Input:  output/MiXeR/results/*.fit.json (combined replicates from mixer_figures.py combine)
# Output: output/MiXeR/mixer_summary.csv  +  venn diagrams
#
# WSL paths; assumes script is run from /mnt/c/Users/Administrator/Desktop/trait pair (A x B)

library(jsonlite)
library(ggplot2)
library(dplyr)

# ---- Detect project root (works whether run from project dir or absolute path) ----
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_arg) > 0) {
    SCRIPT_DIR <- dirname(normalizePath(script_arg))
} else {
    SCRIPT_DIR <- getwd()
}
# analysis/ is one level deep
PROJ_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
MIXER_DIR <- file.path(PROJ_DIR, "output", "MiXeR", "results")
OUT_DIR   <- file.path(PROJ_DIR, "output", "MiXeR")

cat("Project root:", PROJ_DIR, "\n")
cat("MiXeR dir   :", MIXER_DIR, "\n")
cat("Output dir  :", OUT_DIR, "\n\n")

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
# After mixer_figures.py combine, fields are:
#   ci.nc1.mean / std          - number of causal variants for trait 1 (× 1000)
#   ci.nc2.mean / std          - number of causal variants for trait 2
#   ci.nc12.mean / std         - number of SHARED causal variants
#   ci.dice.mean / std         - Dice coefficient
#   ci.rho_beta.mean / std     - genetic-effect correlation (per MiXeR bivariate)
parse_bivar_json <- function(json_file, label) {
  if (!file.exists(json_file)) {
    message("Missing: ", json_file)
    return(NULL)
  }
  j <- fromJSON(json_file)
  if (!"ci" %in% names(j)) {
    message("No 'ci' field in ", json_file)
    return(NULL)
  }

  get_ci <- function(field) {
    if (field %in% names(j$ci)) j$ci[[field]] else list(mean=NA, std=NA, median=NA)
  }

  nc1   <- get_ci("nc1")
  nc2   <- get_ci("nc2")
  nc12  <- get_ci("nc12")
  dice  <- get_ci("dice")
  rho_beta <- get_ci("rho_beta")
  rg    <- get_ci("rg")

  data.frame(
    label   = label,
    nc1     = nc1$mean,    nc1_se    = nc1$std,
    nc2     = nc2$mean,    nc2_se    = nc2$std,
    nc12    = nc12$mean,   nc12_se   = nc12$std,
    dice    = dice$mean,   dice_se   = dice$std,
    rho_ge  = rho_beta$mean, rho_ge_se = rho_beta$std,  # genetic-effect correlation
    rg      = rg$mean,     rg_se     = rg$std,
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
cat("\n=== MiXeR Bivariate Results (combined across 20 reps) ===\n")
print(df[, c("label", "nc12", "nc12_se", "dice", "dice_se", "rg", "rg_se", "rho_ge", "rho_ge_se")],
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
# nc1, nc2, nc12 are in units of thousands of causal variants
venn_data <- df %>%
  mutate(
    unique_t1  = nc1 - nc12,
    unique_t2  = nc2 - nc12,
    shared     = nc12,
    total_t1   = nc1,
    total_t2   = nc2
  ) %>%
  select(label, unique_t1, shared, unique_t2, nc12, dice)

cat("\n=== Shared vs unique causal variants (in thousands) ===\n")
print(venn_data, row.names = FALSE)
write.csv(venn_data, file.path(OUT_DIR, "mixer_venn_summary.csv"), row.names = FALSE)

# ── Stacked bar: proportion shared ──────────────────────────────────────────
if (nrow(df) > 0) {
  long_df <- df %>%
    mutate(trait2_only = nc2 - nc12, asd_only = nc1 - nc12) %>%
    select(label, asd_only, nc12, trait2_only) %>%
    tidyr::pivot_longer(-label, names_to = "component", values_to = "count") %>%
    mutate(component = factor(component,
      levels = c("asd_only", "nc12", "trait2_only"),
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
