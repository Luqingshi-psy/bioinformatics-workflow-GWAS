# ============================================================
# 03_LAVA_trait pair (A x B).R
# trait pair (A x B) Local genetic correlation analysis（LAVA）
#
#
# required：
#   - devtools::install_github("josefin-werme/LAVA")
#   - 1000G EUR Reference panel：${PROJECT_ROOT},bim,fam}
#   - loci partition file：${PROJECT_ROOT}
#   - TRAIT_B extracted file：${PROJECT_ROOT}
# ============================================================

library(LAVA)
library(data.table)

# ── path ─────────────────────────────────────────────────────
DATA_DIR   <- "${PROJECT_ROOT}"
REF_PREFIX <- "${PROJECT_ROOT}"
LOC_FILE   <- "${PROJECT_ROOT}"
ASD_FILE   <- "${PROJECT_ROOT}"
OUT_DIR    <- "${PROJECT_ROOT}"
PREP_DIR   <- "${PROJECT_ROOT}"

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PREP_DIR, showWarnings = FALSE, recursive = TRUE)

# ── TRAIT_A Sample size（Diekstra et al. 2021, Nat Commun）───────────
# verify against original：IBS_cases + IBS_controls = N_total
IBS_CASES    <- 53400
IBS_CONTROLS <- 433201    # 53400 + 433201 = 486601

ibs_prepped_file <- file.path(PREP_DIR, "IBS_GWAS_lava.txt")

if (!file.exists(ibs_prepped_file)) {
  cat(Prepare TRAIT_A GWAS msg（msg）...\n)

  ibs_raw <- fread(file.path(DATA_DIR, "gwas_TRAIT_A.txt"), showProgress = FALSE)
  ibs_raw[, Z := BETA / SE]
  ibs_raw[, P := 2 * pnorm(-abs(Z))]

  bim <- fread(paste0(REF_PREFIX, ".bim"), header = FALSE,
               col.names = c("CHR", "SNP", "CM", "BP", "A1_ref", "A2_ref"),
               showProgress = FALSE)
  bim_coord <- bim[, .(SNP, position = BP)]

  ibs_m <- merge(ibs_raw, bim_coord, by = "SNP", all.x = FALSE)
  # LAVA format：SNP, CHR, position, A1, A2, BETA, SE, P, N
  ibs_out <- ibs_m[, .(SNP, CHR, position, A1, A2, BETA, SE, P, N)]
  ibs_out <- ibs_out[!is.na(P) & is.finite(BETA) & is.finite(SE)]

  fwrite(ibs_out, ibs_prepped_file, sep = "\t")
  cat(sprintf(  TRAIT_A retains %d SNPs（msg %d，msg %d）\n,
              nrow(ibs_out), nrow(ibs_raw), nrow(ibs_m)))
} else {
  cat(IBS GWAS msg，msg。\n)
}

# ── Step 2：generate input.info.txt ──────────────────────────
info_file <- file.path(PREP_DIR, "input.info.txt")

info_df <- data.frame(
  phenotype = c("TRAIT_A_GWAS", "TRAIT_B"),
  cases     = c(IBS_CASES, 18382L),
  controls  = c(IBS_CONTROLS, 27969L),
  filename  = c(ibs_prepped_file, ASD_FILE),
  stringsAsFactors = FALSE
)
write.table(info_df, info_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(input.info.txt msg：\n); print(info_df)

cat(\nmsg LAVA msg...\n)
input <- process.input(
  input.info.file     = info_file,
  sample.overlap.file = NULL,
  ref.prefix          = REF_PREFIX,
  phenos              = c("TRAIT_A_GWAS", "TRAIT_B")
)
cat(msg。phenotypes:, paste(names(input$sum.stats), collapse=", "), "\n")

loci <- read.loci(LOC_FILE)
n_loci <- nrow(loci)
cat(sprintf(msg：%d\n, n_loci))

cat(\nmsg LAVA msg（msg：ASD）...\n)

univ_thresh    <- 0.05 / n_loci
bivar_results  <- list()
n_tested       <- 0L

for (i in seq_len(n_loci)) {
  if (i %% 250 == 0)
    cat(sprintf(  msg：%d / %d  (%.1f%%)  msg：%d\n,
                i, n_loci, 100 * i / n_loci, n_tested))

  locus <- tryCatch(
    process.locus(loci[i, ], input),
    error = function(e) {
      message(process.locus msg locus , i, ": ", e$message)
      NULL
    }
  )
  if (is.null(locus)) next
  if (!"TRAIT_B" %in% locus$phenos) next

  bivar <- tryCatch(
    run.bivar(locus, target = "TRAIT_B"),
    error = function(e) NULL
  )
  if (is.null(bivar)) next

  n_tested <- n_tested + nrow(bivar)
  bivar_results[[i]] <- cbind(
    locus_id = i,
    chr      = locus$chr,
    start    = locus$start,
    stop     = locus$stop,
    n_snps   = locus$n.snps,
    bivar
  )
}

bivar_df <- do.call(rbind, bivar_results[!sapply(bivar_results, is.null)])
bivar_df  <- bivar_df[!is.na(bivar_df$p), ]

bivar_df$p_adj_fdr  <- p.adjust(bivar_df$p, method = "BH")
bivar_df$p_adj_bonf <- p.adjust(bivar_df$p, method = "bonferroni")

cat(sprintf(\nmsg。msg %d msg。\n, nrow(bivar_df)))
cat(sprintf(  msg (p < 0.05):         %d\n, sum(bivar_df$p < 0.05)))
cat(sprintf(  msg (FDR < 0.05):       %d\n, sum(bivar_df$p_adj_fdr < 0.05)))
cat(sprintf(  msg (Bonferroni < 0.05):%d\n, sum(bivar_df$p_adj_bonf < 0.05)))

out_all  <- file.path(OUT_DIR, "LAVA_TRAIT_PAIR_all.txt")
out_sig  <- file.path(OUT_DIR, "LAVA_TRAIT_PAIR_significant.txt")

write.table(bivar_df, out_all, sep = "\t", row.names = FALSE, quote = FALSE)
cat(\nmsg：, out_all, "\n")

if (sum(bivar_df$p < 0.05) > 0) {
  sig <- subset(bivar_df, p < 0.05)
  sig <- sig[order(sig$p), ]
  write.table(sig, out_sig, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(msg：, out_sig, "\n")

  cat(\n=== msg（p < 0.05）===\n)
  print(sig[, c("phen1", "phen2", "chr", "start", "stop",
                "rho", "rho.lower", "rho.upper", "p", "p_adj_fdr")],
        row.names = FALSE)

  # FDR significantlocus
  sig_fdr <- subset(bivar_df, p_adj_fdr < 0.05)
  if (nrow(sig_fdr) > 0) {
    cat(sprintf(\n=== FDR msg（%d msg）===\n, nrow(sig_fdr)))
    print(sig_fdr[order(sig_fdr$p),
                  c("phen1", "phen2", "chr", "start", "stop",
                    "rho", "p", "p_adj_fdr")],
          row.names = FALSE)
  }
}

cat(\nLAVA TRAIT_A x TRAIT_B msg。\n)
