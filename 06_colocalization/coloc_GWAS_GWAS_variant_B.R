#!/usr/bin/env Rscript
# ============================================================
# coloc_GWAS_GWAS_Anxiety.R
# GWAS-GWAS colocalization analysis (coloc package)
#
# Candidate region sources:
#   A) HDL.L significant regions from local genetic correlation
#      (output of 03/04_HDL_L_local_rg_*.R — run those first)
#      * Fill in the significant region file path below before running.
#   B) Trait A GWAS genome-wide significant loci (p < 5e-8, +/-500 kb)
#   C) Trait B GWAS genome-wide significant loci (p < 5e-8, +/-500 kb)
# ============================================================

library(coloc)
library(data.table)

DATA_DIR <- "${PROJECT_ROOT}"
BIM_FILE <- "${PROJECT_ROOT}"
OUT_DIR  <- "${PROJECT_ROOT}"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Sample size ──────────────────────────────────────────────────
IBS_N     <- 486601
IBS_CASES <- 53400
IBS_S     <- IBS_CASES / IBS_N   # ≈ 0.110

ANX_N     <- 418399
ANX_S     <- 0.1434               # Nca / (Nca + Nco) ≈ 0.1434

# ── Read GWAS data ────────────────────────────────────────
cat(msg GWAS msg...\n)

ibs_raw <- fread(file.path(DATA_DIR, "gwas_TRAIT_A.txt"), showProgress = FALSE)
ibs_raw[, Z := BETA / SE]
ibs_raw[, P := 2 * pnorm(-abs(Z))]
ibs_raw[, A1 := toupper(A1)][, A2 := toupper(A2)]

bim <- fread(BIM_FILE, header = FALSE,
             col.names = c("CHR_ref","SNP","CM","BP","A1_ref","A2_ref"),
             showProgress = FALSE)
bim[, CHR_ref := as.integer(CHR_ref)]

ibs_dt <- merge(ibs_raw[, .(SNP, CHR, A1, A2, BETA, SE, P, N)],
                bim[, .(SNP, CHR_ref, BP)], by = "SNP", all.x = FALSE)
ibs_dt <- ibs_dt[CHR == CHR_ref & !is.na(P) & is.finite(BETA)]
ibs_dt[, CHR := as.integer(CHR)]
cat(sprintf(  IBS：%d SNPs（msg）\n, nrow(ibs_dt)))

anx_dt <- fread(file.path(DATA_DIR, "gwas_anxiety_fmt.txt"), showProgress = FALSE)
anx_dt[, A1 := toupper(A1)][, A2 := toupper(A2)]
anx_dt <- anx_dt[!is.na(P) & is.finite(BETA)]
anx_dt[, CHR := as.integer(CHR)]
setnames(anx_dt, "position", "BP")
cat(sprintf(  msg：%d SNPs\n, nrow(anx_dt)))


# A) trait pair (A x B) HDL.L significant regions
hdll_regions <- data.frame(
  label = character(0), chr = integer(0),
  start = integer(0),   stop = integer(0),
  stringsAsFactors = FALSE
)
# hdll_regions <- data.frame(
#   label = c("chrX_HDL-L"),
#   chr   = c(X),
#   start = c(XXXXXXX),
#   stop  = c(XXXXXXX),
#   stringsAsFactors = FALSE
# )

FLANK <- 500000L

# B) TRAIT_A GWAS genome-wide significant loci
sig_ibs   <- ibs_dt[P < 5e-8][order(CHR, P)]
if (nrow(sig_ibs) > 0) {
  sig_ibs[, window := BP %/% 1000000L]
  lead_ibs <- sig_ibs[, .SD[which.min(P)], by = .(CHR, window)]
  ibs_regions <- data.frame(
    label = paste0("chr", lead_ibs$CHR, "_IBS_p",
                   formatC(lead_ibs$P, format = "e", digits = 1)),
    chr   = lead_ibs$CHR,
    start = pmax(1L, lead_ibs$BP - FLANK),
    stop  = lead_ibs$BP + FLANK,
    stringsAsFactors = FALSE
  )
  cat(sprintf(  TRAIT_A significant locimsg：%d msg\n, nrow(ibs_regions)))
} else {
  ibs_regions <- data.frame(label=character(0), chr=integer(0),
                             start=integer(0), stop=integer(0))
  cat(  TRAIT_A msg（p<5e-8）\n)
}

# C) disease_anxiety GWAS genome-wide significant loci
sig_anx <- anx_dt[P < 5e-8][order(CHR, P)]
if (nrow(sig_anx) > 0) {
  sig_anx[, window := BP %/% 1000000L]
  lead_anx <- sig_anx[, .SD[which.min(P)], by = .(CHR, window)]
  anx_regions <- data.frame(
    label = paste0("chr", lead_anx$CHR, "_Anx_p",
                   formatC(lead_anx$P, format = "e", digits = 1)),
    chr   = lead_anx$CHR,
    start = pmax(1L, lead_anx$BP - FLANK),
    stop  = lead_anx$BP + FLANK,
    stringsAsFactors = FALSE
  )
  cat(sprintf(  disease_anxiety significant locimsg：%d msg\n, nrow(anx_regions)))
} else {
  anx_regions <- data.frame(label=character(0), chr=integer(0),
                             start=integer(0), stop=integer(0))
  cat(  disease_anxiety has no genome-wide significant locus（p<5e-8），msg p<1e-5\n)
  sig_anx2 <- anx_dt[P < 1e-5][order(CHR, P)]
  if (nrow(sig_anx2) > 0) {
    sig_anx2[, window := BP %/% 1000000L]
    lead_anx2 <- sig_anx2[, .SD[which.min(P)], by = .(CHR, window)]
    anx_regions <- data.frame(
      label = paste0("chr", lead_anx2$CHR, "_Anx_p",
                     formatC(lead_anx2$P, format = "e", digits = 1)),
      chr   = lead_anx2$CHR,
      start = pmax(1L, lead_anx2$BP - FLANK),
      stop  = lead_anx2$BP + FLANK,
      stringsAsFactors = FALSE
    )
    cat(sprintf(  msg p<1e-5：%d msg\n, nrow(anx_regions)))
  }
}

all_regions <- rbind(hdll_regions, ibs_regions, anx_regions)
cat(sprintf(\nmsg：%d（HDL.L:%d  IBS:%d  disease_anxiety:%d）\n,
            nrow(all_regions), nrow(hdll_regions),
            nrow(ibs_regions), nrow(anx_regions)))

run_coloc_gwas <- function(gwas1, gwas2, n1, s1, n2, s2, label, chr, start, stop) {
  g1 <- gwas1[CHR == chr & BP >= start & BP <= stop & !is.na(BETA) & !is.na(SE)]
  g2 <- gwas2[CHR == chr & BP >= start & BP <= stop & !is.na(BETA) & !is.na(SE)]

  common <- intersect(g1$SNP, g2$SNP)
  if (length(common) < 50) {
    return(data.table(label=label, chr=chr, start=start, stop=stop,
                      nSNP=length(common), PP.H0=NA, PP.H1=NA, PP.H2=NA,
                      PP.H3=NA, PP.H4=NA, note="insufficient_SNPs"))
  }
  g1 <- g1[SNP %in% common]; setkey(g1, SNP)
  g2 <- g2[SNP %in% common]; setkey(g2, SNP)
  g2 <- g2[g1$SNP]

  d1 <- list(beta=g1$BETA, varbeta=g1$SE^2, snp=g1$SNP, position=g1$BP,
             type="cc", N=n1, s=s1)
  d2 <- list(beta=g2$BETA, varbeta=g2$SE^2, snp=g2$SNP, position=g2$BP,
             type="cc", N=n2, s=s2)

  res <- tryCatch(coloc.abf(d1, d2), error = function(e) {
    message(coloc msg [, label, "]: ", e$message); NULL
  })
  if (is.null(res))
    return(data.table(label=label, chr=chr, start=start, stop=stop,
                      nSNP=length(common), PP.H0=NA, PP.H1=NA,
                      PP.H2=NA, PP.H3=NA, PP.H4=NA, note="coloc_error"))

  pp <- res$summary
  top_snp <- res$results[which.max(res$results$SNP.PP.H4), "snp"]
  cat(sprintf("  %-45s  SNPs=%d  H4=%.3f  top=%s\n",
              label, pp["nsnps"], pp["PP.H4.abf"],
              ifelse(length(top_snp) > 0, top_snp, "NA")))

  data.table(label=label, chr=chr, start=start, stop=stop,
             nSNP=pp["nsnps"], PP.H0=round(pp["PP.H0.abf"],4),
             PP.H1=round(pp["PP.H1.abf"],4), PP.H2=round(pp["PP.H2.abf"],4),
             PP.H3=round(pp["PP.H3.abf"],4), PP.H4=round(pp["PP.H4.abf"],4),
             note="ok")
}

cat(\nmsg coloc msg（TRAIT_A x disease_anxiety）...\n)
cat(sprintf("%-45s  %s\n", "Region", "Result"))
cat(strrep("-", 70), "\n")

coloc_results <- vector("list", nrow(all_regions))
for (i in seq_len(nrow(all_regions))) {
  r <- all_regions[i, ]
  coloc_results[[i]] <- run_coloc_gwas(
    gwas1 = ibs_dt, gwas2 = anx_dt,
    n1    = IBS_N,  s1    = IBS_S,
    n2    = ANX_N,  s2    = ANX_S,
    label = r$label, chr = r$chr, start = r$start, stop = r$stop
  )
}

# ── summary ─────────────────────────────────────────────────
if (length(coloc_results) > 0 && any(!sapply(coloc_results, is.null))) {
  summary_dt <- rbindlist(coloc_results)
  setorder(summary_dt, -PP.H4)
  out_file <- file.path(OUT_DIR, "coloc_IBS_Anxiety_all.txt")
  fwrite(summary_dt, out_file, sep = "\t")

  cat("\n\n", strrep("=", 70), "\n", sep = "")
  cat(coloc TRAIT_A x disease_anxiety msg（msg H4 msg）\n)
  cat(strrep("=", 70), "\n")
  print(summary_dt[note == "ok"][order(-PP.H4)], row.names = FALSE)

  coloc_pos <- summary_dt[!is.na(PP.H4) & PP.H4 >= 0.5]
  cat(sprintf(\nmsg（H4 ≥ 0.5）：%d msg\n, nrow(coloc_pos)))
  if (nrow(coloc_pos) > 0) {
    print(coloc_pos[, .(label, chr, start, stop, nSNP, PP.H3, PP.H4)], row.names = FALSE)
    fwrite(coloc_pos, file.path(OUT_DIR, "coloc_IBS_Anxiety_H4pos.txt"), sep = "\t")
  }
  cat(\nmsg：, out_file, "\n")
} else {
  cat(\n⚠️  msg SNP msg，msg HDL.L msg hdll_regions\n)
}

cat(coloc TRAIT_A x disease_anxietymsg。\n)
