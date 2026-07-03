# ============================================================
# coloc_GWAS_GWAS.R
# GWAS-GWAS colocalization analysis (coloc package)
#
# Candidate region sources (3 classes):
#   A) Local genetic correlation (HDL.L) reliable-significant regions
#      reproduced in >=2 datasets, CI excludes 0
#   B) Trait A GWAS genome-wide significant loci (p < 5e-8, +/-500 kb)
#   C) Trait B GWAS genome-wide significant loci (p < 5e-8, +/-500 kb)
#
# Output: posterior probabilities H0-H4 per region;
#         H4 >= 0.5 = considered colocalised.
# ============================================================

library(coloc)
library(data.table)

# ── path ─────────────────────────────────────────────────────
DATA_DIR   <- "${PROJECT_ROOT}"
BIM_FILE   <- "${PROJECT_ROOT}"
OUT_DIR    <- "${PROJECT_ROOT}"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Sample size（verify against original）────────────────────────────────────
IBS_N     <- 486601
IBS_CASES <- 53400
IBS_S     <- IBS_CASES / IBS_N
ASD_N     <- 46351
ASD_CASES <- 18382
ASD_S     <- ASD_CASES / ASD_N

# ── Read GWAS data ────────────────────────────────────────
cat(msg GWAS msg...\n)

ibs_raw <- fread(file.path(DATA_DIR, "gwas_TRAIT_A.txt"), showProgress = FALSE)
ibs_raw[, Z := BETA / SE]
ibs_raw[, P := 2 * pnorm(-abs(Z))]
ibs_raw[, A1 := toupper(A1)][, A2 := toupper(A2)]

bim <- fread(BIM_FILE, header = FALSE,
             col.names = c("CHR_ref", "SNP", "CM", "BP", "A1_ref", "A2_ref"),
             showProgress = FALSE)
bim[, CHR_ref := as.integer(CHR_ref)]

ibs_dt <- merge(ibs_raw[, .(SNP, CHR, A1, A2, BETA, SE, P, N)],
                bim[, .(SNP, CHR_ref, BP)],
                by = "SNP", all.x = FALSE)
ibs_dt <- ibs_dt[CHR == CHR_ref & !is.na(P) & is.finite(BETA)]
ibs_dt[, CHR := as.integer(CHR)]
cat(sprintf(  IBS：%d SNPs（msg）\n, nrow(ibs_dt)))

asd_dt <- fread(file.path(DATA_DIR, "gwas_TRAIT_B.txt"), showProgress = FALSE)
asd_dt[, A1 := toupper(A1)][, A2 := toupper(A2)]
asd_dt <- asd_dt[!is.na(P) & is.finite(BETA)]
asd_dt[, CHR := as.integer(CHR)]
setnames(asd_dt, "position", "BP")
cat(sprintf("  ASD：%d SNPs\n", nrow(asd_dt)))


ibd_asd_regions <- data.frame(
  label   = c("chr6p116_IBD-ASD", "chr1p175_IBD-ASD",
              "chr10p98_IBD-ASD", "chr1p159_IBD-ASD",
              "chr10p72_IBD-ASD", "chr16p59_IBD-ASD"),
  chr     = c(6, 1, 10, 1, 10, 16),
  start   = c(129850179, 224943338, 107877788, 205917549, 77144962,  80559686),
  stop    = c(130550137, 226534711, 108725997, 208162951, 78665481,  81437861),
  direction = c("pos", "neg", "pos", "pos", "pos", "pos"),
  stringsAsFactors = FALSE
)

lava_regions <- data.frame(
  label   = c("chr14L1963_LAVA", "chr6L928_LAVA",  "chr12L1835_LAVA",
              "chr3L431_IL17RC_LAVA", "chr8L1282_LAVA", "chr14L1998_LAVA",
              "chr9L1466_LAVA",  "chr14L2013_LAVA", "chr15L2089_LAVA",
              "chr20L2419_LAVA", "chr8L1273_LAVA",  "chr3L435_LAVA",
              "chr10L1602_LAVA", "chr12L1858_LAVA"),
  chr     = c(14L,  6L, 12L,  3L,  8L, 14L,  9L, 14L, 15L, 20L,  8L,  3L, 10L, 12L),
  start   = c(29029225L,  2746532L, 105079101L,  8664893L, 40179578L,
              67891840L, 121008530L, 85562937L,  89385188L, 54743695L,
              27406512L, 12859210L, 129134739L, 129887419L),
  stop    = c(30831154L,  3964072L, 106153787L,  9970731L, 40974254L,
              68976912L, 122677719L, 86653343L,  90632718L, 55448075L,
              28344176L, 14312007L, 129831969L, 130629840L),
  direction = c("pos","neg","pos","neg","pos","pos","pos","pos","pos",
                "neg","pos","neg","pos","neg"),
  stringsAsFactors = FALSE
)

FLANK <- 500000L
sig_ibs[, window := BP %/% 1000000L]
lead_ibs <- sig_ibs[, .SD[which.min(P)], by = .(CHR, window)]

ibs_regions <- data.frame(
  label = paste0("chr", lead_ibs$CHR, "_IBS_p", 
                 formatC(lead_ibs$P, format = "e", digits = 1)),
  chr   = lead_ibs$CHR,
  start = pmax(1L, lead_ibs$BP - FLANK),
  stop  = lead_ibs$BP + FLANK,
  direction = NA_character_,
  stringsAsFactors = FALSE
)

# C) TRAIT_B GWAS genome-wide significant loci
sig_asd <- asd_dt[P < 5e-8][order(CHR, P)]
sig_asd[, window := BP %/% 1000000L]
lead_asd <- sig_asd[, .SD[which.min(P)], by = .(CHR, window)]

asd_regions <- data.frame(
  label = paste0("chr", lead_asd$CHR, "_ASD_p", 
                 formatC(lead_asd$P, format = "e", digits = 1)),
  chr   = lead_asd$CHR,
  start = pmax(1L, lead_asd$BP - FLANK),
  stop  = lead_asd$BP + FLANK,
  direction = NA_character_,
  stringsAsFactors = FALSE
)

all_regions <- rbind(ibd_asd_regions, lava_regions, ibs_regions, asd_regions)
cat(sprintf(\nmsg：%d\n  A) IBD×ASD HDL.L：%d\n  D) LAVA TRAIT_A x TRAIT_B：%d\n  B) TRAIT_A significant loci：%d\n  C) TRAIT_B significant loci：%d\n,
            nrow(all_regions), nrow(ibd_asd_regions), nrow(lava_regions),
            nrow(ibs_regions), nrow(asd_regions)))

run_coloc_gwas <- function(gwas1, gwas2, n1, s1, n2, s2, label,
                           chr, start, stop) {
  g1 <- gwas1[CHR == chr & BP >= start & BP <= stop & !is.na(BETA) & !is.na(SE)]
  g2 <- gwas2[CHR == chr & BP >= start & BP <= stop & !is.na(BETA) & !is.na(SE)]

  common <- intersect(g1$SNP, g2$SNP)
  if (length(common) < 50) {
    return(data.table(
      label    = label, chr = chr, start = start, stop = stop,
      nSNP     = length(common),
      PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
      note  = "insufficient_SNPs"
    ))
  }

  g1 <- g1[SNP %in% common]; setkey(g1, SNP)
  g2 <- g2[SNP %in% common]; setkey(g2, SNP)
  g2 <- g2[g1$SNP]

  d1 <- list(
    beta     = g1$BETA,
    varbeta  = g1$SE^2,
    snp      = g1$SNP,
    position = g1$BP,
    type     = "cc",
    N        = n1,
    s        = s1
  )
  d2 <- list(
    beta     = g2$BETA,
    varbeta  = g2$SE^2,
    snp      = g2$SNP,
    position = g2$BP,
    type     = "cc",
    N        = n2,
    s        = s2
  )

  res <- tryCatch(coloc.abf(d1, d2), error = function(e) {
    message(coloc msg [, label, "]: ", e$message)
    NULL
  })
  if (is.null(res)) {
    return(data.table(label=label, chr=chr, start=start, stop=stop,
                      nSNP=length(common), PP.H0=NA, PP.H1=NA,
                      PP.H2=NA, PP.H3=NA, PP.H4=NA, note="coloc_error"))
  }

  pp <- res$summary
  out <- data.table(
    label    = label,
    chr      = chr,
    start    = start,
    stop     = stop,
    nSNP     = pp["nsnps"],
    PP.H0    = round(pp["PP.H0.abf"], 4),
    PP.H1    = round(pp["PP.H1.abf"], 4),
    PP.H2    = round(pp["PP.H2.abf"], 4),
    PP.H3    = round(pp["PP.H3.abf"], 4),
    PP.H4    = round(pp["PP.H4.abf"], 4),
    note     = "ok"
  )

  top_snp <- res$results[which.max(res$results$SNP.PP.H4), "snp"]
  cat(sprintf("  %-45s  SNPs=%d  H4=%.3f  top=%s\n",
              label, pp["nsnps"], pp["PP.H4.abf"],
              ifelse(length(top_snp) > 0, top_snp, "NA")))
  out
}

cat(\nmsg coloc msg（TRAIT_PAIR）...\n)
cat(sprintf("%-45s  %s\n", "Region", "Result"))
cat(strrep("-", 70), "\n")

coloc_results <- vector("list", nrow(all_regions))

for (i in seq_len(nrow(all_regions))) {
  r <- all_regions[i, ]
  coloc_results[[i]] <- run_coloc_gwas(
    gwas1 = ibs_dt,   gwas2 = asd_dt,
    n1    = IBS_N,    s1    = IBS_S,
    n2    = ASD_N,    s2    = ASD_S,
    label = r$label,
    chr   = r$chr,    start = r$start,   stop = r$stop
  )
}

# ── summary ─────────────────────────────────────────────────
summary_dt <- rbindlist(coloc_results)
setorder(summary_dt, -PP.H4)

out_file <- file.path(OUT_DIR, "coloc_TRAIT_PAIR_all.txt")
fwrite(summary_dt, out_file, sep = "\t")

cat("\n\n", strrep("=", 70), "\n", sep = "")
cat(coloc TRAIT_A x TRAIT_B msg（msg H4 msg）\n)
cat(strrep("=", 70), "\n")
print(summary_dt[note == "ok"][order(-PP.H4)], row.names = FALSE)

coloc_pos <- summary_dt[!is.na(PP.H4) & PP.H4 >= 0.5]
cat(sprintf(\nmsg（H4 ≥ 0.5）：%d msg\n, nrow(coloc_pos)))
if (nrow(coloc_pos) > 0) {
  print(coloc_pos[, .(label, chr, start, stop, nSNP, PP.H3, PP.H4)],
        row.names = FALSE)
  fwrite(coloc_pos, file.path(OUT_DIR, "coloc_TRAIT_PAIR_H4pos.txt"), sep = "\t")
}

cat(\nmsg：, out_file, "\n")
cat(coloc TRAIT_A x TRAIT_B msg。\n)
