# ============================================================
# 03_gene_annotation.R
# HDL.L significant regions → gene annotation → intersect with shared-gene list
#
# Analysis pipeline：
#
# ============================================================

library(data.table)
library(dplyr)

# ── pathsetting ──────────────────────────────────────────────────
LD_PATH    <- "${PROJECT_ROOT}"
HDLL_DIR   <- "${PROJECT_ROOT}"
OUT_DIR    <- "${PROJECT_ROOT}"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load(file.path(LD_PATH, "HDLL_LOC_snps.RData"))
cat(sprintf(NEWLOC：msg %d msg（hg19 msg）\n, nrow(NEWLOC)))

read_hdll <- function(tag) {
  f <- file.path(HDLL_DIR, paste0("res_HDLL_", tag, ".csv"))
  if (!file.exists(f)) {
    cat(msg：msg，msg：, f, "\n")
    return(NULL)
  }
  d <- fread(f)
  d$pair <- tag
  d
}

res_ibs  <- read_hdll("IBS_GWAS_x_ASD")
res_finn <- read_hdll("IBS_FinnGen_x_ASD")

if (is.null(res_ibs) || is.null(res_finn)) {
  stop(HDL.L msg，msg 02_HDL_L_TRAIT_PAIR.R)
}

merge_with_loc <- function(res) {
  left_join(res, NEWLOC, by = c("chr" = "CHR", "piece" = "piece"))
}
res_ibs  <- merge_with_loc(res_ibs)
res_finn <- merge_with_loc(res_finn)

cat(sprintf(TRAIT_A_GWAS × TRAIT_B  msg：%d\n, nrow(res_ibs)))
cat(sprintf(TRAIT_A_FinnGen × TRAIT_B msg：%d\n, nrow(res_finn)))

# ── Step 3：Filtersignificant regions ─────────────────────────────────────
sig_ibs  <- res_ibs[!is.na(res_ibs$p) & res_ibs$p < 0.05, ]
sig_finn <- res_finn[!is.na(res_finn$p) & res_finn$p < 0.05, ]
cat(sprintf(\nmsg（p<0.05）：\n  IBS_GWAS:    %d msg\n  IBS_FinnGen: %d msg\n,
            nrow(sig_ibs), nrow(sig_finn)))

if ("p" %in% names(res_ibs)) {
  res_ibs$fdr  <- p.adjust(res_ibs$p,  method = "BH")
  res_finn$fdr <- p.adjust(res_finn$p, method = "BH")
  fdr_ibs  <- res_ibs[!is.na(res_ibs$fdr)  & res_ibs$fdr < 0.2, ]
  fdr_finn <- res_finn[!is.na(res_finn$fdr) & res_finn$fdr < 0.2, ]
  cat(sprintf(FDR<0.2：\n  IBS_GWAS:    %d msg\n  IBS_FinnGen: %d msg\n,
              nrow(fdr_ibs), nrow(fdr_finn)))
}

shared_pieces <- inner_join(
  sig_ibs[, c("chr", "piece", "rg", "p", "START", "STOP")],
  sig_finn[, c("chr", "piece", "rg", "p")],
  by = c("chr", "piece"),
  suffix = c("_GWAS", "_FinnGen")
) %>%
  filter(sign(rg_GWAS) == sign(rg_FinnGen))

cat(sprintf(\nmsg：%d\n, nrow(shared_pieces)))

candidate_pieces <- res_ibs %>%
  select(chr, piece, rg_GWAS = rg, p_GWAS = p, START, STOP) %>%
  left_join(
    res_finn %>% select(chr, piece, rg_FinnGen = rg, p_FinnGen = p),
    by = c("chr", "piece")
  ) %>%
  filter((!is.na(p_GWAS) & p_GWAS < 0.05) | (!is.na(p_FinnGen) & p_FinnGen < 0.05)) %>%
  mutate(
    direction_consistent = !is.na(rg_GWAS) & !is.na(rg_FinnGen) &
                           sign(rg_GWAS) == sign(rg_FinnGen),
    min_p = pmin(p_GWAS, p_FinnGen, na.rm = TRUE)
  ) %>%
  arrange(min_p)

cat(sprintf(msg（msgp<0.05）：%d\n, nrow(candidate_pieces)))
if (nrow(candidate_pieces) > 0) {
  print(candidate_pieces[1:min(10, nrow(candidate_pieces)), ], row.names = FALSE)
}

write.csv(candidate_pieces,
          file.path(OUT_DIR, "candidate_pieces.csv"),
          row.names = FALSE)

if (nrow(candidate_pieces) == 0) {
  cat(msg，msg。\n)
  stop(msg HDL.L msg。)
}

map_pieces_to_genes <- function(pieces_df, flank_bp = 0) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(msg biomaRt：BiocManager::install('biomaRt'))
  }
  library(biomaRt)
  cat(\nmsg Ensembl GRCh37 msg...\n)
  mart <- useEnsembl(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     GRCh = 37)  # hg19

  gene_list <- lapply(1:nrow(pieces_df), function(i) {
    row   <- pieces_df[i, ]
    start <- max(1, row$START - flank_bp)
    stop  <- row$STOP + flank_bp
    tryCatch({
      genes <- getBM(
        attributes = c("hgnc_symbol", "ensembl_gene_id",
                       "chromosome_name", "start_position", "end_position",
                       "gene_biotype"),
        filters    = c("chromosome_name", "start", "end"),
        values     = list(as.character(row$chr), start, stop),
        mart       = mart
      )
      if (nrow(genes) > 0) {
        genes$piece_chr   <- row$chr
        genes$piece_num   <- row$piece
        genes$piece_start <- row$START
        genes$piece_stop  <- row$STOP
        genes$p_GWAS      <- row$p_GWAS
        genes$p_FinnGen   <- row$p_FinnGen
        genes$rg_GWAS     <- row$rg_GWAS
        genes$rg_FinnGen  <- row$rg_FinnGen
      }
      genes
    }, error = function(e) {
      cat(sprintf(  chr%s piece%s msg：%s\n, row$chr, row$piece, e$message))
      data.frame()
    })
  })

  gene_df <- do.call(rbind, Filter(nrow, gene_list))
  gene_df <- gene_df[gene_df$gene_biotype == "protein_coding" &
                     gene_df$hgnc_symbol != "", ]
  gene_df
}

cat(\nmsg（biomaRt，hg19）...\n)
all_genes_in_pieces <- tryCatch(
  map_pieces_to_genes(candidate_pieces, flank_bp = 0),
  error = function(e) {
    cat(biomaRt msg：, conditionMessage(e), "\n")
    cat(msg，msg（msg）。\n)
    NULL
  }
)

if (!is.null(all_genes_in_pieces) && nrow(all_genes_in_pieces) > 0) {
  write.csv(all_genes_in_pieces,
            file.path(OUT_DIR, "genes_in_candidate_pieces.csv"),
            row.names = FALSE)
  cat(sprintf(msg %d msg\n,
              length(unique(all_genes_in_pieces$hgnc_symbol))))
} else {
  cat(msg，msg。\n)
}

# ── Step 5：intersect with shared-gene list ───────────────
#
#
tesfaye_70_genes <- c(
  "GENE1", "GENE2", "GENE3"
)
# ─────────────────────────────────────────────────────────────

if (!is.null(all_genes_in_pieces) && nrow(all_genes_in_pieces) > 0 &&
    !all(tesfaye_70_genes == c("GENE1", "GENE2", "GENE3"))) {

  genes_in_pieces <- unique(all_genes_in_pieces$hgnc_symbol)
  overlap_genes   <- intersect(genes_in_pieces, tesfaye_70_genes)

  cat(sprintf(\n external_cohort_2023 70 msg，msg %d msg HDL.L msg：\n,
              length(overlap_genes)))

  if (length(overlap_genes) > 0) {
    cat(paste0("  ", overlap_genes, collapse = "\n"), "\n")

    target_genes_df <- all_genes_in_pieces %>%
      filter(hgnc_symbol %in% overlap_genes) %>%
      select(Gene = hgnc_symbol, Ensembl = ensembl_gene_id,
             Chr = piece_chr, Piece = piece_num,
             Gene_Start = start_position, Gene_End = end_position,
             Piece_Start = piece_start, Piece_Stop = piece_stop,
             rg_TRAIT_A_GWAS = rg_GWAS, p_TRAIT_A_GWAS = p_GWAS,
             rg_FinnGen = rg_FinnGen, p_FinnGen = p_FinnGen) %>%
      distinct(Gene, .keep_all = TRUE) %>%
      arrange(p_IBS_GWAS)

    print(target_genes_df, row.names = FALSE)
    write.csv(target_genes_df,
              file.path(OUT_DIR, "target_genes_Tesfaye_intersection.csv"),
              row.names = FALSE)
    cat(\nmsg:,
        file.path(OUT_DIR, "target_genes_Tesfaye_intersection.csv"), "\n")
  }
} else if (all(tesfaye_70_genes == c("GENE1", "GENE2", "GENE3"))) {
  cat(\n[msg] msg Step 5 msg external_cohort_2023 msg 70 msg。\n)
  cat(msg。\n)
}

if (!is.null(all_genes_in_pieces) && nrow(all_genes_in_pieces) > 0) {
  summary_genes <- all_genes_in_pieces %>%
    group_by(hgnc_symbol, ensembl_gene_id, chromosome_name,
             start_position, end_position) %>%
    summarise(
      n_pieces_overlap = n(),
      min_p_GWAS    = min(p_GWAS,    na.rm = TRUE),
      min_p_FinnGen = min(p_FinnGen, na.rm = TRUE),
      rg_GWAS_range = paste0(round(min(rg_GWAS, na.rm = TRUE), 3),
                             " ~ ",
                             round(max(rg_GWAS, na.rm = TRUE), 3)),
      .groups = "drop"
    ) %>%
    arrange(min_p_GWAS)

  write.csv(summary_genes,
            file.path(OUT_DIR, "all_genes_in_sig_pieces.csv"),
            row.names = FALSE)
  cat(sprintf(\nmsg（%d msg）：%s\n,
              nrow(summary_genes),
              file.path(OUT_DIR, "all_genes_in_sig_pieces.csv")))
}

cat(\nmsg。\n)

# ═══════════════════════════════════════════════════════════════
# ───────────────────────────────────────────────────────────────
# library(GenomicFeatures)
# txdb <- makeTxDbFromGFF("/path/to/Homo_sapiens.GRCh37.87.gtf")
# genes_gr <- genes(txdb)
#
# library(GenomicRanges)
# pieces_gr <- GRanges(
#   seqnames = paste0("chr", candidate_pieces$chr),
#   ranges   = IRanges(start = candidate_pieces$START,
#                      end   = candidate_pieces$STOP)
# )
# overlaps <- findOverlaps(pieces_gr, genes_gr)
# ═══════════════════════════════════════════════════════════════
