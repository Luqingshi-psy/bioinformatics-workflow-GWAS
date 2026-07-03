#!/usr/bin/env Rscript
# ============================================================
# 03_PLACO_trait pair (A x B).R
# ============================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(parallel)
})

PLACO_SRC <- "${PROJECT_ROOT}"
IBS_FILE  <- "${PROJECT_ROOT}"
ANX_FILE  <- "${PROJECT_ROOT}"
OUT_DIR   <- "${PROJECT_ROOT}"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

P_THRESH_NULL     <- 1e-4
P_FDR_CUTOFF      <- 0.05
CHECKPOINT_FILE   <- file.path(OUT_DIR, "placo_checkpoint.rds")
CHECKPOINT_INTERVAL <- 100000
N_CORES <- min(2, parallel::detectCores() - 1)
cat(sprintf(msg %d msg CPU msg\n, N_CORES))

suppressMessages(source(PLACO_SRC))
cat(PLACO msg\n)

cat("Read TRAIT_A GWAS...\n")
ibs <- fread(IBS_FILE, showProgress = FALSE)
setnames(ibs, names(ibs), tolower(names(ibs)))
ibs[, z_ibs := as.numeric(beta) / as.numeric(se)]
ibs[, p_ibs := 2 * pnorm(-abs(z_ibs))]
ibs[, a1 := toupper(a1)]; ibs[, a2 := toupper(a2)]
cat(sprintf(  TRAIT_A total SNP msg：%d\n, nrow(ibs)))

cat("Read disease_anxiety GWAS...\n")
anx <- fread(ANX_FILE, showProgress = FALSE)
setnames(anx, names(anx), tolower(names(anx)))
if ("position" %in% names(anx) && !"bp" %in% names(anx))
  setnames(anx, "position", "bp")
anx[, z_anx := as.numeric(beta) / as.numeric(se)]
anx[, a1 := toupper(a1)]; anx[, a2 := toupper(a2)]
cat(sprintf(  disease_anxiety total SNP msg：%d\n, nrow(anx)))

# ── 3. combine ─────────────────────────────────────────────────
cat(\nMerge TRAIT_A x disease_anxietymsg...\n)
ibs_sub <- ibs[, .(snp, a1_ibs=a1, a2_ibs=a2, z_ibs, p_ibs)]
anx_sub <- anx[, .(snp, a1_anx=a1, a2_anx=a2, z_anx, p_anx=p)]
merged  <- merge(ibs_sub, anx_sub, by = "snp")
cat(sprintf(  msg SNP msg：%d\n, nrow(merged)))

cat(msg...\n)
is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT","TA","CG","GC")
merged[, palindromic := is_palindromic(a1_ibs, a2_ibs)]
cat(sprintf(  msg SNP（msg）：%d\n, sum(merged$palindromic)))
merged <- merged[palindromic == FALSE]

merged[, align := fcase(
  a1_ibs == a1_anx & a2_ibs == a2_anx, "same",
  a1_ibs == a2_anx & a2_ibs == a1_anx, "flip",
  default = "mismatch"
)]
cat(sprintf(  msg：%d，msg：%d，msg（msg）：%d\n,
            sum(merged$align == "same"), sum(merged$align == "flip"),
            sum(merged$align == "mismatch")))
merged <- merged[align != "mismatch"]
merged[align == "flip", z_anx := -z_anx]
cat(sprintf(  msg SNP msg：%d\n, nrow(merged)))

cat(\nmsg PLACO msg...\n)
Z.matrix <- as.matrix(merged[, .(z_ibs, z_anx)])
P.matrix <- as.matrix(merged[, .(p_ibs, p_anx)])
colnames(Z.matrix) <- c("z1", "z2")
colnames(P.matrix) <- c("p1", "p2")
cat(sprintf(  msg：%d SNP × 2 msg\n, nrow(Z.matrix)))

cat(sprintf(msg VarZ（p_thresh = %g）...\n, P_THRESH_NULL))
VarZ <- tryCatch(
  var.placo(Z.matrix, P.matrix, p.threshold = P_THRESH_NULL),
  error = function(e) { cat(sprintf(  ⚠️  var.placo msg：%s\n, e$message)); NULL }
)
if (is.null(VarZ)) stop(msg VarZ)
cat(sprintf(  VarZ: TRAIT_A=%.4f，msg=%.4f\n, VarZ[1], VarZ[2]))

cat(msg CorZ...\n)
CorZ <- tryCatch(
  cor.pearson(Z.matrix, P.matrix, p.threshold = P_THRESH_NULL, returnMatrix = FALSE),
  error = function(e) { cat(sprintf(  ⚠️  cor.pearson msg：%s，msg 0\n, e$message)); 0 }
)
cat(sprintf("  CorZ：%.4f\n", CorZ))

process_batch <- function(idx_vec) {
  results <- vector("list", length(idx_vec))
  for (j in seq_along(idx_vec)) {
    i   <- idx_vec[j]
    res <- tryCatch(
      placo.plus(as.numeric(Z.matrix[i, ]), VarZ, CorZ),
      error = function(e) list(T.placo.plus = NA_real_, p.placo.plus = NA_real_)
    )
    results[[j]] <- c(score = res$T.placo.plus, pvalue = res$p.placo.plus)
  }
  results
}

total_snps <- nrow(Z.matrix)
res_list   <- list()
done       <- 0L
goto_final <- FALSE

if (file.exists(CHECKPOINT_FILE)) {
  cp <- readRDS(CHECKPOINT_FILE)
  if (is.list(cp) && all(c("done","results") %in% names(cp))) {
    done     <- cp$done
    res_list <- cp$results
    cat(sprintf(msg：msg %d / %d msg SNP\n, done, total_snps))
    if (done >= total_snps) {
      cat(msg SNP msg，msg。\n)
      placo_results <- res_list[1:total_snps]
      goto_final    <- TRUE
    }
  }
}

if (!goto_final) {
  test_res <- tryCatch(
    placo.plus(as.numeric(Z.matrix[1, ]), VarZ, CorZ),
    error = function(e) { cat(sprintf(⚠️  placo.plus msg：%s\n, e$message)); NULL }
  )
  if (is.null(test_res)) stop(placo.plus msg)
  cat(sprintf(placo.plus msg：T=%.4f, p=%.4f\n,
              test_res$T.placo.plus, test_res$p.placo.plus))

  remaining_idx <- (done + 1):total_snps
  n_remaining   <- length(remaining_idx)
  batch_starts  <- seq(1, n_remaining, by = CHECKPOINT_INTERVAL)
  n_batches     <- length(batch_starts)
  cat(sprintf(msg %d msg SNP msg %d msg batch\n, n_remaining, n_batches))

  cl <- makeCluster(N_CORES)
  clusterExport(cl, c("Z.matrix", "VarZ", "CorZ", "PLACO_SRC", "process_batch"))
  clusterEvalQ(cl, suppressMessages(source(PLACO_SRC)))

  t_start <- proc.time()
  for (b in seq_along(batch_starts)) {
    start_in_batch      <- batch_starts[b]
    end_in_batch        <- min(start_in_batch + CHECKPOINT_INTERVAL - 1, n_remaining)
    batch_orig_idx      <- remaining_idx[start_in_batch:end_in_batch]
    cat(sprintf(msg batch %d / %d：SNP %d～%d（msg %d msg）...\n,
                b, n_batches, min(batch_orig_idx), max(batch_orig_idx), length(batch_orig_idx)))
    chunk_size <- ceiling(length(batch_orig_idx) / N_CORES)
    chunks     <- split(batch_orig_idx, ceiling(seq_along(batch_orig_idx) / chunk_size))
    batch_results <- unlist(parLapply(cl, chunks, process_batch), recursive = FALSE)
    res_list   <- c(res_list, batch_results)
    done       <- length(res_list)
    saveRDS(list(done = done, results = res_list), CHECKPOINT_FILE)
    cat(sprintf(  msg %d / %d SNP\n, done, total_snps))
  }
  elapsed <- (proc.time() - t_start)[3]
  cat(sprintf(PLACO+ msg，msg %.1f msg\n, elapsed))
  stopCluster(cl)
  placo_results <- res_list[1:total_snps]
  file.remove(CHECKPOINT_FILE)
  cat(msg。\n)
}

cat(msg...\n)
placo_dt <- as.data.table(do.call(rbind, placo_results))
setnames(placo_dt, c("placo_score", "placo_p"))

result <- cbind(
  merged[, .(snp, a1 = a1_ibs, a2 = a2_ibs, z_ibs, p_ibs, z_anx, p_anx, align)],
  placo_dt
)
result <- result[!is.na(placo_p)]
result[, p_fdr  := p.adjust(placo_p, method = "fdr")]
result[, p_bonf := p.adjust(placo_p, method = "bonferroni")]
setorder(result, placo_p)

cat(sprintf(  PLACO+ msg，msg SNP：%d\n, nrow(result)))
cat(sprintf(  FDR<0.05：%d msg SNP\n, sum(result$p_fdr < P_FDR_CUTOFF, na.rm = TRUE)))
cat(sprintf(  msg p<0.05：%d msg SNP\n, sum(result$placo_p < 0.05, na.rm = TRUE)))

fwrite(result, file.path(OUT_DIR, "01_PLACO_full_results.txt.gz"), sep = "\t")
sig_fdr <- result[p_fdr < P_FDR_CUTOFF]
sig_nom <- result[placo_p < 0.05]
fwrite(sig_fdr, file.path(OUT_DIR, "02_PLACO_FDR05.txt"),  sep = "\t")
fwrite(sig_nom, file.path(OUT_DIR, "03_PLACO_nom05.txt"),   sep = "\t")
cat(sprintf(  FDR<0.05 msg：%d msg\n  msg：%d msg\n, nrow(sig_fdr), nrow(sig_nom)))

cat(\nmsg/msg...\n)
pos_col <- if ("position" %in% names(anx)) "position" else "bp"
anx_pos <- anx[, .(snp, chr = as.integer(chr), bp = as.integer(get(pos_col)))]
result2  <- merge(result, anx_pos, by = "snp", all.x = TRUE)
result2  <- result2[!is.na(chr) & !is.na(bp)]
setorder(result2, chr, bp)
cat(sprintf(  msg：%d msg SNP\n, nrow(result2)))

cat(msg Manhattan msg...\n)
if (nrow(result2) > 0) {
  set.seed(42)
  man_data <- rbind(
    result2[placo_p <= 0.1],
    result2[placo_p  > 0.1][sample(.N, min(.N, ceiling(.N * 0.10)))]
  )
  setorder(man_data, chr, bp)
  man_data[, chr := as.integer(chr)]
  chr_info <- man_data[, .(chr_len = max(bp, na.rm = TRUE)), by = chr]
  setorder(chr_info, chr)
  chr_info[, chr_start := c(0, cumsum(as.numeric(chr_len))[-nrow(chr_info)])]
  man_data <- merge(man_data, chr_info[, .(chr, chr_start)], by = "chr")
  man_data[, pos_cum := bp + chr_start]
  axis_df  <- man_data[, .(center = mean(pos_cum)), by = chr]
  man_data[, sig_label := fcase(
    p_fdr   < 0.05, "FDR<0.05",
    placo_p < 0.05, "p<0.05",
    default = "NS"
  )]

  p_manhattan <- ggplot(man_data[sig_label == "NS"],
                        aes(x = pos_cum, y = -log10(placo_p))) +
    geom_point(aes(color = factor(chr %% 2)), size = 0.3, alpha = 0.4, show.legend = FALSE) +
    scale_color_manual(values = c("0" = "#4E79A7", "1" = "#A0CBE8")) +
    geom_point(data = man_data[sig_label == "p<0.05"],
               aes(x = pos_cum, y = -log10(placo_p)), color = "#F4A582", size = 0.8) +
    geom_point(data = man_data[sig_label == "FDR<0.05"],
               aes(x = pos_cum, y = -log10(placo_p)), color = "#D6604D", size = 1.2) +
    geom_hline(yintercept = -log10(0.05 / nrow(result2)),
               linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "orange", linewidth = 0.3) +
    scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor    = element_blank(),
          axis.text.x         = element_text(size = 6)) +
    labs(title    = "PLACO+: TRAIT_A x disease_anxiety Shared Pleiotropic SNPs",
         subtitle = sprintf("VarZ=(%.3f,%.3f)  CorZ=%.3f  FDR<0.05: %d SNP",
                            VarZ[1], VarZ[2], CorZ, nrow(sig_fdr)),
         x = "Chromosome", y = expression(-log[10](p["PLACO+"])))

  ggsave(file.path(OUT_DIR, "04_manhattan.pdf"), p_manhattan, width = 12, height = 5)
  cat(  Manhattan msg\n)
} else {
  cat(  msg，msg\n)
}

cat(msg QQ msg...\n)
if (nrow(result2) > 0) {
  set.seed(42)
  qq_data <- rbind(
    result2[placo_p <= 0.1, .(placo_p)],
    result2[placo_p  > 0.1][sample(.N, min(.N, 200000)), .(placo_p)]
  )
  qq_data <- qq_data[!is.na(placo_p)]
  setorder(qq_data, placo_p)
  n_qq <- nrow(qq_data)
  qq_data[, expected := -log10((seq_len(.N) - 0.5) / n_qq)]
  qq_data[, observed := -log10(placo_p)]

  lambda_gc <- round(median(qchisq(result2$placo_p, df = 1, lower.tail = FALSE), na.rm = TRUE) /
                     qchisq(0.5, df = 1, lower.tail = FALSE), 3)
  cat(sprintf("  lambda_GC（PLACO p）= %.3f\n", lambda_gc))

  p_qq <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_point(size = 0.3, color = "#4E79A7", alpha = 0.5) +
    theme_bw(base_size = 10) +
    labs(title = sprintf("QQ Plot (PLACO+, lambda_GC=%.3f)", lambda_gc),
         x = expression("Expected" ~ -log[10](p)),
         y = expression("Observed" ~ -log[10](p["PLACO+"])))

  ggsave(file.path(OUT_DIR, "05_QQ.pdf"), p_qq, width = 6, height = 6)
  cat(  QQ msg\n)
}

# ── 12. summary ─────────────────────────────────────────────────
cat(\n========== PLACO+ msg（TRAIT_A x disease_anxiety）==========\n)
cat(sprintf(msg SNP msg（msg/msg）：%d\n, nrow(result)))
cat(sprintf(VarZ: TRAIT_A=%.4f，msg=%.4f\n, VarZ[1], VarZ[2]))
cat(sprintf("CorZ：%.4f\n", CorZ))
cat(sprintf(msg（p<0.05）：%d\n, nrow(sig_nom)))
cat(sprintf(FDR msg（FDR<0.05）：%d\n, nrow(sig_fdr)))

if (nrow(sig_fdr) > 0) {
  cat(\nFDR<0.05 msg SNP（msg 20 msg）：\n)
  print(sig_fdr[1:min(20, nrow(sig_fdr)),
                .(snp, a1, a2,
                  z_ibs = round(z_ibs, 3),
                  p_ibs = formatC(p_ibs, format = "e", digits = 2),
                  z_anx = round(z_anx, 3),
                  p_anx = formatC(p_anx, format = "e", digits = 2),
                  placo_p = formatC(placo_p, format = "e", digits = 2),
                  p_fdr   = formatC(p_fdr,   format = "e", digits = 2))],
        row.names = FALSE)
}

cat(\n=== msg ===\nmsg：, OUT_DIR, "\n")
