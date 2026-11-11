#!/usr/bin/env Rscript
## _fig4_deg.R
## Compute T-cell + NK DEG tables for Fig4 (scheme A: T main + NK companion).
## Method correction vs naive approach:
##   - Real UMI counts live in the SCT assay "counts" layer (RNA assay data slot is empty).
##   - We build a log-normalized assay (log1p(CPM*1e4)) from SCT counts -> meaningful avg_log2FC.
##   - FindMarkers (wilcox) gives significance on that normalized data.
##   - "APRR-restored" genes = significant in PRRSV-vs-Control AND reversed & significant in APRR-vs-PRRSV.
##   - KEGG enrichment on the restored union (human-symbol proxy, consistent with SCENIC pipeline).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

ROOT <- "data/scRNA-seq/reanalysis_v2"
T_OBJ <- file.path(ROOT, "playbook_outputs/T_cell_subcluster/T_cell_subcluster_3grp.rds")
CANON <- file.path(ROOT, "data/derived/seurat_reanalysis_v2_3grp.rds")
OUT  <- file.path(ROOT, "fig4_data")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat(sprintf("[%s] loading objects...\n", Sys.time()))

## ---- helper: build log1p(CPM*1e4) matrix from SCT counts ----
compute_ln <- function(obj, assay = "SCT") {
  counts <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || sum(counts) == 0) stop("SCT counts layer empty")
  tot <- colSums(counts); tot[tot == 0] <- 1
  log1p(counts / tot * 1e4)
}

## ---- helper: manual Wilcoxon DEG on logNorm matrix ----
## Avoids Seurat's SCT FindMarkers path entirely (which forces PrepSCTFindMarkers
## on Pearson residuals). Gives identical wilcox statistics + BH-adjusted p.
manual_deg <- function(ln, group, ident1, ident2) {
  stopifnot(all(c(ident1, ident2) %in% group))
  c1 <- colnames(ln)[group == ident1]
  c2 <- colnames(ln)[group == ident2]
  m1 <- as.matrix(ln[, c1, drop = FALSE])
  m2 <- as.matrix(ln[, c2, drop = FALSE])
  mean1 <- rowMeans(m1); mean2 <- rowMeans(m2)
  avg_log2FC <- mean1 - mean2
  pct1 <- rowMeans(ln[, c1, drop = FALSE] > 0)
  pct2 <- rowMeans(ln[, c2, drop = FALSE] > 0)
  p <- sapply(seq_len(nrow(ln)), function(i) {
    x <- m1[i, ]; y <- m2[i, ]
    if (var(c(x, y)) == 0) return(1)
    tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA)
  })
  p_adj <- p.adjust(p, method = "BH")
  data.frame(gene = rownames(ln), p_val = p, avg_log2FC = avg_log2FC,
             pct.1 = pct1, pct.2 = pct2, p_val_adj = p_adj,
             stringsAsFactors = FALSE)
}

## ---- helper: run 3 contrasts + restore set ----
run_deg <- function(obj, label, group_col = "group") {
  ln <- compute_ln(obj, "SCT")
  group <- obj@meta.data[[group_col]]
  names(group) <- rownames(obj@meta.data)
  stopifnot(all(colnames(ln) == names(group)))

  contrasts <- list(
    PvC = list(i1 = "PRRSV", i2 = "Control"),
    AvP = list(i1 = "APRR",  i2 = "PRRSV"),
    AvC = list(i1 = "APRR",  i2 = "Control")
  )
  ctag <- c(PvC = "PRRSV_vs_Control", AvP = "APRR_vs_PRRSV", AvC = "APRR_vs_Control")
  res <- list()
  for (nm in names(contrasts)) {
    df <- manual_deg(ln, group, contrasts[[nm]]$i1, contrasts[[nm]]$i2)
    write.csv(df, file.path(OUT, sprintf("%s_DEG_%s.csv", label, ctag[[nm]])), row.names = FALSE)
    res[[nm]] <- df
    cat(sprintf("  %s %s: %d genes, %d sig(adj<0.05)\n", label, nm, nrow(df),
                sum(df$p_val_adj < 0.05, na.rm = TRUE)))
  }

  ## restored set: significant in PvC and reversed & significant in AvP
  pvc <- res$PvC; avp <- res$AvP
  m <- merge(pvc[, c("gene", "avg_log2FC", "p_val_adj")],
             avp[, c("gene", "avg_log2FC", "p_val_adj")],
             by = "gene", suffixes = c(".PvC", ".AvP"))
  m$restored <- (abs(m$avg_log2FC.PvC) >= 0.25 & m$p_val_adj.PvC < 0.05 &
                 sign(m$avg_log2FC.AvP) == -sign(m$avg_log2FC.PvC) &
                 abs(m$avg_log2FC.AvP) >= 0.25 & m$p_val_adj.AvP < 0.05)
  restored <- m$gene[m$restored]
  write.csv(m, file.path(OUT, sprintf("%s_restored_matrix.csv", label)), row.names = FALSE)
  writeLines(restored, file.path(OUT, sprintf("%s_restored_genes.txt", label)))
  cat(sprintf("  %s restored genes (APRR reverses PRRSV): %d\n", label, length(restored)))

  list(res = res, restored = restored)
}

## ---- helper: manual logNorm means per group for curated functional genes ----
functional_logmean <- function(counts, genes, group) {
  avail <- intersect(genes, rownames(counts))
  if (length(avail) < length(genes)) {
    cat(sprintf("  [warn] functional_logmean dropped %d absent genes: %s\n",
                length(genes) - length(avail),
                paste(setdiff(genes, rownames(counts)), collapse = ",")))
  }
  genes <- avail
  tot <- colSums(counts); tot[tot == 0] <- 1
  ln <- log1p(counts / tot * 1e4)
  out <- data.frame(gene = genes)
  for (g in c("Control", "PRRSV", "APRR")) {
    cells <- colnames(ln)[group == g]
    if (length(cells) == 0) { out[[g]] <- NA; next }
    sub <- ln[genes, cells, drop = FALSE]
    out[[g]] <- rowMeans(sub)
  }
  out
}

## ============ T CELLS ============
cat(sprintf("[%s] === T cells ===\n", Sys.time()))
tobj <- readRDS(T_OBJ)
cat(sprintf("  T obj cells=%d  subcluster levels=%d  group levels=%s\n",
            ncol(tobj), length(unique(tobj$subcluster)),
            paste(sort(unique(tobj$group)), collapse = "/")))
tdeg <- run_deg(tobj, "Tcell")

## T UMAP coords + subcluster + group
temb <- tryCatch(Embeddings(tobj, reduction = "umap"),
                 error = function(e) Embeddings(tobj, reduction = "tsne"))
tcoord <- data.frame(
  cell = colnames(tobj),
  umap1 = temb[, 1], umap2 = temb[, 2],
  subcluster = tobj$subcluster,
  group = tobj$group,
  stringsAsFactors = FALSE
)
write.csv(tcoord, file.path(OUT, "Tcell_umap_coords.csv"), row.names = FALSE)

## T functional genes (cytotoxicity + exhaustion + checkpoint, with directionality)
t_func <- c("IFNG", "GNLY", "NKG7", "GZMB", "PRF1", "GZMK",
            "PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2", "TGFB1",
            "IL7R", "TCF7", "CCR7", "CD69", "CD44")
t_counts <- GetAssayData(tobj, assay = "SCT", layer = "counts")
tfunc_mean <- functional_logmean(t_counts, t_func, tobj$group)
write.csv(tfunc_mean, file.path(OUT, "Tcell_functional_logmean.csv"), row.names = FALSE)

## ============ NK CELLS ============
cat(sprintf("[%s] === NK cells ===\n", Sys.time()))
canon <- readRDS(CANON)
nk_cells <- colnames(canon)[canon$cell_type == "NK"]
cat(sprintf("  NK cells in canonical: %d\n", length(nk_cells)))
nkobj <- subset(canon, cells = nk_cells)
cat(sprintf("  NK group sizes: %s\n",
            paste(table(nkobj$group), collapse = "/")))
ndeg <- run_deg(nkobj, "NK")

## NK UMAP coords (from canonical umap)
nemb <- tryCatch(Embeddings(canon, reduction = "umap"),
                 error = function(e) Embeddings(canon, reduction = "tsne"))
ncoord <- data.frame(
  cell = nk_cells,
  umap1 = nemb[nk_cells, 1], umap2 = nemb[nk_cells, 2],
  group = nkobj$group,
  stringsAsFactors = FALSE
)
write.csv(ncoord, file.path(OUT, "NK_umap_coords.csv"), row.names = FALSE)

nk_func <- c("GNLY", "NKG7", "GZMB", "PRF1", "GZMH",
             "KLRK1", "NCR1", "FCGR3A", "KLRC1", "KLRD1",
             "IFNG", "TBX21", "EOMES", "IL2RB", "XCL1")
nk_counts <- GetAssayData(nkobj, assay = "SCT", layer = "counts")
nkfunc_mean <- functional_logmean(nk_counts, nk_func, nkobj$group)
write.csv(nkfunc_mean, file.path(OUT, "NK_functional_logmean.csv"), row.names = FALSE)

## ============ KEGG enrichment on restored unions ============
cat(sprintf("[%s] === KEGG enrichment (restored unions) ===\n", Sys.time()))
all_t_syms <- rownames(GetAssayData(tobj, assay = "SCT", layer = "counts"))
all_nk_syms <- rownames(nk_counts)

run_kegg <- function(syms_restored, all_syms, tag) {
  if (length(syms_restored) < 5) {
    cat(sprintf("  %s: too few restored genes (%d) for KEGG\n", tag, length(syms_restored)))
    return(NULL)
  }
  eg <- bitr(syms_restored, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  uni <- bitr(all_syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  ek <- tryCatch(enrichKEGG(gene = eg$ENTREZID, universe = uni$ENTREZID,
                             organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                             keyType = "ncbi-geneid"),
                 error = function(e) { cat("  KEGG error:", conditionMessage(e), "\n"); NULL })
  if (is.null(ek) || nrow(ek@result) == 0) {
    cat(sprintf("  %s: no KEGG pathways enriched\n", tag))
    return(NULL)
  }
  ek <- setReadable(ek, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  res <- ek@result
  res <- res[res$pvalue < 0.05, ]
  write.csv(res, file.path(OUT, sprintf("%s_KEGG_restored.csv", tag)), row.names = FALSE)
  cat(sprintf("  %s: %d KEGG pathways enriched\n", tag, nrow(res)))
  res
}

## T restored, NK restored, combined union
comb_restored <- union(tdeg$restored, ndeg$restored)
writeLines(comb_restored, file.path(OUT, "Combined_restored_genes.txt"))
t_kegg  <- run_kegg(tdeg$restored, all_t_syms, "Tcell")
nk_kegg <- run_kegg(ndeg$restored, all_nk_syms, "NK")
comb_kegg <- run_kegg(comb_restored, union(all_t_syms, all_nk_syms), "Combined")

## ---- summary ----
summary_lines <- c(
  sprintf("Tcell DEG PRRSV-vs-Control sig: %d", sum(tdeg$res$PvC$p_val_adj < 0.05, na.rm = TRUE)),
  sprintf("Tcell DEG APRR-vs-PRRSV  sig: %d", sum(tdeg$res$AvP$p_val_adj < 0.05, na.rm = TRUE)),
  sprintf("Tcell restored genes: %d", length(tdeg$restored)),
  sprintf("NK DEG PRRSV-vs-Control sig: %d", sum(ndeg$res$PvC$p_val_adj < 0.05, na.rm = TRUE)),
  sprintf("NK DEG APRR-vs-PRRSV  sig: %d", sum(ndeg$res$AvP$p_val_adj < 0.05, na.rm = TRUE)),
  sprintf("NK restored genes: %d", length(ndeg$restored)),
  sprintf("Combined restored genes: %d", length(comb_restored))
)
writeLines(summary_lines, file.path(OUT, "deg_summary.txt"))
cat(sprintf("[%s] DONE. Outputs in %s\n", Sys.time(), OUT))
cat(paste(summary_lines, collapse = "\n"), "\n")
