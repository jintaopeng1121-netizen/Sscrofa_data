## ============================================================================
## _ssgsea_reproduce_and_fix.R
## 目的：查清 tables/ssGSEA_pathway_activity.csv（旧表）到底怎么算出来的，
##       并在确认后用正确口径重算。
##
## 线索：scripts/pipeline_r.R [9] 段是唯一的写入点：
##   rnk <- apply(sp, 2, rank, ties.method="average")   # 每细胞基因排名, 1=最高表达
##   N   <- nrow(rnk)
##   gsva_mat <- sapply(keep_sets, function(gs) {
##       idx <- which(rownames(sp) %in% gs)
##       colMeans(rnk[idx, , drop=FALSE]) / N            # 0..1, 越高表示该基因集越富集
##   })
##   sp 只含"分析基因宇宙"（markers + LR + pathway_sets + checkpoint），几百个基因
##
## 疑点：真正的 ssGSEA/GSVA 中，基因集内基因排名越靠前(1) → 富集分越高。
##       上式用 colMeans(rank)/N，基因集内基因排名越靠后 → 分数越高，**方向相反**。
##
## 输出：
##   tables/ssGSEA_pathway_activity_reproduced.csv  原算法复现（应与旧表逐值一致）
##   tables/ssGSEA_pathway_activity_gsva_curated.csv 真 GSVA ssGSEA + 同一套手写基因集
##   tables/ssGSEA_audit.txt                        三方对比
## ============================================================================
suppressMessages({ library(Seurat); library(GSVA); library(reshape2) })
set.seed(42)
OUT <- "data/scRNA-seq/reanalysis_v2"

## ---------- 与 pipeline_r.R 完全一致的定义 ----------
markers <- list(
  Epithelial   = c("EPCAM","KRT8","KRT18","KRT19","SFTPA1","SFTPC"),
  Neutrophil   = c("CXCL8","S100A8","S100A9","CXCR2"),
  Macrophage   = c("CD68","CD163","LYZ","MARCO","MSR1","FOLR2"),
  Monocyte     = c("LYZ","CD14","FCN1","VCAN"),
  DC           = c("FCER1A","CD1C","CLEC10A","IRF8"),
  Tcell_CD8    = c("CD3D","CD3E","CD8A","CD8B","GZMK","GZMB"),
  Tcell_CD4    = c("CD3D","CD3E","CD4","IL7R","FOXP3"),
  Bcell        = c("MS4A1","CD79A","CD79B","IGHM"),
  NK           = c("NCAM1","KLRF1","NKG7","GNLY"),
  Endothelial  = c("PECAM1","VWF","CLDN5","FLT1"),
  Fibroblast   = c("COL1A1","COL1A2","PDGFRB","ACTA2"),
  Ciliated     = c("FOXJ1","PIFO","DYDC2"))
LR <- list(
  PDCD1_PDCD1LG1 = c("PDCD1","CD274"), CTLA4_CD80 = c("CTLA4","CD80"),
  HAVCR2_LGALS9  = c("HAVCR2","LGALS9"), CXCL8_CXCR1 = c("CXCL8","CXCR1"),
  CCL2_CCR2      = c("CCL2","CCR2"),   CXCL10_CXCR3 = c("CXCL10","CXCR3"),
  IFNG_IFNGR1    = c("IFNG","IFNGR1"), IL1B_IL1R1 = c("IL1B","IL1R1"),
  TNF_TNFRSF1A   = c("TNF","TNFRSF1A"), TGFB1_TGFBR1 = c("TGFB1","TGFBR1"),
  VEGFA_KDR      = c("VEGFA","KDR"),  IL10_IL10RA = c("IL10","IL10RA"))
pathway_sets <- list(
  IFN_alpha_response = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","MX2","OAS1","OAS2","OASL","RSAD2","ISG20","IRF7","IRF9","STAT1","STAT2"),
  IFN_gamma_response = c("GZMB","CXCL9","CXCL10","CIITA","NLRC5","IRF1","STAT1","GBP1","GBP5"),
  IL1_NFKB_signaling = c("NFKB1","NFKBIA","RELA","IL1B","CXCL8","TNF","CCL2","PTGS2"),
  TNF_signaling = c("TNFAIP3","NFKBIA","CXCL2","CCL20","JUNB","NFKB2"),
  Apoptosis = c("CASP3","CASP8","BAX","BCL2","BIRC3","TNFRSF10B"),
  Hypoxia = c("HIF1A","VEGFA","LDHA","SLC2A1","PGK1"),
  EMT = c("CDH2","VIM","SNAI1","TWIST1","FN1","ACTA2"),
  T_cell_activation = c("CD3D","CD8A","GZMK","GZMB","IL2RG","LCK"),
  Neutrophil_degranulation = c("MPO","ELANE","CXCL8","S100A8","S100A9","DEFA1","LCN2"),
  Macrophage_alternative = c("CD163","MRC1","MSR1","FOLR2","CD209"),
  Complement = c("C1QA","C1QB","C1QC","C3","C4B"),
  Antiviral = c("OAS1","IFIH1","DDX58","MAVS","IRF3","IRF7"))
checkpt_genes <- c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","GZMB","TNFRSF9")

## ---------- 读当前 3 组对象 ----------
cat(">>> loading 3grp object ...\n")
obj <- readRDS(file.path(OUT, "data/derived/seurat_reanalysis_v2_3grp.rds"))
obj$group <- factor(obj$group, levels = c("Control","PRRSV","APRR"))

sct_genes <- rownames(GetAssayData(obj, assay = "SCT", layer = "data"))
need_genes <- unique(c(unlist(markers), unlist(LR), unlist(pathway_sets), checkpt_genes))
need_genes <- need_genes[need_genes %in% sct_genes]
sp <- as.matrix(GetAssayData(obj, assay = "SCT", layer = "data")[need_genes, , drop = FALSE])
cat(sprintf("    分析基因宇宙: %d genes x %d cells\n", nrow(sp), ncol(sp)))

## ---------- A. 原算法复现 ----------
keep_sets <- lapply(pathway_sets, function(g) g[g %in% rownames(sp)])
keep_sets <- keep_sets[sapply(keep_sets, length) >= 3]
cat("    通路基因在宇宙内命中:", sapply(keep_sets, length), "\n")
t0 <- Sys.time()
rnk <- apply(sp, 2, rank, ties.method = "average")
N   <- nrow(rnk)
gsva_mat <- sapply(keep_sets, function(gs) {
  idx <- which(rownames(sp) %in% gs)
  colMeans(rnk[idx, , drop = FALSE]) / N
})
cat(sprintf("    原算法完成 %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
df <- as.data.frame(gsva_mat)
df$cell_type <- obj$cell_type
df$group     <- obj$group
pa_rep <- aggregate(. ~ cell_type + group, data = df, mean)
pa_rep$cell_type <- as.character(pa_rep$cell_type)
pa_rep$group     <- as.character(pa_rep$group)
pa_rep <- pa_rep[order(pa_rep$cell_type, pa_rep$group), ]
write.csv(pa_rep, file.path(OUT, "tables/ssGSEA_pathway_activity_reproduced.csv"), row.names = FALSE)

## ---------- B. 真 GSVA ssGSEA + 同一套手写基因集（全基因空间） ----------
cat(">>> GSVA ssGSEA with the same curated gene sets ...\n")
mat_full <- GetAssayData(obj, assay = "SCT", layer = "data")
gs_use <- lapply(pathway_sets, function(g) intersect(g, rownames(mat_full)))
gs_use <- gs_use[lengths(gs_use) >= 3]
cat("    全基因空间命中:\n"); for (nm in names(gs_use)) cat(sprintf("      %-26s %2d\n", nm, length(gs_use[[nm]])))
cat("    被丢弃:", paste(setdiff(names(pathway_sets), names(gs_use)), collapse=", "), "\n")
blocks <- split(seq_len(ncol(mat_full)), ceiling(seq_len(ncol(mat_full)) / 6000))
res <- NULL
for (b in blocks) {
  p <- ssgseaParam(mat_full[, b, drop = FALSE], gs_use, minSize = 3, maxSize = 1000, normalize = TRUE)
  rk <- gsva(p, verbose = FALSE)
  res <- if (is.null(res)) rk else cbind(res, rk)
}
res <- res[, colnames(mat_full), drop = FALSE]
d2 <- as.data.frame(t(res))
d2$cell_type <- obj$cell_type[match(rownames(d2), colnames(obj))]
d2$group     <- obj$group[match(rownames(d2), colnames(obj))]
pa_new <- aggregate(. ~ cell_type + group, data = d2, mean)
pa_new$cell_type <- as.character(pa_new$cell_type)
pa_new$group     <- as.character(pa_new$group)
pa_new <- pa_new[order(pa_new$cell_type, pa_new$group), ]
write.csv(pa_new, file.path(OUT, "tables/ssGSEA_pathway_activity_gsva_curated.csv"), row.names = FALSE)

## ---------- C. 三方对比 ----------
old <- read.csv(file.path(OUT, "tables/ssGSEA_pathway_activity.csv"), check.names = FALSE)
old$cell_type <- as.character(old$cell_type); old$group <- as.character(old$group)
key <- function(d) paste(d$cell_type, d$group, sep = "|")
m <- merge(old, pa_rep, by = c("cell_type","group"), suffixes = c(".old",".rep"))
m2 <- merge(m, pa_new, by = c("cell_type","group"))
qc <- c()
qc <- c(qc, sprintf("rows: old=%d reproduced=%d gsva=%d ; merged=%d",
                     nrow(old), nrow(pa_rep), nrow(pa_new), nrow(m2)))
qc <- c(qc, sprintf("cells=%d  gene_universe=%d", ncol(sp), nrow(sp)))
qc <- c(qc, "")
qc <- c(qc, sprintf("%-24s %10s %10s %10s", "pathway", "repro_vs_old", "gsva_vs_old", "gsva_vs_repro"))
for (nm in names(pathway_sets)) {
  a <- m2[[paste0(nm, ".old")]]; b <- m2[[paste0(nm, ".rep")]]; c2 <- m2[[paste0(nm, ".new")]]
  if (is.null(a) || is.null(b) || is.null(c2)) next
  r1 <- suppressWarnings(cor(a, b, method = "spearman"))
  r2 <- suppressWarnings(cor(a, c2, method = "spearman"))
  r3 <- suppressWarnings(cor(b, c2, method = "spearman"))
  qc <- c(qc, sprintf("%-24s %10.3f %10.3f %10.3f", nm, r1, r2, r3))
}
## 复现的最大绝对差（验证是否逐值一致）
maxd <- 0
for (nm in names(pathway_sets)) {
  a <- m[[paste0(nm,".old")]]; b <- m[[paste0(nm,".rep")]]
  if (!is.null(a) && !is.null(b)) maxd <- max(maxd, max(abs(a-b)))
}
qc <- c(qc, "", sprintf("复现表 vs 旧表 最大绝对差 = %.8f", maxd))
qc <- c(qc, sprintf("旧表数值范围 %.3f ~ %.3f ; 复现表 %.3f ~ %.3f ; GSVA表 %.3f ~ %.3f",
                     min(as.matrix(old[, -c(1,2)])), max(as.matrix(old[, -c(1,2)])),
                     min(as.matrix(pa_rep[, -c(1,2)])), max(as.matrix(pa_rep[, -c(1,2)])),
                     min(as.matrix(pa_new[, -c(1,2)])), max(as.matrix(pa_new[, -c(1,2)]))))
writeLines(qc, file.path(OUT, "tables/ssGSEA_audit.txt"))
cat("\n"); cat(paste(qc, collapse = "\n"), "\n")
cat(">>> audit DONE\n")
