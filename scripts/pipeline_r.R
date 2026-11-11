# ---------------------------------------------------------------------------
# NOTE (repository release): the pathway-activity block in this script is the
# LEGACY rank-mean calculation. It is retained for provenance only and was NOT
# used for the manuscript. The published ssGSEA scores come from
# _ssgsea_reproduce_and_fix.R (GSVA 2.0.7 + 12 curated gene sets). See
# docs/METHODS_NOTES.md before re-running any analysis from this file.
# ---------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(reshape2)
})
options(future.globals.maxSize = 10 * 1024^3)
set.seed(42)

BASE   <- "data/scRNA-seq"
SRC    <- file.path(BASE, "raw_4samples")
OUT    <- "data/scRNA-seq/reanalysis_v2"
FIG    <- file.path(OUT, "figures")
TAB    <- file.path(OUT, "tables")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)

SAMPLES <- c("S14-control"="Control", "S5-PRRSV"="PRRSV",
             "S23-APPR4"="APRR")
sample_ids <- names(SAMPLES)
group_lab  <- SAMPLES

# ---------- gene reference tables (defined once) ----------
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
  Ciliated     = c("FOXJ1","PIFO","DYDC2")
)
LR <- list(
  PDCD1_PDCD1LG1 = c("PDCD1","CD274"),
  CTLA4_CD80     = c("CTLA4","CD80"),
  HAVCR2_LGALS9  = c("HAVCR2","LGALS9"),
  CXCL8_CXCR1    = c("CXCL8","CXCR1"),
  CCL2_CCR2      = c("CCL2","CCR2"),
  CXCL10_CXCR3   = c("CXCL10","CXCR3"),
  IFNG_IFNGR1    = c("IFNG","IFNGR1"),
  IL1B_IL1R1     = c("IL1B","IL1R1"),
  TNF_TNFRSF1A   = c("TNF","TNFRSF1A"),
  TGFB1_TGFBR1   = c("TGFB1","TGFBR1"),
  VEGFA_KDR      = c("VEGFA","KDR"),
  IL10_IL10RA    = c("IL10","IL10RA")
)
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
  Antiviral = c("OAS1","IFIH1","DDX58","MAVS","IRF3","IRF7")
)
checkpt_genes <- c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","GZMB","TNFRSF9")

# ---------- [1] load ----------
cat(">>> [1] load 4 samples (Read10X)\n")
# 稳健读取：部分样本(如 S23-APPR4)的 matrix.mtx 被转置存成 细胞×基因，
# 这里按 features/barcodes 行数自动判方向，反了就 t()，避免 ReadMtx 维度报错
read_counts_robust <- function(mx, bc, ft, gz) {
  feats <- read.delim(if (gz) gzfile(ft) else ft, header = FALSE, stringsAsFactors = FALSE)
  gene_names <- make.unique(as.character(feats[[min(2, ncol(feats))]]))
  cells <- readLines(if (gz) gzfile(bc) else bc)
  m <- Matrix::readMM(if (gz) gzfile(mx) else mx)
  nf <- length(gene_names); nc <- length(cells)
  if (nrow(m) == nf && ncol(m) == nc) {
    # 标准 基因×细胞
  } else if (nrow(m) == nc && ncol(m) == nf) {
    m <- Matrix::t(m)   # 转置样本：细胞×基因 -> 基因×细胞
    cat("    [transposed matrix detected -> t()]\n")
  } else {
    stop(sprintf("dim mismatch: mtx %dx%d, features=%d, barcodes=%d", nrow(m), ncol(m), nf, nc))
  }
  rownames(m) <- gene_names
  colnames(m) <- cells
  as(m, "CsparseMatrix")
}
obj.list <- list()
for (s in sample_ids) {
  mat_dir <- file.path(SRC, s, "filtered_feature_bc_matrix")
  # 兼容未压缩(.tsv/.mtx)与压缩(.gz)两种 10x 格式
  gz <- file.exists(file.path(mat_dir, "matrix.mtx.gz"))
  bc <- file.path(mat_dir, if (gz) "barcodes.tsv.gz" else "barcodes.tsv")
  ft <- file.path(mat_dir, if (gz) "features.tsv.gz" else "features.tsv")
  mx <- file.path(mat_dir, if (gz) "matrix.mtx.gz"   else "matrix.mtx")
  counts <- read_counts_robust(mx, bc, ft, gz)
  so <- CreateSeuratObject(counts = counts, project = s,
                           min.cells = 3, min.features = 200)
  so$sample <- s
  # 必须 as.character 去掉 names：group_lab[s] 是具名向量，Seurat 会把 name 当条码匹配 → "No cell overlap"
  so$group  <- as.character(group_lab[s])
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-|^mt-")
  # Seurat v5 兼容：先给每个对象手动加样本前缀，避免合并时 barcode 撞车
  so <- RenameCells(so, new.names = paste0(s, "_", colnames(so)))
  obj.list[[s]] <- so
  cat(sprintf("  %s: %d cells x %d genes\n", s, ncol(so), nrow(so)))
}

# ---------- [2] merge raw -> SCTransform (带缓存，加速重跑) ----------
# 先合并原始 RNA 对象（不含 SCT assay，merge 稳定），再一次性 SCTransform，
# 规避 Seurat v5 直接 merge SCT 对象报 "No cell overlap" 的坑。
# SCTransform 最耗时，产物缓存为 seurat_sct.rds，重跑时直接读缓存跳过 ~25min SCT。
sct_rds <- file.path(OUT, "data/derived/seurat_sct.rds")
if (file.exists(sct_rds)) {
  cat(">>> [2] load cached SCT+QC object\n")
  obj <- readRDS(sct_rds)
} else {
  cat(">>> [2] merge raw -> SCTransform(glmGamPoi)\n")
  obj <- Reduce(function(x, y) merge(x, y), obj.list)
  cat(sprintf("  merged: %d cells x %d genes\n", ncol(obj), nrow(obj)))
  obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  obj$group <- factor(obj$group, levels = c("Control","PRRSV","APRR"))
  obj <- subset(obj, subset = percent.mt < 15 & nFeature_RNA > 300 & nFeature_RNA < 6000)
  # ---------- [doublet] scDblFinder 双细胞剔除（SCT 之前，作用于原始 RNA 对象）----------
  cat(">>> [doublet] scDblFinder\n")
  obj <- JoinLayers(obj, assay = "RNA")
  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment(assays = list(counts = counts_mat))
  colData(sce)$sample <- obj$sample
  sce <- scDblFinder(sce, samples = "sample")
  obj <- obj[, sce$scDblFinder.class != "doublet"]
  cat(sprintf("  doublets removed: %d, remaining %d cells\n",
              sum(sce$scDblFinder.class == "doublet"), ncol(obj)))
  saveRDS(obj, file.path(OUT, "data/derived/seurat_raw_doubletclean.rds"))   # 备份：双细胞剔除后原始对象
  saveRDS(obj, sct_rds)
}
cat(sprintf("  after QC: %d cells x %d genes\n", ncol(obj), nrow(obj)))

# ---------- [3] integration ----------
cat(">>> [3] PCA(SCT) -> Harmony\n")
obj <- RunPCA(obj, npcs = 30, assay = "SCT", verbose = FALSE)
# harmony 1.2.4 的 RunHarmony.Seurat 用 reduction.use / reduction.save（老参数名），
# 写 Seurat 风格 reduction/reduction.name 会被部分匹配到多个形参 → "argument N matches multiple formal arguments"
obj <- RunHarmony(obj, group.by.vars = "sample",
                  reduction.use = "pca", dims.use = 1:30,
                  reduction.save = "harmony", assay.use = "SCT", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
# algorithm=1 (Louvain, 纯 R/igraph) 替代 algorithm=4 (leiden, 需 Python/numpy；
# 本机 pypi 不可达，reticulate 装不上 numpy，leiden 会崩)
obj <- FindClusters(obj, resolution = 0.6, algorithm = 1, verbose = FALSE)

# ---------- [4] annotation ----------
cat(">>> [4] annotation via canonical markers\n")
sct_genes <- rownames(GetAssayData(obj, assay = "SCT", layer = "data"))
# 仅取"分析基因宇宙"(marker + LR + pathway + checkpoint，几百个基因)成密集矩阵，
# 避免对全转录本(2万基因)建巨大密集矩阵撑爆内存；SCT data 本身是稀疏矩阵，取子集很省。
need_genes <- unique(c(unlist(markers), unlist(LR), unlist(pathway_sets), checkpt_genes))
need_genes <- need_genes[need_genes %in% sct_genes]
sp <- as.matrix(GetAssayData(obj, assay = "SCT", layer = "data")[need_genes, , drop = FALSE])
EX <- sp
# Seurat v5 的 AddModuleScore 输出列名带尾随数字、不可控，改为手动均值打分
for (k in names(markers)) {
  gs <- markers[[k]][markers[[k]] %in% need_genes]
  obj[[paste0("mod_", k)]] <- colMeans(sp[gs, , drop = FALSE], na.rm = TRUE)
}
mod_cols <- paste0("mod_", names(markers))
mmat <- as.matrix(obj@meta.data[, mod_cols, drop = FALSE]); colnames(mmat) <- names(markers)
obj$cell_type <- colnames(mmat)[max.col(mmat)]
cd8_expr <- sp["CD8A", , drop = FALSE]
cd4_expr <- sp["CD4", , drop = FALSE]
is_t <- obj$cell_type %in% c("Tcell_CD8","Tcell_CD4")
obj$cell_type[is_t] <- ifelse(cd8_expr[1, is_t, drop = TRUE] > cd4_expr[1, is_t, drop = TRUE], "Tcell_CD8", "Tcell_CD4")
cat("  cell types:\n"); print(table(obj$cell_type))
cat(sprintf("  analysis universe: %d genes x %d cells (dense)\n", nrow(sp), ncol(sp)))

# ---------- [5] DEG ----------
cat(">>> [5] DEG (4 comparisons)\n")
deg_list <- list()
comp <- list(
  PRRSV_vs_Control   = c("PRRSV","Control"),
  APRR4_vs_PRRSV     = c("APRR","PRRSV")
)
# 断点续跑：若 4 组 DEG csv 已存在则直接读回，跳过耗时的 PrepSCTFindMarkers + FindMarkers
all_cached <- all(file.exists(file.path(TAB, paste0("DEG_", names(comp), ".csv"))))
if (!all_cached) {
  # SCT 多样本合并后每样本各有独立模型，FindMarkers(SCT) 前必须先 PrepSCTFindMarkers（只需一次）
  obj <- PrepSCTFindMarkers(obj)
}
Idents(obj) <- obj$group
for (nm in names(comp)) {
  fp <- file.path(TAB, paste0("DEG_", nm, ".csv"))
  if (file.exists(fp)) {
    ct <- read.csv(fp, check.names = FALSE); rownames(ct) <- ct$gene
    deg_list[[nm]] <- ct
    cat(sprintf("  %s: cached (%d rows)\n", nm, nrow(ct)))
    next
  }
  ct <- FindMarkers(obj, ident.1 = comp[[nm]][1], ident.2 = comp[[nm]][2],
                    assay = "SCT", test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1,
                    recorrect_umi = FALSE)
  ct$gene <- rownames(ct); ct$comparison <- nm
  deg_list[[nm]] <- ct
  write.csv(ct, fp, row.names=FALSE)
  cat(sprintf("  %s: %d DEGs total; padj<0.05: %d\n",
              nm, nrow(ct), sum(ct$p_val_adj<0.05, na.rm=TRUE)))
}

# ---------- [6] proportions ----------
cat(">>> [6] cell proportions\n")
prop <- prop.table(table(obj$cell_type, obj$group), margin=2)
prop_df <- as.data.frame(as.table(prop))
names(prop_df) <- c("cell_type", "group", "Freq")   # as.table 默认列名 Var1/Var2/Freq, 显式改名供 aes 引用
write.csv(prop_df, file.path(TAB,"proportions.csv"), row.names=FALSE)
p_prop <- ggplot(prop_df,
                 aes(x=group, y=Freq*100, fill=cell_type)) +
  geom_bar(stat="identity", position="stack") +
  labs(x="", y="% of cells", fill="Cell type", title="Cell-type composition by group") +
  theme(axis.text.x=element_text(angle=30, hjust=1)) + scale_fill_brewer(palette="Set3")

# ---------- [7] myeloid pseudotime (PCA-based; 免去 slingshot/BioC 依赖) ----------
cat(">>> [7] myeloid pseudotime (PCA-based)\n")
mye <- subset(obj, cell_type %in% c("Macrophage","Monocyte","DC"))
mye <- RunPCA(mye, assay = "SCT", npcs = 10, verbose = FALSE)
pt <- as.numeric(Embeddings(mye, "pca")[, 1])
# 方向校正：让 Monocyte 起始端 pseudotime 最小
mon_mean <- mean(pt[mye$cell_type == "Monocyte"], na.rm = TRUE)
mac_mean <- mean(pt[mye$cell_type == "Macrophage"], na.rm = TRUE)
if (is.finite(mon_mean) && is.finite(mac_mean) && mon_mean > mac_mean) pt <- -pt
mye$pseudotime <- pt
um <- Embeddings(mye, "umap")
p_traj <- ggplot(data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], PT = pt),
                 aes(UMAP1, UMAP2, color = PT)) +
  geom_point(size = 0.5, alpha = 0.6) + scale_color_viridis_c() +
  labs(title = "Myeloid pseudotime (PCA-based trajectory)", color = "pseudotime") + theme_void()

# ---------- [8] predicted LR ----------
cat(">>> [8] predicted ligand-receptor interactions\n")
Idents(obj) <- obj$cell_type
ctypes <- unique(obj$cell_type)
lr_score <- matrix(NA, nrow=length(LR), ncol=length(ctypes), dimnames=list(names(LR), ctypes))
for (lr in names(LR)) {
  L <- LR[[lr]][1]; R <- LR[[lr]][2]
  for (ct in ctypes) {
    cells <- WhichCells(obj, idents=ct)
    if (length(cells)==0) next
    if (L %in% need_genes && R %in% need_genes)
      lr_score[lr, ct] <- mean(EX[L, cells]) * mean(EX[R, cells])
  }
}
lr_score <- lr_score / max(lr_score, na.rm=TRUE)
lr_df <- as.data.frame(as.table(lr_score)); names(lr_df) <- c("LR_pair","CellType","score")
write.csv(lr_df, file.path(TAB,"LR_predicted.csv"), row.names=FALSE)
p_lr <- ggplot(lr_df, aes(x=CellType, y=LR_pair, fill=score)) +
  geom_tile() + scale_fill_gradient(low="white", high="firebrick") +
  labs(title="Predicted ligand-receptor activity (SCT expr product, normalized)",
       x="Cell type", y="LR pair") +
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size=8))

# ---------- [9] ssGSEA pathway (rank-based, pure R; 免去 GSVA/BioC 依赖) ----------
cat(">>> [9] ssGSEA pathway activity (rank-based)\n")
keep_sets <- lapply(pathway_sets, function(g) g[g %in% rownames(sp)])
keep_sets <- keep_sets[sapply(keep_sets, length) >= 3]
cat("  ranking genes per cell (", nrow(sp), "genes x ", ncol(sp), "cells)...\n", sep="")
rnk <- apply(sp, 2, rank, ties.method = "average")   # 每细胞基因排名, 1=最高表达
N <- nrow(rnk)
gsva_mat <- sapply(keep_sets, function(gs) {
  idx <- which(rownames(sp) %in% gs)
  colMeans(rnk[idx, , drop = FALSE]) / N               # 0..1, 越高表示该基因集越富集
})
gsva_df <- as.data.frame(gsva_mat)                    # cells x genesets (sapply 已按列绑定, 无需转置)
gsva_df$cell_type <- obj$cell_type
gsva_df$group <- obj$group
pa <- aggregate(. ~ cell_type + group, data = gsva_df, mean)
write.csv(pa, file.path(TAB, "ssGSEA_pathway_activity.csv"), row.names = FALSE)
pa_long <- melt(pa, id.vars = c("cell_type", "group"))
p_path <- ggplot(pa_long, aes(x = group, y = value, fill = group)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free") +
  labs(title = "ssGSEA pathway activity (rank-based)", y = "enrichment score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------- [10] Phase 2 feature axes ----------
cat(">>> [10] Phase 2: wet-result guided axes\n")
# Neutrophil states
neu <- subset(obj, cell_type=="Neutrophil")
neu_d <- GetAssayData(neu, layer="data", assay="SCT")
neu$MPO_s   <- if ("MPO"   %in% rownames(neu_d)) as.numeric(neu_d["MPO", ])   else 0
neu$CXCL8_s <- if ("CXCL8" %in% rownames(neu_d)) as.numeric(neu_d["CXCL8", ]) else 0
neu$ISG15_s <- if ("ISG15" %in% rownames(neu_d)) as.numeric(neu_d["ISG15", ]) else 0
p_neu <- wrap_plots(
  DimPlot(neu, reduction="umap", group.by="group") + ggtitle("Neutrophils by group"),
  FeaturePlot(neu, reduction="umap", features=c("MPO","CXCL8","ISG15")), ncol=2) + plot_layout(ncol=1)
write.csv(aggregate(cbind(MPO_s,CXCL8_s,ISG15_s)~group, data=neu@meta.data, mean),
          file.path(TAB,"neutrophil_state_scores.csv"), row.names=FALSE)

# Myeloid landscape
mac <- subset(obj, cell_type=="Macrophage")
mac_d <- GetAssayData(mac, layer="data", assay="SCT")
mac$CD163_s <- if ("CD163" %in% rownames(mac_d)) as.numeric(mac_d["CD163", ]) else 0
mac$CD68_s  <- if ("CD68"  %in% rownames(mac_d)) as.numeric(mac_d["CD68", ])  else 0
p_mac <- wrap_plots(
  DimPlot(mac, reduction="umap", group.by="group") + ggtitle("Macrophages by group"),
  FeaturePlot(mac, reduction="umap", features=c("CD163","CD68")), ncol=2) + plot_layout(ncol=1)
write.csv(aggregate(cbind(CD163_s,CD68_s)~group, data=mac@meta.data, mean),
          file.path(TAB,"macrophage_scores.csv"), row.names=FALSE)

# T-cell checkpoint
chk <- checkpt_genes[checkpt_genes %in% rownames(sp)]
tc <- subset(obj, cell_type %in% c("Tcell_CD8","Tcell_CD4"))
p_tcheck <- FeaturePlot(tc, reduction="umap", features=chk, ncol=4) +
  ggtitle("T-cell checkpoint molecules")
tc_m <- as.data.frame(t(GetAssayData(tc, layer="data", assay="SCT")[chk, , drop=FALSE]))
tc_m$group <- tc$group
write.csv(aggregate(. ~ group, data=tc_m, mean),
          file.path(TAB,"Tcell_checkpoint_scores.csv"), row.names=FALSE)

# ISG / antiviral program
isg_genes <- keep_sets[["Antiviral"]]
obj$ISG_score <- colMeans(sp[isg_genes, , drop = FALSE])
isg_stat <- aggregate(ISG_score ~ group, data=obj@meta.data, mean)
write.csv(isg_stat, file.path(TAB,"ISG_response_scores.csv"), row.names=FALSE)
p_isg <- ggplot(isg_stat, aes(x=group, y=ISG_score, fill=group)) + geom_bar(stat="identity") +
  labs(title="Mean ISG/antiviral program score by group", y="ISG score") +
  theme(axis.text.x=element_text(angle=30, hjust=1))

# ---------- [11] figures ----------
cat(">>> [11] render figures\n")
render <- function(name, plot, w=10, h=7) {
  pdf(file.path(FIG, paste0(name,".pdf")), width=w, height=h); print(plot); dev.off()
  png(file.path(FIG, paste0(name,".png")), width=w*120, height=h*120); print(plot); dev.off()
}
render("Fig1_overview",
       wrap_plots(DimPlot(obj, reduction="umap", group.by="sample") + ggtitle("by sample"),
                   DimPlot(obj, reduction="umap", group.by="group") + ggtitle("by group"),
                   FeaturePlot(obj, reduction="umap", features=c("percent.mt","nFeature_RNA")), ncol=2), 11, 8)
render("Fig2_clusters_annotation",
       wrap_plots(DimPlot(obj, reduction="umap", label=TRUE) + ggtitle("Louvain clusters"),
                   DimPlot(obj, reduction="umap", group.by="cell_type") + ggtitle("annotated cell type"), ncol=2), 12, 6)
top_genes <- unique(na.omit(unlist(lapply(deg_list, function(d) head(rownames(d[order(d$p_val_adj), ]), 15)))))
render("Fig3_DEG_heatmap",
       DoHeatmap(obj, features=top_genes, group.by="group", size=3) + NoLegend(), 10, 14)
render("Fig4_proportions", p_prop, 9, 6)
render("Fig5_myeloid_trajectory", p_traj, 8, 7)
render("Fig6_LR_interactions", p_lr, 10, 8)
render("Fig7_pathway_activity", p_path, 12, 10)
render("Fig8_neutrophil_states", p_neu, 11, 9)
render("Fig9_myeloid_landscape", p_mac, 11, 9)
render("Fig10_Tcell_checkpoint", p_tcheck, 14, 8)
render("Fig11_ISG_response", p_isg, 8, 6)

# ---------- summary ----------
summary_lines <- c(
  paste("Total cells:", ncol(obj)),
  paste("Cell types:", paste(table(obj$cell_type), collapse=", ")),
  "Groups: Control(S14) / PRRSV(S5) / APRR(S23)",
  "Phase1: QC+SCTransform+Harmony+Louvain(cluster)+marker-annot+DEG(3)+proportions+myeloid-pseudotime(PCA)+LR+ssGSEA(rank-based)",
  "Phase2: neutrophil(MPO/CXCL8/ISG15), macrophage(CD163/CD68), T-checkpoint(PDCD1/CD274..), ISG"
)
writeLines(summary_lines, file.path(OUT, "notes/summary.txt"))
saveRDS(obj, file.path(OUT, "data/derived/seurat_reanalysis_v2.rds"))
cat(">>> DONE.\n", paste(summary_lines, collapse="\n"), "\n")
