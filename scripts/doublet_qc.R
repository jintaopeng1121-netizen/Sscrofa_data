## ============================================================================
## doublet_qc.R — 实测双细胞率（独立于 pipeline，便于快速验证）
## 框架：仅 3 组 Control/PRRSV/APRR4()，不含高剂量。
## 流程：读 3 个原始 10x → merge → QC(subset) → DoubletFinder → 报率 + 存干净原始对象
## ============================================================================
suppressMessages({ library(Seurat); library(scDblFinder); library(SingleCellExperiment) })
BASE <- "data/scRNA-seq"
SRC  <- file.path(BASE, "raw_4samples")
OUT  <- "data/scRNA-seq/reanalysis_v2"
set.seed(42)

## 仅 3 组（无高剂量）
SAMPLES <- c("S14-control"="Control", "S5-PRRSV"="PRRSV", "S23-APPR4"="APRR")

read_counts_robust <- function(mx, bc, ft, gz) {
  feats <- read.delim(if (gz) gzfile(ft) else ft, header = FALSE, stringsAsFactors = FALSE)
  gene_names <- make.unique(as.character(feats[[min(2, ncol(feats))]]))
  cells <- readLines(if (gz) gzfile(bc) else bc)
  m <- Matrix::readMM(if (gz) gzfile(mx) else mx)
  nf <- length(gene_names); nc <- length(cells)
  if (nrow(m) == nf && ncol(m) == nc) {} 
  else if (nrow(m) == nc && ncol(m) == nf) { m <- Matrix::t(m) } 
  else stop(sprintf("dim mismatch: mtx %dx%d, features=%d, barcodes=%d", nrow(m), ncol(m), nf, nc))
  rownames(m) <- gene_names; colnames(m) <- cells
  as(m, "CsparseMatrix")
}

obj.list <- list()
for (s in names(SAMPLES)) {
  mat_dir <- file.path(SRC, s, "filtered_feature_bc_matrix")
  gz <- file.exists(file.path(mat_dir, "matrix.mtx.gz"))
  bc <- file.path(mat_dir, if (gz) "barcodes.tsv.gz" else "barcodes.tsv")
  ft <- file.path(mat_dir, if (gz) "features.tsv.gz" else "features.tsv")
  mx <- file.path(mat_dir, if (gz) "matrix.mtx.gz"   else "matrix.mtx")
  counts <- read_counts_robust(mx, bc, ft, gz)
  so <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
  so$sample <- s
  so$group  <- as.character(SAMPLES[s])
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-|^mt-")
  so <- RenameCells(so, new.names = paste0(s, "_", colnames(so)))
  obj.list[[s]] <- so
  cat(sprintf("  %s: %d cells x %d genes\n", s, ncol(so), nrow(so)))
}
obj <- Reduce(function(x, y) merge(x, y), obj.list)
cat(sprintf("  merged raw: %d cells x %d genes\n", ncol(obj), nrow(obj)))
obj <- subset(obj, subset = percent.mt < 15 & nFeature_RNA > 300 & nFeature_RNA < 6000)
cat(sprintf("  after QC: %d cells\n", ncol(obj)))

## scDblFinder（DoubletFinder 已从 CRAN 下架，scDblFinder 为更现代替代）
cat(">>> running scDblFinder ...\n")
# Seurat v5 合并对象后 counts 被拆成 per-sample 子层(counts.S14-control / counts.S5-PRRSV /
# counts.S23-APPR4)，[[ "RNA" ]]$counts 与 LayerData 都只返回首样本。先 JoinLayers 合并为单一层。
obj <- JoinLayers(obj, assay = "RNA")
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
cat(sprintf("  counts matrix: %d genes x %d cells\n", nrow(counts_mat), ncol(counts_mat)))
sce <- SingleCellExperiment(assays = list(counts = counts_mat))
colData(sce)$sample <- obj$sample
colData(sce)$group  <- obj$group
sce <- scDblFinder(sce, samples = "sample")       # 按 sample 分样本检测，避免批次混淆
obj$is_doublet <- sce$scDblFinder.class == "doublet"
n_dbl <- sum(obj$is_doublet)
cat(sprintf("  DOUBLETS: %d / %d (%.2f%%)\n", n_dbl, ncol(obj), 100 * n_dbl / ncol(obj)))
cat("  per-group doublet counts:\n")
print(table(obj$group, obj$is_doublet))
write.csv(obj@meta.data[, c("orig.ident", "group", "sample", "is_doublet")],
          file.path(OUT, "data/derived/doublet_calls.csv"), row.names = TRUE)
obj_clean <- obj[, !obj$is_doublet]
saveRDS(obj_clean, file.path(OUT, "data/derived/seurat_raw_doubletclean.rds"))
cat(sprintf("  saved doublet-clean raw object: %d cells -> seurat_raw_doubletclean.rds\n", ncol(obj_clean)))
cat(">>> doublet_qc DONE\n")
