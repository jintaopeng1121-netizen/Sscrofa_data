## ============================================================================
## ssgsea_pathway_activity.R
## 重新生成 tables/ssGSEA_pathway_activity.csv（原生成脚本已遗失，本脚本为重建版）
##
## 输入：data/derived/seurat_reanalysis_v2_3grp.rds（SCT assay, data layer = log-normalized）
## 方法：GSVA::gsva(method = "ssgsea", ssgsea.norm = TRUE)，对每个细胞单独算 ES 并归一化到 [0,1]，
##       再按 cell_type × group 取均值（与旧表结构一致：12 cell types × 3 groups × 12 pathways）
##
## 基因集来源（全部在脚本内显式声明，可追溯）：
##   - MSigDB Hallmark（msigdbr collection = "H"）：IFN-α / IFN-γ / TNF-NF-κB / Apoptosis /
##     Hypoxia / EMT / Complement
##   - MSigDB C5 GO:BP：IL-1 mediated signaling、T cell activation、neutrophil degranulation
##   - 自定义文献集：Macrophage_alternative（M2 标志物）、Antiviral（ISG 效应集）
##   —— 自定义集基因列表直接写在本文件里，来源见注释
##
## 输出：
##   tables/ssGSEA_pathway_activity_rebuilt.csv
##   tables/ssGSEA_rebuild_qc.txt
## 用法：Rscript ssgsea_pathway_activity.R [n_cells_for_test]
## ============================================================================
suppressMessages({
  library(Seurat); library(GSVA); library(msigdbr); library(Matrix)
})
set.seed(42)

OUT <- "data/scRNA-seq/reanalysis_v2"
args <- commandArgs(trailingOnly = TRUE)
N_TEST <- if (length(args) >= 1) as.integer(args[1]) else NA_integer_

## ---------- 1. 读对象 ----------
cat(">>> loading Seurat object ...\n")
obj <- readRDS(file.path(OUT, "data/derived/seurat_reanalysis_v2_3grp.rds"))
obj$group     <- factor(obj$group, levels = c("Control", "PRRSV", "APRR"))
obj$cell_type <- factor(obj$cell_type)

mat <- GetAssayData(obj, assay = "SCT", layer = "data")
cat(sprintf("    matrix: %d genes x %d cells | class %s\n",
            nrow(mat), ncol(mat), class(mat)[1]))

if (!is.na(N_TEST) && N_TEST < ncol(mat)) {
  idx <- sample(seq_len(ncol(mat)), N_TEST)
  mat <- mat[, idx]
  cat(sprintf("    [TEST MODE] subsampled to %d cells\n", ncol(mat)))
  obj_test <- subset(obj, cells = colnames(mat))
} else {
  obj_test <- obj
}

## ---------- 2. 基因集 ----------
cat(">>> building gene sets ...\n")
msig <- function(...) {
  tryCatch(
    msigdbr(species = "Homo sapiens", ...),
    error = function(e) {
      ## 旧版 msigdbr 用 category / subcategory
      a <- list(...)
      if (!is.null(a$collection)) a$category <- a$collection
      if (!is.null(a$subcollection)) a$subcategory <- a$subcollection
      a$collection <- NULL; a$subcollection <- NULL
      do.call(msigdbr, c(list(species = "Homo sapiens"), a))
    })
}
take <- function(df, name) unique(df$gene_symbol[df$gs_name == name])

H  <- msig(collection = "H")
GO <- msig(collection = "C5", subcollection = "GO:BP")
RE <- msig(collection = "C2", subcollection = "CP:REACTOME")

## 自定义 M2 巨噬细胞集：MRC1/CD163/MSR1/MARCO/STAB1/LYVE1/F13A1/CCL18/IL10/TGFB1/
##   CCL17/CCL22/CD209/CLEC7A/MERTK/AXL/TGM2/SELENOP/SLC40A1/HMOX1/IL4R/PPARG/NR1H3
##   （Mantovani et al. 2004 Trends Immunol; Murray et al. 2014 Immunity M2 共识标志物）
M2 <- c("MRC1","CD163","MSR1","MARCO","STAB1","LYVE1","F13A1","CCL18","IL10","TGFB1",
        "CCL17","CCL22","CD209","CLEC7A","MERTK","AXL","TGM2","SELENOP","SLC40A1",
        "HMOX1","IL4R","PPARG","NR1H3","ARG1","RETNLA","CLEC10A","CD200R1")

## 自定义抗病毒 ISG 效应集（Schneider et al. 2014 Annu Rev Immunol; Schoggins 2019）
ISG <- c("ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1","IFIT2","IFIT3","IFIT5",
         "DDX58","IFIH1","STAT1","STAT2","IRF7","IRF9","RSAD2","BST2","IFI6","IFI44",
         "IFI44L","XAF1","EIF2AK2","ADAR","SAMD9","HERC5","USP18","ISG20","TRIM25",
         "PARP9","PARP12","PARP14","DDX60","IFITM3","IFITM1","LY6E","SP100","TREX1")

gs <- list(
  IFN_alpha_response       = take(H,  "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  IFN_gamma_response       = take(H,  "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  IL1_NFKB_signaling       = unique(c(
                               take(RE, "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING"),
                               "IL1R1","IL1B","IL1RN","NFKB1","NFKB2","RELA","RELB","NFKBIA",
                               "NFKBIZ","TNFAIP3","MAP3K7","IRAK1","IRAK4","MYD88","TRAF6")),
  TNF_signaling            = take(H,  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  Apoptosis                = take(H,  "HALLMARK_APOPTOSIS"),
  Hypoxia                  = take(H,  "HALLMARK_HYPOXIA"),
  EMT                      = take(H,  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  T_cell_activation        = take(GO, "GOBP_T_CELL_ACTIVATION"),
  Neutrophil_degranulation = take(RE, "REACTOME_NEUTROPHIL_DEGRANULATION"),
  Macrophage_alternative   = M2,
  Complement               = take(H,  "HALLMARK_COMPLEMENT"),
  Antiviral                = ISG
)

## 报告基因集规模与在猪注释里的匹配率
allg <- rownames(mat)
gs_use <- list()
cat("\n    pathway                        n_set  matched(seen in object)\n")
for (nm in names(gs)) {
  g0 <- unique(gs[[nm]])
  hit <- intersect(g0, allg)
  gs_use[[nm]] <- hit
  cat(sprintf("    %-30s %5d  %5d\n", nm, length(g0), length(hit)))
}
drop <- names(gs_use)[lengths(gs_use) < 10]
if (length(drop)) {
  cat("    !! 丢弃基因数 <10 的通路:", paste(drop, collapse = ", "), "\n")
  gs_use <- gs_use[!names(gs_use) %in% drop]
}

## ---------- 3. ssGSEA ----------
cat("\n>>> running GSVA ssGSEA ...\n")
t0 <- Sys.time()
## GSVA >= 1.52 / 2.x：必须用参数对象，旧 gsva(..., method=) 已 defunct
## 全量 as.matrix 约 3.4 GB，故按细胞分块（ssGSEA 逐细胞独立，块间无依赖）
CHUNK <- if (length(args) >= 2) as.integer(args[2]) else 6000L
run_ssgsea_block <- function(m) {
  p <- ssgseaParam(m, gs_use, minSize = 10, maxSize = 1000, normalize = TRUE)
  gsva(p, verbose = FALSE)
}
blocks <- split(seq_len(ncol(mat)), ceiling(seq_len(ncol(mat)) / CHUNK))
res <- NULL
for (k in seq_along(blocks)) {
  sub <- mat[, blocks[[k]], drop = FALSE]
  rk <- run_ssgsea_block(sub)
  res <- if (is.null(res)) rk else cbind(res, rk)
  cat(sprintf("    block %d/%d: %d cells done (elapsed %.1f min)\n",
              k, length(blocks), ncol(sub),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- res[, colnames(mat), drop = FALSE]
cat(sprintf("    done in %.1f min | output %d x %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res), ncol(res)))

## ---------- 4. 汇总 cell_type × group ----------
df <- as.data.frame(t(res))
df$cell_type <- obj_test$cell_type[match(rownames(df), colnames(obj_test))]
df$group     <- obj_test$group[match(rownames(df), colnames(obj_test))]
df <- df[!is.na(df$cell_type) & !is.na(df$group), ]

agg <- aggregate(df[, names(gs_use)], by = list(cell_type = df$cell_type, group = df$group),
                 FUN = mean)
agg <- agg[order(agg$cell_type, factor(agg$group, levels = c("Control","PRRSV","APRR"))), ]

write.csv(agg, file.path(OUT, "tables/ssGSEA_pathway_activity_rebuilt.csv"), row.names = FALSE)

## ---------- 5. 与旧表比对 ----------
qc <- c()
qc <- c(qc, sprintf("run_time_min: %.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
qc <- c(qc, sprintf("cells: %d", ncol(mat)))
qc <- c(qc, sprintf("pathways: %d", length(gs_use)))
old_p <- file.path(OUT, "tables/ssGSEA_pathway_activity.csv")
if (file.exists(old_p)) {
  old <- read.csv(old_p, check.names = FALSE)
  key <- function(d) paste(d$cell_type, d$group, sep = "|")
  m <- merge(old, agg, by = c("cell_type","group"), suffixes = c(".old",".new"))
  qc <- c(qc, sprintf("matched rows old vs new: %d", nrow(m)))
  for (nm in names(gs_use)) {
    a <- m[[paste0(nm, ".old")]]; b <- m[[paste0(nm, ".new")]]
    if (is.null(a) || is.null(b)) next
    r <- suppressWarnings(cor(a, b, method = "spearman"))
    rc <- suppressWarnings(cor(a, b, method = "pearson"))
    qc <- c(qc, sprintf("%-24s spearman=%.3f pearson=%.3f maxabsdiff=%.3f meanabsdiff=%.3f",
                        nm, r, rc, max(abs(a-b)), mean(abs(a-b))))
  }
}
writeLines(qc, file.path(OUT, "tables/ssGSEA_rebuild_qc.txt"))
cat("\n>>> QC\n"); cat(paste(qc, collapse = "\n"), "\n")
cat(">>> DONE\n")
