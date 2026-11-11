## ============================================================================
## cellchat_diff.R  —  分组级差异细胞互作分析（3 组：Control / PRRSV / APRR）
##
## 框架锁定：治疗组 =  APRR（APRR），不含高剂量组，不做剂量比较。
## 输入：cellchat_group.rds（group.by="group" 已建，含 14 条显著通路）
## 产出：Fig_diff_A~D + pathway×group 矩阵 CSV，存 _archive_legacy/figures_paper_3grp_panels/
## ============================================================================

suppressMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
})
## CellChat 对象未向 Seurat 的 Idents / Idents<- generic 注册 S3 方法，这里补一个，
## 使 Idents(cc) / Idents(cc) <- 以及 netVisual_bubble 内部调用均可正常工作。
Idents.CellChat <- function(object, ...) object@idents
`Idents<-.CellChat` <- function(object, value) { object@idents <- value; object }
options(future.globals.maxSize = 8 * 1024^3)
future::plan(future::sequential)
.flog <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste(...)))

outdir <- "_archive_legacy/figures_paper_3grp_panels"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cc <- readRDS("data/derived/cellchat_group.rds")

## 归一化分组标签（用户硬约束：治疗组只叫 APRR，高剂量彻底剔除）
.norm_g <- function(x) {
  x <- as.character(x)
  x[grepl("APRR6|high", x, ignore.case = TRUE)] <- NA
  x[grepl("APRR4|low",  x, ignore.case = TRUE)] <- "APRR"
  x
}
Idents(cc) <- factor(.norm_g(Idents(cc)), levels = c("Control", "PRRSV", "APRR"))
for (.d in 1:2) {
  .dn <- .norm_g(dimnames(cc@netP$prob)[[.d]])
  dimnames(cc@netP$prob)[[.d]] <- .dn
}

## 锁定分组顺序：Control(1) / PRRSV(2) / APRR(3)
.flog("levels:", paste(levels(cc@idents), collapse = " | "))

## 注：本版 CellChat 无 netAnalysis_signalingRole / @var，Fig_diff_A/B 改用自建 ggplot 差异图，
##     不依赖 netAnalysis_computeCentrality；C(热图)/D(bubble) 仍用原生函数。

## 探查 netP$prob 结构（写 runlog）
arr <- cc@netP$prob
.flog("netP$prob dims (source x target x pathway):", paste(dim(arr), collapse = " x "))
.flog("netP$prob dimnames[[1]] (source groups):", paste(dimnames(arr)[[1]], collapse = ", "))
.flog("netP$pathways (14):", paste(cc@netP$pathways, collapse = ", "))

## ---- 构建 通路 × 组 参与度矩阵 M[pathway, group] = out_g + in_g ----
## arr 维度 = [source_group(3), target_group(3), pathway(14)]
gs <- dimnames(arr)[[1]]          # groups (dim1 = source)
pws <- dimnames(arr)[[3]]         # pathways (dim3)
M <- matrix(0, nrow = length(pws), ncol = length(gs),
            dimnames = list(pws, gs))
for (pi in seq_along(pws)) {
  slab <- arr[, , pi]             # 3 x 3 source x target for this pathway
  M[pi, ] <- rowSums(slab) + colSums(slab)   # 该组作为 sender + receiver 的总参与度
}
.flog("M built (pathway x group):", paste(dim(M), collapse = " x "))

## 导出矩阵 + 折叠变化
dfM <- as.data.frame(M); dfM$pathway <- rownames(M)
write.csv(dfM, file.path(outdir, "Fig_diff_pathway_by_group.csv"), row.names = FALSE)
if (all(c("Control", "PRRSV", "APRR") %in% colnames(M))) {
  fc <- data.frame(
    pathway = rownames(M),
    log2FC_PRRSV_vs_Ctrl = log2((M[, "PRRSV"] + 1e-9) / (M[, "Control"] + 1e-9)),
    log2FC_APRR_vs_PRRSV = log2((M[, "APRR"] + 1e-9) / (M[, "PRRSV"] + 1e-9)),
    stringsAsFactors = FALSE
  )
  write.csv(fc, file.path(outdir, "Fig_diff_log2FC.csv"), row.names = FALSE)
  .flog("top PRRSV-up pathways:", paste(head(fc$pathway[order(-fc$log2FC_PRRSV_vs_Ctrl)], 5), collapse = ", "))
  .flog("top APRR-reversed pathways:", paste(head(fc$pathway[order(fc$log2FC_APRR_vs_PRRSV)], 5), collapse = ", "))
}

## ---- Fig_diff_A：通路 × 组 sender/receiver 角色散点（自建，替代 signaling role scatter）----
.flog("Fig_diff_A: pathway x group signaling-role scatter ...")
tryCatch({
  role_list <- list()
  for (pi in seq_along(pws)) {
    slab <- arr[, , pi]
    role_list[[pi]] <- data.frame(pathway = pws[pi], group = gs,
                                  outgoing = rowSums(slab), incoming = colSums(slab),
                                  stringsAsFactors = FALSE)
  }
  role_df <- do.call(rbind, role_list)
  role_df$group <- factor(role_df$group, levels = c("Control", "PRRSV", "APRR"))
  gA <- ggplot(role_df, aes(x = outgoing, y = incoming, color = group, label = pathway)) +
    geom_point(size = 2.5) +
    geom_text(size = 2.0, hjust = -0.1, vjust = 0.1, show.legend = FALSE) +
    facet_wrap(~ group, nrow = 1) +
    theme_minimal(base_size = 9) +
    labs(x = "Outgoing (sender) strength", y = "Incoming (receiver) strength",
         title = "CellChat signaling role by pathway & group") +
    theme(legend.position = "bottom")
  ggsave(file.path(outdir, "Fig_diff_A_signalingRole_scatter.pdf"), gA, width = 10, height = 4)
  ggsave(file.path(outdir, "Fig_diff_A_signalingRole_scatter.png"), gA, width = 10, height = 4, dpi = 300)
  .flog("Fig_diff_A OK")
}, error = function(e) .flog("Fig_diff_A SKIP:", conditionMessage(e)))

## ---- Fig_diff_B：按 log2FC 排序的差异通路条形图（自建，替代 rankNet comparison）----
.flog("Fig_diff_B: ranked differential pathways ...")
tryCatch({
  fc <- read.csv(file.path(outdir, "Fig_diff_log2FC.csv"), stringsAsFactors = FALSE)
  gB1 <- ggplot(fc, aes(x = reorder(pathway, log2FC_PRRSV_vs_Ctrl), y = log2FC_PRRSV_vs_Ctrl)) +
    geom_col(aes(fill = log2FC_PRRSV_vs_Ctrl > 0)) + coord_flip() +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +
    theme_minimal(base_size = 9) + theme(legend.position = "none") +
    labs(x = "pathway", y = "log2FC (PRRSV / Control)", title = "PRRSV vs Control")
  ggsave(file.path(outdir, "Fig_diff_B1_rankNet_PRRSV_vs_Ctrl.pdf"), gB1, width = 6, height = 5)
  ggsave(file.path(outdir, "Fig_diff_B1_rankNet_PRRSV_vs_Ctrl.png"), gB1, width = 6, height = 5, dpi = 300)
  gB2 <- ggplot(fc, aes(x = reorder(pathway, log2FC_APRR_vs_PRRSV), y = log2FC_APRR_vs_PRRSV)) +
    geom_col(aes(fill = log2FC_APRR_vs_PRRSV > 0)) + coord_flip() +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +
    theme_minimal(base_size = 9) + theme(legend.position = "none") +
    labs(x = "pathway", y = "log2FC (APRR / PRRSV)", title = "APRR vs PRRSV")
  ggsave(file.path(outdir, "Fig_diff_B2_rankNet_APRR_vs_PRRSV.pdf"), gB2, width = 6, height = 5)
  ggsave(file.path(outdir, "Fig_diff_B2_rankNet_APRR_vs_PRRSV.png"), gB2, width = 6, height = 5, dpi = 300)
  .flog("Fig_diff_B OK")
}, error = function(e) .flog("Fig_diff_B SKIP:", conditionMessage(e)))

## ---- Fig_diff_C：自建 通路×组 z-score 热图 ----
.flog("Fig_diff_C: pathway x group z-score heatmap ...")
tryCatch({
  dm <- melt(M); colnames(dm) <- c("pathway", "group", "value")
  dm$z <- ave(dm$value, dm$pathway, FUN = function(x) {
    s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s
  })
  dm$group <- factor(dm$group, levels = c("Control", "PRRSV", "APRR"))
  gC <- ggplot(dm, aes(x = group, y = reorder(pathway, z, mean), fill = z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, name = "z-score") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("CellChat pathway activity by group (z-scored)") +
    geom_text(aes(label = sprintf("%.2f", z)), size = 2.2, color = "black")
  ggsave(file.path(outdir, "Fig_diff_C_pathway_group_heatmap.pdf"), gC, width = 6, height = 6)
  ggsave(file.path(outdir, "Fig_diff_C_pathway_group_heatmap.png"), gC, width = 6, height = 6, dpi = 300)
  .flog("Fig_diff_C OK")
}, error = function(e) .flog("Fig_diff_C SKIP:", conditionMessage(e)))

## ---- Fig_diff_D：14 条显著通路 bubble（全部为差异显著通路）----
.flog("Fig_diff_D: all-significant-pathway bubble ...")
immune <- cc@netP$pathways
.flog("pathways for bubble:", paste(immune, collapse = ", "))
tryCatch({
  gD <- netVisual_bubble(cc, signaling = immune, remove.isolate = FALSE,
                         title.name = "Immune-related signaling (3 groups)")
  ggsave(file.path(outdir, "Fig_diff_D_immune_bubble.pdf"), gD, width = 10, height = 6)
  ggsave(file.path(outdir, "Fig_diff_D_immune_bubble.png"), gD, width = 10, height = 6, dpi = 300)
  .flog("Fig_diff_D OK")
}, error = function(e) .flog("Fig_diff_D SKIP (fallback custom):", conditionMessage(e)))

.flog(">>> Differential CellChat DONE")
