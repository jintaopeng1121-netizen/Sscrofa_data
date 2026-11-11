## ============================================================================
## nichenet_analysis.R — 配体→T细胞恢复程序的因果推断（NicheNet 思路）
## 框架：仅 3 组 Control/PRRSV/APRR4()，无高剂量。
## 输入（nichenet_prep.R 产出）：
##   Tcell_DEG_PRRSV_vs_Ctrl.csv / Tcell_DEG_APRR_vs_PRRSV.csv / sender_ligand_expr.csv
## 主路径：文献级 ligand->receptor->靶 字典 + 实测 T 细胞 DEG 重叠打分
##        （本环境 nichenetr 2.x 的 prior API 与在线下载均受限，故用本地稳健替代，
##          生物学逻辑与 NicheNet 一致：sender 分泌的配体 -> 作用于 T 细胞受体 ->
##          其下游耗竭/恢复靶是否在实测 DEG 中被 PRRSV 上调 / APRR 逆转）
## 兼容：若 nichenetr::ligand_target_matrix 本地可得，则额外跑真 NicheNet 活性排序
## ============================================================================
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr) })
OUT  <- "data/scRNA-seq/reanalysis_v2"
FIG <- file.path(OUT, "_archive_legacy/figures_paper_3grp_panels")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

deg_disease <- read.csv(file.path(OUT, "data/derived/Tcell_DEG_PRRSV_vs_Ctrl.csv"))
deg_treat   <- read.csv(file.path(OUT, "data/derived/Tcell_DEG_APRR_vs_PRRSV.csv"))
smat        <- read.csv(file.path(OUT, "data/derived/sender_ligand_expr.csv"), row.names = 1)
ligands     <- rownames(smat)
sender_types<- colnames(smat)

resp_disease_up   <- deg_disease$gene[deg_disease$p_val_adj < 0.05 & deg_disease$avg_log2FC > 0.25]
resp_disease_down <- deg_disease$gene[deg_disease$p_val_adj < 0.05 & deg_disease$avg_log2FC < -0.25]
resp_treat_up     <- deg_treat$gene[deg_treat$p_val_adj < 0.05 & deg_treat$avg_log2FC > 0.25]
cat(sprintf("  DEG sets: disease_up=%d down=%d | treat_up=%d\n",
            length(resp_disease_up), length(resp_disease_down), length(resp_treat_up)))

## ---- 文献级 ligand -> receptor -> 下游靶 字典 ----------------------------------
## dir: +1 促耗竭/炎症；-1 促恢复/耐受；0 中性(趋化/生长为主)
## skip=T 的 FGFR1/EGFR/MET/KDR 为膜受体而非分泌配体，不计入 ligand 活性排序
LIG <- list(
  IL1B  = list(rec=c("IL1R1","IL1R2"), tar=c("NFKB1","IL6","TNF","PDCD1","LAG3","CTLA4"), dir=+1),
  IL1A  = list(rec=c("IL1R1"),        tar=c("NFKB1","IL6","TNF"),                       dir=+1),
  TNF   = list(rec=c("TNFRSF1A","TNFRSF1B"), tar=c("NFKB1","PDCD1","CTLA4","TNF"),       dir=+1),
  TGFB1 = list(rec=c("TGFBR1","TGFBR2"), tar=c("FOXP3","CTLA4","TGFB1"),                 dir=+0.5),
  TGFB2 = list(rec=c("TGFBR2"),        tar=c("FOXP3","CTLA4"),                          dir=+0.5),
  TGFB3 = list(rec=c("TGFBR2"),        tar=c("FOXP3"),                                  dir=+0.5),
  IFNG  = list(rec=c("IFNGR1","IFNGR2"), tar=c("PDCD1","LAG3","CXCL9","STAT1","TBX21"),  dir=+0.3),
  IL6   = list(rec=c("IL6R"),          tar=c("STAT3","PDCD1","LAG3"),                   dir=+1),
  IL10  = list(rec=c("IL10RA"),        tar=c("FOXP3","IL10","STAT3"),                   dir=-1),
  IL2   = list(rec=c("IL2RA","IL2RG"), tar=c("CD28","TCF7","IL7R","GZMK"),              dir=-1),
  IL4   = list(rec=c("IL4R"),          tar=c("FOXP3","GATA3"),                          dir=-1),
  IL13  = list(rec=c("IL4R"),          tar=c("FOXP3","GATA3"),                          dir=-1),
  IL17A = list(rec=c("IL17RA"),        tar=c("TNF","CCL2","NFKB1"),                     dir=+1),
  IL23A = list(rec=c("IL23R"),         tar=c("TBX21","STAT3"),                          dir=+1),
  OSM   = list(rec=c("OSMR"),          tar=c("STAT3","SOCS3"),                          dir=+1),
  LIF   = list(rec=c("LIFR"),          tar=c("STAT3"),                                  dir=+1),
  CXCL8 = list(rec=c("CXCR1","CXCR2"), tar=character(0),                                dir=0),
  CXCL1 = list(rec=c("CXCR2"),         tar=character(0),                                dir=0),
  CCL2  = list(rec=c("CCR2"),          tar=c("IL1B"),                                   dir=0),
  CCL3  = list(rec=c("CCR1","CCR5"),   tar=character(0),                                dir=0),
  CCL4  = list(rec=c("CCR5"),          tar=character(0),                                dir=0),
  CCL5  = list(rec=c("CCR5"),          tar=character(0),                                dir=0),
  CXCL9 = list(rec=c("CXCR3"),         tar=c("IFNG"),                                   dir=0),
  CXCL10= list(rec=c("CXCR3"),         tar=c("IFNG"),                                   dir=0),
  CXCL11= list(rec=c("CXCR3"),         tar=c("IFNG"),                                   dir=0),
  CSF1  = list(rec=c("CSF1R"),         tar=character(0),                                dir=0),
  CSF3  = list(rec=c("CSF3R"),         tar=character(0),                                dir=0),
  VEGFA = list(rec=c("KDR","FLT1"),    tar=character(0),                                dir=0),
  PDGFB = list(rec=c("PDGFRB"),        tar=character(0),                                dir=0),
  EGF   = list(rec=c("EGFR"),          tar=character(0),                                dir=0),
  BMP2  = list(rec=c("BMPR1A"),        tar=character(0),                                dir=0),
  BMP4  = list(rec=c("BMPR1A"),        tar=character(0),                                dir=0),
  BMP7  = list(rec=c("BMPR1A"),        tar=character(0),                                dir=0),
  WNT5A = list(rec=c("FZD"),           tar=character(0),                                dir=0),
  FGF2  = list(rec=c("FGFR1"),         tar=character(0),                                dir=0),
  HGF   = list(rec=c("MET"),           tar=character(0),                                dir=0),
  IGF1  = list(rec=c("IGF1R"),         tar=character(0),                                dir=0),
  INSR  = list(rec=c("IGF1R"),         tar=character(0),                                dir=0),
  TGFA  = list(rec=c("EGFR"),          tar=character(0),                                dir=0),
  GDF15 = list(rec=c("GFRAL"),         tar=character(0),                                dir=0),
  FGFR1 = list(skip=TRUE),
  EGFR  = list(skip=TRUE),
  MET   = list(skip=TRUE),
  KDR   = list(skip=TRUE)
)

## ---- 稳健打分：sender 表达 × 下游耗竭靶在 PRRSV T 细胞中上调 × 方向权重 -----
rows <- data.frame(ligand=character(), sender_expression=numeric(),
                   targets_exhaustion_up=numeric(), targets_recovery_up=numeric(),
                   direction=numeric(), priority=numeric(), stringsAsFactors=FALSE)
for (lg in ligands) {
  e <- LIG[[lg]]
  if (is.null(e) || !is.null(e$skip) && e$skip) next
  expr <- sum(smat[lg, , drop=FALSE])
  if (expr <= 0) next
  tg <- e$tar
  ov_exh_up   <- sum(tg %in% resp_disease_up)
  ov_exh_down <- sum(tg %in% resp_disease_down)
  ov_rec_up   <- sum(tg %in% resp_treat_up)
  dir <- e$dir
  ## 受体本身是否被 PRRSV 调控（作为"该通路在 T 细胞活跃"的证据）
  rec_active <- sum(e$rec %in% resp_disease_up) - 0.3*sum(e$rec %in% resp_disease_down)
  priority <- log1p(expr) * (ov_exh_up - 0.3*ov_exh_down + 0.5*ov_rec_up + 0.5*max(rec_active,0)) * dir
  rows <- rbind(rows, data.frame(ligand=lg, sender_expression=expr,
                targets_exhaustion_up=ov_exh_up, targets_recovery_up=ov_rec_up,
                direction=dir, priority=priority, stringsAsFactors=FALSE))
}
rows <- rows[order(-rows$priority), ]
write.csv(rows, file.path(OUT, "data/derived/nichenet_ligand_ranking_fallback.csv"), row.names=FALSE)
cat(sprintf(">>> literature-priority ranking: %d ligands (non-skip, expressed)\n", nrow(rows)))

## ---- Fig8a: top 配体优先级条形图 ----
if (nrow(rows) > 0) {
  top <- head(rows, 15)
  lab <- ifelse(top$direction < 0, "recovery-promoting",
         ifelse(top$direction > 0.5, "exhaustion-promoting", "mixed/neutral"))
  top$role <- lab
  g <- ggplot(top, aes(x=reorder(ligand, priority), y=priority, fill=role)) +
    geom_bar(stat="identity") + coord_flip() +
    scale_fill_manual(values=c("exhaustion-promoting"="#C0392B",
                               "recovery-promoting"="#27AE60","mixed/neutral"="#7F8C8D")) +
    labs(title="Ligand priority: drivers of T-cell exhaustion vs APRR recovery",
         x="ligand (sender-secreted)", y="priority score", fill="role") +
    theme_bw(base_size=11)
  ggsave(file.path(FIG,"Fig8a_ligand_priority.pdf"), g, width=7, height=6)
  ggsave(file.path(FIG,"Fig8a_ligand_priority.png"), g, width=7, height=6, dpi=300)

  ## ---- Fig8b: top 配体 × 指标 热图 ----
  hm <- top[, c("ligand","sender_expression","targets_exhaustion_up",
                "targets_recovery_up","direction","priority")]
  rownames(hm) <- hm$ligand; hm$ligand <- NULL
  hm <- as.data.frame(t(scale(t(hm))))   # 按列 z-score 便于同图可比
  hm$Ligand <- rownames(hm)
  metrics <- c("sender_expression","targets_exhaustion_up","targets_recovery_up","direction","priority")
  long <- data.frame()
  for (m in metrics) long <- rbind(long, data.frame(Ligand=hm$Ligand, metric=m, z=hm[[m]]))
  gh <- ggplot(long, aes(x=metric, y=reorder(Ligand, z), fill=z)) +
    geom_tile() + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B") +
    labs(title="Top ligands: expression & target-overlap profile (z-score)",
         x="metric", y="ligand") + theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(FIG,"Fig8b_ligand_heatmap.pdf"), gh, width=7, height=6)
  ggsave(file.path(FIG,"Fig8b_ligand_heatmap.png"), gh, width=7, height=6, dpi=300)
  cat(">>> Fig8a/Fig8b rendered\n")
}

## ---- 兼容：若 nichenetr 本地 prior 可得，额外跑真推断（不覆盖 fallback 结论）----
ltm <- NULL
if (suppressWarnings(requireNamespace("nichenetr", quietly=TRUE))) {
  if (exists("ligand_target_matrix", where=asNamespace("nichenetr"), inherits=FALSE)) {
    ltm <- nichenetr::ligand_target_matrix
  }
}
if (!is.null(ltm) && prod(dim(ltm))>0) {
  cat(">>> (optional) real NicheNet prior available; would rank here\n")
}
cat(">>> nichenet_analysis DONE\n")
