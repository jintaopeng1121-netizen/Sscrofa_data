## Fig 3 C: Neutrophil state score violin (3 pathways x 3 groups)
## 手工算 score（不依赖 AddModuleScore），基因均值法
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(cowplot); library(grid); library(gridExtra); library(patchwork)
  library(tidyr)
})

OUT <- r"(figures_final/Fig3)"
obj <- readRDS(r"(playbook_outputs/neu_subcluster.rds)")

# 3 个 state 的代表基因集（取自旧版 Figure 3 B/C 命名）
gene_sets <- list(
  proinf = c("CXCL8","IL1B","NLRP3","S100A8","S100A9","SAA3","FCN2","THBS1","SERPINB1"),
  isg    = c("ISG15","IFIT1","IFIT2","IFIT3","ISG20","IFI6","OASL","MX1","USP18"),
  pdl1   = c("CD274","IL1R2","UVRAG","THRB")
)

# 配色（Okabe-Ito canonical, colorblind-safe; 与 export_panels_nature.R GRP_COLS 一致）
grp_cols <- c("Control"="#999999", "PRRSV"="#E69F00", "APRR"="#0072B2")
obj$group <- factor(obj$group, levels=c("Control","PRRSV","APRR"))

# 手工算 score：每细胞 = 该 state 基因在 SCT normalized 矩阵的均值
mat <- GetAssayData(obj, assay="SCT", layer="data")
cat("SCT matrix dim:", dim(mat), "\n")

calc_score <- function(genes) {
  common <- intersect(genes, rownames(mat))
  cat(sprintf("  using %d/%d genes for score\n", length(common), length(genes)))
  if (length(common) < 1) return(rep(0, ncol(mat)))
  colMeans(mat[common, , drop=FALSE])
}

scores <- data.frame(
  cell   = colnames(obj),
  group  = obj$group,
  proinf = calc_score(gene_sets$proinf),
  isg    = calc_score(gene_sets$isg),
  pdl1   = calc_score(gene_sets$pdl1)
)

scores_long <- pivot_longer(scores, cols=c("proinf","isg","pdl1"),
                            names_to="state", values_to="score")
scores_long$state <- factor(scores_long$state, levels=c("proinf","pdl1","isg"))

# 配对比较
stat_test <- function(df) {
  v_ctrl <- df$score[df$group=="Control"]
  v_prr  <- df$score[df$group=="PRRSV"]
  v_apr  <- df$score[df$group=="APRR"]
  p1 <- suppressWarnings(wilcox.test(v_prr, v_ctrl)$p.value)
  p2 <- suppressWarnings(wilcox.test(v_apr, v_prr)$p.value)
  list(p1=p1, p2=p2)
}
sig_stars <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  return("ns")
}

plots <- list()
for (st in c("proinf","pdl1","isg")) {
  df <- scores_long[scores_long$state==st, ]
  pval <- stat_test(df)
  s1 <- sig_stars(pval$p1); s2 <- sig_stars(pval$p2)

  ymax <- max(df$score, na.rm=TRUE)
  p <- ggplot(df, aes(x=group, y=score, fill=group)) +
    geom_violin(scale="width", alpha=0.85, width=0.8, color="grey30", linewidth=0.3, trim=TRUE) +
    geom_boxplot(width=0.12, outlier.size=0.4, alpha=0.7, color="grey20") +
    scale_fill_manual(values=grp_cols) +
    guides(fill=guide_none()) +
    ggtitle(st) +
    theme_minimal(base_size=11) +
    theme(plot.title=element_text(face="bold", hjust=0.5, size=12),
          axis.title=element_blank(),
          axis.text.x=element_text(face="bold", size=10),
          panel.grid.major.x=element_blank(),
          plot.margin=margin(5,5,15,5))

  p <- p + annotate("segment", x=1, xend=2, y=ymax*1.05, yend=ymax*1.05, linewidth=0.4) +
           annotate("text",   x=1.5, y=ymax*1.08, label=s1, size=3.5) +
           annotate("segment", x=2, xend=3, y=ymax*1.12, yend=ymax*1.12, linewidth=0.4) +
           annotate("text",   x=2.5, y=ymax*1.15, label=s2, size=3.5)

  plots[[st]] <- p
}

combined <- wrap_plots(plots, ncol=3) +
  plot_annotation(title="Fig 3C · Neutrophil pathway state score (gene-mean)",
                  theme=theme(plot.title=element_text(face="bold", size=14, hjust=0.5)))

for (ext in c("png","pdf","svg")) {
  f <- file.path(OUT, paste0("Fig3C_pathway_score_violin.", ext))
  if (ext=="svg") { svg(f, width=11, height=5); print(combined); dev.off() }
  else if (ext=="pdf") { pdf(f, width=11, height=5); print(combined); dev.off() }
  else { png(f, width=11*120, height=5*120, res=120, bg="white"); print(combined); dev.off() }
  cat("  [saved]", basename(f), "\n")
}
cat("p values:\n")
for (st in c("proinf","pdl1","isg")) {
  df <- scores_long[scores_long$state==st, ]
  pv <- stat_test(df)
  cat(sprintf("  %s: PRRSV_vs_Control p=%.2e  | APRR_vs_PRRSV p=%.2e\n",
              st, pv$p1, pv$p2))
}
