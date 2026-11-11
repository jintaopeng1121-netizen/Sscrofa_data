## Fig 4 E: T 细胞功能/耗竭/活化 marker violin (PD-1/PD-L1/TIM-3/GZMB/PRF1/IFNG/CD69) × 3 组
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(cowplot); library(grid); library(gridExtra); library(patchwork)
  library(tidyr)
})

OUT <- r"(figures_final/Fig4)"
obj <- readRDS(r"(playbook_outputs/T_cell_subcluster/T_cell_subcluster_3grp.rds)")

grp_cols <- c("Control"="#999999", "PRRSV"="#E69F00", "APRR"="#0072B2")  # Okabe-Ito canonical
obj$group <- factor(obj$group, levels=c("Control","PRRSV","APRR"))

# 关键 marker（7 个，2 行排列）
markers <- c("PDCD1","HAVCR2","LAG3","CTLA4","TIGIT",
             "IFNG","GZMB","PRF1","CD69","CD28")
markers <- markers[markers %in% rownames(obj)]
cat("using markers:", paste(markers, collapse=","), "\n")

DefaultAssay(obj) <- "SCT"
# SCT assay v5 already merged, no JoinLayers needed

# 抽取表达矩阵
mat <- GetAssayData(obj, assay="SCT", layer="data")

# 长格式
df_long <- data.frame()
for (g in markers) {
  v <- mat[g, ]
  d <- data.frame(gene=g, group=obj$group, expr=as.numeric(v))
  df_long <- rbind(df_long, d)
}
df_long$gene <- factor(df_long$gene, levels=markers)

# violin plot 网格 (2 行 x 5 列)
plots <- list()
for (g in markers) {
  d <- df_long[df_long$gene==g, ]
  # 三组比较
  p1 <- suppressWarnings(wilcox.test(d$expr[d$group=="PRRSV"], d$expr[d$group=="Control"])$p.value)
  p2 <- suppressWarnings(wilcox.test(d$expr[d$group=="APRR"],  d$expr[d$group=="PRRSV"])$p.value)
  s1 <- ifelse(p1<0.0001,"****",ifelse(p1<0.001,"***",ifelse(p1<0.01,"**",ifelse(p1<0.05,"*","ns"))))
  s2 <- ifelse(p2<0.0001,"****",ifelse(p2<0.001,"***",ifelse(p2<0.01,"**",ifelse(p2<0.05,"*","ns"))))
  ymax <- max(d$expr, na.rm=TRUE)

  p <- ggplot(d, aes(x=group, y=expr, fill=group)) +
    geom_violin(scale="width", alpha=0.85, width=0.8, color="grey30", linewidth=0.3, trim=TRUE) +
    geom_boxplot(width=0.12, outlier.size=0.3, alpha=0.7, color="grey20") +
    scale_fill_manual(values=grp_cols) +
    guides(fill=guide_none()) +
    ggtitle(g) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold.italic", hjust=0.5, size=11),
          axis.title=element_blank(),
          axis.text.x=element_text(face="bold", size=9),
          panel.grid.major.x=element_blank(),
          plot.margin=margin(5,5,15,5)) +
    annotate("segment", x=1, xend=2, y=ymax*1.05, yend=ymax*1.05, linewidth=0.4) +
    annotate("text",   x=1.5, y=ymax*1.08, label=s1, size=3.0) +
    annotate("segment", x=2, xend=3, y=ymax*1.12, yend=ymax*1.12, linewidth=0.4) +
    annotate("text",   x=2.5, y=ymax*1.15, label=s2, size=3.0)

  plots[[g]] <- p
}

# 5 列 2 行
ncol_layout <- 5
nrow_layout <- ceiling(length(plots)/ncol_layout)
combined <- wrap_plots(plots, ncol=ncol_layout) +
  plot_annotation(title="Fig 4E · T cell functional markers (exhaustion / effector / activation)",
                  theme=theme(plot.title=element_text(face="bold", size=14, hjust=0.5)))

h <- 4 + 1.8*nrow_layout
w <- 3*ncol_layout
for (ext in c("png","pdf","svg")) {
  f <- file.path(OUT, paste0("Fig4E_tc_functional_markers.", ext))
  if (ext=="svg") { svg(f, width=w, height=h); print(combined); dev.off() }
  else if (ext=="pdf") { pdf(f, width=w, height=h); print(combined); dev.off() }
  else { png(f, width=w*120, height=h*120, res=120, bg="white"); print(combined); dev.off() }
  cat("  [saved]", basename(f), "\n")
}
