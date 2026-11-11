## FigS9 v2（2026-11-03）— ssGSEA pathway activity heatmap，基于可复现的重建表
## 与 fix_supp_S9_20260930.R 版式一致（12 面 × 12 通路，y 轴只标一次，固定色标）
## 输入：tables/ssGSEA_pathway_activity_gsva_curated.csv（由 _ssgsea_reproduce_and_fix.R 生成）
## 2026-11-03：色标方向改为常规（深色=高分）。旧版用 rev() 使深色=低分，与阅读直觉相反。
## 输出：_submission_assets/supplementary/FigS9_rebuilt.{png,pdf,svg}（不覆盖原文件）
setwd("data/scRNA-seq/reanalysis_v2")
suppressPackageStartupMessages({ library(ggplot2); library(reshape2); library(svglite) })

IN  <- commandArgs(trailingOnly = TRUE)
infile <- if (length(IN) >= 1) IN[1] else "tables/ssGSEA_pathway_activity_rebuilt.csv"
tag    <- if (length(IN) >= 2) IN[2] else "rebuilt"

ss <- read.csv(infile, check.names = FALSE)
stopifnot(nrow(ss) == 36)

pw_order <- c("T_cell_activation","TNF_signaling","Neutrophil_degranulation",
              "Macrophage_alternative","IL1_NFKB_signaling","IFN_gamma_response",
              "IFN_alpha_response","Hypoxia","EMT","Complement","Apoptosis","Antiviral")
pw_lab <- c("T cell activation","TNF signaling","Neutrophil degranulation",
            "Macrophage (alternative)","IL-1/NF-\u03baB signaling","IFN-\u03b3 response",
            "IFN-\u03b1 response","Hypoxia","EMT","Complement","Apoptosis","Antiviral")
names(pw_lab) <- pw_order

ct_lab <- c(Bcell="B cell", Ciliated="Ciliated", DC="DC", Endothelial="Endothelial",
            Epithelial="Epithelial", Fibroblast="Fibroblast", Macrophage="Macrophage",
            Monocyte="Monocyte", NK="NK cell", Neutrophil="Neutrophil",
            Tcell_CD4="CD4\u207a T cell", Tcell_CD8="CD8\u207a T cell")
grp_order <- c("APRR","Control","PRRSV")
grp_lab <- c(APRR="APRR", Control="Control", PRRSV="PRRSV")

m <- melt(ss, id.vars = c("cell_type","group"), variable.name = "pathway", value.name = "score")
m$pathway <- factor(m$pathway, levels = rev(pw_order))
m$cell_type <- factor(m$cell_type, levels = names(ct_lab))
m$group <- factor(m$group, levels = grp_order)
stopifnot(!any(is.na(m$pathway)), !any(is.na(m$cell_type)), !any(is.na(m$group)))

rng <- range(m$score, na.rm = TRUE)
pad <- diff(rng) * 0.02
brk <- pretty(rng, n = 5)

p <- ggplot(m, aes(x = group, y = pathway, fill = score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~cell_type, ncol = 4, scales = "fixed",
             labeller = labeller(cell_type = ct_lab)) +
  scale_x_discrete(labels = grp_lab, expand = c(0, 0)) +
  scale_y_discrete(labels = pw_lab, expand = c(0, 0)) +
  scale_fill_gradientn(colours = c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1",
                                       "#6BAED6","#4292C6","#2171B5","#08519C","#08306B"),
                       limits = c(rng[1] - pad, rng[2] + pad), breaks = brk,
                       name = "ssGSEA\nscore") +
  labs(x = NULL, y = NULL,
       title = "Pathway activity by ssGSEA across cell types and groups") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 6)),
    strip.background = element_rect(fill = "grey92", color = "grey70", linewidth = 0.4),
    strip.text = element_text(size = 9.5, face = "bold", margin = margin(2, 2, 2, 2)),
    axis.text.x = element_text(size = 7.2, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 7.2, color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey55", linewidth = 0.5),
    legend.title = element_text(size = 8, lineheight = 0.9),
    legend.text = element_text(size = 7.5),
    legend.key.height = unit(0.32, "cm"),
    legend.margin = margin(0, 0, 0, 2)
  )

dir.create("_submission_assets/supplementary", showWarnings = FALSE, recursive = TRUE)
base <- sprintf("_submission_assets/supplementary/FigS9_%s", tag)
ggsave(paste0(base, ".png"), p, width = 7.0, height = 5.9, dpi = 300, bg = "white")
ggsave(paste0(base, ".pdf"), p, width = 7.0, height = 5.9, bg = "white")
ggsave(paste0(base, ".svg"), p, width = 7.0, height = 5.9, bg = "white")
cat(sprintf("done: %s.{png,pdf,svg} | fill range %.3f to %.3f\n", base, rng[1], rng[2]))
