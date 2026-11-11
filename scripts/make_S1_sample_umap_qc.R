# ============================================================================
#  S1  Supplement  -  Sample-level UMAP + quality control
#  Premium 6-panel figure (PDF / PNG / SVG)
#  Panels:
#    A  UMAP colored by sample        (batch / integration check)
#    B  UMAP colored by group          (Okabe-Ito canonical)
#    C  UMAP colored by nCount_RNA      (continuous QC metric)
#    D  Violin nFeature_RNA per sample
#    E  Violin nCount_RNA   per sample
#    F  Cell-type composition per sample (100% stacked bar)
# ============================================================================
suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork)
  library(dplyr); library(scales)
})

obj <- readRDS("data/derived/seurat_reanalysis_v2.rds")
meta <- obj@meta.data

## ---- sample colors (deliberately distinct from Okabe-Ito group palette:
##      group = grey/orange/blue  ->  sample = purple/pink/teal) ----
samp_cols <- c("S14-control" = "#7E03A8",   # bold purple
               "S5-PRRSV"    = "#1B9E77",   # teal
               "S23-APPR4"   = "#D81B60")   # magenta / pink
grp_cols  <- c("Control" = "#999999", "PRRSV" = "#E69F00", "APRR" = "#0072B2") # Okabe-Ito

## ---- draw order: largest first so smaller samples show on top ----
samp_order <- names(sort(table(obj@meta.data$sample)))
obj@meta.data$sample <- factor(obj@meta.data$sample, levels = samp_order)

## ---- Panel A: UMAP by sample  (shuffle=TRUE randomizes draw order,
##      alpha<1 lets smaller samples show through) ----
pA <- DimPlot(obj, reduction = "umap", group.by = "sample",
              cols = samp_cols, raster = FALSE, pt.size = 0.5,
              shuffle = TRUE, seed = 42, alpha = 0.72) +
  ggtitle("A  UMAP by sample") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.7, "cm"))

## ---- Panel B: UMAP by group (Okabe-Ito) ----
pB <- DimPlot(obj, reduction = "umap", group.by = "group",
              cols = grp_cols, raster = FALSE, pt.size = 0.5,
              shuffle = TRUE, seed = 42, alpha = 0.72) +
  ggtitle("B  UMAP by treatment group") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.7, "cm"))

## ---- Panel C: UMAP colored by nCount_RNA (full dynamic range) ----
nc_max <- max(obj@meta.data$nCount_RNA)
pC <- FeaturePlot(obj, features = "nCount_RNA", reduction = "umap",
                  raster = FALSE, pt.size = 0.45,
                  cols = c("#F7FCF0", "#41B6C4", "#225EA8"),
                  min.cutoff = 0, max.cutoff = nc_max) +
  ggtitle("C  UMAP by nCount_RNA") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.8, "cm"),
        legend.title.align = 0.5)

## ---- Panel D / E: QC violins per sample (y-limit clipped, no scale conflict) ----
nf_max <- quantile(obj@meta.data$nFeature_RNA, 0.99)
nc_clip <- quantile(obj@meta.data$nCount_RNA, 0.99)

pD <- VlnPlot(obj, features = "nFeature_RNA", group.by = "sample",
              cols = samp_cols, pt.size = 0, raster = FALSE,
              y.max = nf_max) +
  ggtitle("D  nFeature_RNA per sample") +
  xlab("") + ylab("nFeature_RNA") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        legend.position = "none")

pE <- VlnPlot(obj, features = "nCount_RNA", group.by = "sample",
              cols = samp_cols, pt.size = 0, raster = FALSE,
              y.max = nc_clip) +
  ggtitle("E  nCount_RNA per sample") +
  xlab("") + ylab("nCount_RNA") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        legend.position = "none")

## ---- Panel F: cell-type composition per sample (100% stacked) ----
comp <- table(meta$cell_type, meta$sample)
comp_pct <- as.data.frame(prop.table(comp, 2) * 100)
names(comp_pct) <- c("cell_type", "sample", "pct")
comp_pct$cell_type <- factor(comp_pct$cell_type,
                             levels = rev(levels(factor(comp_pct$cell_type))))
comp_pct$sample <- factor(comp_pct$sample, levels = names(samp_cols))

# 12 distinct cell-type colors (colorblind-safe-ish qualitative)
ct_cols <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
             "#E377C2","#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78")
pF <- ggplot(comp_pct, aes(x = sample, y = pct, fill = cell_type)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = ct_cols) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Cell fraction (%)", fill = "Cell type") +
  ggtitle("F  Cell-type composition per sample") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major.x = element_blank())

## ---- compose ----
top <- pA + pB + pC + plot_layout(ncol = 3, widths = c(1, 1, 1))
bot <- pD + pE + pF + plot_layout(ncol = 3, widths = c(1, 1, 1.15))
p_all <- top / bot +
  plot_annotation(
    title = "S1.  Sample-level UMAP embedding and single-cell quality control",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                             margin = margin(b = 12))))

## ---- save 3 formats ----
outdir <- "_for_supplementary_2026-08-05"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ggsave(file.path(outdir, "S1_sample_umap_qc.png"),
       p_all, width = 16, height = 9.2, dpi = 300, bg = "white")
ggsave(file.path(outdir, "S1_sample_umap_qc.pdf"),
       p_all, width = 16, height = 9.2, bg = "white")
ggsave(file.path(outdir, "S1_sample_umap_qc.svg"),
       p_all, width = 16, height = 9.2, bg = "white")
cat("S1 sample UMAP + QC written (png/pdf/svg) to", outdir, "\n")
