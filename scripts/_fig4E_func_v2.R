#!/usr/bin/env Rscript
## _fig4E_func_v2.R  -- Fig 4 panel E (T-cell functional dot plot) -- CORRECTED classification
##
## WHY THIS FILE EXISTS
##   The original builder (_fig4D_func.R, line ~44) classified genes with a 2-way ifelse:
##       ifelse(gene %in% c(IFNG,GNLY,NKG7,GZMB,PRF1,GZMK),
##              "Cytotoxicity / effector", "Exhaustion / checkpoint")
##   Everything not in the 6-gene cytotoxicity list was dumped into the
##   "Exhaustion / checkpoint" facet -- including CCR7, TCF7, CD44, IL7R, CD69,
##   which are progenitor / survival / activation markers, not exhaustion markers.
##   The figure label then propagated into the manuscript text as a biological claim.
##
##   This v2 keeps the SAME data (Tcell_functional_logmean.csv), the SAME dot encoding
##   and the SAME 3-group palette; ONLY the facet assignment and the facet subtitle
##   are corrected. No numbers are touched.
##
##   Original script is left untouched for traceability (_fig4D_func.R).
##
## OUTPUT: figures_final/Fig4E_tc_functional_v2.{svg,pdf,png}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(svglite)
})

ROOT <- "data/scRNA-seq/reanalysis_v2"
OUT  <- file.path(ROOT, "data", "fig4_data")
FIG  <- file.path(ROOT, "figures_final")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

W <- 6.0
H <- 9.6          # was 8.0 for a 2-facet layout; 3 facets need more vertical room
FIGNAME <- "Fig4E_tc_functional_v2"

grp_col <- c(Control = "#5B8DB8", PRRSV = "#D1495B", APRR = "#2A9D8F")

theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "#888", linewidth = 0.4),
          axis.ticks = element_line(color = "#888", linewidth = 0.4),
          legend.key = element_blank(),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "#555",
                                       margin = margin(b = 8)),
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "#f0f0f0", color = NA))
}

## ---- corrected functional grouping (3 classes) ----
GENE_CLASS <- c(
  IFNG = "Cytotoxicity / effector", GNLY = "Cytotoxicity / effector",
  NKG7 = "Cytotoxicity / effector", GZMB = "Cytotoxicity / effector",
  PRF1 = "Cytotoxicity / effector", GZMK = "Cytotoxicity / effector",
  PDCD1 = "Checkpoint / inhibitory", CTLA4 = "Checkpoint / inhibitory",
  LAG3  = "Checkpoint / inhibitory", TIGIT = "Checkpoint / inhibitory",
  HAVCR2 = "Checkpoint / inhibitory", TGFB1 = "Checkpoint / inhibitory",
  IL7R = "Progenitor / activation", CD44 = "Progenitor / activation",
  CD69 = "Progenitor / activation", TCF7 = "Progenitor / activation",
  CCR7 = "Progenitor / activation"
)

build_Fig4E <- function() {
  tf <- read.csv(file.path(OUT, "Tcell_functional_logmean.csv"), stringsAsFactors = FALSE)
  stopifnot(all(tf$gene %in% names(GENE_CLASS)))          # no silent leftovers
  tf$cat <- factor(unname(GENE_CLASS[tf$gene]),
                   levels = c("Cytotoxicity / effector",
                              "Checkpoint / inhibitory",
                              "Progenitor / activation"))

  tf_long <- tf %>% pivot_longer(cols = c("Control", "PRRSV", "APRR"),
                                 names_to = "group", values_to = "val")
  tf_long$group <- factor(tf_long$group, levels = c("Control", "PRRSV", "APRR"))

  ## gene order inside each facet (ggplot draws y from bottom up, so levels are reversed
  ## relative to the printed top-to-bottom order kept from the original figure)
  gene_order <- c("GNLY", "NKG7", "GZMB", "PRF1", "GZMK", "IFNG",
                  "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TGFB1",
                  "IL7R", "CD44", "CD69", "TCF7", "CCR7")
  tf_long$gene <- factor(tf_long$gene, levels = gene_order)

  ggplot(tf_long, aes(x = group, y = gene)) +
    geom_line(aes(group = gene), color = "#D9D9D9", linewidth = 0.6) +
    geom_point(aes(fill = val), shape = 21, size = 5.8,
               color = "white", stroke = 0.4) +
    scale_fill_gradientn(colors = c("#EAF6F3", "#7FCDC2", "#2A9D8F", "#176B63"),
                         name = "log1p(CPMx1e4)") +
    facet_wrap(~ cat, ncol = 1, scales = "free_y") +
    theme_pub() +
    theme(legend.position = "right",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = "T-cell functional gene programmes",
         subtitle = "Mean expression per group (n = 1 animal per group)",
         fill = "Expression")
}

if (sys.nframe() == 0) {
  p <- build_Fig4E()
  svglite(file.path(FIG, paste0(FIGNAME, ".svg")), width = W, height = H)
  print(p); dev.off()
  cairo_pdf(file.path(FIG, paste0(FIGNAME, ".pdf")), width = W, height = H)
  print(p); dev.off()
  png(file.path(FIG, paste0(FIGNAME, ".png")), width = W * 200, height = H * 200,
      res = 200, bg = "white")
  print(p); dev.off()
  cat(sprintf("wrote %s.{svg,pdf,png} (%.1f x %.1f in)\n", FIGNAME, W, H))
}
