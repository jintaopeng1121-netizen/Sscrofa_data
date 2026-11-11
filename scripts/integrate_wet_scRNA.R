suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(readr)
})
options(stringsAsFactors=FALSE, width=120)

BASE <- "data/scRNA-seq/reanalysis_v2"
FIG  <- file.path(BASE, "figures_integration")
TAB  <- file.path(BASE, "tables", "integration")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)

cat(">>> load rds\n")
obj <- readRDS(file.path(BASE, "data/derived/seurat_reanalysis_v2.rds"))
DefaultAssay(obj) <- "SCT"
Idents(obj) <- "cell_type"
obj$cell_type <- factor(as.character(obj$cell_type))
obj$group     <- factor(as.character(obj$group),
                        levels=c("Control","PRRSV","APRR"))
mat <- GetAssayData(obj, layer="data", assay="SCT")
grp <- as.character(obj$group)
message("cells=", ncol(mat), " genes=", nrow(mat))

genesets <- list(
  PRRSV_receptor = c("CD163","LY86","MARCO","MSR1"),
  inflam_myeloid = c("IL1B","CXCL8","TNF","CCL2","CCL3","CCL4","S100A8","S100A9","IL6"),
  neutrophil     = c("MPO","CXCL8","S100A8","S100A9","IL1B","DEFA1","DEFA3"),
  tcell_ifn      = c("IFNG","GZMB","IL2","TBX21","IL7R","CD3D","CD3E"),
  checkpoint     = c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT"),
  ISG            = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","MX2","OAS1","OAS2","OASL","RSAD2","ISG20","IRF7","IRF9","STAT1","STAT2"),
  antigen_pres   = c("THRB","CD74","SLA-DRA","SLA-DRB1","CIITA","HLA-DRA","HLA-DRB1"),
  mhc2           = c("SLA-DRA","SLA-DRB1","CD74","CIITA")
)
genesets <- lapply(genesets, function(g) g[g %in% rownames(mat)])
genesets <- genesets[sapply(genesets, length)>0]

sc_long <- data.frame()
for (nm in names(genesets)) {
  gs <- genesets[[nm]]
  m  <- colMeans(mat[gs, , drop=FALSE])
  sc_long <- rbind(sc_long, data.frame(gene_set=nm, group=grp, score=as.numeric(m)))
}
sc_mean <- sc_long %>% group_by(gene_set, group) %>% summarise(score=mean(score), .groups="drop")

single <- c("CD163","IL1B","CXCL8","TNF","IFNG","CD274","PDCD1","THRB","ISG15","MX1","SLA-DRA")
single <- single[single %in% rownames(mat)]
sc_single <- data.frame()
for (g in single) {
  m <- as.numeric(mat[g, ])
  sc_single <- rbind(sc_single, data.frame(gene=g, group=grp, expr=as.numeric(m)))
}
sc_single_mean <- sc_single %>% group_by(gene, group) %>% summarise(expr=mean(expr), .groups="drop")

write_csv(sc_mean,        file.path(TAB, "scRNA_module_by_group.csv"))
write_csv(sc_single_mean, file.path(TAB, "scRNA_gene_by_group.csv"))
cat(sprintf("  scRNA module rows=%d; single-gene rows=%d\n", nrow(sc_mean), nrow(sc_single_mean)))

## ============ 湿实验读取 + 组名别名对齐 ============
WET <- file.path(BASE, "tables", "wet")
wet_vl   <- read_csv(file.path(WET,"wet_viral_load_by_tissue.csv"))
wet_cyto <- read_csv(file.path(WET,"wet_serum_cytokines.csv"))
wet_phe  <- read_csv(file.path(WET,"wet_phenotype.csv"))
wet_av   <- read_csv(file.path(WET,"wet_antiviral_dose.csv"))

alias <- c("Control"="Control","PRRSV"="PRRSV","APRR4"="APRR")
wet_vl$group   <- alias[as.character(wet_vl$group)]
wet_cyto$group <- alias[as.character(wet_cyto$group)]
wet_phe$group  <- alias[as.character(wet_phe$group)]
wet_groups <- c("Control","PRRSV","APRR")
gcol <- c("Control"="#4C72B0","PRRSV"="#C44E52","APRR"="#8172B2")

## ============ Fig A: 病毒趋向性 ============
safe <- function(expr, name) {
  r <- tryCatch({ force(expr); NULL },
                error=function(e){ cat(sprintf("  [WARN] %s skipped: %s\n", name, conditionMessage(e))); NULL })
  invisible(NULL)
}
pA <- tryCatch({
  vl_sum <- wet_vl %>% group_by(tissue) %>% summarise(mean=mean(value,na.rm=TRUE)) %>% arrange(desc(mean))
  top_t <- head(vl_sum$tissue, 10)
  pA1 <- ggplot(wet_vl %>% filter(tissue %in% top_t),
                aes(x=reorder(tissue,value,mean), y=value, fill=factor(group,levels=wet_groups))) +
    geom_boxplot(outlier.size=0.5) + coord_flip() + scale_fill_manual(values=gcol,name="Group") +
    labs(title="A1. PRRSV viral load by tissue (wet)", x="", y="Viral load") +
    theme_bw(base_size=10) + theme(legend.position="top")
  cd163 <- sc_single_mean %>% filter(gene=="CD163")
  pA2 <- ggplot(cd163, aes(x=group, y=expr, fill=group)) + geom_bar(stat="identity") +
    scale_fill_manual(values=gcol) + labs(title="A2. CD163 (PRRSV receptor) scRNA mean", x="", y="SCT expr") +
    theme_bw(base_size=10) + theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1))
  pA1/pA2
}, error=function(e){ cat("  [WARN] FigA skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(pA)) { ggsave(file.path(FIG,"FigA_viral_tropism.pdf"), pA, width=8, height=9);
                    ggsave(file.path(FIG,"FigA_viral_tropism.png"), pA, width=8, height=9, dpi=150); cat("  saved FigA\n") }

## ============ Fig B: 血清细胞因子 ↔ scRNA 炎症/IFN 基因 ============
pB <- tryCatch({
  cyto_mean <- wet_cyto %>% group_by(cytokine, group) %>% summarise(value=mean(value,na.rm=TRUE), .groups="drop")
  # 模式匹配 (避免特殊字符 β/α/γ 编码不一致)
  map_cyto <- list()
  for (cx in unique(cyto_mean$cytokine)) {
    if (grepl("IL-1", cx, ignore.case=TRUE) && !grepl("IL-10", cx, ignore.case=TRUE)) map_cyto[[cx]] <- "IL1B"
    else if (grepl("TNF", cx, ignore.case=TRUE)) map_cyto[[cx]] <- "TNF"
    else if (grepl("INF|IFN", cx, ignore.case=TRUE)) map_cyto[[cx]] <- "IFNG"
  }
  int_df <- data.frame()
  for (cx in names(map_cyto)) {
    g <- map_cyto[[cx]]
    w <- cyto_mean %>% filter(cytokine==cx)
    s <- sc_single_mean %>% filter(gene==g)
    if (nrow(w)==0 || nrow(s)==0) { cat("  [skip]", cx, "no match\n"); next }
    int_df <- rbind(int_df, data.frame(cytokine=cx,
      group=c(as.character(w$group), as.character(s$group)),
      src=c(rep("wet serum",nrow(w)), rep("scRNA gene",nrow(s))),
      val=c(w$value, s$expr)))
  }
  if (nrow(int_df)==0) stop("int_df empty")
  ggplot(int_df, aes(x=group, y=val, fill=src)) + geom_bar(stat="identity", position="dodge") +
    facet_wrap(~cytokine, scales="free_y") + scale_fill_brewer(palette="Set2") +
    labs(title="B. Serum cytokines (wet) vs scRNA gene expr", x="", y="value") +
    theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=30,hjust=1))
}, error=function(e){ cat("  [WARN] FigB skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(pB)) { ggsave(file.path(FIG,"FigB_serum_inflammation.pdf"), pB, width=9, height=4);
                    ggsave(file.path(FIG,"FigB_serum_inflammation.png"), pB, width=9, height=4, dpi=150); cat("  saved FigB\n") }

## ============ Fig C: 表型 ↔ 疾病严重度模块 ============
pC <- tryCatch({
  phe_metric <- wet_phe %>% group_by(metric, group) %>% summarise(value=mean(value,na.rm=TRUE), .groups="drop")
  sev <- sc_mean %>% filter(gene_set %in% c("inflam_myeloid","neutrophil","ISG")) %>%
    group_by(group) %>% summarise(score=mean(score), .groups="drop")
  pC1 <- ggplot(phe_metric, aes(x=metric, y=value, fill=factor(group,levels=wet_groups))) +
    geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=gcol,name="Group") +
    labs(title="C1. In vivo phenotype (wet)", x="", y="value") +
    theme_bw(base_size=10) + theme(legend.position="top", axis.text.x=element_text(angle=30,hjust=1))
  pC2 <- ggplot(sev, aes(x=group, y=score, fill=group)) + geom_bar(stat="identity") +
    scale_fill_manual(values=gcol) + labs(title="C2. scRNA disease-severity module", x="", y="SCT score") +
    theme_bw(base_size=10) + theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1))
  pC1/pC2
}, error=function(e){ cat("  [WARN] FigC skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(pC)) { ggsave(file.path(FIG,"FigC_phenotype_severity.pdf"), pC, width=9, height=9);
                    ggsave(file.path(FIG,"FigC_phenotype_severity.png"), pC, width=9, height=9, dpi=150); cat("  saved FigC\n") }

## ============ Fig D: 抗病毒剂量 ↔ ISG 模块 ============
pD <- tryCatch({
  dose_mean <- wet_av %>% group_by(assay, dose) %>% summarise(value=mean(value,na.rm=TRUE), .groups="drop")
  isg <- sc_mean %>% filter(gene_set=="ISG")
  pD1 <- ggplot(dose_mean, aes(x=dose, y=value, color=assay)) + geom_line() + geom_point() +
    scale_color_brewer(palette="Dark2") + labs(title="D1. APRR antiviral dose-response (wet)", x="Dose", y="Inhibition") +
    theme_bw(base_size=10)
  pD2 <- ggplot(isg, aes(x=group, y=score, fill=group)) + geom_bar(stat="identity") +
    scale_fill_manual(values=gcol) + labs(title="D2. scRNA ISG module by group", x="", y="SCT ISG score") +
    theme_bw(base_size=10) + theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1))
  pD1/pD2
}, error=function(e){ cat("  [WARN] FigD skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(pD)) { ggsave(file.path(FIG,"FigD_antiviral_ISG.pdf"), pD, width=9, height=9);
                    ggsave(file.path(FIG,"FigD_antiviral_ISG.png"), pD, width=9, height=9, dpi=150); cat("  saved FigD\n") }

## ============ 点对点整合表 ============
cmp <- data.frame()
for (cx in names(map_cyto)) {
  g <- map_cyto[[cx]]
  w <- cyto_mean %>% filter(cytokine==cx) %>% select(group, wet=value)
  s <- sc_single_mean %>% filter(gene==g) %>% select(group, scRNA=expr)
  m <- merge(w, s, by="group", all=TRUE)
  m$cytokine <- cx; m$scRNA_gene <- g
  cmp <- rbind(cmp, m)
}
write_csv(cmp, file.path(TAB, "integration_pointwise.csv"))
cat(sprintf("  integration pointwise rows=%d\n", nrow(cmp)))

cat(">>> DONE. Phase 2a integration figures in ", FIG, "\n")
