## Fig 4 B: T cell 5 亚群 × 3 组比例堆叠条
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(cowplot); library(grid); library(gridExtra); library(patchwork)
})

OUT <- r"(figures_final/Fig4)"
obj <- readRDS(r"(playbook_outputs/T_cell_subcluster/T_cell_subcluster_3grp.rds)")

sub_levels <- c("0","1","2","3","4")
sub_cols   <- c("0"="#E41A1C", "1"="#377EB8", "2"="#4DAF4A", "3"="#984EA3", "4"="#FF7F00")
short_name <- c("0"="T_CD8_eff", "1"="T_CD4", "2"="T_prolif",
                "3"="T_JCHAIN+", "4"="T_innate/MAIT")
obj$subcluster <- factor(obj$subcluster, levels=sub_levels)
obj$group <- factor(obj$group, levels=c("Control","PRRSV","APRR"))

# 算每个 group 内的 subcluster 比例
prop_df <- obj@meta.data %>%
  group_by(group, subcluster) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(group) %>%
  mutate(total=sum(n), prop=n/total) %>%
  ungroup()

# 堆叠条
p <- ggplot(prop_df, aes(x=group, y=prop, fill=subcluster)) +
  geom_bar(stat="identity", position="stack", color="white", linewidth=0.4) +
  scale_fill_manual(values=sub_cols,
                    labels=short_name[sub_levels],
                    name="T cell subset") +
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=expansion(mult=c(0,0.05))) +
  geom_text(aes(label=ifelse(prop>0.04, paste0(round(prop*100, 1), "%"), "")),
            position=position_stack(vjust=0.5), size=3.6, color="white", fontface="bold") +
  labs(title="Fig 4B · T cell subset composition across groups",
       x=NULL, y="Proportion") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14, margin=margin(b=10)),
        legend.position="right",
        axis.text.x=element_text(face="bold", size=12),
        panel.grid.major.x=element_blank())

for (ext in c("png","pdf","svg")) {
  f <- file.path(OUT, paste0("Fig4B_tc_subcluster_prop.", ext))
  if (ext=="svg") { svg(f, width=8, height=5.5); print(p); dev.off() }
  else if (ext=="pdf") { pdf(f, width=8, height=5.5); print(p); dev.off() }
  else { png(f, width=8*120, height=5.5*120, res=120, bg="white"); print(p); dev.off() }
  cat("  [saved]", basename(f), "\n")
}
print(prop_df, n=Inf)
