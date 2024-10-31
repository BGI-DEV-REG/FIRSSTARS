library(viridis)
library(Seurat)
library(tidyverse)

Immu_subset <- readRDS("E:/CodePrograms/Immu/Immu_0524/Immu_scvi_0526.rds")

genes <- c("Lyz2", "Ly6c", "Ccr2", "Lyve1")

for (gene in genes) {
  p <- FeaturePlot(Immu_subset, gene, reduction = "scvi_umap", pt.size = 0.1, cols = viridis(100)) + 
    theme_minimal() +
    ggtitle(gene) +
    theme(
      plot.title = element_text(color = "white",size = 40, hjust = 0.5),
      legend.key.size = unit(1, "cm"),
      plot.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white", size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  ggsave(paste0("Immu_", gene, "_0526.pdf"), p, width = 18, height = 16)
}