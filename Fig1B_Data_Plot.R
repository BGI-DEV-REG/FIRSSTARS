library(tidyverse)

df <- read.csv("E:/CodePrograms/All_Rat/All_RAT_0910/allcell_last_obs_umap.csv")

color_celltype <- read.csv("E:/CodePrograms/plot/plot_0513/Fig1_color_seq_20240512.csv")
all_celltype <- color_celltype$celltype
all_color <- color_celltype$color
names(all_color) <- all_celltype

p0 <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = as.factor(celltype))) +
  geom_point(size = 0.05) +
  labs(color = NULL) +
  theme_minimal() + 
  scale_color_manual(values = all_color) +
  theme(
    legend.position = "left",
    legend.key.size = unit(1, "cm"),
    plot.background = element_rect(fill = "black") ,
    legend.text = element_text(color = "white",size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
ggsave("fig1_allcell_umap.png", plot = p0, width = 23, height = 20)