library(tidyverse)

df <- read.csv("E:/CodePrograms/Immu/Immu_0703/Immu_meta_0703.csv")

df$celltype <- factor(df$celltype, rev(Immu_celltype))
df$Age <- factor(df$Age, Age_factor)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = celltype)) +
  geom_point(size = 0.3) +
  labs(color = NULL) +
  dark_theme_minimal() + 
  scale_color_manual(values = color_select) +
  theme(
    legend.position = "left",  
    legend.key.size = unit(1, "cm"), 
    plot.background = element_rect(fill = "black") ,
    legend.text = element_text(color = "white",size = 20),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),       
    axis.text = element_blank(),         
    axis.ticks = element_blank()
  )+
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 6)))
ggsave("C:/Users/20972/Desktop/Fig4_plot/Immu_celltype_umap_0704.pdf", plot = p0, width = 25, height = 20)