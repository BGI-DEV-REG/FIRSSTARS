library(tidyverse)

df <- read.csv("E:/CodePrograms/Immu/Immu_0703/Immu_meta_0703.csv")

Age_color <- c("#00FF00","#A020F0","#F2AD00","#26c6da","#f9320c","#FFF300")

p_Age <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = Age)) +
  geom_point(size = 0.3) +
  labs(color = NULL) +
  dark_theme_minimal() +
  scale_color_manual(values = Age_color) +
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
ggsave("C:/Users/20972/Desktop/Fig4_plot/Immu_Age_umap_0704.pdf", plot = p_Age, width = 25, height = 20)