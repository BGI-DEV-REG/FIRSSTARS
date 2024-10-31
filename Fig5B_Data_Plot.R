library(tidyverse)
library(viridis)
library(ggdark)

Immu_df <- read.csv("E:/CodePrograms/Immu/Immu_0602/monocel/set_seed/All_Immu_Seed_5.0_3.0_obs_umap.csv")
TPI_factor <- c("Nor", "6H", "1D", "2D", "3D", "5D", "6D", "10D", "19D")
Immu_df$TimePostInjury <- factor(Immu_df$TimePostInjury, TPI_factor)

TPI_color <- c("#00FF00","#A020F0","#F2AD00","#26c6da","#f9320c","#FFF300")

for (age in c("P5","Adult")) {
  
  Immu_df_filter <- Immu_df %>% 
    filter(Age == age)
  
  p0 <- ggplot(Immu_df_filter, aes(x = UMAP1 , y = UMAP2, color = TimePostInjury)) +
    geom_point(size = 0.3) +
    labs(color = NULL) +
    dark_theme_minimal() +
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
    ) +
    scale_color_manual(values = TPI_color)
  
  ggsave(paste0(age,"_Immu_TPI.pdf"), p0, w = 9, h = 8)
  
}