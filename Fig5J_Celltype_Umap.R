library(tidyverse)
library(ggdark)

st_split_df <- read.csv("C:/Users/20972/Desktop/spatial_selection_rds_path.csv")
st_name <- st_split_df$stage

color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- color_df$color 
names(color_select) <- color_df$celltype

# 免疫细胞分群显示
for ( n in st_name) {
  
  metadata <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv"))
  
  for ( m in na.omit(unique(metadata$celltype))) {
    single_celltype_meta <- metadata
    single_celltype_meta$celltype[single_celltype_meta$celltype != m] <- NA
    
    na_data <- single_celltype_meta[is.na(single_celltype_meta$celltype), ]
    non_na_data <- single_celltype_meta[!is.na(single_celltype_meta$celltype), ]
    
    p_st <- ggplot() +
      coord_fixed() +
      geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 0.1, colour = "#303030") +
      geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 0.5) +
      coord_fixed() +
      labs(color = NULL) +
      dark_theme_minimal() +
      scale_color_manual(values = color_select) +
      theme(
        legend.key.size = unit(1, "cm"),
        plot.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(color = "white", size = 20),
        axis.text = element_text(color = "white", size = 15),
        axis.ticks = element_line(color = "white"),
        axis.line = element_line(color = "white")
      ) +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2)))
    
    ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0725/Fig4_split_Immu_celltype_plot/",n,"_",m,".pdf"),p_st, w = 16, h = 20)
    
  }
}