library(tidyverse)
library(ggdark)

st_split_df <- read.csv("C:/Users/20972/Desktop/spatial_selection_rds_path.csv")
st_name <- st_split_df$stage

color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- color_df$color 
names(color_select) <- color_df$celltype

E_st_name <- grep("^E", st_name, value = T) 

for ( n in st_name) {
  
  metadata <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv"))
  na_data <- metadata[is.na(metadata$celltype), ]
  non_na_data <- metadata[!is.na(metadata$celltype), ]
  
  na_num <- nrow(na_data)
  Immu_num <- nrow(non_na_data)
  cell_ratio <- Immu_num / na_num
  
  if (cell_ratio > 0.02 & !(n %in% E_st_name)){
    p_st <- ggplot() +
      dark_theme_minimal() +
      coord_fixed() +
      geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 0.5, colour = "#484848") +
      geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 0.5) +
      coord_fixed() +
      labs(color = NULL) +
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
  }else if(cell_ratio <= 0.02 & !(n %in% E_st_name)){
    p_st <- ggplot() +
      dark_theme_minimal() +
      coord_fixed() +
      geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 0.5, colour = "#484848") +
      geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 1.5) +
      coord_fixed() +
      labs(color = NULL) +
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
  }else if(cell_ratio > 0.02 & n %in% E_st_name) {
    p_st <- ggplot() +
      dark_theme_minimal() +
      coord_fixed() +
      geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 1.5, colour = "#484848") +
      geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 2) +
      coord_fixed() +
      labs(color = NULL) +
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
  }else if(cell_ratio <= 0.02 & n %in% E_st_name){
    p_st <- ggplot() +
      dark_theme_minimal() +
      coord_fixed() +
      geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 1.5, colour = "#484848") +
      geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 3) +
      coord_fixed() +
      labs(color = NULL) +
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
  }
  
  
  if (n %in% E_st_name) {
    ggsave(paste0("./Fig4_plot/st_new/",n,".png"),p_st, w = 12, h = 16)
  }else{
    ggsave(paste0("./Fig4_plot/st_new/",n,".png"),p_st, w = 16, h = 20)
  }
}