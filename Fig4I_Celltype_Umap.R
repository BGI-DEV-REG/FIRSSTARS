library(tidyverse)
library(ggdark)

st_split_df <- read.csv("C:/Users/20972/Desktop/spatial_selection_rds_path.csv")
st_name <- st_split_df$stage

color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- color_df$color 
names(color_select) <- color_df$celltype

Mono_color <- c("#ad5f7d")
names(Mono_color) <- c("cMono")
color_select <- append(color_select, Mono_color)


for ( n in st_name) {
  
  metadata <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv"))
  
  single_celltype_meta <- metadata
  single_celltype_meta$celltype[!(single_celltype_meta$celltype %in% c('Mono_classical'))] <- NA
  single_celltype_meta$celltype[single_celltype_meta$celltype %in% c('Mono_classical')] <- "cMono"
  
  na_data <- single_celltype_meta[is.na(single_celltype_meta$celltype), ]
  non_na_data <- single_celltype_meta[!is.na(single_celltype_meta$celltype), ]
  
  p_st <- ggplot() +
    coord_fixed() +
    geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 0.1, colour = "#303030") +
    geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 0.5) +
    coord_fixed() +
    labs(color = NULL) +
    dark_theme_minimal() +
    scale_color_manual(values = color_select) + #Mono
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
  
  ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0725/Fig4_3_mono/",n,"_3_mono_1color",".pdf"),p_st, w = 16, h = 20)
  
}


FB_color_df <- read.csv("E:/CodePrograms/Immu/Immu_plot_0718/Fig2_color_seq_20240521.csv")
FB_color <- FB_color_df$colors
names(FB_color) <- FB_color_df$X

Myo_color <- c("#cddc39")
names(Myo_color) <- c("Myo_FB")
FB_color <- append(FB_color, Myo_color)
FB_ST_name <- list.files("C:/Users/20972/Desktop/ST_FB_detail_anno/")


for (index in seq_along(st_name)) {
  
  n <- st_name[index]
  i <- FB_ST_name[index]
  
  metadata <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv"))
  
  if (str_detect(i,"Adult")) {
    FB_meta <- read.csv(paste0("C:/Users/20972/Desktop/ST_FB_detail_anno/", i, "/FB_pred_update240524.csv.gz"))
  }else{
    FB_meta <- read.csv(paste0("C:/Users/20972/Desktop/ST_FB_detail_anno/", i, "/pred.csv.gz"))
  }
  
  metadata$celltype <- ifelse(metadata$X %in% FB_meta$X, FB_meta$pred_celltype[match(metadata$X, FB_meta$X)], NA)
  
  metadata$celltype[!(metadata$celltype %in% c("FB_fascia_1_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult", "FB_myo_Adult", 
                                               "FB_fascia_EP", "FB_inflammatory_EP", "FB_oxidative_stress_EP"))] <- NA
  metadata$celltype[metadata$celltype %in% c("FB_fascia_1_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult", 
                                             "FB_myo_Adult", "FB_fascia_EP", "FB_inflammatory_EP", "FB_oxidative_stress_EP")] <- "Myo_FB"
  
  na_data <- metadata[is.na(metadata$celltype), ]
  non_na_data <- metadata[!is.na(metadata$celltype), ]
  
  p_st <- ggplot() +
    coord_fixed() +
    geom_point(data = na_data, aes(x = UMAP_X, y = UMAP_Y), size = 0.1, colour = "#303030") +
    geom_point(data = non_na_data, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 0.5) +
    coord_fixed() +
    labs(color = NULL) +
    dark_theme_minimal() +
    scale_color_manual(values = FB_color) +
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
  
  ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0725/Fig4_only_Myo_FB/",n,"_only_Myo_FB.pdf"),p_st, w = 30, h = 30)
  
}