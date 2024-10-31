library(tidyverse)
library(ggdark)
library(gridExtra)

st_split_df <- read.csv("C:/Users/20972/Desktop/spatial_selection_rds_path.csv")
st_name <- st_split_df$stage

color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- color_df$color 
names(color_select) <- color_df$celltype

Mono_color <- c("#ad5f7d")
names(Mono_color) <- c("Mono")
color_select <- append(color_select, Mono_color)

st_angle_df <- read.csv("E:/CodePrograms/Immu/Immu_plot_0802/spatial_slice_angle.csv")
angle_list <- st_angle_df$angle
width_list <- st_angle_df$width
height_list <- st_angle_df$height
size_list <- st_angle_df$size 

num_count <- 0
for ( n in st_name) {
  num_count <- num_count + 1
  
  metadata_xy <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv")) %>% 
    select(UMAP_X, UMAP_Y)
  
  metadata_celltype <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0703/ST_meta/",n,"_metadata.csv")) %>% 
    select(celltype)
  
  angle <- angle_list[num_count]
  theta <- angle * pi / 180
  rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
  rotated_data <- as.data.frame(as.matrix(metadata_xy) %*% rotation_matrix)
  metadata <- as.data.frame(rotated_data)
  names(metadata) <- c("UMAP_X", "UMAP_Y")
  metadata$celltype <- metadata_celltype
  
  max_X <- max(metadata$UMAP_X)
  min_X <- min(metadata$UMAP_X)
  max_Y <- max(metadata$UMAP_Y)
  min_Y <- min(metadata$UMAP_Y)
  
  number <- 0
  displacement_NA_point_list <- list()
  displacement_celltype_point_list <- list()
  rare_celltype_point_list <- list()
  
  unique_celltype <- unlist(na.omit(unique(metadata$celltype)))
  for ( m in unique_celltype){
    
    number <- number + 1
    single_celltype_meta <- metadata
    
    single_celltype_meta$celltype[single_celltype_meta$celltype != m] <- NA
    
    if (number %% 5 == 1 ) {
      single_celltype_meta$UMAP_X <- single_celltype_meta$UMAP_X
    }else{
      single_celltype_meta$UMAP_X <- single_celltype_meta$UMAP_X+(number - 1)%%5*(max_X-min_X+200)
    }
    
    if ((number-1) %/% 5 == 0) {
      single_celltype_meta$UMAP_Y <- single_celltype_meta$UMAP_Y
    }else{
      single_celltype_meta$UMAP_Y <- single_celltype_meta$UMAP_Y + (number - 1)%/% 5*(max_Y-min_Y+200)
    }
    
    row.names(single_celltype_meta) <- paste(m, seq_len(nrow(single_celltype_meta)), sep = ".")
    
    na_data <- single_celltype_meta[is.na(single_celltype_meta$celltype), ]
    displacement_NA_point_list[[m]] <- na_data
    
    non_na_data <- single_celltype_meta[!is.na(single_celltype_meta$celltype), ]
    
    cell_ratio <- nrow(non_na_data)/nrow(na_data)
    if (cell_ratio > 0.002) {
      displacement_celltype_point_list[[m]] <- non_na_data
    }else{
      rare_celltype_point_list[[m]] <- non_na_data
    }
  }
  
  na_rbind <- bind_rows(displacement_NA_point_list)
  non_na_rbind <- bind_rows(displacement_celltype_point_list)
  rare_rebind <- bind_rows(rare_celltype_point_list)
  
  na_rbind$celltype <- na_rbind$celltype$celltype
  non_na_rbind$celltype <- non_na_rbind$celltype$celltype
  rare_rebind$celltype <- rare_rebind$celltype$celltype
  
  p_st <- ggplot() +
    coord_fixed() +
    geom_point(data = na_rbind, aes(x = UMAP_X, y = UMAP_Y), size = 0.1, colour = "#484848") +
    geom_point(data = non_na_rbind, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 0.5) +
    coord_fixed() +
    labs(color = NULL) +
    dark_theme_void() +
    scale_color_manual(values = color_select) +
    theme(
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(color = "white", size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis.title = element_text(color = "white", size = 20),
      # axis.text = element_text(color = "white", size = 15),
      # axis.ticks = element_line(color = "white"),
      # axis.line = element_line(color = "white"),
      #legend.position = "none",
      plot.background = element_rect(fill = "black"),
      axis.title = element_blank(),
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
  
  if (n != "P5_nor_C02926B1") {
    p_st <- p_st + geom_point(data = rare_rebind, aes(x = UMAP_X, y = UMAP_Y, colour = celltype), size = 3)A
  }
  
  ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0802/Immu_split_one_plot/",n,"_split_OnePlot.png"), 
         p_st, width = width_list[num_count], height = height_list[num_count], limitsize = FALSE)
}