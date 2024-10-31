library(tidyverse)
library(ggdark)

color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- color_df$color 
names(color_select) <- color_df$celltype

df <- read.csv("E:/CodePrograms/Immu/Immu_0703/Immu_meta_0703.csv")
Mono_type <- c("Macr_remodeling","Macr_neovascularization","Mono_noclassical",
               "Mono_intermediate","Mono_classical","Macr_resident")
TPI_factor <- c("Nor", "6H", "1D", "2D", "3D", "Born", "5D", "6D", "10D", "19D")

for (age in c("P5","Adult")) {
  Mono_df <- df %>% 
    filter(Age == age & celltype %in% Mono_type)
  
  TPI_num <- as.data.frame(table(Mono_df$celltype, Mono_df$TimePostInjury))
  names(TPI_num) <- c("celltype", "TPI", "number")
  TPI_num$Age <- sub("_.*","",TPI_num$TPI)
  #计算每种细胞在某Age_TPI的比例
  TPI_ratio <- TPI_num %>%
    group_by(TPI) %>%
    mutate(prop_number = number / sum(number))
  
  TPI_ratio$TPI <- factor(TPI_ratio$TPI, TPI_factor)
  
  p_line <- ggplot(TPI_ratio, aes(x = TPI, y = prop_number, fill = celltype, group = celltype, color = celltype)) +
    geom_smooth(method = "loess", se = FALSE, size = 0.5) +
    ggtitle(age) +
    dark_theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "white"),
      axis.text.y = element_text(size = 12, color = "white"),
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      axis.line.x=element_line(size=0.2),
      axis.line.y=element_line(size=0.2)
    ) + 
    guides(color = FALSE) +
    scale_color_manual(values = color_select)
  
  ggsave(paste0(age,"_mono_line_plot.pdf"), p_line, w = 6, h = 6)
  
}