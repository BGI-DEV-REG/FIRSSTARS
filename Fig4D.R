library(tidyverse)
library(ggdark)

all_cell <- read.csv("E:/CodePrograms/Immu/Immu_0703/all_cell_meta_0703.csv.gz")
Immu_cell <- read.csv("E:/CodePrograms/Immu/Immu_0703/Immu_meta_0703.csv")
setwd("E:/CodePrograms/Immu/Immu_plot_0802/fig4 pieplot/")

Immu_color_df <- read.csv("E:/CodePrograms/Immu/Immu_0703/color_select_last.csv")
color_select <- Immu_color_df$color
names(color_select) <- Immu_color_df$celltype

all_cell$celltype = ifelse(!(all_cell$celltype %in% unique(Immu_cell$celltype)), "other_celltype", all_cell$celltype)
Immu_celltype <- c("other_celltype","Macr_remodeling","Macr_neovascularization","Mono_noclassical","Mono_intermediate",
                   "Mono_classical","Macr_resident","Neutrophils","moDC","mDC","cDC","pDC","Langerhans",
                   "Myeloid_naive","NK_ILC","T_cells","Mast_cells","B_cells")
Immu_celltype_de <- c("Macr_remodeling","Macr_neovascularization","Mono_noclassical","Mono_intermediate","Mono_classical",
                      "Macr_resident","Neutrophils","moDC","mDC","cDC","pDC","Langerhans",
                      "Myeloid_naive","NK_ILC","T_cells","Mast_cells","B_cells")
all_cell$celltype <- factor(all_cell$celltype, rev(Immu_celltype))

for ( age in c("E16.5","E17.5","E18.5","E","P5","Adult")) {
  
  if (age == "E") {
    len_all <- sum(all_cell$Age %in% c("E16.5","E17.5","E18.5"))
    select_cell <- all_cell %>%
      filter( Age %in% c("E16.5","E17.5","E18.5") ) %>%
      mutate(Age = "E") %>% 
      group_by(celltype) %>% 
      mutate(celltype_prop = length(celltype) / len_all) %>% 
      select(Age, celltype, celltype_prop) %>% 
      unique() %>% 
      arrange(celltype)
  }else{
    len_all <- sum(all_cell$Age == age)
    select_cell <- all_cell %>% 
      filter( Age == age ) %>% 
      group_by(celltype) %>% 
      mutate(celltype_prop = length(celltype) / len_all) %>% 
      select(Age, celltype, celltype_prop) %>% 
      unique() %>% 
      arrange(celltype)
  }
  
  p_pie <- ggplot()+
    coord_fixed() +
    dark_theme_void() +
    scale_fill_manual(values = color_select)+
    geom_arc_bar(data = select_cell,
                 stat = "pie",
                 aes(x0=0,y0=0,r0=0,r=2,
                     amount=celltype_prop,
                     fill=celltype,
                     explode=c(rep(0, nrow(select_cell) -1), 0.1))
                 ,color = NA)
  
  print(p_pie)
  ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0802/fig4 pieplot/pie_plot/Pie_", age, ".pdf"), p_pie)
  
  
  if (age == "E") {
    len_all_Immu <- sum(all_cell$Age %in% c("E16.5","E17.5","E18.5") & all_cell$celltype %in% Immu_celltype_de)
    select_cell_Immu <- all_cell %>%
      filter( celltype %in% Immu_celltype_de) %>% 
      filter( Age %in% c("E16.5","E17.5","E18.5") ) %>%
      mutate(Age = "E") %>% 
      group_by(celltype) %>% 
      mutate(celltype_prop = length(celltype) / len_all_Immu) %>% 
      select(Age, celltype, celltype_prop) %>% 
      unique() %>% 
      arrange(celltype)
  }else{
    len_all_Immu <- sum(all_cell$Age == age & all_cell$celltype %in% Immu_celltype_de)
    select_cell_Immu <- all_cell %>% 
      filter( celltype %in% Immu_celltype_de) %>% 
      filter( Age == age ) %>% 
      group_by(celltype) %>% 
      mutate(celltype_prop = length(celltype) / len_all_Immu) %>% 
      select(Age, celltype, celltype_prop) %>% 
      unique() %>% 
      arrange(celltype)
  }
  
  p_pie_Immu <- ggplot()+
    coord_fixed() +
    dark_theme_void() +
    scale_fill_manual(values = color_select)+
    geom_arc_bar(data = select_cell_Immu,
                 stat = "pie",
                 aes(x0=0,y0=0,r0=1,r=2,
                     amount=celltype_prop,
                     fill=celltype), color = NA)
  
  print(p_pie_Immu)
  ggsave(paste0("E:/CodePrograms/Immu/Immu_plot_0802/fig4 pieplot/pie_plot/PieImmu_", age, ".pdf"), p_pie_Immu)
}