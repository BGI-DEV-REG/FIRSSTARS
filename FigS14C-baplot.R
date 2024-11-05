library(ggplot2)
library(tidyverse)

colorlist = read.csv("Fig2_color_seq_20240521.csv")
colors = colorlist$colors
names(colors) = colorlist$X

meta_all = read.csv("meta_all.csv.gz",row.names=1)
AgeTPI_num <- as.data.frame(table(meta_all$celltype_240520, meta_all$Age_TimePostInjury))
names(AgeTPI_num) <- c("celltype", "Age_TPI", "number")
#计算每种细胞在某Age_TPI的比例
AgeTPI_ratio <- AgeTPI_num %>%
  group_by(Age_TPI) %>%
  mutate(prop_number = number / sum(number))
data = AgeTPI_ratio

AgeTPI_ratio = data
unique_age_tpi <- unique(AgeTPI_ratio$Age_TPI)

celltype  <-  c("FB_DP_Adult", "FB_DS_Adult", "FB_adipogenic_Adult", "FB_basement_resident_Adult",  "FB_papi_Adult", "FB_reti_Adult","FB_muscle_Adult", "FB_fascia_2_Adult", "FB_fascia_1_Adult","FB_basement_activated_Adult", "FB_inflammatory_Adult", "FB_exosome_Adult", "FB_apoptotic_Adult", "FB_oxidative_stress_Adult", "FB_myo_Adult","FB_DP_EP", "FB_DS_EP", "FB_adipogenic_EP", "FB_basement_resident_EP",  "FB_papi_EP", "FB_reti_EP","FB_muscle_EP", "FB_fascia_EP","FB_basement_activated_EP", "FB_inflammatory_EP", "FB_exosome_EP", "FB_apoptotic_EP", "FB_oxidative_stress_EP", "FB_myo_EP")
data$celltype <- factor(data$celltype, celltype)
data$Age_TPI <- factor(data$Age_TPI, c("E16.5_Nor", "E16.5_6H", "E16.5_1D", "E16.5_2D", "E16.5_3D", "E16.5_Born",
             "E17.5_Nor", "E17.5_6H", "E17.5_1D", "E17.5_2D", "E17.5_3D", "E17.5_Born",
             "E18.5_Nor", "E18.5_6H", "E18.5_1D", "E18.5_2D", "E18.5_3D", "E18.5_Born",
             "P0_Nor", "P5_Nor", "P5_6H", "P5_1D", "P5_2D", "P5_3D", "P5_6D",
             "Adult_Nor", "Adult_1D", "Adult_2D", "Adult_5D", "Adult_10D", "Adult_19D"))
p2 <-ggplot(data, aes(x = Age_TPI, y = prop_number, fill = celltype)) +
  geom_col(width = 0.8) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "white"),
    axis.text.y = element_text(size = 12, color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) + 
  guides(fill = FALSE) +
  scale_fill_manual(values = colors)
  #guides(fill = guide_legend(ncol = 1))
  
ggsave("AgeTPI_ratio_all_240523_no_tab.png", p2, height =5 , width = 20)
