library(ggplot2)
library(tidyverse)
meta_all = read.csv("meta_all.csv.gz",row.names=1)

AgeTPI_num <- as.data.frame(table(meta_all$celltype_240520, meta_all$Age_TimePostInjury))
names(AgeTPI_num) <- c("celltype", "Age_TPI", "number")
AgeTPI_ratio <- AgeTPI_num %>%
  group_by(Age_TPI) %>%
  mutate(prop_number = number / sum(number))
AgeTPI_ratio$Age_TPI <- factor(AgeTPI_ratio$Age_TPI, c("E16.5_Nor", "E16.5_6H", "E16.5_1D", "E16.5_2D", "E16.5_3D", "E16.5_Born",
             "E17.5_Nor", "E17.5_6H", "E17.5_1D", "E17.5_2D", "E17.5_3D", "E17.5_Born",
             "E18.5_Nor", "E18.5_6H", "E18.5_1D", "E18.5_2D", "E18.5_3D", "E18.5_Born",
             "P0_Nor", "P5_Nor", "P5_6H", "P5_1D", "P5_2D", "P5_3D", "P5_6D",
             "Adult_Nor", "Adult_1D", "Adult_2D", "Adult_5D", "Adult_10D", "Adult_19D"))
  
df = read.csv("Fig2_color_seq_20240521.csv")
colors = df$colors
names(colors) = df$X

########### Adult #####################
da = AgeTPI_ratio[grepl("^Adult",AgeTPI_ratio$Age_TPI),]
da = da[da$celltype %in% c("FB_myo_Adult" ,"FB_fascia_1_Adult", "FB_oxidative_stress_Adult", "FB_inflammatory_Adult"),]
da$celltype = factor(da$celltype,levels=c("FB_fascia_1_Adult","FB_inflammatory_Adult","FB_oxidative_stress_Adult","FB_myo_Adult"))

p <- ggplot(da,aes(x = Age_TPI , y = prop_number, colour = celltype,shape = celltype, group = celltype)) + 
  geom_line(linewidth=1.5)+  
  geom_point(size=3)+ 
  scale_shape_manual(values = c(10,5,2,7)) +
  scale_color_manual(values= colors)+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, color = "white"),
    axis.text.x = element_text(size = 8, color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    )
ggsave(plot=p,"Adult_fascia_myo_time_line_plot.pdf",w=8,h=4)

########### EP #####################
for (stage in c("E16.5","E17.5","E18.5","P5")){
da = AgeTPI_ratio[grepl(stage,AgeTPI_ratio$Age_TPI),]
da = da[da$celltype %in% c("FB_fascia_EP","FB_inflammatory_EP", "FB_oxidative_stress_EP"),]

p <- ggplot(da,aes(x = Age_TPI , y = prop_number, colour = celltype,shape = celltype, group = celltype)) + 
  geom_line(linewidth=1.5)+  
  geom_point(size=3)+ 
  scale_shape_manual(values = c(10,5,2)) +
  scale_color_manual(values= colors)+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, color = "white"),
    axis.text.x = element_text(size = 8, color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    )
ggsave(plot=p,paste0(stage,"_fascia_myo_time_line_plot.pdf"),w=7,h=4)
}


