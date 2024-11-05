library(ggplot2)
library(tidyverse)
library(Seurat)
library(viridis)

setwd("../Fascia")
obj <- readRDS("fascia_Reduc.rds")

##################### UMAP ##########################################
color_df = read.csv("Fig2_color_seq_20240521.csv", row.names=1)
colors = color_df$colors
names(colors) = rownames(color_df)
p <- DimPlot(obj, reduction='fascia', group.by="celltype", cols=colors, pt.size=0.5) + theme(panel.background = element_rect(fill = "black"),
plot.background = element_rect(fill = "black"),
# text = element_text(color = "white", size=1),
plot.title = element_text(color = "white"),
panel.grid=element_blank(),#strip.text = element_blank(),
panel.border = element_rect(fill=NA,color="white", size=1, linetype="solid")) + NoLegend()
ggsave(plot=p,"fascia_umap_celltype_black.pdf")

##################### Age 着色 
colors = inferno(6, alpha = 1, begin = 0.2, end = 0.8, direction = 1)
names(colors) = c("E16.5","E17.5","E18.5","P0","P5","Adult")
barplot(1:6, col = colors)
dev.off()
obj$Age= factor(obj$Age,levels=c("E16.5","E17.5","E18.5","P0","P5","Adult"))
p <- DimPlot(obj, reduction='fascia', group.by="Age", cols=colors, pt.size=0.2) + theme(panel.background = element_rect(fill = "black"),
plot.background = element_rect(fill = "black"),
text = element_text(color = "white", size=12),
plot.title = element_text(color = "white"),
panel.grid=element_blank(),#strip.text = element_blank(),
panel.border = element_rect(fill=NA,color="white", size=1, linetype="solid"))
ggsave("fascia_umap_age_black.pdf")

##################### Dotplot ##########################################
obj$celltype= factor(obj$celltype,levels=c("FB_fascia_EP","FB_fascia_1_Adult",  "FB_fascia_2_Adult"))
p <- DotPlot(obj, features=c("Mest","Dlk1", "Cd55","Plac8", "Pi16", "Procr", "Irf4","Smarcd2","Clec2g","Nova1"), group.by="celltype", dot.scale=20) +theme_classic() + 
theme( 
axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "white"), 
axis.text.y = element_text(size = 12, color = "white"), 
plot.background = element_rect(fill = "black"), 
panel.background = element_rect(fill = "black"), 
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), 
)+scale_color_viridis()+coord_flip()
ggsave("dot_markers_black.png",w=6)


##################### Lineplot ##########################################
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
  
da = AgeTPI_ratio[AgeTPI_ratio$celltype == "FB_fascia_1_Adult",]
da = da[grepl("^Adult",da$Age_TPI),]

p <- ggplot(da,aes(x = Age_TPI , y = prop_number, colour = celltype,shape = celltype, group = 1)) + 
  geom_line(linewidth=1.5)+  
  geom_point(size=3)+ 
  scale_shape_manual(values = 10) +
  scale_color_manual(values= "#db2400")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, color = "white"),
    axis.text.x = element_text(size = 8, color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    )

ggsave(plot=p,"FB_fascia_1_Adult_time_line_plot.png",w=8,h=2)

