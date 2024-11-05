library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

######## arguments ##########
parser = argparse::ArgumentParser(description="celltype plot")
parser$add_argument('-S','--stage', help='section stage')
args = parser$parse_args()

stage = args$stage
outpath=paste0("./",stage)
if(!(dir.exists(outpath))){
    dir.create(outpath)
}

setwd(outpath)

######## 从全部坐标提取空间组坐标 ########
coord = read.csv("All_coord_adjusted_20240519.csv.gz",row.names=1)

table(coord$stage)
head(coord)

coord_sub = coord[coord$stage == stage,]

######## 添加注释 ########
df = read.csv('spatial_selection_rds_path_fig2_meta.txt',sep='\t')
df = df[c(1:6,25:30),]
coord_sub[,"celltype"]="others"
meta_path = df[df$stage == stage,"anno"]
meta = read.csv(meta_path,row.names=1)
rownames(meta) = paste0(stage,".",rownames(meta))
print(dim(meta))
print(dim(coord_sub))
coord_sub[rownames(meta),"celltype"]=meta[,"pred_celltype"]
print(table(coord_sub$celltype))

######## 设置颜色 ########
df = read.csv("Fig2_color_seq_20240521.csv")
colors = df$colors
names(colors) = df$X

######## 可视化 ########
for( ct in c("FB_fascia_1_Adult","FB_oxidative_stress_Adult","FB_inflammatory_Adult","FB_myo_Adult")){
print(ct)
p <- ggplot(coord_sub, aes(x = x_adjust, y = y_adjust)) +
  geom_point(data = coord_sub, color = "#9E9E9E",size = 0.5) +
  geom_point(data = coord_sub[coord_sub$celltype == ct,], color = colors[ct], size = 1) +
  coord_fixed() +
  labs(title = ct) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    text = element_text(color = "white"),
    plot.title = element_text(hjust = 0.5, color = "white", size = 18)
  ) +
  guides(size = FALSE)
ggsave(plot = p, filename = paste0(stage,"_",ct,"_Spatial.png"), height = 10, width = 10,limitsize = FALSE)
}

