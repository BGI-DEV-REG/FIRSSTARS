# Co-localization analysis visulization
## take FigS29A as an instance

############################ Bubble Plot #############################################
library(ggplot2)
library(dplyr)
library(viridis)
path <- "../19.colocation/"
setwd(path)

df_list = list()
for (st in c("Adult_1D_C01830C1D1","Adult_2D_C01830E4F4","Adult_5D_C01830E6F6","Adult_10D_A02988C5","Adult_19D_B03424F4")){
    csvpath = paste0(path,st,"_20_colocation.csv")
    df = read.csv(csvpath,row.names=1)
    fn = sub(".csv", "", basename(csvpath))
    prefix = sub(strsplit(st,"_")[[1]][3],"",st)
    df$From = paste0(prefix,df$From)
    df$Log2FC[df$Log2FC <= -2] = -2
    df$neglogP = -log10(df$Pvalue)
    df$neglogP[df$neglogP >= 20] = 20
    df_list[[st]]=df
}

da = do.call(rbind,df_list)

ct_select = c("Adult_1D_FB_fascia_1_Adult","Adult_1D_FB_inflammatory_Adult","Adult_2D_FB_inflammatory_Adult","Adult_2D_FB_oxidative_stress_Adult","Adult_5D_FB_oxidative_stress_Adult","Adult_5D_FB_myo_Adult","Adult_10D_FB_myo_Adult","Adult_19D_FB_myo_Adult")
subdf = da[da$From %in% ct_select,]
subdf$From = factor(subdf$From,levels=rev(ct_select))

white_theme <- function() {
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black",color="black"),
        legend.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        axis.line = element_line( linewidth = 0.5,color="white"),
        axis.text = element_text( size = 15,color="white"),
        axis.title = element_text( size = 14,color="white"),
        legend.title = element_text( size = 12,color="white"),
        legend.text = element_text( size = 12,color="white"),
        legend.position = "right",
        #legend.key = element_rect(fill = "black")
       )
}

p <- ggplot(subdf, aes(x=To, y=From, color=Log2FC, size=neglogP)) + geom_point()+scale_color_viridis(option="B",begin=0.2)+theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=5),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size=6)
     )+white_theme()
ggsave(plot=p,paste0("./fig3_supp/colocation/All_stage_FB_fascia_myo_tra_Colocation_legend.pdf"),h=3.5,w=5)

############################ Contour Plot #############################################
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

df = read.csv("../16.Table/spatial_slice_angle.csv")

stages = c("Adult_2d_C01830E4F4","Adult_5d1_C01830E6F6")
celltypes_1 = list("Adult_2d_C01830E4F4" = c("FB_inflammatory_Adult"),"Adult_5d1_C01830E6F6"=c("FB_oxidative_stress_Adult"))
celltypes_2 = list("Adult_2d_C01830E4F4" = c("FB_oxidative_stress_Adult"),"Adult_5d1_C01830E6F6"=c("FB_myo_Adult"))

for (stage in stages){
### 提取空间坐标
df1 = read.csv('./00.data/spatial_selection_rds_path.txt',sep='\t')
objpath = df1[df1$stage == stage,"path"]
print(objpath)
obj = readRDS(objpath)
obj$coor_y = obj@images[[names(obj@images)]]@coordinates$row
obj$coor_x = obj@images[[names(obj@images)]]@coordinates$col
coord=data.frame(x = obj$coor_x, y = obj$coor_y)

angle <- df[df$stage == stage,"angle"]+90
angle <- -angle * (pi / 180)
rotation_matrix <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2)
rotated_coords <- as.data.frame(t(rotation_matrix %*% t(as.matrix(coord))))
colnames(rotated_coords) <- c("x", "y")
coord=rotated_coords

name = df[df$stage == stage,"name"]
meta = read.csv(paste0("./00.data/",name,"/tacco_pred.csv.gz"),row.names=1)
coord$celltype = meta$pred_celltype_lowres
print(table(coord$celltype))

co = read.csv("./00.data/color_for_cellchat_0529.csv")
colors = co$color
names(colors)=co$celltype

white_theme <- function() {
  theme(panel.background = element_rect(fill = "black",color="black"),
        plot.background = element_rect(fill = "black",color="black"),
        legend.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        #panel.border = element_rect(fill = NA, color = "black", size = 0, linetype = "solid"),
        axis.line = element_line( linewidth = 0.5),
        axis.text = element_text( size = 15),
        axis.title = element_text( size = 14),
        legend.title = element_text( size = 12,color="white"),
        legend.text = element_text( size = 12,color="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "black"),
        axis.text.x = element_text(margin = margin(t = 10)),
        axis.title.x = element_text(margin = margin(t = 10))
       )
}

setwd(paste0("./fig3_supp/colocation"))
ct1 = celltypes_1[[stage]]
p1 = ggplot(coord, aes(x = x, y = y)) +
      geom_point(data = coord, color = "#C4C4C4",size = 0.5) +
      geom_point(data = coord[coord$celltype %in% c(ct1),], aes(color = celltype),size = 0.1) + scale_color_manual(values = colors) +
      geom_density_2d_filled(data = coord[coord$celltype == ct1,c("x","y")],alpha = 0.8,bins=8) +
      geom_density_2d(data = coord[coord$celltype == ct1,c("x","y")],linewidth=0.25,bins=8)+
      coord_fixed() +
      theme_bw() +
      xlim(min(coord$x),max(coord$x)) +
      ylim(min(coord$y),max(coord$y)) +
      white_theme()

ct2 = celltypes_2[[stage]]
p2 = ggplot(coord, aes(x = x, y = y)) +
      geom_point(data = coord, color = "#C4C4C4",size = 0.5) +
      geom_point(data = coord[coord$celltype %in% c(ct2),], aes(color = celltype),size = 0.1) + scale_color_manual(values = colors) +
      geom_density_2d_filled(data = coord[coord$celltype == ct2,c("x","y")],alpha = 0.8,bins=8) + scale_fill_viridis_d(option = "plasma")+
      geom_density_2d(data = coord[coord$celltype == ct2,c("x","y")],linewidth=0.25,bins=8,color="red")+
      coord_fixed() +
      theme_bw() +
      xlim(min(coord$x),max(coord$x)) +
      ylim(min(coord$y),max(coord$y)) +
      white_theme()
g = grid.arrange(p1,p2,nrow=1)
p <- cowplot::ggdraw(g) + 
  theme(plot.background = element_rect(fill="black", color = NA))

ggsave(plot=p,paste0(stage,"_",ct1,"_",ct2,"_density_contour.pdf"),w=20,h=7)

}



