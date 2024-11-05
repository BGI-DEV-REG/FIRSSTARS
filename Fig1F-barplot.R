############################# get zoom-in region cellid ###########################
### E16.5_6d 根据坐标直接截取
coord = read.csv("All_coord_adjusted_20240519.csv.gz",row.names=1) 
df = coord[coord$stage == "E16.5_6d_C02938F4",]
meta = read.csv("./E16.5_6D_C02938F4_2/pred.csv.gz",row.names=1)
rownames(meta)=paste0("E16.5_6d_C02938F4.",rownames(meta))
df$celltype = meta$pred_celltype2

co = read.csv("Fig1_color_seq_20240512.csv")
colors = co$color
names(colors)=co$celltype
colors = colors[names(colors) %in% unique(df$celltype)]
output="/04.Fig1_reg_vs_fibro/"

p = ggplot() +
    geom_point(data = df,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.1) +
    geom_vline(xintercept = c(186500, 190200), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(179400,182700), linetype = "dashed", color = "black") +
    #geom_hline(yintercept = 3100, linetype = "dashed", color = "black")+
    scale_color_manual(values = colors) +  # Define colors for each cluster
     theme_bw() +  # Use theme_void to remove most plot elements
    coord_fixed(ratio = 1) 
ggsave(plot=p,paste0(output,"E16.5_6d.png"),w=10,h=5)

df$cell_id = sub("E16.5_6d_C02938F4.","",rownames(df))

no_Wound = df[(df$x_adjust >186500 & df$x_adjust < 190200),]
p = ggplot() +
    geom_point(data = no_Wound,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.5) + coord_fixed()+
        theme_void()+
            theme(panel.background = element_rect(fill = "black"),
                  plot.background = element_rect(fill = "black"),
                  text = element_text(color = "white"),
                  plot.title = element_text(color = "white"))+
            scale_color_manual(values = colors)
ggsave(plot=p,paste0(output,"E16.5_6d_no_Wound.png"),w=5,h=3)


no_Wound_cellid = no_Wound$cell_id
write.table(no_Wound_cellid,"E16.5_6d_no_Wound_select.txt",quote=F,row.names=F)

Wound = df[(df$x_adjust >179400 & df$x_adjust < 182700),]
p = ggplot() +
    geom_point(data = Wound,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.5) +theme_void()+ coord_fixed()+
            theme(panel.background = element_rect(fill = "black"),
                  plot.background = element_rect(fill = "black"),
                  text = element_text(color = "white"),
                  plot.title = element_text(color = "white"))+
    scale_color_manual(values = colors)   # Define colors for each cluster
      
ggsave(plot=p,paste0(output,"E16.5_6d_Wound.png"),w=5,h=3)
Wound_cellid = Wound$cell_id
write.table(Wound_cellid,"E16.5_6d_Wound_select.txt",quote=F,row.names=F)

### Adult_19d 根据坐标直接截取
coord = read.csv("All_coord_adjusted_20240519.csv.gz",row.names=1) 
df = coord[coord$stage == "Adult_19d_B03424F4",]
meta = read.csv("./Adult_19D_B03424F4_2/pred_forFig1.csv.gz",row.names=1)
rownames(meta)=paste0("Adult_19d_B03424F4.",rownames(meta))
df$celltype = meta$pred_celltype2

p = ggplot() +
    geom_point(data = df,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.1) +
    geom_vline(xintercept = c(185750, 188750), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(178300,181300), linetype = "dashed", color = "black") +
    geom_hline(yintercept = 3100, linetype = "dashed", color = "black")+
    scale_color_manual(values = colors) +  # Define colors for each cluster
     theme_bw() +  # Use theme_void to remove most plot elements
    coord_fixed(ratio = 1) 
ggsave(plot=p,paste0(output,"Adult_19d.png"),w=10,h=5)

df$cell_id = sub("Adult_19d_B03424F4.","",rownames(df))
no_Wound = df[(df$x_adjust >185750 & df$x_adjust < 188750),]
p = ggplot() +
    geom_point(data = no_Wound,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.1) + coord_fixed()+
        theme_void()+
            theme(panel.background = element_rect(fill = "black"),
                  plot.background = element_rect(fill = "black"),
                  text = element_text(color = "white"),
                  plot.title = element_text(color = "white"))+
            scale_color_manual(values = colors)
ggsave(plot=p,paste0(output,"Adult_19d_no_Wound.png"),w=6,h=3)


no_Wound_cellid = no_Wound$cell_id
write.table(no_Wound_cellid,"Adult_19d_no_Wound_select.txt",quote=F,row.names=F)

Wound = df[(df$x_adjust >178300 & df$x_adjust < 181300 & df$y_adjust > 3100),]
p = ggplot() +
    geom_point(data = Wound,
               aes(x = x_adjust, y = y_adjust, color = as.character(celltype)),
               size = 0.1) +
    scale_color_manual(values = colors) +  # Define colors for each cluster
     theme_bw() +  # Use theme_void to remove most plot elements
    coord_fixed(ratio = 1) 
ggsave(plot=p,paste0(output,"Adult_19d_Wound.png"),w=6,h=3)
Wound_cellid = Wound$cell_id
write.table(Wound_cellid,"Adult_19d_Wound_select.txt",quote=F,row.names=F)

############################# Bar plot ###########################
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(readxl)
library(reshape2)
output="/04.Fig1_reg_vs_fibro/"
setwd(output)

### Adult_19d
meta = read.csv("./Adult_19D_B03424F4_2/pred_forFig1.csv.gz",row.names=1)

Wound = read.csv("Adult_19d_Wound_select.txt")
colnames(Wound)="cellid"
no_Wound = read.csv("Adult_19d_no_Wound_select.txt")
colnames(no_Wound)="cellid"
Wound$type = "Wound"
no_Wound$type = "non_Wound"
Wound$celltype = meta[rownames(meta) %in% Wound$cellid,"pred_celltype2"]
no_Wound$celltype = meta[rownames(meta) %in% no_Wound$cellid,"pred_celltype2"]

data = rbind(Wound,no_Wound)
rownames(data)=data$cellid

Cellratio <- prop.table(table(data$type,data$celltype), margin = 1) # margin = 1
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("type", "celltype", "Percentage")

df1 = Cellratio[grepl("FB",Cellratio$celltype),]
df2 = Cellratio[grepl("Immu",Cellratio$celltype),]
df3 = rbind(df1,df2)
df3$celltype = factor(df3$celltype,levels=c("FB_Injury","FB_Fascia","FB_Papillary","FB_Reticular","FB_DP","Immu_Myeloid","Immu_T_NK","Immu_Neu","Immu_pDC_B"))
df3$type = factor(df3$type,levels=c("Wound","non_Wound"))
p = ggplot(df3,aes(x=celltype,y=Percentage,fill=type))+
    geom_bar(stat="identity",position="dodge")+
    theme_classic()+
    theme( axis.text.x = element_text(angle = 45, hjust = 1, color = "white"), 
    axis.text.y = element_text(size = 12, color = "white"), 
    plot.background = element_rect(fill = "black"), 
    panel.background = element_rect(fill = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"))
ggsave(plot=p,paste0(output,"Adult_19d_FB_Immu_percentage.pdf"),w=6,h=4)

### E16.5_6d
meta = read.csv("./E16.5_6D_C02938F4_2/pred.csv.gz",row.names=1)

Wound = read.csv("E16.5_6d_Wound_select.txt")
colnames(Wound)="cellid"
no_Wound = read.csv("E16.5_6d_no_Wound_select.txt")
colnames(no_Wound)="cellid"
Wound$type = "Wound"
no_Wound$type = "non_Wound"
Wound$celltype = meta[rownames(meta) %in% Wound$cellid,"pred_celltype2"]
no_Wound$celltype = meta[rownames(meta) %in% no_Wound$cellid,"pred_celltype2"]

data = rbind(Wound,no_Wound)
rownames(data)=data$cellid

Cellratio <- prop.table(table(data$type,data$celltype), margin = 1) # margin = 1
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("type", "celltype", "Percentage")

df1 = Cellratio[grepl("FB",Cellratio$celltype),]
df2 = Cellratio[grepl("Immu",Cellratio$celltype),]
df3 = rbind(df1,df2)
df3$celltype = factor(df3$celltype,levels=c("FB_Injury","FB_Fascia","FB_Papillary","FB_Reticular","FB_DP","FB_DS","Immu_Myeloid","Immu_Neu","Immu_pDC_B","Immu_T_NK"))
df3$type = factor(df3$type,levels=c("Wound","non_Wound"))

p = ggplot(df3,aes(x=celltype,y=Percentage,fill=type))+
    geom_bar(stat="identity",position="dodge")+
    theme_classic()+
    theme( axis.text.x = element_text(angle = 45, hjust = 1, color = "white"), 
    axis.text.y = element_text(size = 12, color = "white"), 
    plot.background = element_rect(fill = "black"), 
    panel.background = element_rect(fill = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"))
ggsave(plot=p,paste0(output,"E16.5_6d_FB_Immu_percentage.pdf"),w=6,h=4)

