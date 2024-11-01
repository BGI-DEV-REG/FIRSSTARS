library(stringr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(CellChat)
library(Seurat)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(reshape)

ump = read.csv("./table/All_Immu_Seed_5.0_3.0_obs_umap.csv")
ump$day = gsub("D","",ump$TimePostInjury)
ump[ump$day=="Nor","day"] = 0
ump[ump$day=="6H","day"] = 0.25
ump$day = as.numeric(ump$day)

ump[ump$celltype=="Macr_Cd93","celltype"] = "Macr_Mmp14"
ump[ump$celltype=="Macr_Il10","celltype"] = "Macr_Ms4a7"

acmat = table(ump[ump$Age=="Adult","celltype"],ump[ump$Age=="Adult","day"])%>%as.data.frame()%>%melt()
acmat = data.frame(num = acmat$value/sum(acmat$value),cell =paste(acmat$Var1,acmat$Var2,sep = "@") )
aflow = read.csv("./table/0814_Adult_wot_flow.csv")
aflow$cell = paste(aflow$celltype,aflow$day,sep = "@")
flow = melt(cflow[,c(3:8)])
flow$day1 = str_split(flow$cell,"@",simplify = T)[,2]
flow$cell1 = str_split(flow$cell,"@",simplify = T)[,1]
flow <- flow %>%
  mutate(day2 = recode(day1, '1' = 2, '2' = 5, '5' = 10,'10' =19))
colnames(flow) = c("cell","cell2","flow","day1","cell1","day2")

agrid = read.csv("./table/grid_xy.csv")
rownames(agrid) = paste(agrid$celltype,agrid$day,sep = "@")
flow$gridX1<-agrid$gridX[match(flow$cell,rownames(agrid))]
flow$gridY1<-agrid$gridY[match(flow$cell,rownames(agrid))]
flow$gridX2<-agrid$gridX[match(paste(flow$cell2,flow$day2,sep = "@"),rownames(agrid))]
flow$gridY2<-agrid$gridY[match(paste(flow$cell2,flow$day2,sep = "@"),rownames(agrid))]
flow$num_weight =acmat$num[match(flow$cell,acmat$cell)]

pcmat = table(ump[ump$Age=="P5","celltype"],ump[ump$Age=="P5","day"])%>%as.data.frame()%>%melt()
pcmat = data.frame(num = pcmat$value/sum(pcmat$value),cell =paste(pcmat$Var1,pcmat$Var2,sep = "@") )
pflow = read.csv("./table/0814_P5_wot_flow.csv")
pflow$cell = paste(pflow$celltype,pflow$day,sep = "@")
flow = melt(pflow[,c(3:8)])
flow$day1 = str_split(flow$cell,"@",simplify = T)[,2]
flow$cell1 = str_split(flow$cell,"@",simplify = T)[,1]
flow <- flow %>%
  mutate(day2 = recode(day1, '0.25' = 1, '1' = 2, '2' = 3,'3' =6))
colnames(flow) = c("cell","cell2","flow","day1","cell1","day2")

pgrid = read.csv("./table/P5_grid_xy.csv")
rownames(pgrid) = paste(pgrid$celltype,pgrid$day,sep = "@")
flow$gridX1<-pgrid$gridX[match(flow$cell,rownames(pgrid))]
flow$gridY1<-pgrid$gridY[match(flow$cell,rownames(pgrid))]
flow$gridX2<-pgrid$gridX[match(paste(flow$cell2,flow$day2,sep = "@"),rownames(pgrid))]
flow$gridY2<-pgrid$gridY[match(paste(flow$cell2,flow$day2,sep = "@"),rownames(pgrid))]
flow$num_weight =pcmat$num[match(flow$cell,pcmat$cell)]

# ump <- ump %>%
#   mutate(day = recode(day, '0.25' = 1, '1' = 2, '2' = 3,'3' =4,'5'=5,'6'=6,'10'=7,'19'=8))
#ump_pt = read.csv("./table/ump_point.csv")


# 


GetPY = function(flow,i,px=0.5,py=0.0,size_factor = 2,weight_factor = 5){
  x = mean(c(flow[i,"gridX1"],flow[i,"gridX2"]))
  y = mean(c(flow[i,"gridY1"],flow[i,"gridY2"]))
  #c = mean(c(ump_pt[r1,"color"],ump_pt[r2,"color"]))
  pt = data.frame(gridX=c(flow[i,"gridX1"],x+px,x-px,flow[i,"gridX2"]),
                  gridY=c(flow[i,"gridY1"],y+py,y-py,flow[i,"gridY2"]),
                  size = size_factor*flow$flow[i],weight = weight_factor*flow$num_weight[i],
                  group=i)
  return(pt)
}


py  =data.frame()
for(i in 1:nrow(flow)){
  py = rbind(py,GetPY(flow,i))
}

py$gridX <- as.numeric(as.character(py$gridX))

p1 = ggplot() + 
  labs(title ="P5",x = 'Day',y = "celltype")+
  theme_black()+
  geom_bezier2(aes(x = gridX, y = gridY, group = group, linewidth = weight, alpha = size / 2, color = gridY),
               data = py, show.legend = TRUE) +  # 确保显示图例
  scale_color_gradientn(colors = c("#f9320c","#FF1493", "#A020F0","#26c6da","#FFF300"))+
  new_scale("color")+
  geom_point(pgrid, mapping = aes(x = gridX, y = gridY, color = celltype, size = num * 20), alpha = 1)+
  scale_color_manual(values = col1)+
  scale_size_continuous(range = c(1, 10))+  
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(size = guide_legend(title = "Cell number",override.aes = list(color = "white", shape = 16)))

p2 = ggplot() + 
  labs(title ="P5",x = 'Day',y = "celltype")+
  theme_black()+
  geom_bezier2(aes(x = gridX, y = gridY, group = group, linewidth = weight, alpha = size / 2, color = gridY),
               data = py, show.legend = TRUE) +  # 确保显示图例
  scale_color_gradientn(colors = c("#f9320c","#FF1493", "#A020F0","#26c6da","#FFF300"))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(linewidth = guide_legend(title = "weighted prob sum",override.aes = list(color = "white")),      # 为 linewidth 添加图例
         alpha = guide_legend(title = "Cell percent",override.aes = list(color = "white")))    # 为 alpha 添加图例

ggsave(filename = "./figures/P5_flow.pdf",plot = p1)
ggsave(filename = "./figures/P5_flow_legend.pdf",plot = p2)
