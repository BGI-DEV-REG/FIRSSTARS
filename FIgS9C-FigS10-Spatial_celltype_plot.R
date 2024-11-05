library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

################################ combined all spatial sections coordinates ################################

df = read.csv('spatial_selection_rds_path.txt',sep='\t')
df = df[1:30,]

######## 合并数据,调整角度 ########

rot <- list("E16.5_nor_A02999D1"=130,"E16.5_6h_C02938E1" = -135,"E16.5_1d_A03091D1" = 60, "E16.5_2d_A03091D4" = -115, "E16.5_3d_D03056D5" = 70,  "E16.5_6d_C02938F4" = -82,"E17.5_nor_C02845B1"=45,"E17.5_6h_B03415D4" = 60,"E17.5_1d_A02988A6" = 180, "E17.5_2d_C02938E3" = 90, "E17.5_3d_C02926E1" = 0,  "E17.5_5d_C02926A6" = -175,"E18.5_nor_C02926A1"=0,"E18.5_6h_A02885C5"=10,"E18.5_1d_C03049D4" = 100,"E18.5_2d_A02988D2" = 0, "E18.5_3d1_C02926E5" = 0, "E18.5_5d_C02926D4" = -85,"P5_nor_C02926B1"=175,"P5_6h_C02926C6" =  180,"P5_1d_C02926D3" = 170, "P5_2d_A02998D6" = 0, "P5_3d_B02621C2" = -100,"P5_6d_A03000C5"=190,"Adult_nor_A02988C1"=-95,"Adult_1d_C01830C1D1" = 180,"Adult_2d_C01830E4F4" = -90, "Adult_5d1_C01830E6F6" = 90, "Adult_10d_A02988C5" = 90,"Adult_19d_B03424F4"=45)

coord_list=list()
for (path in df$path){
    objpath=paste0(path)
    st = df[df$path == path,"stage"]
    print(st)
    obj = readRDS(objpath)
    obj$coor_y = obj@images[[names(obj@images)]]@coordinates$row
    obj$coor_x = obj@images[[names(obj@images)]]@coordinates$col
    coord=data.frame(x = obj$coor_x, y = obj$coor_y)
    
    print(rot[[st]])
    angle <- rot[[st]]
    angle <- -angle * (pi / 180)
    rotation_matrix <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2)
    rotated_coords <- as.data.frame(t(rotation_matrix %*% t(as.matrix(coord))))
    colnames(rotated_coords) <- c("x", "y")
    
    cellid = read.csv(paste0("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/02.ST/12.DEG/00.Data/",st,"_Injury.txt"),header=F)
    obj$type = "non-Wounds"
    obj$type[cellid$V1]="Wounds"
    obj$stage = st
    print(table(obj$type))
    meta = obj@meta.data[,c("type","stage")]
  
    data = cbind(rotated_coords,meta)
    coord_list[[st]]=data
}

coord_merge  = do.call(rbind, coord_list)
write.csv(coord_merge,gzfile("coord_angle_adjust.csv.gz"))

######## 分组缩放样本 ########
group = list(nor = c("Adult_nor_A02988C1","P5_nor_C02926B1","E18.5_nor_C02926A1","E17.5_nor_C02845B1","E16.5_nor_A02999D1"),
            inflam1 = c("Adult_1d_C01830C1D1","P5_6h_C02926C6","E18.5_6h_A02885C5","E17.5_6h_B03415D4","E16.5_6h_C02938E1"),
            inflam2 = c("Adult_2d_C01830E4F4","P5_1d_C02926D3","E18.5_1d_C03049D4","E17.5_1d_A02988A6","E16.5_1d_A03091D1"),
            pro1 = c("Adult_5d1_C01830E6F6","P5_2d_A02998D6","E18.5_2d_A02988D2","E17.5_2d_C02938E3","E16.5_2d_A03091D4"),
            pro2 = c("Adult_10d_A02988C5" ,"P5_3d_B02621C2" ,"E18.5_3d1_C02926E5","E17.5_3d_C02926E1","E16.5_3d_D03056D5"),
            re1 = c("Adult_19d_B03424F4","P5_6d_A03000C5","E18.5_5d_C02926D4","E17.5_5d_C02926A6","E16.5_6d_C02938F4"))

coord_merge$x_adjust = coord_merge$x
coord_merge$y_adjust = coord_merge$y

######################## case1:计算缩放因子并保存为csv 

scale_factor_list=data.frame(stage = df$stage)
scale_factor_list$scale_factor=1
for (i in 1:length(names(group))){
    n = names(group)[i]
    print(n)
    stages = group[[n]]
    sub_coord = coord_merge[coord_merge$stage %in% stages,]
    range_x = sub_coord %>% 
    group_by(stage) %>%
    summarise(x_range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    max_x = max(range_x$x_range) 
    for (st in stages){
    print(st)
    x_len = max(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE) - min(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE)
    scale_factor = max_x / x_len
    print(scale_factor)
    scale_factor_list[scale_factor_list$stage == st,"scale_factor"]=scale_factor
    coord_merge[coord_merge$stage == st,]$x_adjust = coord_merge[coord_merge$stage == st,]$x_adjust*scale_factor
    coord_merge[coord_merge$stage == st,]$y_adjust = coord_merge[coord_merge$stage == st,]$y_adjust*scale_factor
    print(max(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE) - min(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE))
    }
}
write.csv(scale_factor_list,"scale_factor_list_all.csv")

######################## case2: 读入缩放因子表格，缩放坐标
scale_factor_list= read.csv("scale_factor_list_all.csv",row.names=1)
for (i in 1:length(names(group))){
    n = names(group)[i]
    print(n)
    stages = group[[n]]
    sub_coord = coord_merge[coord_merge$stage %in% stages,]
    range_x = sub_coord %>% 
    group_by(stage) %>%
    summarise(x_range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    max_x = max(range_x$x_range) 
    for (st in stages){
    print(st)
    scale_factor = scale_factor_list[scale_factor_list$stage == st,"scale_factor"]
    print(scale_factor)
    coord_merge[coord_merge$stage == st,]$x_adjust = coord_merge[coord_merge$stage == st,]$x_adjust*scale_factor
    coord_merge[coord_merge$stage == st,]$y_adjust = coord_merge[coord_merge$stage == st,]$y_adjust*scale_factor
    print(max(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE) - min(coord_merge[coord_merge$stage == st,]$x_adjust, na.rm = TRUE))
    }
}

######## 调整X轴坐标 ########

ord = c(
"Adult_nor_A02988C1","Adult_1d_C01830C1D1" ,"Adult_2d_C01830E4F4" , "Adult_5d1_C01830E6F6" , "Adult_10d_A02988C5" ,"Adult_19d_B03424F4",
"P5_nor_C02926B1","P5_6h_C02926C6","P5_1d_C02926D3" , "P5_2d_A02998D6", "P5_3d_B02621C2" ,"P5_6d_A03000C5",
"E18.5_nor_C02926A1","E18.5_6h_A02885C5","E18.5_1d_C03049D4" ,"E18.5_2d_A02988D2", "E18.5_3d1_C02926E5", "E18.5_5d_C02926D4","E17.5_nor_C02845B1","E17.5_6h_B03415D4","E17.5_1d_A02988A6", "E17.5_2d_C02938E3", "E17.5_3d_C02926E1",  "E17.5_5d_C02926A6","E16.5_nor_A02999D1","E16.5_6h_C02938E1","E16.5_1d_A03091D1", "E16.5_2d_A03091D4", "E16.5_3d_D03056D5", "E16.5_6d_C02938F4"
)

coord_merge$stage = factor(coord_merge$stage, levels = ord)
coord_merge = coord_merge %>% arrange(coord_merge$stage)
data = coord_merge

range_x = lapply(group,function(l){
    sub_coord = coord_merge[coord_merge$stage %in% l,]
    range = sub_coord %>% 
    group_by(stage) %>%
    summarise(x_range = max(x_adjust, na.rm = TRUE) - min(x_adjust, na.rm = TRUE)) 
    return(max(range$x_range))})

for (i in 1:length(names(group))){
    n = names(group)[i]
    print(n)
    if (i > 1){
    x_base=sum(unlist(range_x[1:i-1]))
    }else{
    x_base=0
    }
    x_sep = x_base+(i-1)*10000
    stages = group[[n]]
    for (st in stages){
    print(st)
    min_x =  min(coord_merge[coord_merge$stage == st,]$x_adjust)
    coord_merge[coord_merge$stage == st,]$x_adjust = coord_merge[coord_merge$stage == st,]$x_adjust - min_x +x_sep
    print(min(coord_merge[coord_merge$stage == st,]$x_adjust))
    }
}

######## 调整Y轴坐标 ########
group_y = list(Adult=c("Adult_nor_A02988C1","Adult_1d_C01830C1D1" ,"Adult_2d_C01830E4F4" , "Adult_5d1_C01830E6F6" , "Adult_10d_A02988C5" ,"Adult_19d_B03424F4"),P5=c("P5_nor_C02926B1","P5_6h_C02926C6","P5_1d_C02926D3" , "P5_2d_A02998D6", "P5_3d_B02621C2" ,"P5_6d_A03000C5"),E18.5 = c("E18.5_nor_C02926A1","E18.5_6h_A02885C5","E18.5_1d_C03049D4" ,"E18.5_2d_A02988D2", "E18.5_3d1_C02926E5", "E18.5_5d_C02926D4"),E17.5=c("E17.5_nor_C02845B1","E17.5_6h_B03415D4","E17.5_1d_A02988A6", "E17.5_2d_C02938E3", "E17.5_3d_C02926E1",  "E17.5_5d_C02926A6"),E16.5=c("E16.5_nor_A02999D1","E16.5_6h_C02938E1","E16.5_1d_A03091D1", "E16.5_2d_A03091D4", "E16.5_3d_D03056D5", "E16.5_6d_C02938F4")
)


range_y = lapply(group_y,function(l){
    sub_coord = coord_merge[coord_merge$stage %in% l,]
    range = sub_coord %>% 
    group_by(stage) %>%
    summarise(y_range = max(y_adjust, na.rm = TRUE) - min(y_adjust, na.rm = TRUE)) 
    return(max(range$y_range))})

for (i in 1:length(names(group_y))){
    n = names(group_y)[i]
    print(n)
    if (i > 1){
    y_base=sum(unlist(range_y[1:i-1]))
    }else{
    y_base=0
    }
    y_sep = y_base+(i-1)*10000
    stages = group_y[[n]]
    for (st in stages){
    print(st)
    min_y =  min(coord_merge[coord_merge$stage == st,]$y_adjust)
    coord_merge[coord_merge$stage == st,]$y_adjust = coord_merge[coord_merge$stage == st,]$y_adjust - min_y +y_sep
    print(min(coord_merge[coord_merge$stage == st,]$y_adjust))
    }
}

##### 单独处理E16.5——2d样本
st="E16.5_2d_A03091D4"
mean_x = mean(coord_merge[coord_merge$stage == st,]$x_adjust)
mean_group = mean(coord_merge[coord_merge$stage == "E17.5_2d_C02938E3",]$x_adjust)
x_diff = mean_group-mean_x
coord_merge[coord_merge$stage == st,]$x_adjust = coord_merge[coord_merge$stage == st,]$x_adjust + x_diff

data = coord_merge[,c("stage","x_adjust","y_adjust")]
write.csv(data,gzfile("All_coord_adjusted_20240519.csv.gz"))

################################ Major celltype visualization ################################
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

######## 读取全部切片坐标 ########
coord = read.csv("All_coord_adjusted_20240519.csv.gz",row.names=1)
df = read.csv('/spatial_selection_rds_path_fig1_meta.txt',sep='\t')
df = df[1:30,]
coord[,"celltype"]=""
for (path in df$anno){
    meta_path=paste0(path)
    st = df[df$anno == path,"stage"]
    print(st)
    meta = read.csv(meta_path,row.names=1)
    rownames(meta) = paste0(st,".",rownames(meta))
coord[coord$stage == st,"celltype"]=meta[rownames(coord[coord$stage == st,]),"pred_celltype"]
}
coord$celltype[coord$celltype == "FB_basement"]="FB"


df = read.csv("/data/work/Figure1/spatial_plot/Fig1_color_seq_supp_20240519.csv")
colors = df$color
names(colors)=df$celltype

scale_factor_list= read.csv("/data/work/Figure1/spatial_plot/scale_factor_list_all.csv",row.names=1)
stages = scale_factor_list$stage
xlist=c()
ylist=c()
xendlist=c()
yendlist=c()

for (i in 1:length(stages)){
    st=stages[i]
    #print(st)
   scale_factor = scale_factor_list[scale_factor_list$stage == st,"scale_factor"]
   sub_coord = coord[coord$stage == st,]
  # print(sub_coord[1:5,1:3])
    xlist[i]=floor(min(sub_coord$x_adjust)+500)
    ylist[i]=floor(max(sub_coord$y_adjust)+1000)
    xendlist[i]=floor(min(sub_coord$x_adjust)+500+2000*scale_factor)
    yendlist[i]=floor(max(sub_coord$y_adjust)+1000)
}

p <-ggplot(coord,aes(x=x_adjust,y=y_adjust,color = celltype))+
            geom_point(aes(x=x_adjust,y=y_adjust,color = celltype),size=0.1)+coord_fixed()+
            theme_void()+
            theme(panel.background = element_rect(fill = "black"),
                  plot.background = element_rect(fill = "black"),
                  text = element_text(color = "white"),
                  plot.title = element_text(color = "white"))+scale_color_manual(values = colors)
segments_df <- data.frame(
  xlist = xlist,
  ylist = ylist,
  xendlist = xendlist,
  yendlist = yendlist
)
p = p+
geom_segment(
    data = segments_df,
    aes(x = xlist, y = ylist, xend = xendlist, yend = yendlist), color = "white",inherit.aes = FALSE,linewidth = 2
  )
ggsave(plot=p,paste0("All_Spatial_Fig1_celltype.png"),height = 25,width = 42

################################ Wounds region visualization ################################
coord = read.csv("All_coord_adjusted_20240519.csv.gz",row.names=1)
df = read.csv('/Fig1/spatial_selection_rds_path_fig1_meta.txt',sep='\t')
df = df[1:30,]
coord[,"type"]="non-Wounds"
for (st in df$stage){
    print(st)
    cellid = read.csv(paste0(st,"_Injury.txt"),header=F)
    cellid$V1 = paste0(st,".",cellid$V1)
    coord[cellid$V1,"type"]="Wounds"
}

coord_merge=coord
coord$type = factor(coord$type,levels=c("non-Wounds","Wounds"))
colors = c("red","lightgrey")
names(colors)=c("Wounds","non-Wounds")
outpath="/Fig1"
p <-ggplot(coord,aes(x=x_adjust,y=y_adjust,color = type))+
            geom_point(aes(x=x_adjust,y=y_adjust,color = type),size=0.1)+coord_fixed()+
            theme_void()+
            theme(panel.background = element_rect(fill = "black"),
                  plot.background = element_rect(fill = "black"),
                  text = element_text(color = "white"),
                  plot.title = element_text(color = "white"))+
            scale_color_manual(values = colors)
ggsave(plot=p,paste0(outpath,"/","All_Spatial_Fig1_region.pdf"),height = 25,width = 42)



