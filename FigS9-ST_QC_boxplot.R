library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

coord = read.csv("All_coord_adjusted.csv.gz",row.names=1)
df = read.csv("spatial_selection_rds_path.txt",sep='\t')
df = df[1:30,]
data_list <- lapply(df$path, function(path) {
    objpath=paste0(path)
    st = df[df$path == path,"stage"]
    print(st)
    obj = readRDS(objpath)
    meta = obj@meta.data[,c("nCount_Spatial","nFeature_Spatial")]
    rownames(meta) = paste0(st,".",rownames(meta))
    return(meta)
    })
head(data_list[[1]])
max_cols_df <- data_list[[which.max(sapply(data_list, ncol))]]
    max_cols_names <- colnames(max_cols_df)

    # 遍历列表中的每个数据框，添加缺失的列并用0填充
    adjusted_data_frames <- lapply(data_list, function(df) {
     missing_cols <- setdiff(max_cols_names, colnames(df))
    for (col in missing_cols) {
        df[[col]] <- 0  # 添加缺失的列并用0填充
    }
    # 确保列的顺序与列数最多的数据框一致
    df <- df[, max_cols_names]
    return(df)
})

combined_data = do.call(rbind, adjusted_data_frames)

combined_data$RowName <- row.names(combined_data)
coord$RowName <- row.names(coord)

#dim(coord)
#dim(combined_data)
data_merged <- merge(combined_data, coord, by = "RowName")
row.names(data_merged) <- data_merged$RowName

library(ggplot2)
ord=c("E16.5_nor_A02999D1","E16.5_6h_C02938E1","E16.5_1d_A03091D1", "E16.5_2d_A03091D4" , "E16.5_3d_D03056D5",  "E16.5_6d_C02938F4" ,"E17.5_nor_C02845B1","E17.5_6h_B03415D4","E17.5_1d_A02988A6", "E17.5_2d_C02938E3", "E17.5_3d_C02926E1",  "E17.5_5d_C02926A6" ,"E18.5_nor_C02926A1","E18.5_6h_A02885C5","E18.5_1d_C03049D4" ,"E18.5_2d_A02988D2", "E18.5_3d1_C02926E5", "E18.5_5d_C02926D4" ,"P5_nor_C02926B1","P5_6h_C02926C6","P5_1d_C02926D3", "P5_2d_A02998D6", "P5_3d_B02621C2" ,"P5_6d_A03000C5","Adult_nor_A02988C1","Adult_1d_C01830C1D1","Adult_2d_C01830E4F4", "Adult_5d1_C01830E6F6", "Adult_10d_A02988C5","Adult_19d_B03424F4")
data_merged$stage = factor(data_merged$stage,levels=ord)

colorlist = c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
              "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
              "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
              "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
              "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
              "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
              "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
              "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
              "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
              "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
              "#191970", "#000080")
              

p0 <- ggplot(data_merged, aes(x = stage, y = log10(nCount_Spatial), fill = stage)) +
  geom_boxplot(color = "black") + # 设置箱线图的填充颜色和边框颜色
  scale_fill_manual(values = colorlist) +
  labs(y = "UMI count", title = "UMIs count", x = "") +# 设置 x 和 y 轴的标签以及图表标题
  #ylim(c(5000,30000))+ #修改y轴刻度范围
  theme_classic() +
  guides(fill = 'none') +
  theme(
    panel.grid.major = element_blank(),  # 隐藏主要网格线
    panel.grid.minor = element_blank(),# 隐藏次要网格线
    axis.text.x = element_text(angle = 55, hjust = 1, size = 15, face = "bold"), 
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, vjust = 1, size = 30),#改变title位置
    axis.title.y = element_text(size = 30, margin = margin(0, 20, 0, 0))
    #axis.line = element_line(linewidth = 1, colour = "black"), # 设置轴线粗细和颜色
  )
ggsave("All_Spatial_nCountSpatial_boxplot.png", plot = p0, width = 30, height = 10)

p0 <- ggplot(data_merged, aes(x = stage, y = log10(nFeature_Spatial), fill = stage)) +
  geom_boxplot(color = "black") + # 设置箱线图的填充颜色和边框颜色
  scale_fill_manual(values = colorlist) +
  labs(y = "Gene count", title = "Gene count", x = "") +# 设置 x 和 y 轴的标签以及图表标题
  #ylim(c(5000,30000))+ #修改y轴刻度范围
  theme_classic() +
  guides(fill = 'none') +
  theme(
    panel.grid.major = element_blank(),  # 隐藏主要网格线
    panel.grid.minor = element_blank(),# 隐藏次要网格线
    axis.text.x = element_text(angle = 55, hjust = 1, size = 15, face = "bold"), 
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, vjust = 1, size = 30),#改变title位置
    axis.title.y = element_text(size = 30, margin = margin(0, 20, 0, 0))
    #axis.line = element_line(linewidth = 1, colour = "black"), # 设置轴线粗细和颜色
  )
ggsave("All_Spatial_nFeatureSpatial_boxplot.png", plot = p0, width = 30, height = 10)