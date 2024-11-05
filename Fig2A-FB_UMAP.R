library(Seurat)
library(ggplot2)

########### EP
obj = readRDS("EP_FB_fig1_merge_harmony.rds")\
args = list()
args$new_annotation = "celltype_0519"
colorlist = read.csv("Fig2_color_seq_20240521.csv")
 
colors = colorlist$colors
names(colors) = colorlist$X

p <- DimPlot(obj, group.by = args$new_annotation, cols = colors, raster=F)  +
     theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "black"),   #将整个图形背景颜色改为黑色
        panel.background = element_rect(fill = "black"),  #将绘图区域背景颜色改为黑色
        #legend.position = "none" #去除图例
        legend.text = element_text(color = "white") #图例标签为白色
        ) 
ggsave("EP_UMAP_black_20240530.png", p, width = 10, height = 7)

########### Adult
obj = readRDS("Adult_FB_fig1_merge_harmony.rds")
args = list()
args$new_annotation = "celltype_240516"
 
colors = colorlist$colors
names(colors) = colorlist$X

p <- DimPlot(obj, group.by = args$new_annotation, pt.size=1,cols = colors, raster=F)  +
     theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "black"),   #将整个图形背景颜色改为黑色
        panel.background = element_rect(fill = "black"),  #将绘图区域背景颜色改为黑色
        #legend.position = "none" #去除图例
        legend.text = element_text(color = "white") #图例标签为白色
        ) 
ggsave("UMAP_black_20240523.png", p, width = 10, height = 7)

