library(Seurat)

obj = readRDS("Adult_Nor_fig1_harmony.rds")
meta = read.csv("anno_metadata_2024051.csv",row.names=1)

obj@meta.data = meta

color_df = read.csv('Fig1_color_seq_20240512.csv',row.names=1)
color_df = color_df[color_df$celltype %in% unique(obj@meta.data[["celltype"]]),]
prefix = "Adult_Nor_fig1_harmony"
colors = color_df$color
names(colors) = color_df$celltype
p <- DimPlot(obj, group.by="celltype", pt.size=0.1, cols=colors, label=F)
ggsave(paste0(prefix, "_", "celltype_20240520", ".pdf"),width=9,height=6)

