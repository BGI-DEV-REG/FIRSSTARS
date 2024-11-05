# segmented data to Seurat SpatailObject
parser = argparse::ArgumentParser(description="Script segmented data to Seurat SpatailObject")
parser$add_argument('-I','--input', help='input directory')
parser$add_argument('-D','--id', help='tissue ID')
parser$add_argument('-O','--output', help='out directory')
parser$add_argument('-N','--nCount_Spatial', help='nCount_RNA')
args = parser$parse_args()

library(dplyr)
library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(SeuratDisk)

source('/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/fish/02.Clustering/01.script/LoadBGI_Spatial_fish.R')

id = args$id
wd = paste0(args$input)

setwd(wd)

#### Creat obj 
obj = LoadBGI_Spatial(paste0(id,'_scgem.csv.gz'),
                             outdir = getwd(),
                             bin_data = F,
                             bin_size = 50,
                             cell_mask = T,
                             area_mask = F,
                             save_as = "rds",
                             pro_name = id,
                             UMI_GreyScale_Image = F,
                             assay = "Spatial",
                             slice = id,
                             delete_bg = T,csv=F)



#### QC and selecting cells for further analysis
obj <- obj[,obj$nCount_Spatial > as.numeric(args$nCount_Spatial)] # nCount_RNA can be adjusted

print(summary(obj@meta.data))
#### SCT
DefaultAssay(obj) = 'Spatial'
obj <- SCTransform(obj, verbose = FALSE,assay = 'Spatial')
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE) ## or assay = "RNA"
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

saveRDS(obj,paste0(args$output,'/',args$id,'.SCT.rds'))

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
             "#191970", "#000080",
             "#6495ED", "#1E90FF", "#00BFFF", "#00FFFF", "#FF1493",
             "#FF00FF", "#A020F0", "#63B8FF", "#008B8B", "#54FF9F",
             "#00FF00", "#76EE00", "#FFF68F")

p <- DimPlot(obj,label = T,reduction = "umap", cols = colorlist)
ggplot2::ggsave(plot = p ,file = paste0(args$output,'/',args$id,'.UMAP.SCT.png'),width = 8,height = 8)

#p <- SpatialDimPlot(obj, label = TRUE, label.size = 1 , cols = colorlist,pt.size.factor = 5, stroke = 0)

p <- SpatialDimPlot(obj, label = TRUE, label.size = 1, pt.size.factor = 5, stroke = 0)
ggplot2::ggsave(plot = p ,file = paste0(args$output,'/',args$id,'.Spatial.SCT.png'),width = 3.5,height = 3.5)



#### Meta.data 
meta = as.data.frame(obj@meta.data)
print(summary(meta))
write.csv(meta,paste0(args$output,'/',args$id,'.Seurat.metadata.SCT.csv'),quote=F)

color.list <- colorRampPalette(c("blue", "yellow", "red"))(length(unique(obj$nFeature_Spatial)))
names(color.list) = unique(obj$nFeature_Spatial)

p <- SpatialFeaturePlot(obj, features = "nFeature_SCT",
                        #min.cutoff = q0,
                       # max.cutoff = quantile(obj$nFeature_SCT, 0.95),
                        pt.size.factor = 120,
                        stroke = 0) +
  scale_color_gradientn(colours  = color.list) +
  theme(legend.position = "right") +
  ggtitle("nFeature_SCT")

 ggsave(p, file = paste0(args$output,'/', args$id,'_nfeature_dimplot.SCT.png'), width = 8, height = 5)


#### FindAllMarkers
DefaultAssay(obj) = 'SCT'
markers <- FindAllMarkers(obj,only.pos=T)

write.csv(markers,paste0(args$output,'/',args$id,'.markers.SCT.csv'),quote=F)

#### Convert as h5ad file
library(reticulate)
use_python("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/USER/liuyang9/soft/MINIconda/miniconda3/bin/python")

genes <- as.data.frame(rownames(obj), row.names = rownames(obj))
names(genes) <- "Gene"

cells <- as.data.frame(colnames(obj), row.names = colnames(obj))
names(cells) <- "CellID"

row <- obj@images[[1]]@coordinates$row
col <- obj@images[[1]]@coordinates$col
coordinates <- list(matrix(c(row, col), ncol = 2))
names(coordinates) <- "spatial"

ann <- import("anndata")
ad <- ann$AnnData(X = obj@assays$SCT@data, obs = genes, var = cells, varm = coordinates,layers = list(counts = obj@assays$SCT@counts))

ad <- ad$T
ad$write_h5ad(file.path(args$output, paste0(args$id, "_SCT.h5ad")))
