library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(patchwork)

### Get the parameters
parser <- argparse::ArgumentParser(description=cat("scRNA-seq data processing"))
parser$add_argument('-I','--input', help='input path of project')
parser$add_argument('-s','--samples',help='sample name')
parser$add_argument('-O','--out', help='out directory')

args <- parser$parse_args()
print(args)
dir.create(args$out)
setwd(args$out)

if (endsWith(args$input, '.gz')){
    obj <- read.table(args$input)
} else {
    obj <- tryCatch(Read10X(args$input), error = function(err) Read10X(args$input, gene.column=1))
}
obj <- CreateSeuratObject(obj, project=args$samples)
obj <- RenameCells(obj, args$samples)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-")

p = VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
ggsave(paste0(args$samples, "_QC_vln_original.png"), p)

saveRDS(obj, paste0(args$samples, "_raw.rds"))

obj <- subset(obj, subset = nFeature_RNA > 800 & 
        nFeature_RNA < 7500 & 
        nCount_RNA > 3000 & 
        nCount_RNA < 25000 &
        percent.mt < 7.5)

p = VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
ggsave(paste0(args$samples, "_QC_vln_1.png"), p)      

Find_doublet <- function(data){
        sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        nExp_poi <- round(0.05*ncol(data))
        p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
        data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
        colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
        data
}      

obj <- obj %>% 
    NormalizeData %>%  
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>% 
    RunUMAP(1:30) %>%
    Find_doublet      

obj$doublet_score = obj@meta.data[,5]                      
p1 = DimPlot(obj, group.by='doublet_info') + NoLegend()
p2 = FeaturePlot(obj, "doublet_score")
ggsave(paste0(args$samples, "_QC_doublet.png"), p1 + p2, width=14)

p = VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
ggsave(paste0(args$samples, "_QC_vln_2.png"), p)   

obj <- obj %>%
    subset(doublet_info=="Singlet") %>%
    NormalizeData %>%  
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>% 
    RunUMAP(1:30) %>%
    FindNeighbors %>%
    FindClusters

p1 = DimPlot(obj, label=T) + NoLegend()
p2 = FeaturePlot(obj, "percent.mt")
p3 = FeaturePlot(obj, "nFeature_RNA")
p4 = FeaturePlot(obj, "nCount_RNA")
ggsave(paste0(args$samples, "QC_umap.png"), (p1 | p2)/(p3 | p4))

saveRDS(obj, paste0(args$samples, ".rds"))       

#markers <- FindAllMarkers(obj, only.pos=T)
#write.csv(markers, paste0(args$samples, "_markers.csv"))   
#print("Done")                               
