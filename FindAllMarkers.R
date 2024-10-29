# findallmarkers

### Get the parameters
parser <- argparse::ArgumentParser(description=cat("findallmarkers\n"))
parser$add_argument('-P','--path',help='full path of integrated object')
parser$add_argument('-G','--group',default="seurat_clusters", help='group by [default %(default)s]')
parser$add_argument('-L','--logfc',type="double",default=0.5,help='logfc.threshold [default %(default)s]')

args <- parser$parse_args()
print(args)

library(Seurat)
library(dplyr)
library(future)
options(future.globals.maxSize = 500 * 1024 ^ 3)
plan('multisession', workers=4)

obj <- readRDS(args$path)
DefaultAssay(obj) = 'RNA'
Idents(obj) <- obj@meta.data[[args$group]]

markers <- FindAllMarkers(obj, only.pos=T, logfc.threshold = args$logfc)
write.csv(markers, sub(".rds", paste(args$group, args$logfc, "markers.csv", sep="_"), args$path))
