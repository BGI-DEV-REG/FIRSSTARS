library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)

################# E/P ######################
df_EP = read.csv("../Myo_trac_EP_All_trajectory_genes_MI_05.csv")
cds = readRDS("../Myo_trac_EP_All_monocle3_cds.rds")
genes=df_EP$gene_short_name
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
max_col_indices <- max.col(pt.matrix, ties.method = "first")
# 根据最大值所在列的索引对行进行排序
sorted_mat <- pt.matrix[order(max_col_indices), ]

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  sorted_mat,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #clustering_method_rows = "ward.D2",
  #clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  #use_raster                   = FALSE,
  cluster_columns              = FALSE
)

pdf("plot_pseudotime_heatmap_EP.pdf",w=20,h=30)
draw(hthc,background="black")
dev.off()

################# Adult ######################

df = read.csv("../Myo_trac_Adult_All_trajectory_genes_MI_05.csv")
cds = readRDS("../Myo_trac_Adult_All_monocle3_cds.rds")
genes=df$gene_short_name
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
max_col_indices <- max.col(pt.matrix, ties.method = "first")
# 根据最大值所在列的索引对行进行排序
sorted_mat <- pt.matrix[order(max_col_indices), ]


hthc <- Heatmap(
  sorted_mat,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #clustering_method_rows = "ward.D2",
  #clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  #use_raster                   = FALSE,
  cluster_columns              = FALSE
)

pdf("plot_pseudotime_heatmap_Adult.pdf",w=20,h=30)
draw(hthc,background="black")
dev.off()
