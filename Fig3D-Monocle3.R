library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(igraph)
library(viridis)
library(SeuratWrappers)

################# E/P ######################
obj = readRDS("../00.data/EP_FB_fascia_myo_trajectory.rds")
#Idents(obj)=obj$celltype_0519
data <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation) 
# 预处理：标准化和PCA降维                          
cds <- preprocess_cds(cds, num_dim = 50) # preprocess_cds函数相当于seurat中NormalizeData + ScaleData + RunPCA
plot_pc_variance_explained(cds)
# 可视化
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds)                                
cds <- cluster_cells(cds, reduction_method = "UMAP")

fData(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))

head(fData(cds))
head(counts(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
head(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <-obj@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-obj@reductions$umap@cell.embeddings
data = cds
cds = data
cds <- cluster_cells(cds, reduction_method = "UMAP",resolution = 0.0005,k=50) #
cds <- learn_graph(cds, use_partition = T) #,learn_graph_control = list(minimal_branch_len = 20,nn.k=20)

## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin=c('FB_fascia_1_Adult')){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "celltype_0519"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}

nodes_vec <- c(get_earliest_principal_node(cds,"FB_fascia_EP"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")

p1 = plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = FALSE,
  group_cells_by = 'cluster',cell_size = 1,label_leaves=FALSE) + scale_color_viridis(option = "D")

p1 = p1+theme_void()+
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        plot.title = element_text(color = "white"))
#+theme(legend.position="top",legend.key.size=unit(1,'mm'),legend.key.width=unit(1,'mm'))
p1
ggsave(plot=p1,"Myo_trac_EP_All_pseudotime.pdf",w=4,h=3.5)

df_color = read.csv("../Fig2_color_seq_20240521.csv",row.names=1)
color = df_color$colors
names(color)=rownames(df_color)
p2 = plot_cells(cds,
                color_cells_by = "celltype_0519",
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,cell_size = 1, label_cell_groups=FALSE)+theme_void()+
  scale_color_manual(values = color)+
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        plot.title = element_text(color = "white"))
p2
ggsave(plot=p2,"Myo_trac_EP_All_celltype.pdf",w=5,h=3.5)

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8)

Track_genes <- Track_genes[,c(2,3,4,5,6)] %>% filter(q_value < 1e-3) %>% filter(abs(morans_I) > 0.05)
Track_genes <- Track_genes[!(grepl("^Mt",rownames(Track_genes))),]
Track_genes <- Track_genes[!(grepl("^Hb",rownames(Track_genes))),]
Track_genes <- Track_genes[!(grepl("^ENSR",rownames(Track_genes))),]

write.csv(Track_genes, "Myo_trac_EP_All_trajectory_genes_MI_05.csv", row.names = F)

################# Adult ######################
obj = readRDS("../Adult_FB_fascia_myo_trajectory.rds")
#DimPlot(obj,group.by = "celltype_240516")

data <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation) 
# 预处理：标准化和PCA降维                          
cds <- preprocess_cds(cds, num_dim = 50) # preprocess_cds函数相当于seurat中NormalizeData + ScaleData + RunPCA
plot_pc_variance_explained(cds)
# 可视化
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds)                                
cds <- cluster_cells(cds, reduction_method = "UMAP")

fData(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))

head(fData(cds))
head(counts(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
head(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <-obj@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-obj@reductions$umap@cell.embeddings

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin=c('FB_fascia_1_Adult')){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "celltype_240516"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}

nodes_vec <- c(get_earliest_principal_node(cds,"FB_fascia_1_Adult"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")

p1 = plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = FALSE,
  group_cells_by = 'cluster',cell_size = 1,label_leaves=FALSE) + scale_color_viridis(option = "D")

p1 = p1+theme_void()+
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        plot.title = element_text(color = "white"))
#+theme(legend.position="top",legend.key.size=unit(1,'mm'),legend.key.width=unit(1,'mm'))
p1
ggsave(plot=p1,"Myo_trac_Adult_All_pseudotime.pdf",w=4,h=3.5)
df_color = read.csv("../Fig2_color_seq_20240521.csv",row.names=1)
color = df_color$colors
names(color)=rownames(df_color)
p2 = plot_cells(cds,
                color_cells_by = "celltype_240516",
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,show_trajectory_graph = FALSE,
                label_branch_points=FALSE,cell_size = 1, label_cell_groups=FALSE)+theme_void()+
  scale_color_manual(values = color)+
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        plot.title = element_text(color = "white"))
p2
ggsave(plot=p2,"Myo_trac_Adult_All_celltype.pdf",w=5,h=3.5)

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8)

Track_genes <- Track_genes[,c(2,3,4,5,6)] %>% filter(q_value < 1e-3) %>% filter(abs(morans_I) > 0.05)
Track_genes <- Track_genes[!(grepl("^Mt",rownames(Track_genes))),]
Track_genes <- Track_genes[!(grepl("^Hb",rownames(Track_genes))),]
Track_genes <- Track_genes[!(grepl("^ENSR",rownames(Track_genes))),]

write.csv(Track_genes, "Myo_trac_Adult_All_trajectory_genes_MI_05.csv", row.names = F)
