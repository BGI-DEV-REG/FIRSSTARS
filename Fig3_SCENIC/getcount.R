# /hwfssz1/ST_SUPERCELLS/P21Z10200N0090/lijinxiu/04.software/anaconda3/envs/r4-base/bin/R
library(Seurat)

obj_adult = readRDS("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/Adult_FB_240414/240516_anno/Adult_FB_fascia_myo_trajectory.rds")
obj_ep = readRDS("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/EP_FB_240414/FB_subset/240519_anno/EP_FB_fascia_myo_trajectory.rds")

df = read.csv("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/Adult_FB_240414/240516_anno/Monocle3_trac/Myo_trac_EP_All_240526_trajectory_genes_MI_05.csv")
genes = df$gene_short_name
df = read.csv("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/Adult_FB_240414/240516_anno/Monocle3_trac/Myo_trac_Adult_All_240526_trajectory_genes_MI_05.csv")
genes = c(genes,df$gene_short_name)

genes = unique(genes)


obj = merge(obj_adult[genes,], obj_ep[genes,])
assay="RNA"
write.csv(t(as.matrix(obj@assays[[assay]]@counts)),file = gzfile('myo_trace_genes_count.csv.gz'),quote=F)

obj$celltype = c(obj_adult$celltype_240516, obj_ep$celltype_0519)
write.table(obj@meta.data,file = 'myo_trace_metadata.xls',sep='\t',quote=F)
