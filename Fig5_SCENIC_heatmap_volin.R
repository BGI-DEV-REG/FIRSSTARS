library(ggplot2)
library(Matrix)
library(Seurat)
library(stringr)
library(ggsignif)
library(ggsci)
library(latex2exp)
library(gghalves)
library(cowplot)
library(Matrix)
library(ggpubr)
library(minerva)
library(energy)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
setwd("/Users/zhaolab3/Documents/skin_immune")

wd = "./Rdata/SCENIC_All/"
celltype <- "celltype"
assay <- "RNA"
loom <- open_loom(paste0(wd,"aucell.loom"))

setwd("/Users/zhaolab3/Documents/skin_immune")
meta = read.table("./table/Adult&P5_MNP_Track_meta.csv",sep = "\t")
# cellinfo <- meta[,c(celltype,"nFeature_RNA" ,"nCount_RNA")]
# colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
# cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype')g)
# selectedResolution <- "celltype"

#celltype <- ct.col
cellsPerGroup <- split(rownames(meta),meta[,celltype])
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

library(Matrix)
sparse.gbm <- Matrix(getAUC(regulonAUC), sparse = T )
write(x = sparse.gbm@Dimnames[[1]], file = "features.tsv")
writeMM(obj = sparse.gbm, file="./Rdata/SCENIC_All/tf_matrix.mtx")


tf =c("Runx3", "Zfp36l1", "Ahr", "Ppard", "Nfkb1", "Aebp2", 
     "Zfp60", "Nr4a3", "Hmga1", "Jund", "Zbtb1", "Tcf4", 
     "Ikzf4", "Rbpj", "Erf", "Crem", "Nme2", "Etv5", "Zfp707", 
     "Cebpb", "Giot1", "Arid5a", "Meis3", "Tsc22d3", "Klf10", 
     "Zfp35", "Mbd2", "Irf7", "Bhlhe41", "Ddit3")

genes_df <-read.csv("./table/Adult_Mmp14_wot_gene.csv")
fare_df = read.csv("./table/Adult_Mmp14_wot_2_5_fate.csv",row.names = 'idx')
realtime = fare_df$fate
names(realtime) = rownames(fare_df)

rownames(regulonAUC) =  gsub("\\([+]\\)","",rownames(regulonAUC))

genes<-c(head(genes_df$x,300),tail(genes_df$x,300))
diff_tf = read.csv("./Rdata/SCENIC_all/tf_diff_res.csv")
tfAUC = regulonAUC[intersect(diff_tf$X,rownames(regulonAUC)),]
library(monocle)
cds =readRDS("./Rdata/Adult&P5_Immu_cds_important_37th.rds")
pt.matrix <-tfAUC[,rownames(fare_df)[order(realtime)]]
pt.matrix = pt.matrix@assays@data$AUC

idx = c()
for(i in 1:dim(pt.matrix)[2]){
  idx = c(idx,i%/%200)
}

library(Matrix.utils)
pt.matrix <- aggregate.Matrix(t(pt.matrix), 
                              groupings = idx, fun = "means")%>%t()
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- gsub("\\([+]\\)","",rownames(tfAUC@assays@data$AUC))
quantile(pt.matrix,na.rm=T)
#pt.matrix[which(is.na(pt.matrix))]=-2.8243446
tmp=pt.matrix%>%na.omit()

h=pheatmap::pheatmap(tmp,#[rnl_order,],
                     show_rownames = T,
                     show_colnames = F,
                     cluster_cols=F,
                     cluster_rows=F,
                     fontsize=5,
                     display_numbers=F,
                     fontsize_row=5,
                     border_color =NA,
                     na_col = "grey",
                     main ="Modul scores",
                     color = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
                     #color = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
                     bg='black')
dev.off()


no_Cd93_score<- c("Bach1", "Irf5",  
                 "Atf4", "Stat2", "Ehf",  
                 "Stat4", "Klf4", "Elf2",  
                 "Jun", "Pbx1",  "Fos",  
                 "Zfp287", "Irf7",  "Runx1", "Hoxd8", "Zfp143", "Smad4", "Sp1",  
                 "Hsf1",  
                 "Fli1",  
                 "Irf6",  
                 "Lhx2",  
                 "Zeb1",  
                 "Irf2",  
                 "Nfatc4", "Sox18", "Zbtb25",  
                 "Hoxa3",  
                 "Hoxa5",  
                 "Msx2",  
                 "Tfap2c", "Max",  
                 "Nr3c1",  
                 "Ikzf1",  
                 "Gabpa", "Zfp148", "Pparg",  
                 "Bhlhe41")

Cd93_score <- c("Gata3",  
                 "Jund",  
                 "Sox12",  
                 "Etv1",  
                 "E2f3",  
                 "Pou2f1", "Mitf",  
                 "Etv5",  
                 "Zbtb44", "Rxra",  
                 "Gata2",  
                 "Gabpb1", "Hoxb7",  
                 "Smad1", "Pura",  
                 "Rfx5",  
                 "Irf3",  
                 "Shox2", "Rfxank",  
                 "Msx1",  
                 "Mzf1",  
                 "Prrx1",  
                 "Gata6", "Klf12",  
                 "Sox5",  
                 "Zfp655", "Zfp182", "Nfatc1", "Foxj2", "Dlx1",  
                 "Maf",  
                 "Rarg", "Nfic",  
                 "Dlx2", "Fosb", "Vezf1", "Irf4",  
                 "Mafg",  
                 "Bcl6b", "Yy1",  
                 "Stat5a", "Ets1")  
Cd93_score = rownames(tmp)[1:40]
non_Cd93_score = rownames(tmp)[76:116]

setwd("/Users/zhaolab3/Documents/skin_immune")
meta = read.table("./table/Adult&P5_MNP_Track_meta.csv",sep = "\t")
cellsPerGroup <- split(rownames(meta),meta[,celltype])
meta = meta[meta$TimePostInjury%in%c("6H","1D","2D"),]
AgePerGroup <- split(rownames(meta),meta[,"Age"])



regulonActivity_byAge <- colSums(getAUC(tfAUC)[Cd93_score ,c(rownames(meta)[meta$Age=="Adult"],rownames(meta)[meta$Age=="P5"])])
regulonActivity_byAge <- colSums(getAUC(tfAUC)[no_Cd93_score ,c(rownames(meta)[meta$Age=="Adult"],rownames(meta)[meta$Age=="P5"])])
wilcox.test(regulonActivity_byAge[rownames(meta)[meta$Age=="Adult"]],
            regulonActivity_byAge[rownames(meta)[meta$Age=="P5"]])
quantile(regulonActivity_byAge[rownames(meta)[meta$Age=="Adult"]])
quantile(regulonActivity_byAge[rownames(meta)[meta$Age=="P5"]])

regulonActivity_byAge <- sapply(AgePerGroup,
                                function(cells) rowMeans(getAUC(tfAUC)[,cells]))
wilcox.test(regulonActivity_byAge[Cd93_score,"Adult"],regulonActivity_byAge[Cd93_score,"P5"],paired = T)

Cd93_score[(regulonActivity_byAge[Cd93_score,"Adult"]-regulonActivity_byAge[Cd93_score,"P5"])<(0)]
no_Cd93_score[(regulonActivity_byAge[no_Cd93_score,"Adult"]-regulonActivity_byAge[no_Cd93_score,"P5"])>(0)]
df = data.frame(gene = c(Cd93_score,Cd93_score,no_Cd93_score,no_Cd93_score),
                Age = c(rep("Adult",length(Cd93_score)),rep("P5",length(Cd93_score)),
                        rep("Adult",length(no_Cd93_score)),rep("P5",length(no_Cd93_score))),
                exp = c(regulonActivity_byAge[Cd93_score,"Adult"],regulonActivity_byAge[Cd93_score,"P5"],
                        regulonActivity_byAge[no_Cd93_score,"Adult"],regulonActivity_byAge[no_Cd93_score,"P5"]))

df$Age = factor(df$Age,levels = c("P5","Adult"))
pdf(file="./figures/Cd93_average_regulon_Activity.pdf",width=3,height=4)
ggplot(df[1:84,],aes(x = Age,y = exp)) +
  geom_violin(aes(color = Age),fill=NA,width=0.5)+
  geom_point(aes(color = Age),size = 1,shape=21) +  #绘制散点
  geom_line(aes(group = gene), color = 'lightgray', lwd = 0.5,alpha=0.3)+
  xlab('') + ylab("Regulon Activity")+ylim(c(0,0.075))+
  theme_black()+scale_color_manual(values = c("#26c6da","#f9320c"))+
  stat_compare_means(method = "wilcox.test",paired = T, 
                     comparisons=list(c("Adult",
                                        "P5")),color="white")
dev.off()

pdf(file="./figures/no_Cd93_average_regulon_Activity.pdf",width=3,height=4)
ggplot(df[85:nrow(df),],aes(x = Age,y = exp)) +
  geom_violin(aes(color = Age),fill=NA,width=0.5)+
  geom_point(aes(color = Age),size = 1,shape=21) +  #绘制散点
  geom_line(aes(group = gene), color = 'lightgray', lwd = 0.5,alpha=0.3)+
  xlab('') + ylab("Regulon Activity")+ylim(c(0,0.065))+
  theme_black()+scale_color_manual(values = c("#26c6da","#f9320c"))+
  stat_compare_means(method = "wilcox.test",paired = T, 
                     comparisons=list(c("Adult",
                                        "P5")),color="white")
dev.off()