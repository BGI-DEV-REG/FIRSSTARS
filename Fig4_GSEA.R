rm(list = ls())
library(magrittr)
library(plyr)
library(stringr)
species ="Rattus norvegicus"
setwd('/Users/zhaolab3/Documents/skin_immune/')
immune_meta = read.csv("./re-analysis/immune_meta_umap_0517.csv",row.names = 'X')
immune_cols_df = read.csv("./re-analysis/Immu_colors.csv",row.names = 'X')
table_dir = '/Users/zhaolab3/Documents/skin_immune/All/'

runGSEA_go<-function(pbmc.genes,clusters,species){
  cluster.genes<- pbmc.genes %>% dplyr::filter(group == clusters) %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
  ranks<- deframe(cluster.genes)
  library(msigdbr)
  library(fgsea)
  m_df<- msigdbr(species = species, category = "C5",subcategory = 'GO:BP')
  head(m_df)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
  fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
  fgseaResTidy$clusters<-rep(clusters,nrow(fgseaResTidy))
  fgseaResTidy
  return(fgseaResTidy)
}

getGOname = function(fgseaResTidy,num){
  fgseaResTidy = fgseaResTidy[order(fgseaResTidy$padj,decreasing =F ),]
  return(fgseaResTidy$pathway[1:num])
}
FoldFunction<-function(results){
  library(stringr)
  gr1 <- as.numeric(str_split(results$GeneRatio,"/",simplify = T)[,1])
  gr2 <- as.numeric(str_split(results$GeneRatio,"/",simplify = T)[,2])
  bg1 <- as.numeric(str_split(results$BgRatio,"/",simplify = T)[,1])
  bg2 <- as.numeric(str_split(results$BgRatio,"/",simplify = T)[,2])
  results$fold <- (gr1/gr2)/(bg1/bg2)
  results$GeneRatio <- (gr1/gr2)
  return(results)
}

adata = readRDS("./Rdata/Immune_Rat2Mouse_MNP_raw_counts.rds")
# homolo = read.csv("./re-analysis/Rat2Mouse_homologene.csv",row.names = 1)
# rat_genes = rownames(adata@assays$RNA@counts)
# homolo=homolo[homolo$Var_names%in%rat_genes,]
# rownames(homolo)=homolo$Var_names
# homolo = homolo[rownames(adata),]
# du_gene = homolo[duplicated(homolo$new_symbol),]
# homolo$new_symbol[homolo$new_symbol%in%du_gene$new_symbol]=homolo$Var_names[homolo$new_symbol%in%du_gene$new_symbol]


immune_meta=immune_meta[rownames(adata@meta.data),]
#counts = adata@assays$RNA$counts
# tail(rownames(counts))
# rownames(counts) = homolo$new_symbol
# length(unique(homolo$new_symbol))
# adata <- CreateSeuratObject(counts = counts, project = "MNP")

adata@meta.data=immune_meta

adata$P5vsAdult = paste(adata$Age,adata$celltype,sep = "_")
adata =subset(adata,Age%in%c("Adult","P5"))
adata =subset(adata,celltype%in%c('Macr_Ms4a7', 'Macr_Trem2','Macr_Mmp14','Mono_Spp1',  'Mono_Ly6c'))

clusters<-adata@meta.data$P5vsAdult%>%as.factor()%>%levels()

head(adata@meta.data)
#pbmc = subset(adata, nFeature_RNA > 1000)
pbmc = Seurat::NormalizeData(adata)
library(presto)
library(tibble)
pbmc.genes <- wilcoxauc(pbmc,assay = "data", group_by='P5vsAdult')
head(pbmc.genes)
write.csv(pbmc.genes,file = "./re-analysis//P5vsAdult_MNP_wilcoxauc.csv",quote = F)

go.res1<-runGSEA_go(pbmc.genes,clusters[1],species)
go.res2<-runGSEA_go(pbmc.genes,clusters[2],species)
go.res3<-runGSEA_go(pbmc.genes,clusters[3],species)
go.res4<-runGSEA_go(pbmc.genes,clusters[4],species)
go.res5<-runGSEA_go(pbmc.genes,clusters[5],species)
go.res6<-runGSEA_go(pbmc.genes,clusters[6],species)

#kegg.res8<-runGSEA_kegg(pbmc.genes,clusters[8],species)

go_list = list(go.res1,go.res2,go.res3,go.res4,go.res5,go.res6)
go_names = Reduce(union,lapply(go_list, function(fgseaResTidy){
  fgseaResTidy = fgseaResTidy[order(fgseaResTidy$padj,decreasing =F ),]
  return(fgseaResTidy$pathway[1:50])
}))
select_data = read.csv("./table/GSEA_selected.csv",row.names = 'X')

df<-rbind(go.res1,go.res2,go.res3,go.res4,go.res5,go.res6)

data<-df[,c(1,5,9)]%>%as.data.frame()
unique(data$pathway[grep(paste0("*",toupper("angiogenesis"),"*"),data$pathway)])
data = data[data$pathway%in%go_nanems,]
head(data)
matBP<-data
matBP$pathway<-lapply(matBP$pathway,function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()

a<-gsub("In","in",matBP$pathway)
a<-gsub("Of","of",a)
matBP$pathway<-gsub("To","to",a)
matBP = matBP[matBP$pathway%in%rownames(select_data),]
library(reshape2)
matBP<-melt(matBP,id.vars = c("clusters", "pathway"))
matBP<-dcast(matBP, pathway ~ clusters, na.rm = TRUE)
rownames(matBP)<-matBP$pathway
matBP<-matBP[,-1]
matBP<-as.matrix(matBP)
write.csv(data,file = "./re-analysis/GSEA_NES_top50.csv",quote = F)
colnames(matBP) = c("Macr_resident","Macr_neovascularization",
                    "Mono_noclassical","Macr_remodeling",
                    "Mono_classical:", "Mono_intermediate")
data<-df[,c(1,3,9)]%>%as.data.frame()
data = data[data$pathway%in%go_nanems,]
matpval<-data#[1:59,]
matpval$pathway<-lapply(matpval$pathway,function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()
a<-gsub("In","in",matpval$pathway)
a<-gsub("Of","of",a)
matpval$pathway<-gsub("To","to",a)
matpval = matpval[matpval$pathway%in%rownames(select_data),]
library(reshape2)
matpval<-melt(matpval,id.vars = c("clusters", "pathway"))
matpval<-dcast(matpval, pathway ~ clusters, na.rm = TRUE)
rownames(matpval)<-matpval$pathway
matpval<-matpval[,-1]
matpval<-as.matrix(matpval)
colnames(matpval) = colnames(matBP)

library(pheatmap)
h=pheatmap::pheatmap(matBP,scale = "row",
                     number_color="gray",
                     display_numbers = matrix(ifelse(matpval < 0.001, "***", ifelse(
                       matpval < 0.01,"**",ifelse(
                         matpval < 0.05,"*","-"
                       )
                     )),nrow(matBP)),
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = T, cluster_col = T,
                     border_color =NA,#'black',
                     angle_col=45,fontsize=10,fontsize_col=10,fontsize_row=6,
                     color = viridis::inferno(100),
                     bg='black',
                     main ="GSEA")
dev.off()