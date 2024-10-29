### Get the parameters
parser <- argparse::ArgumentParser(description=cat("scRNA-seq data integration\n"))
parser$add_argument('-P','--path', help='Full path of object')
parser$add_argument('-N','--name',help='The prefix of integrated object')
parser$add_argument('-O','--out', help='Output directory')
parser$add_argument('-R','--resolution',type="double",default=0.8,help='clustering resolution [default %(default)s]')
parser$add_argument('-G','--group',default="AnalysisName",help='group by [default %(default)s]')

args <- parser$parse_args()
print(args)
print(Sys.time())
dir.create(args$out)
setwd(args$out)

library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(patchwork)

merge_obj = readRDS(args$path)

# cellcycle
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)



merge_obj = merge_obj %>%
    NormalizeData %>%
    FindVariableFeatures(nfeatures=3000) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>%
    RunPCA %>%
    RunHarmony(group.by.vars = c(args$group, "Phase")) %>%
    RunUMAP(dims=1:30, reduction = "harmony") %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution=args$resolution)

saveRDS(merge_obj, paste0(args$name, "_harmony.rds"))

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



plots <- function(obj, draw_obj, group_name, prefix, ...){
    colors = colorlist[1:length(unique(obj@meta.data[[group_name]]))]
    names(colors) = unique(obj@meta.data[[group_name]])
    p <- DimPlot(obj, group.by=group_name, pt.size=0.1, cols=colors, ...) + NoLegend()
    ggsave(paste0(prefix, "_", group_name, "_Nolegend.pdf"),width=14,height=14)

    p <- DimPlot(obj, group.by=group_name, pt.size=0.1, ...)
    ggsave(paste0(prefix, "_", group_name, ".pdf"),width=14,height=14)

    pdf(paste0(prefix, "_", group_name, "_split.pdf"))
    for (i in sort(unique(draw_obj@meta.data[[group_name]]))){
        cells = Cells(draw_obj)[draw_obj@meta.data[[group_name]]==i]
        p <- DimPlot(draw_obj, group.by=group_name, cells.highlight=cells)
        print(p+ggtitle(i))
    }
    dev.off()
}

if (length(Cells(merge_obj))>50000) {
    draw_obj = subset(merge_obj, cells=sample(Cells(merge_obj), 50000))
} else {
    draw_obj = merge_obj
}

plots(merge_obj, draw_obj, "seurat_clusters", paste0(args$name, "_harmony"), label=T)
plots(merge_obj, draw_obj, "Age", paste0(args$name, "_harmony"))
plots(merge_obj, draw_obj, "Age_TimePostInjury", paste0(args$name, "_harmony"))
plots(merge_obj, draw_obj, "Phase", paste0(args$name, "_harmony"))
plots(merge_obj, draw_obj, "SampleIndex", paste0(args$name, "_harmony"))

print("Done!")
print(Sys.time())
